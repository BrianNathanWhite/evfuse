#' Fit nonstationary GEV with covariate (Stage 1)
#'
#' NOAA sites get nonstationary GEV: mu(t) = mu0 + mu1 * (cov(t) - ref_value),
#' sigma and xi stationary. All four parameters (mu0, mu1, log_sigma, xi)
#' are returned for spatial modeling. ADCIRC sites get the usual stationary
#' fit with three parameters (mu0, log_sigma, xi).
#'
#' @param dat An \code{evfuse_data} object.
#' @param df The raw data frame containing the covariate column.
#' @param covariate String column name in \code{df} (e.g., \code{"sst_regional_warm"}).
#' @param ref_value Centering value for the covariate. Default: mean across
#'   all NOAA site-years. Centering improves numerical stability and makes mu0
#'   interpretable as the location at the reference covariate value.
#' @seealso \code{\link{fit_gev_detrended}} for the time-trend version that
#'   discards mu1 before spatial modeling.
#' @return An \code{evfuse_stage1} object with 4-column \code{theta_hat}.
#'   NOAA rows contain (mu0, mu1, log_sigma, xi); ADCIRC rows contain
#'   (mu0, log_sigma, xi, NA). The \code{vcov_list} entries are 4x4 for NOAA
#'   and 3x3 for ADCIRC.
#' @export
fit_gev_ns <- function(dat, df, covariate, ref_value = NULL) {
  n <- dat$n_sites

  # Default ref_value: mean of covariate across all NOAA site-years
  if (is.null(ref_value)) {
    noaa_locs <- dat$sites$location[dat$sites$data_source == "NOAA"]
    noaa_rows <- as.character(df$location) %in% noaa_locs
    ref_value <- mean(df[[covariate]][noaa_rows], na.rm = TRUE)
  }

  theta_hat <- matrix(NA_real_, nrow = n, ncol = 4)
  colnames(theta_hat) <- c("mu0", "mu1", "log_sigma", "xi")
  vcov_list <- vector("list", n)
  fits <- vector("list", n)
  converged <- logical(n)

  for (i in seq_len(n)) {
    loc_name <- dat$sites$location[i]
    is_noaa <- dat$sites$data_source[i] == "NOAA"

    if (is_noaa) {
      # Nonstationary GEV with covariate in mu
      site_df <- df[as.character(df$location) == loc_name, ]
      site_df <- site_df[order(site_df$year), ]
      site_df$cov_c <- site_df[[covariate]] - ref_value

      fit <- tryCatch(
        extRemes::fevd(site_df$max_sea_level, data = site_df,
                       type = "GEV", location.fun = ~cov_c),
        error = function(e) NULL
      )

      if (is.null(fit) || fit$results$convergence != 0) {
        warning(sprintf("Nonstationary GEV failed at NOAA site '%s'.", loc_name))
        converged[i] <- FALSE
        fits[[i]] <- fit
        next
      }

      converged[i] <- TRUE
      fits[[i]] <- fit

      pars <- fit$results$par
      mu0   <- pars["mu0"]
      mu1   <- pars["mu1"]
      sigma <- pars["scale"]
      xi    <- pars["shape"]

      # Vcov via delta method: (mu0, mu1, sigma, xi) -> (mu0, mu1, log_sigma, xi)
      # Jacobian J (4x4): identity except J[3,3] = 1/sigma
      V_raw <- tryCatch(solve(fit$results$hessian), error = function(e) NULL)
      if (is.null(V_raw) || any(diag(V_raw) < 0)) {
        warning(sprintf("Hessian inversion failed at NOAA site '%s'.", loc_name))
        converged[i] <- FALSE
        next
      }

      J <- diag(4)
      J[3, 3] <- 1 / sigma   # d(log_sigma)/d(sigma)
      V_transformed <- J %*% V_raw %*% t(J)

      theta_hat[i, ] <- c(mu0, mu1, log(sigma), xi)
      vcov_list[[i]] <- V_transformed

    } else {
      # ADCIRC: stationary GEV (same as fit_gev_all)
      x <- dat$maxima[[i]]
      fit <- tryCatch(
        extRemes::fevd(x, type = "GEV", method = "MLE"),
        error = function(e) NULL
      )

      if (is.null(fit) || fit$results$convergence != 0) {
        warning(sprintf("GEV fit failed at ADCIRC site '%s'.", loc_name))
        converged[i] <- FALSE
        fits[[i]] <- fit
        next
      }

      converged[i] <- TRUE
      fits[[i]] <- fit

      pars <- fit$results$par
      mu    <- pars["location"]
      sigma <- pars["scale"]
      xi    <- pars["shape"]

      V_raw <- tryCatch(solve(fit$results$hessian), error = function(e) NULL)
      if (is.null(V_raw) || any(diag(V_raw) < 0)) {
        warning(sprintf("Hessian inversion failed at ADCIRC site '%s'.", loc_name))
        converged[i] <- FALSE
        next
      }

      J <- diag(c(1, 1 / sigma, 1))
      # Store in first 3 columns; column 4 (mu1) stays NA
      theta_hat[i, 1:3] <- c(mu, log(sigma), xi)
      vcov_list[[i]] <- J %*% V_raw %*% t(J)
    }
  }

  names(fits) <- dat$sites$location
  names(vcov_list) <- dat$sites$location

  if (any(!converged)) {
    warning(sprintf("%d of %d sites failed to converge.", sum(!converged), n))
  }

  structure(
    list(
      theta_hat = theta_hat,
      vcov_list = vcov_list,
      fits = fits,
      converged = converged,
      log_scale = TRUE,
      nonstationary = TRUE,
      covariate = covariate,
      ref_value = ref_value
    ),
    class = "evfuse_stage1"
  )
}

#' Bootstrap W with nonstationary covariate
#'
#' Like \code{\link{bootstrap_W_detrended}} but uses an arbitrary covariate
#' (e.g., SST) instead of a linear time trend, and retains the covariate
#' sensitivity mu1 in the bootstrap parameter vector. NOAA sites are fitted
#' with nonstationary GEV; ADCIRC sites use stationary GEV.
#'
#' The bootstrap Gamma matrix has \code{L * 4} columns in parameter-major
#' ordering: (mu0, mu1, log_sigma, xi) at all sites. ADCIRC mu1 entries are
#' \code{NA}; the resulting \code{W_bs} has \code{NA} in those positions.
#' This is by design: \code{\link{embed_W}} skips unobserved entries.
#'
#' @param dat An \code{evfuse_data} object.
#' @param df The raw data frame with covariate column.
#' @param covariate String column name in \code{df}.
#' @param B Number of bootstrap replications.
#' @param ref_value Centering value (default: mean across NOAA site-years).
#' @param seed Random seed.
#' @seealso \code{\link{bootstrap_W_detrended}} for the time-trend version.
#' @return A list with \code{W_bs} (4L x 4L), \code{Gamma} (B x 4L),
#'   and \code{n_failures}.
#' @note A small fraction of bootstrap resamples may produce degenerate data
#'   that causes warnings or convergence failures in the underlying GEV fitting
#'   routine (\code{extRemes::fevd}). This is expected behavior: failed fits are
#'   recorded in \code{n_failures} and handled via pairwise-complete covariance
#'   estimation. Typical failure rates are well below 1 percent.
#' @export
bootstrap_W_ns <- function(dat, df, covariate, B = 500,
                            ref_value = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  n <- dat$n_sites

  # Default ref_value
  if (is.null(ref_value)) {
    noaa_locs <- dat$sites$location[dat$sites$data_source == "NOAA"]
    noaa_rows <- as.character(df$location) %in% noaa_locs
    ref_value <- mean(df[[covariate]][noaa_rows], na.rm = TRUE)
  }

  # Build per-site data
  site_years <- vector("list", n)
  site_maxima <- vector("list", n)
  site_covariate <- vector("list", n)

  for (i in seq_len(n)) {
    loc <- dat$sites$location[i]
    site_df <- df[as.character(df$location) == loc, ]
    site_df <- site_df[order(site_df$year), ]
    site_years[[i]] <- site_df$year
    site_maxima[[i]] <- site_df$max_sea_level
    if (covariate %in% names(site_df)) {
      site_covariate[[i]] <- site_df[[covariate]]
    }
  }

  # Per-source year counts
  noaa_idx <- which(dat$sites$data_source == "NOAA")
  adcirc_idx <- which(dat$sites$data_source == "ADCIRC")

  T_noaa <- unique(lengths(site_maxima[noaa_idx]))
  T_adcirc <- unique(lengths(site_maxima[adcirc_idx]))
  if (length(T_noaa) > 1) {
    warning("NOAA sites have differing year counts; using minimum.")
    T_noaa <- min(lengths(site_maxima[noaa_idx]))
  }
  if (length(T_adcirc) > 1) {
    warning("ADCIRC sites have differing year counts; using minimum.")
    T_adcirc <- min(lengths(site_maxima[adcirc_idx]))
  }

  # 4 params per site in parameter-major ordering:
  # positions 1..L = mu0, L+1..2L = mu1, 2L+1..3L = log_sigma, 3L+1..4L = xi
  p_bs <- 4
  Lp <- n * p_bs
  Gamma <- matrix(NA_real_, nrow = B, ncol = Lp)
  n_failures <- 0

  message(sprintf("Running %d bootstrap replicates (%d total GEV fits). ",
                  B, B * n),
          "Occasional warnings from GEV fitting are expected and handled gracefully.")

  for (b in seq_len(B)) {
    noaa_yidx <- sample.int(T_noaa, replace = TRUE)
    adcirc_yidx <- sample.int(T_adcirc, replace = TRUE)
    theta_b <- rep(NA_real_, Lp)

    for (i in seq_len(n)) {
      is_noaa <- dat$sites$data_source[i] == "NOAA"
      yidx <- if (is_noaa) noaa_yidx else adcirc_yidx
      x_boot <- site_maxima[[i]][yidx]

      if (is_noaa) {
        # Resample covariate alongside maxima (preserving temporal pairing)
        cov_boot <- site_covariate[[i]][yidx]
        d_boot <- data.frame(x = x_boot, cov_c = cov_boot - ref_value)
        fit <- tryCatch(
          extRemes::fevd(d_boot$x, data = d_boot, type = "GEV",
                         location.fun = ~cov_c),
          error = function(e) NULL
        )
        if (is.null(fit) || fit$results$convergence != 0) {
          n_failures <- n_failures + 1
          next
        }
        pars <- fit$results$par
        theta_b[i]         <- pars["mu0"]
        theta_b[n + i]     <- pars["mu1"]
        theta_b[2 * n + i] <- log(pars["scale"])
        theta_b[3 * n + i] <- pars["shape"]
      } else {
        # ADCIRC: stationary; mu1 position (n + i) stays NA
        fit <- tryCatch(
          extRemes::fevd(x_boot, type = "GEV", method = "MLE"),
          error = function(e) NULL
        )
        if (is.null(fit) || fit$results$convergence != 0) {
          n_failures <- n_failures + 1
          next
        }
        pars <- fit$results$par
        theta_b[i]         <- pars["location"]
        theta_b[2 * n + i] <- log(pars["scale"])
        theta_b[3 * n + i] <- pars["shape"]
      }
    }

    Gamma[b, ] <- theta_b
    if (b %% 100 == 0) cat(sprintf("  Bootstrap %d / %d\n", b, B))
  }

  W_bs <- cov(Gamma, use = "pairwise.complete.obs")

  if (n_failures > 0) {
    message(sprintf("Bootstrap: %d failures out of %d total fits (%.1f%%)",
                    n_failures, B * n, 100 * n_failures / (B * n)))
  }

  list(W_bs = W_bs, Gamma = Gamma, n_failures = n_failures)
}
