#' Fit GEV at all sites with NOAA detrending (Stage 1, nonstationary)
#'
#' NOAA sites get nonstationary GEV: mu(t) = mu0 + mu1*(t - ref_year),
#' sigma and xi stationary. The detrended triplet (mu0, log_sigma, xi)
#' is returned. ADCIRC sites get the usual stationary fit.
#' This is the primary Stage 1 approach used in the manuscript.
#'
#' @param dat An \code{evfuse_data} object.
#' @param df The raw data frame with year column.
#' @param ref_year Reference year for centering (default 2000).
#' @seealso \code{\link{fit_gev_all}} for the stationary version.
#' @return An \code{evfuse_stage1} object with detrended NOAA parameters.
#' @export
fit_gev_detrended <- function(dat, df, ref_year = 2000) {
  n <- dat$n_sites
  theta_hat <- matrix(NA_real_, nrow = n, ncol = 3)
  colnames(theta_hat) <- c("mu", "log_sigma", "xi")
  vcov_list <- vector("list", n)
  fits <- vector("list", n)
  converged <- logical(n)

  for (i in seq_len(n)) {
    loc_name <- dat$sites$location[i]
    is_noaa <- dat$sites$data_source[i] == "NOAA"

    if (is_noaa) {
      # Nonstationary fit with linear trend in mu
      site_df <- df[as.character(df$location) == loc_name, ]
      site_df <- site_df[order(site_df$year), ]
      site_df$year_c <- site_df$year - ref_year

      fit <- tryCatch(
        extRemes::fevd(site_df$max_sea_level, data = site_df,
                       type = "GEV", location.fun = ~year_c),
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

      # Parameters: mu0, mu1, scale, shape
      pars <- fit$results$par
      mu0   <- pars["mu0"]
      sigma <- pars["scale"]
      xi    <- pars["shape"]

      # Vcov via delta method: project (mu0, mu1, sigma, xi) -> (mu0, log_sigma, xi)
      # Jacobian J (3x4):
      #   row 1 (mu0):       [1, 0, 0, 0]
      #   row 2 (log_sigma): [0, 0, 1/sigma, 0]
      #   row 3 (xi):        [0, 0, 0, 1]
      V_raw <- tryCatch(solve(fit$results$hessian), error = function(e) NULL)
      if (is.null(V_raw) || any(diag(V_raw) < 0)) {
        warning(sprintf("Hessian inversion failed at NOAA site '%s'.", loc_name))
        converged[i] <- FALSE
        next
      }

      J <- matrix(0, 3, 4)
      J[1, 1] <- 1           # d(mu0)/d(mu0)
      J[2, 3] <- 1 / sigma   # d(log_sigma)/d(sigma)
      J[3, 4] <- 1           # d(xi)/d(xi)
      V_transformed <- J %*% V_raw %*% t(J)

      theta_hat[i, ] <- c(mu0, log(sigma), xi)
      vcov_list[[i]] <- V_transformed

    } else {
      # ADCIRC: standard stationary fit (same as fit_gev_all)
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
      theta_hat[i, ] <- c(mu, log(sigma), xi)
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
      detrended = TRUE,
      ref_year = ref_year
    ),
    class = "evfuse_stage1"
  )
}

#' Bootstrap W with NOAA detrending
#'
#' Like \code{\link{bootstrap_W}} but NOAA sites are fitted with nonstationary
#' GEV (linear mu trend), extracting (mu0, log_sigma, xi) per replicate.
#' ADCIRC sites use the standard stationary GEV.
#' This is the primary bootstrap approach used in the manuscript.
#'
#' @param dat An \code{evfuse_data} object.
#' @param df The raw data frame with year column.
#' @param B Number of bootstrap replications.
#' @param ref_year Reference year for centering.
#' @param seed Random seed.
#' @seealso \code{\link{bootstrap_W}} for the stationary version.
#' @return A list with W_bs, Gamma, and n_failures.
#' @note Some bootstrap resamples will produce degenerate data that fails
#'   GEV fitting (convergence warnings from \code{extRemes::fevd}). These
#'   are recorded in \code{n_failures} and excluded via pairwise-complete
#'   covariance estimation. Failure rates are typically well under 1 percent.
#' @export
bootstrap_W_detrended <- function(dat, df, B = 500,
                                   ref_year = 2000, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  n <- dat$n_sites
  noaa_idx <- which(dat$sites$data_source == "NOAA")
  adcirc_idx <- which(dat$sites$data_source == "ADCIRC")

  # Build per-site year vectors and maxima aligned to dat ordering
  site_years <- vector("list", n)
  site_maxima <- vector("list", n)
  for (i in seq_len(n)) {
    loc <- dat$sites$location[i]
    site_df <- df[as.character(df$location) == loc, ]
    site_df <- site_df[order(site_df$year), ]
    site_years[[i]] <- site_df$year
    site_maxima[[i]] <- site_df$max_sea_level
  }

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

  Lp <- n * 3
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
        yr_boot <- site_years[[i]][yidx]
        d_boot <- data.frame(x = x_boot, year_c = yr_boot - ref_year)
        fit <- tryCatch(
          extRemes::fevd(d_boot$x, data = d_boot, type = "GEV",
                         location.fun = ~year_c),
          error = function(e) NULL
        )
        if (is.null(fit) || fit$results$convergence != 0) {
          n_failures <- n_failures + 1
          next
        }
        pars <- fit$results$par
        theta_i <- c(pars["mu0"], log(pars["scale"]), pars["shape"])
      } else {
        fit <- tryCatch(
          extRemes::fevd(x_boot, type = "GEV", method = "MLE"),
          error = function(e) NULL
        )
        if (is.null(fit) || fit$results$convergence != 0) {
          n_failures <- n_failures + 1
          next
        }
        pars <- fit$results$par
        theta_i <- c(pars["location"], log(pars["scale"]), pars["shape"])
      }

      for (j in 1:3) {
        theta_b[(j - 1) * n + i] <- theta_i[j]
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
