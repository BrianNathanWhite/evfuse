#' Estimate W via nonparametric block bootstrap (stationary)
#'
#' Implements the bootstrap procedure from Section 2.4 of Russell et al. (2020).
#' Resamples years (preserving spatial dependence) and re-fits stationary GEV at
#' each site, then estimates the covariance of the stage-1 estimators. For a
#' nonstationary alternative that accounts for NOAA sea level trends, see
#' \code{\link{bootstrap_W_detrended}}.
#'
#' @param dat An \code{evfuse_data} object.
#' @param B Number of bootstrap replications (default 500).
#' @param log_scale Logical. If TRUE, work with log(sigma) (default TRUE).
#' @param seed Random seed for reproducibility.
#' @seealso \code{\link{bootstrap_W_detrended}} for the nonstationary version.
#' @return A list with components:
#'   \describe{
#'     \item{W_bs}{The raw bootstrap covariance matrix (Lp x Lp).}
#'     \item{Gamma}{The bootstrap matrix of MLEs (B x Lp).}
#'     \item{n_failures}{Number of (site, bootstrap) combinations that failed.}
#'   }
#' @note A small fraction of bootstrap resamples may produce degenerate data
#'   that causes warnings or convergence failures in \code{extRemes::fevd}.
#'   Failed fits are recorded in \code{n_failures} and excluded via
#'   pairwise-complete covariance estimation (failure rates are typically
#'   well below 1 percent).
#' @export
bootstrap_W <- function(dat, B = 500, log_scale = TRUE, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  n <- dat$n_sites
  p <- 3  # GEV parameters per source; total columns will be n * p_total
  # But we work per-source: each site contributes 3 parameters
  # Theta ordering: (theta_1(s1), ..., theta_1(sL), theta_2(s1), ..., theta_3(sL))
  # i.e., parameter-major ordering

  # Resample years within each source to preserve spatial dependence.
  # All sites of the same source must have the same number of years.
  # Get number of years per site
  T_per_site <- lengths(dat$maxima)

  # Per-source site indices and year counts
  source_names <- if (!is.null(dat$source_params)) names(dat$source_params) else c("NOAA", "ADCIRC")
  source_idx <- lapply(source_names, function(s) which(dat$sites$data_source == s))
  names(source_idx) <- source_names

  T_source <- list()
  for (s in source_names) {
    T_s <- unique(T_per_site[source_idx[[s]]])
    if (length(T_s) > 1) {
      warning(s, " sites have differing numbers of years. ",
              "Using minimum for bootstrap resampling.")
      T_s <- min(T_per_site[source_idx[[s]]])
    }
    T_source[[s]] <- T_s
  }

  # Total parameter vector length: n_sites * 3 (mu, log_sigma, xi per site)
  Lp <- n * 3
  Gamma <- matrix(NA_real_, nrow = B, ncol = Lp)
  n_failures <- 0

  message(sprintf("Running %d bootstrap replicates (%d total GEV fits). ",
                  B, B * n),
          "Occasional warnings from GEV fitting are expected and will not affect results.")

  for (b in seq_len(B)) {
    # Resample years separately per source (they may have different time spans)
    year_idx <- lapply(source_names, function(s) sample.int(T_source[[s]], replace = TRUE))
    names(year_idx) <- source_names

    theta_b <- rep(NA_real_, Lp)

    for (i in seq_len(n)) {
      x <- dat$maxima[[i]]
      src <- dat$sites$data_source[i]
      x_boot <- x[year_idx[[src]]]

      fit <- tryCatch(
        extRemes::fevd(x_boot, type = "GEV", method = "MLE"),
        error = function(e) NULL
      )

      if (is.null(fit) || fit$results$convergence != 0) {
        n_failures <- n_failures + 1
        next
      }

      pars <- fit$results$par
      mu <- pars["location"]
      sigma <- pars["scale"]
      xi <- pars["shape"]

      if (log_scale) {
        theta_i <- c(mu, log(sigma), xi)
      } else {
        theta_i <- c(mu, sigma, xi)
      }

      # Parameter-major ordering: theta_j at all sites, then theta_{j+1}, ...
      # Position for parameter j at site i: (j-1)*n + i
      for (j in 1:3) {
        theta_b[(j - 1) * n + i] <- theta_i[j]
      }
    }

    Gamma[b, ] <- theta_b
  }

  # Compute sample covariance, handling NAs via pairwise complete obs
  W_bs <- cov(Gamma, use = "pairwise.complete.obs")

  if (n_failures > 0) {
    message(sprintf("Bootstrap: %d failures out of %d total fits (%.1f%%)",
                    n_failures, B * n, 100 * n_failures / (B * n)))
  }

  list(W_bs = W_bs, Gamma = Gamma, n_failures = n_failures)
}

#' Apply covariance tapering to W
#'
#' Computes W_tap = W_bs * T_tap (Hadamard product) where T_tap is built
#' from the Wendland 2 function with range lambda.
#'
#' @param W_bs Raw bootstrap covariance matrix (Lp x Lp).
#' @param D Distance matrix (L x L).
#' @param lambda Taper range in km.
#' @param p Number of GEV parameters per site (default 3).
#' @return Tapered covariance matrix W_tap (Lp x Lp).
#' @export
taper_W <- function(W_bs, D, lambda, p = 3) {
  T_tap <- build_taper(D, lambda, p = p)
  W_bs * T_tap
}

#' Embed W_tap into full parameter space
#'
#' The bootstrap produces a covariance matrix in observed-parameter space
#' (L sites x p_bs params per site). This function embeds it into the full
#' joint parameter space (L sites x p_full params) required by the stage 2 model.
#' Each source's observed params are mapped to their corresponding positions
#' in the joint model via \code{source_params}.
#'
#' When sources observe different numbers of parameters (e.g., NOAA has 4 and
#' ADCIRC has 3 in the nonstationary case), the bootstrap uses a union
#' parameter ordering with \code{NA} entries for unobserved positions.
#' This function detects unobserved positions from \code{NA} on the diagonal
#' of \code{W_tap} and maps only the observed entries.
#'
#' @param W_tap Tapered bootstrap covariance (L*p_bs x L*p_bs).
#' @param dat Data object with site information and \code{source_params}.
#' @param source_params Named list mapping source labels to full-model param
#'   indices. Default: \code{list(NOAA = 1:3, ADCIRC = 4:6)}.
#' @return Embedded covariance matrix (L*p_full x L*p_full).
#' @export
embed_W <- function(W_tap, dat, source_params = NULL) {
  L <- dat$n_sites
  if (is.null(source_params)) {
    source_params <- if (!is.null(dat$source_params)) dat$source_params
                     else list(NOAA = 1:3, ADCIRC = 4:6)
  }
  p_full <- max(unlist(source_params))
  p_bs <- nrow(W_tap) %/% L
  diag_W <- diag(W_tap)

  # For each source, determine which bootstrap param positions are observed
  # by checking the W_tap diagonal for NaN (NaN = no data for this param).
  w_indices <- c()
  f_indices <- c()

  for (src in names(source_params)) {
    sp <- source_params[[src]]
    src_sites <- which(dat$sites$data_source == src)

    # Check first site of this source to find observed bootstrap params
    # Unobserved entries are NA (from cov with all-NA columns)
    i0 <- src_sites[1]
    diag_entries <- diag_W[seq(i0, by = L, length.out = p_bs)]
    obs_bs <- which(!is.na(diag_entries))

    if (length(obs_bs) != length(sp)) {
      # Fallback: sequential mapping (backward compatible for uniform p_obs)
      obs_bs <- seq_along(sp)
    }

    for (j in seq_along(sp)) {
      bs_j <- obs_bs[j]
      full_j <- sp[j]
      for (i in src_sites) {
        w_indices <- c(w_indices, (bs_j - 1) * L + i)
        f_indices <- c(f_indices, (full_j - 1) * L + i)
      }
    }
  }

  W_full <- matrix(0, nrow = L * p_full, ncol = L * p_full)
  W_full[f_indices, f_indices] <- W_tap[w_indices, w_indices]
  W_full
}
