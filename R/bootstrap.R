#' Estimate W via nonparametric block bootstrap
#'
#' Implements the bootstrap procedure from Section 2.4 of Russell et al. (2019).
#' Resamples years (preserving spatial dependence) and re-fits GEV at each site,
#' then estimates the covariance of the stage-1 estimators.
#'
#' @param dat An \code{evfuse_data} object.
#' @param B Number of bootstrap replications (default 500).
#' @param log_scale Logical. If TRUE, work with log(sigma) (default TRUE).
#' @param seed Random seed for reproducibility.
#' @return A list with components:
#'   \describe{
#'     \item{W_bs}{The raw bootstrap covariance matrix (Lp x Lp).}
#'     \item{Gamma}{The bootstrap matrix of MLEs (B x Lp).}
#'     \item{n_failures}{Number of (site, bootstrap) combinations that failed.}
#'   }
#' @export
bootstrap_W <- function(dat, B = 500, log_scale = TRUE, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  n <- dat$n_sites
  p <- 3  # GEV parameters per source; total columns will be n * p_total
  # But we work per-source: each site contributes 3 parameters
  # Theta ordering: (theta_1(s1), ..., theta_1(sL), theta_2(s1), ..., theta_3(sL))
  # i.e., parameter-major ordering

  # All sites share the same year structure? Not necessarily.
  # We need to identify years available per site and resample from the
  # intersection or handle per-site. For simplicity, we assume the data
  # has already been structured so each site has observations keyed by year.
  # We resample years and use whatever data is available.

  # First, identify all unique years across all sites
  # We need the raw data organized by (location, year) -> max_sea_level
  # Since maxima is just a vector per site, we need to know which years they correspond to.
  # For bootstrap, the key is resampling the SAME years at ALL sites.

  # We'll assume the maxima vectors are aligned by index (year 1, year 2, ...)
  # and all sites of the same source have the same number of years.
  # This is the standard setup from the paper.

  # Get number of years per site
  T_per_site <- lengths(dat$maxima)

  # For the block bootstrap to preserve spatial dependence, we need to resample
  # the same time indices at all sites. If sites have different numbers of years,
  # we need to handle this carefully.
  # For now, we require all sites of the same source to have the same number of years.
  noaa_idx <- which(dat$sites$data_source == "NOAA")
  adcirc_idx <- which(dat$sites$data_source == "ADCIRC")

  T_noaa <- unique(T_per_site[noaa_idx])
  T_adcirc <- unique(T_per_site[adcirc_idx])

  if (length(T_noaa) > 1) {
    warning("NOAA sites have differing numbers of years. ",
            "Using minimum for bootstrap resampling.")
    T_noaa <- min(T_per_site[noaa_idx])
  }
  if (length(T_adcirc) > 1) {
    warning("ADCIRC sites have differing numbers of years. ",
            "Using minimum for bootstrap resampling.")
    T_adcirc <- min(T_per_site[adcirc_idx])
  }

  # Total parameter vector length: n_sites * 3 (mu, log_sigma, xi per site)
  Lp <- n * 3
  Gamma <- matrix(NA_real_, nrow = B, ncol = Lp)
  n_failures <- 0

  for (b in seq_len(B)) {
    # Resample years separately for NOAA and ADCIRC
    # (they may have different time spans)
    noaa_years <- sample.int(T_noaa, replace = TRUE)
    adcirc_years <- sample.int(T_adcirc, replace = TRUE)

    theta_b <- rep(NA_real_, Lp)

    for (i in seq_len(n)) {
      x <- dat$maxima[[i]]
      if (dat$sites$data_source[i] == "NOAA") {
        x_boot <- x[noaa_years]
      } else {
        x_boot <- x[adcirc_years]
      }

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

#' Embed W_tap into full 6-dimensional parameter space
#'
#' The bootstrap produces a 387x387 matrix (L sites x 3 observed params).
#' This function embeds it into the 774x774 space (L sites x 6 params)
#' required by the stage 2 model. NOAA sites map to dims 1-3, ADCIRC sites
#' map to dims 4-6.
#'
#' @param W_tap Tapered bootstrap covariance (L*3 x L*3).
#' @param dat Data object with site information.
#' @return Embedded covariance matrix (L*6 x L*6).
#' @export
embed_W <- function(W_tap, dat) {
  L <- dat$n_sites
  p_obs <- 3  # observed params per site
  p_full <- 6  # full params in model

  # Build mapping from W_tap indices to full 774-dim indices
  # W_tap ordering: param-major with 3 params
  #   indices 1:L = param 1 at all sites
  #   indices (L+1):(2L) = param 2 at all sites
  #   indices (2L+1):(3L) = param 3 at all sites
  #
  # Full 774-dim ordering: param-major with 6 params
  #   NOAA sites: params 1-3 stay at same relative positions

  #   ADCIRC sites: params 1-3 map to params 4-6 (shift by 3L)

  noaa_idx <- which(dat$sites$data_source == "NOAA")
  adcirc_idx <- which(dat$sites$data_source == "ADCIRC")

  # Create index mapping: W_tap index -> full index
  idx_map <- integer(L * p_obs)

  for (j in seq_len(p_obs)) {
    # W_tap indices for param j: ((j-1)*L + 1) : (j*L)
    w_start <- (j - 1) * L

    for (i in seq_len(L)) {
      w_idx <- w_start + i

      if (dat$sites$data_source[i] == "NOAA") {
        # NOAA: param j stays as param j
        # Full index: (j-1)*L + i
        idx_map[w_idx] <- (j - 1) * L + i
      } else {
        # ADCIRC: param j becomes param (j+3)
        # Full index: (j+3-1)*L + i = (j+2)*L + i
        idx_map[w_idx] <- (j + 2) * L + i
      }
    }
  }

  # Build embedded matrix
  W_full <- matrix(0, nrow = L * p_full, ncol = L * p_full)

  for (i in seq_len(L * p_obs)) {
    for (k in seq_len(L * p_obs)) {
      W_full[idx_map[i], idx_map[k]] <- W_tap[i, k]
    }
  }

  W_full
}
