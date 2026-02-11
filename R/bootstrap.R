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
#' (L sites x p_obs params per site). This function embeds it into the full
#' joint parameter space (L sites x p_full params) required by the stage 2 model.
#' Each source's observed params are mapped to their corresponding positions
#' in the joint model via \code{source_params}.
#'
#' @param W_tap Tapered bootstrap covariance (L*p_obs x L*p_obs).
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
  p_obs <- length(source_params[[1]])  # observed params per site (all sources same)
  p_full <- max(unlist(source_params))

  # Build mapping from W_tap indices to full-space indices
  # W_tap ordering: param-major with p_obs params
  # Full ordering: param-major with p_full params
  idx_map <- integer(L * p_obs)

  for (j in seq_len(p_obs)) {
    w_start <- (j - 1) * L
    for (i in seq_len(L)) {
      w_idx <- w_start + i
      src <- dat$sites$data_source[i]
      full_param <- source_params[[src]][j]
      idx_map[w_idx] <- (full_param - 1) * L + i
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
