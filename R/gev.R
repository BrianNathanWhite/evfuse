#' Fit GEV at all sites (Stage 1)
#'
#' Fits a stationary GEV distribution to annual block maxima at each site
#' using \code{extRemes::fevd}. Returns pointwise MLEs and their asymptotic
#' covariance matrices.
#'
#' @param dat An \code{evfuse_data} object from \code{load_data}.
#' @param log_scale Logical. If TRUE (default), return log(sigma) instead of
#'   sigma. This is standard for the stage 2 model where we need unconstrained
#'   parameters.
#' @return A list with components:
#'   \describe{
#'     \item{theta_hat}{Matrix of dimension (n_sites x 3). Columns are
#'       mu, log_sigma (or sigma), xi. Rows ordered to match dat$sites.}
#'     \item{vcov_list}{List of n_sites 3x3 asymptotic covariance matrices
#'       for (mu, log_sigma, xi) at each site.}
#'     \item{fits}{List of raw fevd fit objects for diagnostics.}
#'     \item{converged}{Logical vector indicating convergence at each site.}
#'   }
#' @export
fit_gev_all <- function(dat, log_scale = TRUE) {
  n <- dat$n_sites
  theta_hat <- matrix(NA_real_, nrow = n, ncol = 3)
  colnames(theta_hat) <- c("mu", "log_sigma", "xi")
  vcov_list <- vector("list", n)
  fits <- vector("list", n)
  converged <- logical(n)

  for (i in seq_len(n)) {
    loc_name <- dat$sites$location[i]
    x <- dat$maxima[[i]]

    fit <- tryCatch(
      extRemes::fevd(x, type = "GEV", method = "MLE"),
      error = function(e) NULL
    )

    if (is.null(fit) || !fit$results$convergence == 0) {
      warning(sprintf("GEV fit failed or did not converge at site '%s'.", loc_name))
      converged[i] <- FALSE
      fits[[i]] <- fit
      next
    }

    converged[i] <- TRUE
    fits[[i]] <- fit

    # fevd returns parameters as (location, scale, shape)
    pars <- fit$results$par
    mu <- pars["location"]
    sigma <- pars["scale"]
    xi <- pars["shape"]

    # Get the asymptotic covariance matrix from fevd
    # fevd parameterizes as (location, scale, shape)
    V_raw <- vcov(fit)

    if (log_scale) {
      # Transform scale -> log(scale) via delta method
      # Jacobian: d(mu, log_sigma, xi) / d(mu, sigma, xi) = diag(1, 1/sigma, 1)
      J <- diag(c(1, 1 / sigma, 1))
      V_transformed <- J %*% V_raw %*% t(J)

      theta_hat[i, ] <- c(mu, log(sigma), xi)
      vcov_list[[i]] <- V_transformed
    } else {
      colnames(theta_hat) <- c("mu", "sigma", "xi")
      theta_hat[i, ] <- c(mu, sigma, xi)
      vcov_list[[i]] <- V_raw
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
      log_scale = log_scale
    ),
    class = "evfuse_stage1"
  )
}

#' GEV quantile function (return level)
#'
#' Computes the r-year return level from GEV parameters.
#'
#' @param mu Location parameter.
#' @param sigma Scale parameter (NOT log scale).
#' @param xi Shape parameter.
#' @param r Return period in years.
#' @return The r-year return level.
#' @export
gev_return_level <- function(mu, sigma, xi, r) {
  p <- 1 / r
  yp <- -log(1 - p)
  if (abs(xi) < 1e-8) {
    mu - sigma * log(yp)
  } else {
    mu - (sigma / xi) * (1 - yp^(-xi))
  }
}

#' GEV exceedance probability
#'
#' Computes P(Z > z0) for Z ~ GEV(mu, sigma, xi).
#'
#' @param z0 Threshold value.
#' @param mu Location parameter.
#' @param sigma Scale parameter (NOT log scale).
#' @param xi Shape parameter.
#' @return Exceedance probability.
#' @export
gev_exceedance_prob <- function(z0, mu, sigma, xi) {
  if (abs(xi) < 1e-8) {
    1 - exp(-exp(-(z0 - mu) / sigma))
  } else {
    t_val <- 1 + xi * (z0 - mu) / sigma
    if (t_val <= 0) {
      if (xi > 0) return(1)
      else return(0)
    }
    1 - exp(-t_val^(-1 / xi))
  }
}
