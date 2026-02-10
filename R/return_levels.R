#' Compute return levels with confidence intervals
#'
#' Takes kriging predictions (NOAA GEV parameters and their covariance)
#' and computes r-year return levels with 95% confidence intervals via
#' the delta method and/or simulation from the predictive distribution.
#'
#' @param predictions An \code{evfuse_predictions} object from \code{predict_krig}.
#' @param r Return period in years (default 100).
#' @param alpha Confidence level for CIs (default 0.05 for 95% CIs).
#' @param method One of "delta", "simulation", or "both" (default "both").
#' @param n_sim Number of simulations for the simulation method (default 2500).
#' @param seed Random seed for simulation method.
#' @return A data frame with columns:
#'   \describe{
#'     \item{lon, lat}{Coordinates of prediction sites.}
#'     \item{return_level}{Point estimate of the r-year return level.}
#'     \item{se_delta}{Standard error from the delta method (if computed).}
#'     \item{ci_lower_delta, ci_upper_delta}{Delta method CI bounds.}
#'     \item{ci_lower_sim, ci_upper_sim}{Simulation-based CI bounds.}
#'   }
#' @export
compute_return_levels <- function(predictions, r = 100, alpha = 0.05,
                                   method = "both", n_sim = 2500, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  n_new <- nrow(predictions$noaa_mean)
  z_alpha <- qnorm(1 - alpha / 2)
  p <- 1 / r  # exceedance probability
  yp <- -log(1 - p)

  results <- data.frame(
    lon = predictions$new_sites$lon,
    lat = predictions$new_sites$lat,
    return_level = NA_real_,
    se_delta = NA_real_,
    ci_lower_delta = NA_real_,
    ci_upper_delta = NA_real_,
    ci_lower_sim = NA_real_,
    ci_upper_sim = NA_real_
  )

  for (k in seq_len(n_new)) {
    mu <- predictions$noaa_mean[k, "mu"]
    log_sigma <- predictions$noaa_mean[k, "log_sigma"]
    xi <- predictions$noaa_mean[k, "xi"]
    sigma <- exp(log_sigma)
    V <- predictions$noaa_cov[[k]]

    # Point estimate
    rl <- gev_return_level(mu, sigma, xi, r)
    results$return_level[k] <- rl

    # Delta method
    if (method %in% c("delta", "both")) {
      grad <- rl_gradient(mu, log_sigma, xi, r)
      var_rl <- as.numeric(t(grad) %*% V %*% grad)
      se <- sqrt(max(var_rl, 0))
      results$se_delta[k] <- se
      results$ci_lower_delta[k] <- rl - z_alpha * se
      results$ci_upper_delta[k] <- rl + z_alpha * se
    }

    # Simulation
    if (method %in% c("simulation", "both")) {
      # Ensure V is positive definite
      V_reg <- V + diag(1e-8, 3)
      chol_V <- tryCatch(chol(V_reg), error = function(e) NULL)
      if (!is.null(chol_V)) {
        sims <- matrix(rnorm(n_sim * 3), n_sim, 3) %*% chol_V
        sims <- sweep(sims, 2, c(mu, log_sigma, xi), "+")

        rl_sims <- vapply(seq_len(n_sim), function(s) {
          gev_return_level(sims[s, 1], exp(sims[s, 2]), sims[s, 3], r)
        }, numeric(1))

        # Remove any Inf/NaN from extreme xi draws
        rl_sims <- rl_sims[is.finite(rl_sims)]
        if (length(rl_sims) > 10) {
          results$ci_lower_sim[k] <- quantile(rl_sims, alpha / 2)
          results$ci_upper_sim[k] <- quantile(rl_sims, 1 - alpha / 2)
        }
      }
    }
  }

  results
}

#' Gradient of return level w.r.t. (mu, log_sigma, xi)
#'
#' For the delta method. Note: parameterized in terms of log(sigma),
#' not sigma, to match the kriging output.
#'
#' @param mu Location parameter.
#' @param log_sigma Log of scale parameter.
#' @param xi Shape parameter.
#' @param r Return period.
#' @return Numeric vector of length 3.
#' @keywords internal
rl_gradient <- function(mu, log_sigma, xi, r) {
  sigma <- exp(log_sigma)
  p <- 1 / r
  yp <- -log(1 - p)

  # RL = mu - (sigma/xi)(1 - yp^{-xi}) for xi != 0
  # Derivatives w.r.t. (mu, log_sigma, xi):
  if (abs(xi) < 1e-8) {
    # RL = mu - sigma * log(yp)
    # d/d_mu = 1
    # d/d_log_sigma = -sigma * log(yp)  (chain rule: d/d_log_sigma = sigma * d/d_sigma)
    # d/d_xi: use limit, approximately sigma * (log(yp))^2 / 2
    dmu <- 1
    dlog_sigma <- -sigma * log(yp)
    dxi <- sigma * (log(yp))^2 / 2
  } else {
    t_val <- yp^(-xi)
    # d(RL)/d(mu) = 1
    dmu <- 1

    # d(RL)/d(sigma) = -(1/xi)(1 - t_val)
    # d(RL)/d(log_sigma) = sigma * d(RL)/d(sigma)
    dlog_sigma <- sigma * (-(1 / xi) * (1 - t_val))

    # d(RL)/d(xi) = (sigma/xi^2)(1 - t_val) - (sigma/xi) * t_val * log(yp)
    dxi <- (sigma / xi^2) * (1 - t_val) - (sigma / xi) * t_val * log(yp)
  }

  c(dmu, dlog_sigma, dxi)
}
