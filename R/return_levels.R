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
#'     \item{se_sim}{Standard error from simulation (if computed).}
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
    se_sim = NA_real_,
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
          results$se_sim[k] <- sd(rl_sims)
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

#' Compute nonstationary return levels with covariate
#'
#' Takes kriging predictions of 4 GEV parameters (mu0, mu1, log_sigma, xi)
#' and a future covariate value, and computes r-year return levels with
#' confidence intervals via the delta method and/or simulation.
#'
#' The effective location is \code{mu = mu0 + mu1 * covariate_value}, where
#' \code{covariate_value} should be centered relative to the same reference
#' used in \code{\link{fit_gev_ns}}.
#'
#' @param predictions An \code{evfuse_predictions} object from
#'   \code{\link{predict_krig}} with 4-column \code{noaa_mean} and 4x4
#'   covariance matrices.
#' @param covariate_value Scalar or vector (length n_new) of centered
#'   covariate values at prediction sites.
#' @param r Return period in years (default 100).
#' @param alpha Confidence level for CIs (default 0.05 for 95% CIs).
#' @param method One of \code{"delta"}, \code{"simulation"}, or
#'   \code{"both"} (default).
#' @param n_sim Number of simulations (default 2500).
#' @param seed Random seed for simulation method.
#' @return A data frame with columns: lon, lat, return_level, se_delta,
#'   ci_lower_delta, ci_upper_delta, se_sim, ci_lower_sim, ci_upper_sim.
#' @export
compute_return_levels_ns <- function(predictions, covariate_value,
                                      r = 100, alpha = 0.05,
                                      method = "both", n_sim = 2500,
                                      seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  n_new <- nrow(predictions$noaa_mean)
  z_alpha <- qnorm(1 - alpha / 2)

  if (length(covariate_value) == 1) {
    covariate_value <- rep(covariate_value, n_new)
  }

  results <- data.frame(
    lon = predictions$new_sites$lon,
    lat = predictions$new_sites$lat,
    return_level = NA_real_,
    se_delta = NA_real_,
    ci_lower_delta = NA_real_,
    ci_upper_delta = NA_real_,
    se_sim = NA_real_,
    ci_lower_sim = NA_real_,
    ci_upper_sim = NA_real_
  )

  for (k in seq_len(n_new)) {
    mu0 <- predictions$noaa_mean[k, 1]
    mu1 <- predictions$noaa_mean[k, 2]
    log_sigma <- predictions$noaa_mean[k, 3]
    xi <- predictions$noaa_mean[k, 4]
    sigma <- exp(log_sigma)
    sst <- covariate_value[k]
    V <- predictions$noaa_cov[[k]]  # 4x4

    # Effective location parameter
    mu <- mu0 + mu1 * sst

    # Point estimate
    rl <- gev_return_level(mu, sigma, xi, r)
    results$return_level[k] <- rl

    # Delta method
    if (method %in% c("delta", "both")) {
      # 3-param gradient w.r.t. (mu, log_sigma, xi)
      grad3 <- rl_gradient(mu, log_sigma, xi, r)

      # 4-param gradient w.r.t. (mu0, mu1, log_sigma, xi) via chain rule
      grad4 <- c(grad3[1], grad3[1] * sst, grad3[2], grad3[3])

      var_rl <- as.numeric(t(grad4) %*% V %*% grad4)
      se <- sqrt(max(var_rl, 0))
      results$se_delta[k] <- se
      results$ci_lower_delta[k] <- rl - z_alpha * se
      results$ci_upper_delta[k] <- rl + z_alpha * se
    }

    # Simulation
    if (method %in% c("simulation", "both")) {
      V_reg <- V + diag(1e-8, 4)
      chol_V <- tryCatch(chol(V_reg), error = function(e) NULL)
      if (!is.null(chol_V)) {
        sims <- matrix(rnorm(n_sim * 4), n_sim, 4) %*% chol_V
        sims <- sweep(sims, 2, c(mu0, mu1, log_sigma, xi), "+")

        rl_sims <- vapply(seq_len(n_sim), function(s) {
          mu_s <- sims[s, 1] + sims[s, 2] * sst
          sigma_s <- exp(sims[s, 3])
          xi_s <- sims[s, 4]
          gev_return_level(mu_s, sigma_s, xi_s, r)
        }, numeric(1))

        rl_sims <- rl_sims[is.finite(rl_sims)]
        if (length(rl_sims) > 10) {
          results$se_sim[k] <- sd(rl_sims)
          results$ci_lower_sim[k] <- quantile(rl_sims, alpha / 2)
          results$ci_upper_sim[k] <- quantile(rl_sims, 1 - alpha / 2)
        }
      }
    }
  }

  results
}
