#' Leave-one-out cross-validation via Rasmussen & Williams (2006) shortcut
#'
#' Computes closed-form LOO-CV predictions at NOAA sites using the
#' block version of Rasmussen & Williams (2006) eq. 5.12. For site k with
#' observation block B_k, the LOO predictive mean and covariance are
#' computed from the precision matrix without refitting.
#'
#' @param model A fitted model (evfuse_model or evfuse_naive_model).
#' @return An \code{evfuse_loo} object with components:
#'   \describe{
#'     \item{loo_mean}{Matrix (n_sites x 3) of LOO predictive means.}
#'     \item{loo_cov}{List of 3x3 LOO predictive covariance matrices.}
#'     \item{observed}{Matrix (n_sites x 3) of observed Stage 1 MLEs.}
#'     \item{sites}{Data frame of site metadata.}
#'     \item{loo_lpd}{Vector of per-site log predictive densities.}
#'   }
#' @export
loo_cv <- function(model, ...) UseMethod("loo_cv")

#' @rdname loo_cv
#' @param loo_source Which data source to evaluate LOO-CV at (default: first
#'   source in \code{source_params}, typically \code{"NOAA"}).
#' @param ... Additional arguments (currently unused).
#' @export
loo_cv.evfuse_model <- function(model, loo_source = NULL, ...) {
  L <- model$dat$n_sites
  p <- if (!is.null(model$p)) model$p else 6
  source_params <- if (!is.null(model$dat$source_params)) model$dat$source_params
                   else list(NOAA = 1:3, ADCIRC = 4:6)
  if (is.null(loo_source)) loo_source <- names(source_params)[1]
  p_block <- length(source_params[[loo_source]])

  obs <- build_observation_structure(model$stage1, model$dat)
  n_obs <- length(obs$obs_idx)

  obs_param <- ((obs$obs_idx - 1) %/% L) + 1
  obs_site  <- ((obs$obs_idx - 1) %% L) + 1

  V_full <- model$Sigma + model$W_tap
  V_obs <- V_full[obs$obs_idx, obs$obs_idx] + diag(1e-6, n_obs)
  mu_obs <- rep(model$beta, each = L)[obs$obs_idx]

  R_chol <- chol(V_obs)
  V_inv <- chol2inv(R_chol)
  alpha <- V_inv %*% (obs$theta_obs - mu_obs)

  # LOO at each site of loo_source
  target_idx <- which(model$dat$sites$data_source == loo_source)
  n_target <- length(target_idx)

  loo_mean <- matrix(NA_real_, n_target, p_block)
  loo_cov  <- vector("list", n_target)
  loo_lpd  <- numeric(n_target)

  for (k in seq_len(n_target)) {
    site_k <- target_idx[k]
    block_idx <- which(obs_site == site_k)

    V_inv_BB <- V_inv[block_idx, block_idx]
    Sigma_loo <- solve(V_inv_BB)
    alpha_B <- alpha[block_idx]

    loo_mean[k, ] <- obs$theta_obs[block_idx] - Sigma_loo %*% alpha_B
    loo_cov[[k]] <- Sigma_loo

    # Log predictive density: log N(y_B; mu_LOO, Sigma_LOO)
    resid_k <- obs$theta_obs[block_idx] - loo_mean[k, ]
    R_k <- chol(Sigma_loo)
    loo_lpd[k] <- -0.5 * p_block * log(2 * pi) - sum(log(diag(R_k))) -
      0.5 * sum(resid_k * solve(Sigma_loo, resid_k))
  }

  colnames(loo_mean) <- c("mu", "log_sigma", "xi")

  structure(
    list(
      loo_mean = loo_mean,
      loo_cov  = loo_cov,
      observed = model$stage1$theta_hat[target_idx, ],
      sites    = model$dat$sites[target_idx, ],
      loo_lpd  = loo_lpd
    ),
    class = "evfuse_loo"
  )
}

#' @rdname loo_cv
#' @export
loo_cv.evfuse_naive_model <- function(model, ...) {
  p <- 3
  L <- nrow(model$sites)
  n_obs <- L * p

  V_obs <- model$Sigma + model$W_tap + diag(1e-6, n_obs)

  theta_obs <- numeric(n_obs)
  for (j in seq_len(p)) {
    theta_obs[((j - 1) * L + 1):(j * L)] <- model$stage1$theta_hat[, j]
  }
  mu_obs <- rep(model$beta, each = L)

  R_chol <- chol(V_obs)
  V_inv <- chol2inv(R_chol)
  alpha <- V_inv %*% (theta_obs - mu_obs)

  loo_mean <- matrix(NA_real_, L, p)
  loo_cov  <- vector("list", L)
  loo_lpd  <- numeric(L)

  for (k in seq_len(L)) {
    block_idx <- k + (0:(p - 1)) * L

    V_inv_BB <- V_inv[block_idx, block_idx]
    Sigma_loo <- solve(V_inv_BB)
    alpha_B <- alpha[block_idx]

    loo_mean[k, ] <- theta_obs[block_idx] - Sigma_loo %*% alpha_B
    loo_cov[[k]] <- Sigma_loo

    resid_k <- theta_obs[block_idx] - loo_mean[k, ]
    R_k <- chol(Sigma_loo)
    loo_lpd[k] <- -0.5 * p * log(2 * pi) - sum(log(diag(R_k))) -
      0.5 * sum(resid_k * solve(Sigma_loo, resid_k))
  }

  colnames(loo_mean) <- c("mu", "log_sigma", "xi")

  structure(
    list(
      loo_mean = loo_mean,
      loo_cov  = loo_cov,
      observed = model$stage1$theta_hat,
      sites    = model$sites,
      loo_lpd  = loo_lpd
    ),
    class = "evfuse_loo"
  )
}

#' Summarize LOO-CV results
#'
#' Computes RMSE, MAD, and return level metrics from LOO predictions.
#'
#' @param loo An \code{evfuse_loo} object.
#' @param r Return period for return level comparison (default 100).
#' @return A list with parameter-level and return-level accuracy metrics.
#' @export
loo_summary <- function(loo, r = 100) {
  n <- nrow(loo$loo_mean)
  params <- c("mu", "log_sigma", "xi")

  # Per-parameter RMSE and MAD
  param_stats <- data.frame(
    parameter = params,
    rmse = NA_real_,
    mad = NA_real_,
    stringsAsFactors = FALSE
  )
  for (j in seq_along(params)) {
    resid <- loo$loo_mean[, j] - loo$observed[, j]
    param_stats$rmse[j] <- sqrt(mean(resid^2))
    param_stats$mad[j]  <- mean(abs(resid))
  }

  # Return level comparison
  rl_loo <- vapply(seq_len(n), function(k) {
    gev_return_level(loo$loo_mean[k, "mu"],
                     exp(loo$loo_mean[k, "log_sigma"]),
                     loo$loo_mean[k, "xi"], r)
  }, numeric(1))
  rl_obs <- vapply(seq_len(n), function(k) {
    gev_return_level(loo$observed[k, 1],
                     exp(loo$observed[k, 2]),
                     loo$observed[k, 3], r)
  }, numeric(1))

  # Delta-method SE at LOO predictions
  rl_se <- vapply(seq_len(n), function(k) {
    grad <- rl_gradient(loo$loo_mean[k, "mu"],
                        loo$loo_mean[k, "log_sigma"],
                        loo$loo_mean[k, "xi"], r)
    sqrt(max(drop(t(grad) %*% loo$loo_cov[[k]] %*% grad), 0))
  }, numeric(1))

  rl_resid <- rl_loo - rl_obs

  list(
    param_stats = param_stats,
    total_lpd   = sum(loo$loo_lpd),
    mean_lpd    = mean(loo$loo_lpd),
    rl_rmse     = sqrt(mean(rl_resid^2)),
    rl_mad      = mean(abs(rl_resid)),
    rl_mean_se  = mean(rl_se),
    rl_loo      = rl_loo,
    rl_obs      = rl_obs,
    rl_se       = rl_se
  )
}
