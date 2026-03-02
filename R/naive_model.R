#' Subset bootstrap covariance to a single data source
#'
#' Extracts the source-specific block from the full bootstrap covariance
#' matrix W_bs (L*3 x L*3 in parameter-major ordering).
#'
#' @param W_bs Raw bootstrap covariance matrix (L*3 x L*3).
#' @param dat An \code{evfuse_data} object.
#' @param source Character: "NOAA" or "ADCIRC".
#' @return Subsetted covariance matrix (L_sub*3 x L_sub*3).
#' @export
subset_W_bs <- function(W_bs, dat, source) {
  L <- dat$n_sites
  src_idx <- which(dat$sites$data_source == source)
  # W_bs is parameter-major: param j at site i -> index (j-1)*L + i
  w_idx <- c(src_idx, src_idx + L, src_idx + 2 * L)
  W_bs[w_idx, w_idx]
}

#' Fit a naive (single-source) 3-dim spatial model
#'
#' Fits a 3-dimensional GP with coregionalization to data from a single
#' source (NOAA or ADCIRC only). All sites observe all 3 parameters,
#' eliminating the partial-observation complexity of the joint model.
#'
#' @param stage1 An \code{evfuse_stage1} object from \code{\link{fit_gev_all}}
#'   or \code{\link{fit_gev_detrended}}.
#' @param dat An \code{evfuse_data} object.
#' @param W_bs Raw bootstrap covariance matrix (L*3 x L*3).
#' @param D Distance matrix (L x L).
#' @param source Character: "NOAA" or "ADCIRC".
#' @param lambda Wendland taper range in km (default 300).
#' @param n_starts Number of multi-start attempts (default 20).
#' @param control Control list passed to \code{optim}.
#' @return An \code{evfuse_naive_model} object.
#' @export
fit_naive_model <- function(stage1, dat, W_bs, D,
                            source = c("NOAA", "ADCIRC"),
                            lambda = 300, n_starts = 20,
                            control = list(maxit = 2000, trace = 1)) {
  source <- match.arg(source)
  p <- 3

  # ── Subset to source sites ─────────────────────────────────
  src_idx <- which(dat$sites$data_source == source)
  L_sub <- length(src_idx)

  stage1_sub <- list(
    theta_hat = stage1$theta_hat[src_idx, , drop = FALSE],
    converged = stage1$converged[src_idx]
  )
  D_sub <- D[src_idx, src_idx]
  W_bs_sub <- subset_W_bs(W_bs, dat, source)
  W_tap <- taper_W(W_bs_sub, D_sub, lambda, p = p)

  # ── Observation structure (full: all sites observe all params) ──
  n_obs <- L_sub * p
  obs_param <- ((seq_len(n_obs) - 1) %/% L_sub) + 1
  obs_site  <- ((seq_len(n_obs) - 1) %% L_sub) + 1
  param_groups <- lapply(1:p, function(j) which(obs_param == j))
  param_sites  <- lapply(1:p, function(j) obs_site[param_groups[[j]]])

  # theta_obs in parameter-major order
  theta_obs <- numeric(n_obs)
  for (j in seq_len(p)) {
    theta_obs[((j - 1) * L_sub + 1):(j * L_sub)] <- stage1_sub$theta_hat[, j]
  }

  W_tap_obs <- W_tap  # full matrix (no subsetting needed)

  # ── Shared cache between nll and nll_grad ───────────────────
  cache <- new.env(parent = emptyenv())
  cache$par <- NULL

  compute_core <- function(par) {
    if (!is.null(cache$par) && identical(par, cache$par)) {
      return(cache$result)
    }

    params <- unpack_params(par, p = p)
    A <- params$A
    rho <- params$rho
    beta <- params$beta

    Omega <- vector("list", p)
    Omega_dot <- vector("list", p)
    for (i in seq_len(p)) {
      Omega[[i]] <- exp_cov_matrix(D_sub, rho[i])
      Omega_dot[[i]] <- D_sub / rho[i]^2 * Omega[[i]]
    }

    # Build Sigma_obs directly
    Sigma_obs <- matrix(0, n_obs, n_obs)
    for (i in seq_len(p)) {
      a_vec <- A[obs_param, i]
      Sigma_obs <- Sigma_obs + tcrossprod(a_vec) * Omega[[i]][obs_site, obs_site]
    }

    V_obs <- Sigma_obs + W_tap_obs + diag(1e-6, n_obs)
    mu_obs <- rep(beta, each = L_sub)
    resid <- theta_obs - mu_obs

    R_chol <- tryCatch(chol(V_obs), error = function(e) NULL)
    if (is.null(R_chol)) {
      cache$par <- par; cache$result <- NULL; return(NULL)
    }

    alpha <- backsolve(R_chol, backsolve(R_chol, resid, transpose = TRUE))
    log_det <- 2 * sum(log(diag(R_chol)))
    nll_val <- 0.5 * (n_obs * log(2 * pi) + log_det + sum(resid * alpha))

    result <- list(A = A, rho = rho, alpha = alpha, R_chol = R_chol,
                   nll_val = nll_val, Omega = Omega, Omega_dot = Omega_dot)
    cache$par <- par
    cache$result <- result
    result
  }

  nll <- function(par) {
    core <- compute_core(par)
    if (is.null(core)) return(1e20)
    core$nll_val
  }

  nll_grad <- function(par) {
    core <- compute_core(par)
    if (is.null(core)) return(rep(0, length(par)))

    A <- core$A
    alpha <- core$alpha
    R_chol <- core$R_chol
    Omega <- core$Omega
    Omega_dot <- core$Omega_dot

    V_inv <- chol2inv(R_chol)
    Q <- V_inv - tcrossprod(alpha)

    compute_H <- function(Omega_mat) {
      H <- matrix(0, p, p)
      for (j1 in seq_len(p)) {
        idx1 <- param_groups[[j1]]
        s1 <- param_sites[[j1]]
        for (j2 in j1:p) {
          idx2 <- param_groups[[j2]]
          s2 <- param_sites[[j2]]
          val <- sum(Q[idx1, idx2, drop = FALSE] *
                       Omega_mat[s1, s2, drop = FALSE])
          H[j1, j2] <- val
          if (j1 != j2) H[j2, j1] <- val
        }
      }
      H
    }

    H <- lapply(Omega, compute_H)
    H_dot <- lapply(Omega_dot, compute_H)

    # Beta gradient
    grad_beta <- numeric(p)
    for (k in seq_len(p)) {
      grad_beta[k] <- -sum(alpha[param_groups[[k]]])
    }

    # A gradient (lower triangular)
    grad_A <- matrix(0, p, p)
    for (a in seq_len(p)) {
      for (b in seq_len(a)) {
        grad_A[a, b] <- sum(H[[b]][a, ] * A[, b])
      }
    }
    grad_A_vec <- grad_A[lower.tri(grad_A, diag = TRUE)]

    # Rho gradient (log-rho chain rule)
    rho <- core$rho
    grad_rho <- numeric(p)
    for (k in seq_len(p)) {
      Ak <- A[, k]
      grad_rho[k] <- 0.5 * drop(crossprod(Ak, H_dot[[k]] %*% Ak)) * rho[k]
    }

    c(grad_beta, grad_A_vec, grad_rho)
  }

  # ── Starting values ─────────────────────────────────────────
  converged_idx <- which(stage1_sub$converged)
  beta0 <- colMeans(stage1_sub$theta_hat[converged_idx, , drop = FALSE])
  A0 <- diag(0.1, p)
  rho0 <- rep(200, p)
  par0 <- pack_params(beta0, A0, rho0, p = p)

  # ── Initial fit ─────────────────────────────────────────────
  cat(sprintf("Fitting %s-only model (p=%d, L=%d, %d params)...\n",
              source, p, L_sub, length(par0)))

  result <- optim(par0, nll, gr = nll_grad, method = "L-BFGS-B",
                  lower = rep(-Inf, length(par0)),
                  upper = rep(Inf, length(par0)),
                  control = control)

  best_nll <- result$value
  best_result <- result
  cat(sprintf("Initial NLL: %.4f\n", best_nll))

  # ── Multi-start ─────────────────────────────────────────────
  for (i in seq_len(n_starts)) {
    rho_i <- exp(runif(p, log(50), log(5000)))
    par_i <- pack_params(beta0, A0, rho_i, p = p)
    fit_i <- tryCatch(
      suppressWarnings(
        optim(par_i, nll, gr = nll_grad, method = "L-BFGS-B",
              lower = rep(-Inf, length(par_i)),
              upper = rep(Inf, length(par_i)),
              control = list(maxit = 2000, trace = 0))
      ),
      error = function(e) NULL
    )
    if (is.null(fit_i)) next
    cat(sprintf("  Start %2d: NLL = %.4f  rho = %s\n",
                i, fit_i$value,
                paste(round(exp(fit_i$par[(length(par_i) - p + 1):length(par_i)]), 1),
                      collapse = " ")))
    if (fit_i$value < best_nll) {
      best_nll <- fit_i$value
      best_result <- fit_i
    }
  }

  if (best_result$convergence != 0) {
    warning(source, "-only model did not converge. Code: ", best_result$convergence)
  }

  cat(sprintf("Best NLL: %.4f\n", best_nll))

  # ── Build output ────────────────────────────────────────────
  params <- unpack_params(best_result$par, p = p)
  Sigma <- build_sigma(params$A, params$rho, D_sub, p = p)

  structure(
    list(
      beta = params$beta,
      A = params$A,
      rho = params$rho,
      Sigma = Sigma,
      W_tap = W_tap,
      D = D_sub,
      source = source,
      sites = dat$sites[src_idx, ],
      stage1 = stage1_sub,
      optim_result = best_result
    ),
    class = "evfuse_naive_model"
  )
}

#' Predict at new locations from a naive model
#'
#' Universal kriging for the 3-dim naive model. All observed sites
#' contribute all 3 parameters (no partial observations).
#'
#' @param model An \code{evfuse_naive_model} from \code{fit_naive_model}.
#' @param new_sites Data frame with lon and lat columns.
#' @return An \code{evfuse_predictions} object compatible with
#'   \code{compute_return_levels}.
#' @export
predict_krig_naive <- function(model, new_sites) {
  p <- 3
  L <- nrow(model$sites)
  n_new <- nrow(new_sites)

  D_no <- compute_cross_distances(new_sites, model$sites)

  beta <- model$beta
  A <- model$A
  rho <- model$rho

  # Full covariance at observed sites
  V_obs <- model$Sigma + model$W_tap + diag(1e-6, L * p)
  R_obs <- chol(V_obs)

  # Observed theta in parameter-major order
  theta_obs <- numeric(L * p)
  for (j in seq_len(p)) {
    theta_obs[((j - 1) * L + 1):(j * L)] <- model$stage1$theta_hat[, j]
  }
  mu_obs <- rep(beta, each = L)
  resid <- theta_obs - mu_obs

  alpha <- backsolve(R_obs, backsolve(R_obs, resid, transpose = TRUE))

  pred_mean <- matrix(NA_real_, n_new, p)
  pred_cov <- vector("list", n_new)

  for (k in seq_len(n_new)) {
    Sigma_cross <- build_cross_sigma_single(A, rho, D_no[k, ], L, p = p)
    pred_mean[k, ] <- beta + as.vector(Sigma_cross %*% alpha)

    Sigma_new <- build_sigma_single(A, rho, p = p)
    solved <- backsolve(R_obs, backsolve(R_obs, t(Sigma_cross), transpose = TRUE))
    pred_cov[[k]] <- Sigma_new - Sigma_cross %*% solved
  }

  colnames(pred_mean) <- c("mu", "log_sigma", "xi")

  structure(
    list(
      pred_mean = pred_mean,
      pred_cov = pred_cov,
      noaa_mean = pred_mean,
      noaa_cov = pred_cov,
      new_sites = new_sites
    ),
    class = "evfuse_predictions"
  )
}
