#' Fit Stage 2 spatial model via MLE
#'
#' Estimates beta (6x1), A (6x6 lower triangular), and rho (6x1) by maximizing
#' the Gaussian likelihood from Eq. (7) of Russell et al. (2020), extended to
#' handle partial observations. NOAA sites observe components 1-3 and ADCIRC
#' sites observe components 4-6 of the 6-dimensional parameter vector.
#'
#' Uses an analytic gradient derived from the coregionalization structure
#' (see Sections 2, 5, 6 of the gradient notes). The key identity is that
#' dSigma/dA_ab and dSigma/drho_k are rank-1 Kronecker products, enabling
#' O(n_obs^2) gradient computation per parameter instead of O(n_obs^3).
#'
#' @param stage1 An \code{evfuse_stage1} object from \code{\link{fit_gev_all}}
#'   or \code{\link{fit_gev_detrended}}.
#' @param dat An \code{evfuse_data} object.
#' @param W_tap Tapered covariance matrix from \code{taper_W}.
#' @param D Distance matrix from \code{compute_distances}.
#' @param start Optional named list of starting values for beta, A, rho.
#' @param method Optimization method for \code{optim} (default "L-BFGS-B").
#' @param control Control list passed to \code{optim}.
#' @param check_gradient If TRUE, verify analytic gradient against numerical
#'   finite differences at the starting values before optimizing.
#' @return A list with components:
#'   \describe{
#'     \item{beta}{Estimated mean vector (6x1).}
#'     \item{A}{Estimated lower triangular matrix (6x6).}
#'     \item{rho}{Estimated range parameters (6x1).}
#'     \item{Sigma}{Estimated covariance matrix at observed sites.}
#'     \item{optim_result}{Raw output from optim.}
#'     \item{grad_check}{If \code{check_gradient = TRUE}, the output of
#'       \code{verify_gradient}.}
#'   }
#' @export
fit_spatial_model <- function(stage1, dat, W_tap, D,
                               start = NULL, method = "L-BFGS-B",
                               control = list(maxit = 1000, trace = 1),
                               check_gradient = FALSE) {
  L <- dat$n_sites
  p <- if (!is.null(dat$p)) dat$p else 6

  # Embed W_tap from observed-param space into full param space if needed
  p_obs <- length(dat$source_params[[1]])  # params per source
  if (nrow(W_tap) == L * p_obs && p_obs < p) {
    W_tap <- embed_W(W_tap, dat)
  }

  # Build observation structure: which indices of the full 6L vector are observed
  obs <- build_observation_structure(stage1, dat)
  n_obs <- length(obs$obs_idx)

  # Precompute observation metadata for gradient
  # For full-vector index idx: param = ceiling(idx/L), site = ((idx-1) %% L) + 1
  obs_param <- ((obs$obs_idx - 1) %/% L) + 1  # which param (1-6) each obs corresponds to
  obs_site  <- ((obs$obs_idx - 1) %% L) + 1   # which site (1-L) each obs corresponds to
  param_groups <- lapply(1:p, function(j) which(obs_param == j))
  param_sites  <- lapply(1:p, function(j) obs_site[param_groups[[j]]])

  # Precompute W_tap in observed space (constant throughout optimization)
  W_tap_obs <- W_tap[obs$obs_idx, obs$obs_idx]

  # --- Shared cache between nll and nll_grad ---
  # L-BFGS-B calls fn(x) then gr(x) at the same point each iteration.
  # Caching avoids a redundant Cholesky decomposition.
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

    # rho = exp(log_rho) is always positive; no guard needed

    # Precompute Omega_i = exp(-D/rho_i) and Omega_dot_i = (D/rho_i^2) * Omega_i
    Omega <- vector("list", p)
    Omega_dot <- vector("list", p)
    for (i in seq_len(p)) {
      Omega[[i]] <- exp_cov_matrix(D, rho[i])
      Omega_dot[[i]] <- D / rho[i]^2 * Omega[[i]]
    }

    # Build Sigma_obs directly from Kronecker structure, avoiding full 6L x 6L matrix.
    # Sigma_obs[m, n] = sum_i A[param[m], i] * A[param[n], i] * Omega_i[site[m], site[n]]
    Sigma_obs <- matrix(0, n_obs, n_obs)
    for (i in seq_len(p)) {
      a_vec <- A[obs_param, i]
      Sigma_obs <- Sigma_obs + tcrossprod(a_vec) * Omega[[i]][obs_site, obs_site]
    }

    V_obs <- Sigma_obs + W_tap_obs + diag(1e-6, n_obs)
    mu_obs <- rep(beta, each = L)[obs$obs_idx]
    resid <- obs$theta_obs - mu_obs

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

  # --- Negative log-likelihood ---
  nll <- function(par) {
    core <- compute_core(par)
    if (is.null(core)) return(1e20)
    core$nll_val
  }

  # --- Analytic gradient of negative log-likelihood ---
  # Derivation follows Sections 2, 4-6 of the gradient notes.
  # Cost: 1 chol2inv + 12 H-matrix contractions (O(n_obs^2) each) + 33 dot products
  nll_grad <- function(par) {
    core <- compute_core(par)
    if (is.null(core)) return(rep(0, length(par)))

    A <- core$A
    alpha <- core$alpha
    R_chol <- core$R_chol
    Omega <- core$Omega
    Omega_dot <- core$Omega_dot

    # Q = V^{-1} - alpha alpha^T  (computed once, reused for all 33 components)
    V_inv <- chol2inv(R_chol)
    Q <- V_inv - tcrossprod(alpha)

    # Contract Q with each Omega_i into a 6x6 matrix H_i.
    # H_i[j1, j2] = sum_{m: param=j1, n: param=j2} Q[m,n] * Omega_i[site_m, site_n]
    # Then: tr[(V^{-1} - alpha alpha^T) S (u v^T ⊗ Omega) S^T] = u^T H v
    # H is symmetric since Q and Omega are both symmetric.
    compute_H <- function(Omega_mat) {
      H <- matrix(0, p, p)
      for (j1 in seq_len(p)) {
        idx1 <- param_groups[[j1]]
        if (length(idx1) == 0) next
        s1 <- param_sites[[j1]]
        for (j2 in j1:p) {
          idx2 <- param_groups[[j2]]
          if (length(idx2) == 0) next
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

    # --- Beta gradient (Eq. 2): d(nll)/d(beta) = -X^T alpha ---
    # Column k of design matrix X has 1s where obs_param == k, 0 elsewhere
    grad_beta <- numeric(p)
    for (k in seq_len(p)) {
      idx_k <- param_groups[[k]]
      if (length(idx_k) > 0) {
        grad_beta[k] <- -sum(alpha[idx_k])
      }
    }

    # --- A gradient ---
    # dSigma/dA_{ab} = [e_a (Ae_b)^T + (Ae_b) e_a^T] ⊗ Omega_b
    # Only Omega_b appears because only the k=b term in sum_k A_{j1,k} A_{j2,k} Omega_k
    # depends on A_{ab}. Contracting with Q via H_b:
    # d(nll)/dA_{ab} = H_b[a,:] . A[:,b]
    grad_A <- matrix(0, p, p)
    for (a in seq_len(p)) {
      for (b in seq_len(a)) {
        grad_A[a, b] <- sum(H[[b]][a, ] * A[, b])
      }
    }
    grad_A_vec <- grad_A[lower.tri(grad_A, diag = TRUE)]

    # --- Rho gradient (Eq. 8), with log-rho chain rule ---
    # dSigma/drho_k = (A e_k e_k^T A^T) ⊗ Omega_dot_k
    # d(nll)/drho_k = 0.5 * A[:,k]^T H_dot_k A[:,k]
    # Chain rule: d(nll)/d(log_rho_k) = d(nll)/drho_k * drho_k/d(log_rho_k)
    #                                 = d(nll)/drho_k * rho_k
    rho <- core$rho
    grad_rho <- numeric(p)
    for (k in seq_len(p)) {
      Ak <- A[, k]
      grad_rho[k] <- 0.5 * drop(crossprod(Ak, H_dot[[k]] %*% Ak)) * rho[k]
    }

    c(grad_beta, grad_A_vec, grad_rho)
  }

  # Starting values
  if (is.null(start)) {
    start <- default_start(stage1, dat, p = p)
  }
  par0 <- pack_params(start$beta, start$A, start$rho, p = p)

  # Verify gradient before optimizing if requested
  if (check_gradient) {
    grad_check <- verify_gradient(par0, nll, nll_grad)
    message(sprintf("Gradient check: max relative error = %.2e", grad_check$max_rel_err))
    if (grad_check$max_rel_err > 1e-3) {
      warning("Analytic gradient may be incorrect. Max rel error: ",
              sprintf("%.2e", grad_check$max_rel_err))
    }
  }

  # All parameters are unconstrained (rho is optimized in log space)
  n_par <- length(par0)
  lower <- rep(-Inf, n_par)
  upper <- rep(Inf, n_par)

  result <- optim(par0, nll, gr = nll_grad, method = method,
                  lower = lower, upper = upper,
                  control = control)

  if (result$convergence != 0) {
    warning("Stage 2 optimization did not converge. Code: ", result$convergence)
  }

  params <- unpack_params(result$par, p = p)
  Sigma <- build_sigma(params$A, params$rho, D, p = p)

  out <- structure(
    list(
      beta = params$beta,
      A = params$A,
      rho = params$rho,
      Sigma = Sigma,
      W_tap = W_tap,
      D = D,
      dat = dat,
      stage1 = stage1,
      optim_result = result,
      p = p
    ),
    class = "evfuse_model"
  )
  if (check_gradient) out$grad_check <- grad_check
  out
}

#' Verify analytic gradient against numerical finite differences
#'
#' Computes the analytic gradient and a central finite-difference approximation
#' at the same parameter vector, then reports the relative error for each
#' component.
#'
#' @param par Parameter vector to check at.
#' @param fn Objective function (scalar-valued).
#' @param gr Gradient function (returns vector same length as par).
#' @param eps Finite difference step size (default 1e-5).
#' @return List with components:
#'   \describe{
#'     \item{analytic}{Analytic gradient vector.}
#'     \item{numeric}{Numerical gradient vector.}
#'     \item{rel_err}{Relative error for each component.}
#'     \item{max_rel_err}{Maximum relative error across all components.}
#'   }
#' @export
verify_gradient <- function(par, fn, gr, eps = 1e-5) {
  analytic <- gr(par)
  n <- length(par)
  numeric_grad <- numeric(n)

  for (i in seq_len(n)) {
    par_up <- par; par_up[i] <- par[i] + eps
    par_dn <- par; par_dn[i] <- par[i] - eps
    numeric_grad[i] <- (fn(par_up) - fn(par_dn)) / (2 * eps)
  }

  denom <- pmax(abs(analytic), abs(numeric_grad), 1e-8)
  rel_err <- abs(analytic - numeric_grad) / denom

  list(
    analytic = analytic,
    numeric = numeric_grad,
    rel_err = rel_err,
    max_rel_err = max(rel_err)
  )
}

#' Build observation structure for partial observations
#'
#' Maps each site to the indices of its observed components in the
#' full Lp-dimensional parameter vector.
#'
#' @param stage1 Stage 1 fits.
#' @param dat Data object.
#' @return List with obs_idx (indices into the full vector) and theta_obs
#'   (the observed parameter values).
#' @keywords internal
build_observation_structure <- function(stage1, dat) {
  L <- dat$n_sites
  source_params <- if (!is.null(dat$source_params)) dat$source_params
                   else list(NOAA = 1:3, ADCIRC = 4:6)

  obs_idx <- c()
  theta_obs <- c()

  for (i in seq_len(L)) {
    if (!stage1$converged[i]) next

    src <- dat$sites$data_source[i]
    param_indices <- source_params[[src]]

    # In the full Lp vector with parameter-major ordering:
    # param j at site i -> index (j-1)*L + i
    full_indices <- (param_indices - 1) * L + i
    obs_idx <- c(obs_idx, full_indices)
    n_obs_i <- length(param_indices)
    theta_obs <- c(theta_obs, stage1$theta_hat[i, seq_len(n_obs_i)])
  }

  list(obs_idx = obs_idx, theta_obs = theta_obs)
}

#' Pack parameters into a single vector for optim
#'
#' Rho is stored as log(rho) so the optimizer works in unconstrained space.
#'
#' @param beta Mean vector (p).
#' @param A Lower triangular matrix (p x p).
#' @param rho Range parameters (p), in natural (positive) scale.
#' @param p Dimension (default 6).
#' @return Numeric vector.
#' @keywords internal
pack_params <- function(beta, A, rho, p = 6) {
  # A: extract lower triangular elements (including diagonal)
  A_vec <- A[lower.tri(A, diag = TRUE)]
  c(beta, A_vec, log(rho))
}

#' Unpack parameter vector into beta, A, rho
#'
#' The last p entries are log(rho); this function exponentiates to return
#' rho in the natural (positive) scale.
#'
#' @param par Numeric vector from optim.
#' @param p Dimension (default 6).
#' @return Named list with beta, A, rho (positive scale).
#' @keywords internal
unpack_params <- function(par, p = 6) {
  beta <- par[1:p]

  n_A <- p * (p + 1) / 2
  A_vec <- par[(p + 1):(p + n_A)]
  A <- matrix(0, p, p)
  A[lower.tri(A, diag = TRUE)] <- A_vec

  rho <- exp(par[(p + n_A + 1):(p + n_A + p)])

  list(beta = beta, A = A, rho = rho)
}

#' Generate default starting values
#'
#' @param stage1 Stage 1 fits.
#' @param dat Data object.
#' @param p Dimension (default 6).
#' @return Named list with beta, A, rho.
#' @keywords internal
default_start <- function(stage1, dat, p = NULL) {
  source_params <- if (!is.null(dat$source_params)) dat$source_params
                   else list(NOAA = 1:3, ADCIRC = 4:6)
  if (is.null(p)) p <- max(unlist(source_params))

  # beta: mean of observed parameters per source component
  beta <- rep(0, p)
  for (src in names(source_params)) {
    src_idx <- which(dat$sites$data_source == src & stage1$converged)
    if (length(src_idx) > 0) {
      n_obs_src <- length(source_params[[src]])
      beta[source_params[[src]]] <- colMeans(
        stage1$theta_hat[src_idx, seq_len(n_obs_src), drop = FALSE],
        na.rm = TRUE
      )
    }
  }

  # A: start with small diagonal
  A <- diag(0.1, p)

  # rho: start with moderate range (e.g., 200 km)
  rho <- rep(200, p)

  list(beta = beta, A = A, rho = rho)
}
