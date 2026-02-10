#' Fit Stage 2 spatial model via MLE
#'
#' Estimates beta (6x1), A (6x6 lower triangular), and rho (6x1) by maximizing
#' the Gaussian likelihood from Eq. (7) of Russell et al. (2019), extended to
#' handle partial observations. NOAA sites observe components 1-3 and ADCIRC
#' sites observe components 4-6 of the 6-dimensional parameter vector.
#'
#' @param stage1 An \code{evfuse_stage1} object from \code{fit_gev_all}.
#' @param dat An \code{evfuse_data} object.
#' @param W_tap Tapered covariance matrix from \code{taper_W}.
#' @param D Distance matrix from \code{compute_distances}.
#' @param start Optional named list of starting values for beta, A, rho.
#' @param method Optimization method for \code{optim} (default "L-BFGS-B").
#' @param control Control list passed to \code{optim}.
#' @return A list with components:
#'   \describe{
#'     \item{beta}{Estimated mean vector (6x1).}
#'     \item{A}{Estimated lower triangular matrix (6x6).}
#'     \item{rho}{Estimated range parameters (6x1).}
#'     \item{Sigma}{Estimated covariance Sigma_{A,rho} at observed sites.}
#'     \item{optim_result}{Raw output from optim.}
#'   }
#' @export
fit_spatial_model <- function(stage1, dat, W_tap, D,
                               start = NULL, method = "L-BFGS-B",
                               control = list(maxit = 1000, trace = 1)) {
  L <- dat$n_sites
  p <- 6  # total parameters: 3 NOAA + 3 ADCIRC

  # Build the observation vector Theta_hat and selection structure
  # Theta_hat ordering: parameter-major.
  # For the 6-dim model: (param1 at all sites, param2 at all sites, ..., param6 at all sites)
  # But each site only observes 3 of 6 parameters.
  # NOAA sites observe params 1,2,3. ADCIRC sites observe params 4,5,6.
  obs <- build_observation_structure(stage1, dat)

  # Negative log-likelihood
  nll <- function(par) {
    params <- unpack_params(par, p = p)
    beta <- params$beta
    A <- params$A
    rho <- params$rho

    # Ensure rho > 0 (should be handled by box constraints, but just in case)
    if (any(rho <= 0)) return(1e20)

    # Build full covariance: Sigma_{A,rho} + W_tap
    Sigma <- build_sigma(A, rho, D, p = p)
    V_full <- Sigma + W_tap

    # Mean vector: beta kron 1_L
    mu_full <- rep(beta, each = L)

    # Select observed components
    V_obs <- V_full[obs$obs_idx, obs$obs_idx]
    mu_obs <- mu_full[obs$obs_idx]
    theta_obs <- obs$theta_obs

    # Gaussian log-likelihood
    # Add small ridge for numerical stability
    V_obs <- V_obs + diag(1e-6, nrow(V_obs))

    chol_result <- tryCatch(chol(V_obs), error = function(e) NULL)
    if (is.null(chol_result)) return(1e20)

    R <- chol_result
    resid <- theta_obs - mu_obs
    z <- backsolve(R, backsolve(R, resid, transpose = TRUE))

    # -0.5 * (n*log(2*pi) + log|V| + resid^T V^{-1} resid)
    n_obs <- length(theta_obs)
    log_det <- 2 * sum(log(diag(R)))
    nll_val <- 0.5 * (n_obs * log(2 * pi) + log_det + sum(resid * z))

    nll_val
  }

  # Starting values
  if (is.null(start)) {
    start <- default_start(stage1, dat, p = p)
  }
  par0 <- pack_params(start$beta, start$A, start$rho, p = p)

  # Box constraints: rho > 0 (log-transformed internally? No, use lower bounds)
  n_par <- length(par0)
  lower <- rep(-Inf, n_par)
  upper <- rep(Inf, n_par)
  # rho parameters are the last p entries
  rho_start_idx <- n_par - p + 1
  lower[rho_start_idx:n_par] <- 1e-3  # rho > 0

  result <- optim(par0, nll, method = method,
                  lower = lower, upper = upper,
                  control = control)

  if (result$convergence != 0) {
    warning("Stage 2 optimization did not converge. Code: ", result$convergence)
  }

  params <- unpack_params(result$par, p = p)
  Sigma <- build_sigma(params$A, params$rho, D, p = p)

  structure(
    list(
      beta = params$beta,
      A = params$A,
      rho = params$rho,
      Sigma = Sigma,
      W_tap = W_tap,
      D = D,
      dat = dat,
      stage1 = stage1,
      optim_result = result
    ),
    class = "evfuse_model"
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
  p <- 6

  obs_idx <- c()
  theta_obs <- c()

  for (i in seq_len(L)) {
    if (!stage1$converged[i]) next

    if (dat$sites$data_source[i] == "NOAA") {
      # Observe parameters 1, 2, 3 (mu_N, log_sigma_N, xi_N)
      param_indices <- 1:3
    } else {
      # Observe parameters 4, 5, 6 (mu_A, log_sigma_A, xi_A)
      param_indices <- 4:6
    }

    # In the full Lp vector with parameter-major ordering:
    # param j at site i -> index (j-1)*L + i
    full_indices <- (param_indices - 1) * L + i
    obs_idx <- c(obs_idx, full_indices)
    theta_obs <- c(theta_obs, stage1$theta_hat[i, ])
  }

  list(obs_idx = obs_idx, theta_obs = theta_obs)
}

#' Pack parameters into a single vector for optim
#'
#' @param beta Mean vector (p).
#' @param A Lower triangular matrix (p x p).
#' @param rho Range parameters (p).
#' @param p Dimension (default 6).
#' @return Numeric vector.
#' @keywords internal
pack_params <- function(beta, A, rho, p = 6) {
  # A: extract lower triangular elements (including diagonal)
  A_vec <- A[lower.tri(A, diag = TRUE)]
  c(beta, A_vec, rho)
}

#' Unpack parameter vector into beta, A, rho
#'
#' @param par Numeric vector from optim.
#' @param p Dimension (default 6).
#' @return Named list with beta, A, rho.
#' @keywords internal
unpack_params <- function(par, p = 6) {
  beta <- par[1:p]

  n_A <- p * (p + 1) / 2
  A_vec <- par[(p + 1):(p + n_A)]
  A <- matrix(0, p, p)
  A[lower.tri(A, diag = TRUE)] <- A_vec

  rho <- par[(p + n_A + 1):(p + n_A + p)]

  list(beta = beta, A = A, rho = rho)
}

#' Generate default starting values
#'
#' @param stage1 Stage 1 fits.
#' @param dat Data object.
#' @param p Dimension (default 6).
#' @return Named list with beta, A, rho.
#' @keywords internal
default_start <- function(stage1, dat, p = 6) {
  # beta: mean of observed parameters per component
  beta <- rep(0, p)
  noaa_idx <- which(dat$sites$data_source == "NOAA" & stage1$converged)
  adcirc_idx <- which(dat$sites$data_source == "ADCIRC" & stage1$converged)

  if (length(noaa_idx) > 0) {
    beta[1:3] <- colMeans(stage1$theta_hat[noaa_idx, , drop = FALSE])
  }
  if (length(adcirc_idx) > 0) {
    beta[4:6] <- colMeans(stage1$theta_hat[adcirc_idx, , drop = FALSE])
  }

  # A: start with small diagonal
  A <- diag(0.1, p)

  # rho: start with moderate range (e.g., 200 km)
  rho <- rep(200, p)

  list(beta = beta, A = A, rho = rho)
}
