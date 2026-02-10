#' Exponential covariance matrix
#'
#' Computes the correlation matrix Omega_i where
#' Omega_i[j,k] = exp(-D[j,k] / rho) for range parameter rho > 0.
#'
#' @param D Distance matrix.
#' @param rho Range parameter (positive).
#' @return Correlation matrix of same dimension as D.
#' @export
exp_cov_matrix <- function(D, rho) {
  if (rho <= 0) stop("rho must be positive")
  exp(-D / rho)
}

#' Wendland 2 compactly supported covariance function
#'
#' Evaluates the Wendland C2 function: for d in [0, 1],
#' C(d) = (1 - d)^6 * (35d^2 + 18d + 3) / 3.
#' Returns 0 for d > 1.
#'
#' @param d Non-negative distance (scaled by range parameter).
#' @return Covariance value.
#' @export
wendland2 <- function(d) {
  d <- pmax(d, 0)
  ifelse(d >= 1, 0, (1 - d)^6 * (35 * d^2 + 18 * d + 3) / 3)
}

#' Wendland 2 taper correlation matrix
#'
#' Constructs a sparse taper matrix based on the Wendland C2 function
#' with compact support at distance lambda.
#'
#' @param D Distance matrix.
#' @param lambda Taper range (same units as D). Correlations are zero beyond this distance.
#' @return Sparse taper matrix (class "dgCMatrix" from Matrix package).
#' @export
wendland_taper_matrix <- function(D, lambda) {
  if (lambda <= 0) stop("lambda must be positive")
  T_mat <- wendland2(D / lambda)
  Matrix::Matrix(T_mat, sparse = TRUE)
}

#' Build the coregionalization covariance matrix Sigma_{A, rho}
#'
#' Constructs the covariance matrix for the spatial random effects
#' following Eq. (7) in Russell et al. (2019). For p parameters and L sites:
#'
#' Sigma_{A,rho} = (I_L x A) [sum_i (e_i e_i^T x Omega_i(rho))] (I_L x A)^T
#'
#' In our setting, p = 6 (3 NOAA + 3 ADCIRC GEV parameters).
#'
#' @param A Lower triangular matrix (p x p).
#' @param rho Vector of p range parameters.
#' @param D Distance matrix (L x L).
#' @param p Number of parameters (default 6).
#' @return Covariance matrix of dimension (L*p x L*p).
#' @export
build_sigma <- function(A, rho, D, p = 6) {
  L <- nrow(D)

  # Build the inner sum: sum_i (e_i e_i^T kron Omega_i)
  inner <- matrix(0, nrow = L * p, ncol = L * p)
  for (i in seq_len(p)) {
    Omega_i <- exp_cov_matrix(D, rho[i])
    # e_i e_i^T is a p x p matrix with 1 in position (i,i), 0 elsewhere
    E_i <- matrix(0, p, p)
    E_i[i, i] <- 1
    inner <- inner + kronecker(E_i, Omega_i)
  }

  # (I_L kron A) %*% inner %*% (I_L kron A)^T
  I_kron_A <- kronecker(diag(L), A)
  I_kron_A %*% inner %*% t(I_kron_A)
}

#' Build the taper matrix for W
#'
#' Constructs T_tap(lambda) = 1_p 1_p^T kron C_W2(lambda) as in Russell et al.
#' This applies the same Wendland taper across all p parameter blocks,
#' allowing cross-parameter dependence at nearby stations.
#'
#' @param D Distance matrix (L x L).
#' @param lambda Taper range in km.
#' @param p Number of GEV parameters (default 6).
#' @return Taper matrix of dimension (L*p x L*p).
#' @export
build_taper <- function(D, lambda, p = 6) {
  C_W2 <- as.matrix(wendland_taper_matrix(D, lambda))
  ones_p <- matrix(1, p, p)
  kronecker(ones_p, C_W2)
}
