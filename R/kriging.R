#' Predict at new locations via universal kriging
#'
#' Given the fitted stage 2 model, predict the full 6-dimensional parameter
#' vector and its covariance at new spatial locations. Then extract the
#' NOAA components (1-3) and their 3x3 marginal covariance.
#'
#' @param model A fitted \code{evfuse_model} from \code{fit_spatial_model}.
#' @param new_sites Data frame with lon and lat columns for prediction locations.
#' @return A list with components:
#'   \describe{
#'     \item{pred_mean}{Matrix (n_new x 6) of predicted parameter means.}
#'     \item{pred_cov}{List of n_new 6x6 predictive covariance matrices.}
#'     \item{noaa_mean}{Matrix (n_new x 3) of predicted NOAA GEV parameters
#'       (mu, log_sigma, xi).}
#'     \item{noaa_cov}{List of n_new 3x3 marginal covariances for NOAA params.}
#'   }
#' @export
predict_krig <- function(model, new_sites) {
  p <- 6
  L <- model$dat$n_sites
  n_new <- nrow(new_sites)

  # Distances: observed-observed, new-observed
  D_oo <- model$D
  D_no <- compute_cross_distances(new_sites, model$dat$sites)

  beta <- model$beta
  A <- model$A
  rho <- model$rho

  # Build cross-covariance between new and observed sites
  # For each latent GP delta_i, the cross-covariance is exp(-D_no / rho_i)
  # Then the full cross-cov is (I_new_L kron A) [sum_i (e_i e_i^T kron Omega_no_i)] (I_L kron A)^T

  # Full covariance at observed sites
  V_full <- model$Sigma + model$W_tap

  # Observation structure
  obs <- build_observation_structure(model$stage1, model$dat)
  V_obs <- V_full[obs$obs_idx, obs$obs_idx]
  V_obs <- V_obs + diag(1e-6, nrow(V_obs))

  # Cholesky of V_obs for solving
  R_obs <- chol(V_obs)
  mu_obs <- rep(beta, each = L)[obs$obs_idx]
  resid <- obs$theta_obs - mu_obs

  # Precompute alpha = V_obs^{-1} resid (constant across prediction sites)
  alpha <- backsolve(R_obs, backsolve(R_obs, resid, transpose = TRUE))

  # Preallocate outputs
  pred_mean <- matrix(NA_real_, n_new, p)
  pred_cov <- vector("list", n_new)

  for (k in seq_len(n_new)) {
    # Cross-covariance between new site k and all observed sites
    # Sigma_{new_k, obs}: a (p x Lp) matrix, then select observed columns
    Sigma_cross_full <- build_cross_sigma_single(A, rho, D_no[k, ], L, p)
    Sigma_cross <- Sigma_cross_full[, obs$obs_idx, drop = FALSE]

    # Kriging mean: beta + Sigma_cross V_obs^{-1} resid
    pred_mean[k, ] <- beta + as.vector(Sigma_cross %*% alpha)

    # Kriging variance: Sigma_{new,new} - Sigma_cross V_obs^{-1} Sigma_cross^T
    Sigma_new <- build_sigma_single(A, rho, p)
    solved <- backsolve(R_obs, backsolve(R_obs, t(Sigma_cross), transpose = TRUE))
    pred_cov[[k]] <- Sigma_new - Sigma_cross %*% solved
  }

  # Extract NOAA components (indices 1:3)
  noaa_mean <- pred_mean[, 1:3, drop = FALSE]
  colnames(noaa_mean) <- c("mu", "log_sigma", "xi")

  noaa_cov <- lapply(pred_cov, function(V) V[1:3, 1:3])

  structure(
    list(
      pred_mean = pred_mean,
      pred_cov = pred_cov,
      noaa_mean = noaa_mean,
      noaa_cov = noaa_cov,
      new_sites = new_sites
    ),
    class = "evfuse_predictions"
  )
}

#' Cross-distance matrix between two sets of sites
#'
#' @param sites1 Data frame with lon, lat (n1 rows).
#' @param sites2 Data frame with lon, lat (n2 rows).
#' @return Matrix of dimension (n1 x n2) of great-circle distances in km.
#' @export
compute_cross_distances <- function(sites1, sites2) {
  n1 <- nrow(sites1)
  n2 <- nrow(sites2)
  R <- 6371

  lon1 <- sites1$lon * pi / 180
  lat1 <- sites1$lat * pi / 180
  lon2 <- sites2$lon * pi / 180
  lat2 <- sites2$lat * pi / 180

  D <- matrix(0, n1, n2)
  for (i in seq_len(n1)) {
    for (j in seq_len(n2)) {
      dlat <- lat2[j] - lat1[i]
      dlon <- lon2[j] - lon1[i]
      a <- sin(dlat / 2)^2 + cos(lat1[i]) * cos(lat2[j]) * sin(dlon / 2)^2
      D[i, j] <- 2 * R * asin(sqrt(a))
    }
  }
  D
}

#' Cross-covariance for a single new site
#'
#' Builds the p x (L*p) cross-covariance matrix between a single new site
#' and all L observed sites.
#'
#' @param A Lower triangular matrix (p x p).
#' @param rho Range parameters (p).
#' @param d_vec Vector of distances from new site to each of L observed sites.
#' @param L Number of observed sites.
#' @param p Number of parameters (default 6).
#' @return Matrix of dimension (p x Lp).
#' @keywords internal
build_cross_sigma_single <- function(A, rho, d_vec, L, p = 6) {
  # For each GP i, the cross-covariance vector is exp(-d_vec / rho_i)
  # Full cross-cov: A %*% [sum_i e_i (c_i^T kron e_i^T)] ... but simpler:
  # Sigma_cross = (1 kron A) [sum_i (e_i kron omega_i)] (I_L kron A)^T
  # where omega_i = exp(-d_vec / rho_i), a vector of length L

  # Simpler construction: for new site, the cross-covariance with the
  # L*p observed vector is:
  # For parameter components (j_new, j_obs) at sites (new, l):
  # Cov = sum_i A[j_new, i] * A[j_obs, i] * exp(-d_l / rho_i)

  Sigma_cross <- matrix(0, p, L * p)
  for (j_new in seq_len(p)) {
    for (j_obs in seq_len(p)) {
      # Covariance between param j_new at new site and param j_obs at obs sites
      cov_vec <- rep(0, L)
      for (i in seq_len(p)) {
        cov_vec <- cov_vec + A[j_new, i] * A[j_obs, i] * exp(-d_vec / rho[i])
      }
      # These go into columns ((j_obs - 1) * L + 1) : (j_obs * L) of row j_new
      col_idx <- ((j_obs - 1) * L + 1):(j_obs * L)
      Sigma_cross[j_new, col_idx] <- cov_vec
    }
  }
  Sigma_cross
}

#' Marginal covariance at a single site (no spatial separation)
#'
#' At distance 0, Omega_i = 1, so Sigma = A A^T.
#'
#' @param A Lower triangular matrix.
#' @param rho Range parameters (unused at distance 0, included for API consistency).
#' @param p Dimension.
#' @return p x p covariance matrix.
#' @keywords internal
build_sigma_single <- function(A, rho, p = 6) {
  A %*% t(A)
}
