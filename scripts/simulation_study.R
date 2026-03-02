#!/usr/bin/env Rscript
# simulation_study.R
#
# Recover cross-source correlations under non-co-location.
# Generates synthetic data from the fitted model (known A, rho, beta),
# adds measurement error, refits Stage 2, and checks whether the
# cross-source correlations are recovered across 100 replications.
#
# Outputs: printed summary table (Table 6 in manuscript).

library(devtools)
load_all()

set.seed(2026)

# ── Load fitted model ─────────────────────────────────────────────────────────
cat("Loading fitted model...\n")
dat <- evfuse::coast_data
D <- compute_distances(dat$sites)

model_joint <- readRDS("data-raw/model_6dim_ns.rds")

L <- dat$n_sites
noaa_idx <- which(dat$sites$data_source == "NOAA")
adcirc_idx <- which(dat$sites$data_source == "ADCIRC")

# True model parameters
A_true <- model_joint$A
rho_true <- model_joint$rho
beta_true <- model_joint$beta
p <- 6

# Cross-source correlations from AA^T
Sigma_marginal <- A_true %*% t(A_true)
cors_true <- c(
  mu    = Sigma_marginal[1, 4] / sqrt(Sigma_marginal[1, 1] * Sigma_marginal[4, 4]),
  logsig = Sigma_marginal[2, 5] / sqrt(Sigma_marginal[2, 2] * Sigma_marginal[5, 5]),
  xi    = Sigma_marginal[3, 6] / sqrt(Sigma_marginal[3, 3] * Sigma_marginal[6, 6])
)

cat("\n")
cat("============================================================\n")
cat("True cross-source correlations (from fitted model):\n")
cat(sprintf("  mu:      %.4f\n", cors_true["mu"]))
cat(sprintf("  log_sig: %.4f\n", cors_true["logsig"]))
cat(sprintf("  xi:      %.4f\n", cors_true["xi"]))
cat("============================================================\n\n")

# ── Simulation ────────────────────────────────────────────────────────────────
N_sim <- 100
cat(sprintf("Running %d replications...\n", N_sim))

# Precompute Omega matrices (constant across replications)
Omega_list <- lapply(seq_len(p), function(i) exp_cov_matrix(D, rho_true[i]))

# Cholesky factors for generating from each GP
chol_Omega <- lapply(Omega_list, function(Om) {
  chol(Om + diag(1e-10, nrow(Om)))
})

# W_tap in observed space for generating measurement error
obs <- build_observation_structure(model_joint$stage1, dat)
n_obs <- length(obs$obs_idx)
W_tap_obs <- model_joint$W_tap[obs$obs_idx, obs$obs_idx]
W_tap_obs <- W_tap_obs + diag(1e-6, n_obs)
chol_W <- chol(W_tap_obs)

# Map obs indices to (param, site)
obs_param <- ((obs$obs_idx - 1) %/% L) + 1
obs_site  <- ((obs$obs_idx - 1) %% L) + 1

# Storage for recovered correlations
cors_recovered <- matrix(NA_real_, N_sim, 3,
                         dimnames = list(NULL, c("mu", "logsig", "xi")))

for (sim in seq_len(N_sim)) {
  if (sim %% 10 == 0) cat(sprintf("  Replication %d/%d\n", sim, N_sim))

  # Generate latent GPs: delta_i(s) ~ N(0, Omega_i)
  delta <- matrix(0, L, p)
  for (i in seq_len(p)) {
    z <- rnorm(L)
    delta[, i] <- as.vector(crossprod(chol_Omega[[i]], z))
  }

  # theta(s) = beta + A * delta(s) for each site
  theta_full <- matrix(NA_real_, L, p)
  for (l in seq_len(L)) {
    theta_full[l, ] <- beta_true + A_true %*% delta[l, ]
  }

  # Extract observed parameters and add measurement error
  theta_obs_true <- numeric(n_obs)
  for (m in seq_len(n_obs)) {
    theta_obs_true[m] <- theta_full[obs_site[m], obs_param[m]]
  }
  eps <- as.vector(crossprod(chol_W, rnorm(n_obs)))
  theta_obs_sim <- theta_obs_true + eps

  # Build fake stage1 object
  theta_hat_sim <- matrix(NA_real_, L, p)
  for (m in seq_len(n_obs)) {
    theta_hat_sim[obs_site[m], obs_param[m]] <- theta_obs_sim[m]
  }
  # NOAA sites: params 1-3 in cols 1:3; ADCIRC sites: params 4-6 into cols 1:3
  theta_hat_stage1 <- matrix(NA_real_, L, max(3, p))
  for (l in noaa_idx) {
    theta_hat_stage1[l, 1:3] <- theta_hat_sim[l, 1:3]
  }
  for (l in adcirc_idx) {
    theta_hat_stage1[l, 1:3] <- theta_hat_sim[l, 4:6]
  }

  stage1_sim <- list(
    theta_hat = theta_hat_stage1[, 1:3, drop = FALSE],
    converged = rep(TRUE, L)
  )

  # Fit Stage 2 (single start from true values + noise)
  A_start <- A_true + matrix(rnorm(p * p, 0, 0.02), p, p)
  A_start[upper.tri(A_start)] <- 0
  rho_start <- rho_true * exp(rnorm(p, 0, 0.1))
  beta_start <- beta_true + rnorm(p, 0, 0.05)

  fit_sim <- tryCatch({
    suppressWarnings(
      fit_spatial_model(
        stage1_sim, dat, model_joint$W_tap, D,
        start = list(beta = beta_start, A = A_start, rho = rho_start),
        control = list(maxit = 2000, trace = 0)
      )
    )
  }, error = function(e) NULL)

  if (is.null(fit_sim)) {
    cat(sprintf("  Replication %d: fitting failed\n", sim))
    next
  }

  # Extract cross-source correlations
  Sig_hat <- fit_sim$A %*% t(fit_sim$A)
  cors_recovered[sim, "mu"]     <- Sig_hat[1, 4] / sqrt(Sig_hat[1, 1] * Sig_hat[4, 4])
  cors_recovered[sim, "logsig"] <- Sig_hat[2, 5] / sqrt(Sig_hat[2, 2] * Sig_hat[5, 5])
  cors_recovered[sim, "xi"]     <- Sig_hat[3, 6] / sqrt(Sig_hat[3, 3] * Sig_hat[6, 6])
}

# Summarize
n_success <- sum(!is.na(cors_recovered[, 1]))
cat(sprintf("\nSimulation study: %d/%d replications converged\n", n_success, N_sim))
cat("\nRecovered cross-source correlations:\n")
cat(sprintf("%-10s %8s %8s %8s %8s %8s\n",
            "Param", "True", "Mean", "Median", "SD", "IQR"))
for (j in 1:3) {
  nm <- c("mu", "logsig", "xi")[j]
  vals <- cors_recovered[, j]
  vals <- vals[!is.na(vals)]
  cat(sprintf("%-10s %8.4f %8.4f %8.4f %8.4f %8.4f\n",
              nm, cors_true[j], mean(vals), median(vals), sd(vals),
              IQR(vals)))
}
cat("\n")
