#!/usr/bin/env Rscript
# ══════════════════════════════════════════════════════════════
# QC REGRESSION CHECK after source_params refactoring
# ══════════════════════════════════════════════════════════════
# Tests every refactored code path against published results.
# Exits with status 1 on any discrepancy.
devtools::load_all()

pass <- 0L
fail <- 0L
results <- data.frame(Quantity = character(), Expected = character(),
                      Got = character(), Status = character(),
                      stringsAsFactors = FALSE)

check <- function(name, expected, got, tol) {
  status <- if (abs(got - expected) <= tol) "PASS" else "FAIL"
  results <<- rbind(results, data.frame(
    Quantity = name,
    Expected = sprintf("%.4f", expected),
    Got      = sprintf("%.4f", got),
    Status   = status,
    stringsAsFactors = FALSE
  ))
  if (status == "PASS") pass <<- pass + 1L else fail <<- fail + 1L
}

# ══════════════════════════════════════════════════════════════
# 1. Load data with refactored load_data()
# ══════════════════════════════════════════════════════════════
cat("\n== 1. Loading data with refactored load_data() ==\n")
df <- read.csv("data-raw/combined_max.csv")
dat <- load_data(df)

stopifnot(!is.null(dat$source_params))
stopifnot(!is.null(dat$p))
stopifnot(dat$p == 6)
stopifnot(dat$n_sites == 129)
stopifnot(dat$n_noaa == 29)
stopifnot(dat$n_adcirc == 100)
cat("  load_data() structure OK\n")

# ══════════════════════════════════════════════════════════════
# 2. Stationary pipeline: NLL + correlations
# ══════════════════════════════════════════════════════════════
cat("\n== 2. Stationary joint model ==\n")
model <- readRDS("data-raw/model_6dim_best.rds")

# Re-evaluate NLL at saved parameters using refactored code
# This tests: build_observation_structure, embed_W, compute_core
D <- compute_distances(dat$sites)
stage1 <- fit_gev_all(dat)

# Rebuild W_tap
bs <- bootstrap_W(dat, B = 500, seed = 42)
W_tap <- taper_W(bs$W_bs, D, lambda = 300)
W_tap_full <- embed_W(W_tap, dat)

# Re-evaluate NLL using refactored observation structure
obs <- build_observation_structure(stage1, dat)
L <- dat$n_sites
p <- dat$p
n_obs <- length(obs$obs_idx)
obs_param <- ((obs$obs_idx - 1) %/% L) + 1
obs_site  <- ((obs$obs_idx - 1) %% L) + 1

# Build Sigma_obs at saved parameters
Omega <- lapply(seq_len(p), function(i) exp_cov_matrix(D, model$rho[i]))
Sigma_obs <- matrix(0, n_obs, n_obs)
for (i in seq_len(p)) {
  a_vec <- model$A[obs_param, i]
  Sigma_obs <- Sigma_obs + tcrossprod(a_vec) * Omega[[i]][obs_site, obs_site]
}
V_obs <- Sigma_obs + W_tap_full[obs$obs_idx, obs$obs_idx] + diag(1e-6, n_obs)
mu_obs <- rep(model$beta, each = L)[obs$obs_idx]
resid <- obs$theta_obs - mu_obs

R_chol <- chol(V_obs)
alpha <- backsolve(R_chol, backsolve(R_chol, resid, transpose = TRUE))
log_det <- 2 * sum(log(diag(R_chol)))
nll_recomputed <- 0.5 * (n_obs * log(2 * pi) + log_det + sum(resid * alpha))

cat(sprintf("  Saved NLL: %.4f\n", model$optim_result$value))
cat(sprintf("  Recomputed NLL: %.4f\n", nll_recomputed))
check("NLL (stationary)", -103.61, nll_recomputed, 0.5)
check("NLL match saved", model$optim_result$value, nll_recomputed, 0.001)

# Cross-source correlations from A matrix
AAT <- model$A %*% t(model$A)
sds <- sqrt(diag(AAT))
cor_mat <- AAT / outer(sds, sds)
check("Cor(mu_N, mu_A)", 0.992, cor_mat[1, 4], 0.001)
check("Cor(log_sigma_N, log_sigma_A)", 0.521, cor_mat[2, 5], 0.001)
check("Cor(xi_N, xi_A)", 0.952, cor_mat[3, 6], 0.001)

# ══════════════════════════════════════════════════════════════
# 3. Stationary LOO-CV with refactored code
# ══════════════════════════════════════════════════════════════
cat("\n== 3. Stationary LOO-CV (refactored) ==\n")

# Need to re-create a proper model object with refactored dat
model_refactored <- model
model_refactored$dat <- dat
model_refactored$stage1 <- stage1
model_refactored$W_tap <- W_tap_full
model_refactored$p <- p

loo_joint <- loo_cv(model_refactored)
sum_joint <- loo_summary(loo_joint, r = 100)

cat(sprintf("  Total LPD (joint): %.1f\n", sum_joint$total_lpd))

# NOAA-only model LOO-CV
noaa_model <- readRDS("data-raw/model_noaa_only.rds")
noaa_model$stage1 <- list(
  theta_hat = stage1$theta_hat[which(dat$sites$data_source == "NOAA"), ],
  converged = stage1$converged[which(dat$sites$data_source == "NOAA")]
)
loo_noaa <- loo_cv(noaa_model)
sum_noaa <- loo_summary(loo_noaa, r = 100)

cat(sprintf("  Total LPD (NOAA-only): %.1f\n", sum_noaa$total_lpd))

n_joint_wins <- sum(loo_joint$loo_lpd > loo_noaa$loo_lpd)
lpd_gain <- sum_joint$total_lpd - sum_noaa$total_lpd
rl_pct <- 100 * (1 - sum_joint$rl_rmse / sum_noaa$rl_rmse)

check("Joint wins (stat)", 28, n_joint_wins, 0.5)
check("Total LPD NOAA-only (stat)", -2.5, sum_noaa$total_lpd, 0.1)
check("Total LPD Joint (stat)", 30.0, sum_joint$total_lpd, 0.1)
check("LPD gain (stat)", 32.4, lpd_gain, 0.2)
check("RL RMSE reduction % (stat)", 34.9, rl_pct, 0.5)

# ══════════════════════════════════════════════════════════════
# 4. Detrended pipeline
# ══════════════════════════════════════════════════════════════
cat("\n== 4. Detrended pipeline ==\n")

detr <- readRDS("data-raw/detrended_results.rds")
detr_model <- detr$model

# Re-compute detrended cross-source correlations
AAT_dt <- detr_model$A %*% t(detr_model$A)
sds_dt <- sqrt(diag(AAT_dt))
cor_dt <- AAT_dt / outer(sds_dt, sds_dt)
check("Cor(mu_N, mu_A) detr", 0.991, cor_dt[1, 4], 0.001)
check("Cor(log_sigma_N, log_sigma_A) detr", 0.533, cor_dt[2, 5], 0.001)
check("Cor(xi_N, xi_A) detr", 0.808, cor_dt[3, 6], 0.001)

# Re-run detrended LOO-CV with refactored code
cat("  Re-running detrended Stage 1...\n")
stage1_dt <- fit_gev_detrended(dat, df, ref_year = 2000)

# Rebuild detrended model object with refactored dat
detr_model$dat <- dat
detr_model$stage1 <- stage1_dt
detr_model$p <- p

loo_dt <- loo_cv(detr_model)
sum_dt <- loo_summary(loo_dt, r = 100)

cat(sprintf("  Total LPD (detr joint): %.1f\n", sum_dt$total_lpd))

# Detrended NOAA-only
noaa_dt <- readRDS("data-raw/model_noaa_only_detrended.rds")
noaa_dt$stage1 <- list(
  theta_hat = stage1_dt$theta_hat[which(dat$sites$data_source == "NOAA"), ],
  converged = stage1_dt$converged[which(dat$sites$data_source == "NOAA")]
)
loo_noaa_dt <- loo_cv(noaa_dt)
sum_noaa_dt <- loo_summary(loo_noaa_dt, r = 100)

cat(sprintf("  Total LPD (detr NOAA-only): %.1f\n", sum_noaa_dt$total_lpd))

n_dt_wins <- sum(loo_dt$loo_lpd > loo_noaa_dt$loo_lpd)
lpd_gain_dt <- sum_dt$total_lpd - sum_noaa_dt$total_lpd
rl_pct_dt <- 100 * (1 - sum_dt$rl_rmse / sum_noaa_dt$rl_rmse)

check("Joint wins (detr)", 28, n_dt_wins, 0.5)
check("Total LPD NOAA-only (detr)", -8.5, sum_noaa_dt$total_lpd, 0.1)
check("Total LPD Joint (detr)", 23.0, sum_dt$total_lpd, 0.1)
check("LPD gain (detr)", 31.5, lpd_gain_dt, 0.2)
check("RL RMSE reduction % (detr)", 34.7, rl_pct_dt, 0.5)

# ══════════════════════════════════════════════════════════════
# RESULTS TABLE
# ══════════════════════════════════════════════════════════════
cat("\n")
cat("══════════════════════════════════════════════════════════════\n")
cat("REGRESSION CHECK RESULTS\n")
cat("══════════════════════════════════════════════════════════════\n\n")

# Pretty print
max_name <- max(nchar(results$Quantity))
fmt <- sprintf("%%-%ds  %%10s  %%10s  %%s\n", max_name)
cat(sprintf(fmt, "Quantity", "Expected", "Got", "Status"))
cat(sprintf(fmt, paste(rep("-", max_name), collapse = ""),
            "----------", "----------", "------"))
for (i in seq_len(nrow(results))) {
  status_str <- if (results$Status[i] == "PASS") "PASS" else "** FAIL **"
  cat(sprintf(fmt, results$Quantity[i], results$Expected[i],
              results$Got[i], status_str))
}

cat(sprintf("\nTotal: %d PASS, %d FAIL\n", pass, fail))

if (fail > 0) {
  cat("\n** REGRESSION FAILURES DETECTED — NOT regenerating outputs **\n")
  quit(status = 1)
} else {
  cat("\nAll checks passed. Regenerating figures and tables...\n")
  source("scripts/plot_results.R", local = TRUE)
  source("scripts/plot_comparison.R", local = TRUE)
  source("scripts/run_trends.R", local = TRUE)
  cat("\nAll outputs regenerated successfully.\n")
}
