#!/usr/bin/env Rscript
# baseline_comparisons.R
#
# Compare the joint LMC model against simpler fusion baselines:
#   - Nearest-ADCIRC + global bias correction
#   - Nearest-ADCIRC + local bias correction (K=5 neighbors)
#   - LOO regression on nearest ADCIRC
# All evaluated via LOO-CV on 100-yr return levels.
#
# Outputs: printed comparison table (Section 5.1 in manuscript).

library(devtools)
load_all()

# ── Load fitted models ────────────────────────────────────────────────────────
cat("Loading models...\n")
dat <- evfuse::coast_data
D <- compute_distances(dat$sites)

model_joint <- readRDS("data-raw/model_6dim_ns.rds")
model_noaa  <- readRDS("data-raw/model_noaa_only_ns.rds")

L <- dat$n_sites
noaa_idx <- which(dat$sites$data_source == "NOAA")
adcirc_idx <- which(dat$sites$data_source == "ADCIRC")
L_N <- length(noaa_idx)

rmse <- function(x) sqrt(mean(x^2))

# ── Reference: LMC LOO-CV return levels ───────────────────────────────────────
loo_j <- loo_cv(model_joint)
loo_n <- loo_cv(model_noaa)
n_noaa <- nrow(loo_j$loo_mean)

rl_obs <- rl_joint <- rl_noaa <- numeric(n_noaa)
for (k in seq_len(n_noaa)) {
  rl_obs[k]   <- gev_return_level(loo_j$observed[k, 1], exp(loo_j$observed[k, 2]), loo_j$observed[k, 3], 100)
  rl_joint[k] <- gev_return_level(loo_j$loo_mean[k, 1], exp(loo_j$loo_mean[k, 2]), loo_j$loo_mean[k, 3], 100)
  rl_noaa[k]  <- gev_return_level(loo_n$loo_mean[k, 1], exp(loo_n$loo_mean[k, 2]), loo_n$loo_mean[k, 3], 100)
}

rmse_noaa  <- rmse(rl_noaa - rl_obs)
rmse_joint <- rmse(rl_joint - rl_obs)

# ── Baseline setup ────────────────────────────────────────────────────────────
noaa_sites <- dat$sites[noaa_idx, ]
adcirc_sites <- dat$sites[adcirc_idx, ]
D_na <- compute_cross_distances(noaa_sites, adcirc_sites)

nearest_adcirc <- apply(D_na, 1, which.min)
nearest_dist <- apply(D_na, 1, min)

cat(sprintf("Nearest ADCIRC to each NOAA site: median %.1f km, range %.1f--%.1f km\n",
            median(nearest_dist), min(nearest_dist), max(nearest_dist)))

stage1 <- model_joint$stage1
theta_adcirc <- stage1$theta_hat[adcirc_idx, 1:3]
theta_noaa   <- stage1$theta_hat[noaa_idx, 1:3]

# Global bias
bias <- colMeans(theta_noaa) - colMeans(theta_adcirc)
cat(sprintf("Global bias (NOAA - ADCIRC mean):\n"))
cat(sprintf("  mu:      %+.3f m\n", bias[1]))
cat(sprintf("  log_sig: %+.3f\n", bias[2]))
cat(sprintf("  xi:      %+.3f\n", bias[3]))
cat("\n")

# ── Method A: Nearest ADCIRC + global bias ────────────────────────────────────
rl_nearest <- numeric(L_N)
for (k in seq_len(L_N)) {
  bias_loo <- colMeans(theta_noaa[-k, ]) - colMeans(theta_adcirc)
  pred_k <- theta_adcirc[nearest_adcirc[k], ] + bias_loo
  rl_nearest[k] <- gev_return_level(pred_k[1], exp(pred_k[2]), pred_k[3], 100)
}
rmse_nearest <- rmse(rl_nearest - rl_obs)

# ── Method B: Nearest ADCIRC + local bias (K=5) ──────────────────────────────
D_nn <- compute_cross_distances(noaa_sites, noaa_sites)
K <- 5

rl_local_bias <- numeric(L_N)
for (k in seq_len(L_N)) {
  dists_k <- D_nn[k, ]
  dists_k[k] <- Inf
  nn_idx <- order(dists_k)[1:min(K, L_N - 1)]

  local_biases <- matrix(NA_real_, length(nn_idx), 3)
  for (j in seq_along(nn_idx)) {
    local_biases[j, ] <- theta_noaa[nn_idx[j], ] - theta_adcirc[nearest_adcirc[nn_idx[j]], ]
  }
  bias_local <- colMeans(local_biases)

  pred_k <- theta_adcirc[nearest_adcirc[k], ] + bias_local
  rl_local_bias[k] <- gev_return_level(pred_k[1], exp(pred_k[2]), pred_k[3], 100)
}
rmse_local <- rmse(rl_local_bias - rl_obs)

# ── Method C: LOO regression on nearest ADCIRC ───────────────────────────────
rl_regkrig <- numeric(L_N)
for (k in seq_len(L_N)) {
  train_noaa <- theta_noaa[-k, ]
  train_adcirc_nn <- theta_adcirc[nearest_adcirc[-k], ]

  pred_k <- numeric(3)
  for (j in 1:3) {
    fit_lm <- lm(train_noaa[, j] ~ train_adcirc_nn[, j])
    pred_k[j] <- predict(fit_lm, newdata = data.frame(
      `train_adcirc_nn[, j]` = theta_adcirc[nearest_adcirc[k], j]
    ))
  }
  rl_regkrig[k] <- gev_return_level(pred_k[1], exp(pred_k[2]), pred_k[3], 100)
}
rmse_regkrig <- rmse(rl_regkrig - rl_obs)

# ── Summary ───────────────────────────────────────────────────────────────────
cat("100-yr RL LOO-CV RMSE comparison:\n")
cat(sprintf("%-40s %8s %10s\n", "Method", "RMSE (m)", "Reduction"))
cat(sprintf("%-40s %8.3f %10s\n", "NOAA-only LMC (baseline)", rmse_noaa, "--"))
cat(sprintf("%-40s %8.3f %9.1f%%\n", "Nearest ADCIRC + global bias", rmse_nearest,
            100 * (1 - rmse_nearest / rmse_noaa)))
cat(sprintf("%-40s %8.3f %9.1f%%\n",
            sprintf("Nearest ADCIRC + local bias (K=%d)", K), rmse_local,
            100 * (1 - rmse_local / rmse_noaa)))
cat(sprintf("%-40s %8.3f %9.1f%%\n", "LOO regression on nearest ADCIRC", rmse_regkrig,
            100 * (1 - rmse_regkrig / rmse_noaa)))
cat(sprintf("%-40s %8.3f %9.1f%%\n", "Joint 6-dim LMC (paper)", rmse_joint,
            100 * (1 - rmse_joint / rmse_noaa)))
cat("\n")

# Per-site wins
wins_nearest <- sum(abs(rl_nearest - rl_obs) < abs(rl_noaa - rl_obs))
wins_regkrig <- sum(abs(rl_regkrig - rl_obs) < abs(rl_noaa - rl_obs))
wins_joint   <- sum(abs(rl_joint - rl_obs) < abs(rl_noaa - rl_obs))

cat("Sites beating NOAA-only (RL absolute error):\n")
cat(sprintf("  Nearest ADCIRC + global bias: %d/29\n", wins_nearest))
cat(sprintf("  LOO regression:               %d/29\n", wins_regkrig))
cat(sprintf("  Joint LMC:                    %d/29\n", wins_joint))
cat("\n")
