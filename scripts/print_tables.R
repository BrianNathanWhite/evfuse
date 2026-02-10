#!/usr/bin/env Rscript
# Print summary tables from fitted models
devtools::load_all()

# ── Load all model outputs ────────────────────────────────────
joint    <- readRDS("data-raw/model_6dim_best.rds")
noaa_m   <- readRDS("data-raw/model_noaa_only.rds")
adcirc_m <- readRDS("data-raw/model_adcirc_only.rds")

rl_joint  <- readRDS("data-raw/return_levels_100yr.rds")
rl_noaa   <- readRDS("data-raw/rl_noaa_only.rds")
rl_adcirc <- readRDS("data-raw/rl_adcirc_only.rds")

grid   <- read.csv("data-raw/prediction_grid.csv")
sites  <- joint$dat$sites
stage1 <- joint$stage1

# ══════════════════════════════════════════════════════════════
# TABLE 1: Model Fit Summary
# ══════════════════════════════════════════════════════════════
cat("\n")
cat("================================================================================\n")
cat("TABLE 1: Model Fit Summary\n")
cat("================================================================================\n")
cat(sprintf("%-15s %5s %6s %10s   %-40s\n", "Model", "Sites", "Params", "NLL", "rho"))
cat(sprintf("%-15s %5s %6s %10s   %-40s\n", "---------------", "-----", "------", "----------", "----------------------------------------"))
cat(sprintf("%-15s %5d %6d %10.4f   %s\n", "Joint (6-dim)", 129, 33,
    joint$optim_result$value, paste(round(joint$rho, 1), collapse = ", ")))
cat(sprintf("%-15s %5d %6d %10.4f   %s\n", "NOAA-only", 29, 12,
    noaa_m$optim_result$value, paste(round(noaa_m$rho, 1), collapse = ", ")))
cat(sprintf("%-15s %5d %6d %10.4f   %s\n", "ADCIRC-only", 100, 12,
    adcirc_m$optim_result$value, paste(round(adcirc_m$rho, 1), collapse = ", ")))
cat("\n")

# ══════════════════════════════════════════════════════════════
# TABLE 2: NOAA Site Comparison
# ══════════════════════════════════════════════════════════════
noaa_idx <- which(sites$data_source == "NOAA")
noaa_sites <- sites[noaa_idx, ]

cat("Kriging at 29 NOAA sites...\n")
preds_joint_noaa <- predict_krig(joint, noaa_sites)
preds_noaa_noaa  <- predict_krig_naive(noaa_m, noaa_sites)

theta_noaa <- stage1$theta_hat[noaa_idx, ]

# Return levels at NOAA sites
rl_joint_noaa <- vapply(seq_len(nrow(noaa_sites)), function(k) {
  gev_return_level(preds_joint_noaa$noaa_mean[k, "mu"],
                   exp(preds_joint_noaa$noaa_mean[k, "log_sigma"]),
                   preds_joint_noaa$noaa_mean[k, "xi"], r = 100)
}, numeric(1))
rl_noaa_noaa <- vapply(seq_len(nrow(noaa_sites)), function(k) {
  gev_return_level(preds_noaa_noaa$noaa_mean[k, "mu"],
                   exp(preds_noaa_noaa$noaa_mean[k, "log_sigma"]),
                   preds_noaa_noaa$noaa_mean[k, "xi"], r = 100)
}, numeric(1))
rl_obs <- vapply(seq_along(noaa_idx), function(k) {
  gev_return_level(theta_noaa[k, "mu"], exp(theta_noaa[k, "log_sigma"]),
                   theta_noaa[k, "xi"], r = 100)
}, numeric(1))

cat("\n")
cat("================================================================================\n")
cat("TABLE 2: Prediction Accuracy at 29 NOAA Sites (kriging vs Stage 1 MLE)\n")
cat("================================================================================\n")
cat(sprintf("%-12s | %12s %15s | %12s %15s\n",
    "Parameter", "Joint RMSE", "NOAA-only RMSE", "Joint MAD", "NOAA-only MAD"))
cat(sprintf("%-12s | %12s %15s | %12s %15s\n",
    "------------", "------------", "---------------", "------------", "---------------"))

params <- c("mu", "log_sigma", "xi")
for (param in params) {
  j <- match(param, params)
  obs <- theta_noaa[, j]
  pred_j <- preds_joint_noaa$noaa_mean[, param]
  pred_n <- preds_noaa_noaa$noaa_mean[, param]

  cat(sprintf("%-12s | %12.4f %15.4f | %12.4f %15.4f\n", param,
      sqrt(mean((pred_j - obs)^2)), sqrt(mean((pred_n - obs)^2)),
      mean(abs(pred_j - obs)), mean(abs(pred_n - obs))))
}
cat(sprintf("%-12s | %12.4f %15.4f | %12.4f %15.4f\n", "100yr RL",
    sqrt(mean((rl_joint_noaa - rl_obs)^2)), sqrt(mean((rl_noaa_noaa - rl_obs)^2)),
    mean(abs(rl_joint_noaa - rl_obs)), mean(abs(rl_noaa_noaa - rl_obs))))
cat("\n")

# ══════════════════════════════════════════════════════════════
# TABLE 3: Mean SE by Region
# ══════════════════════════════════════════════════════════════
atlantic <- grid$lon > -82
gulf     <- grid$lon <= -82

cat("================================================================================\n")
cat("TABLE 3: Mean (Median) SE of 100-Year Return Levels by Region\n")
cat("================================================================================\n")
cat(sprintf("%-12s | %20s | %20s | %20s\n",
    "Region", "Joint", "NOAA-only", "ADCIRC-only"))
cat(sprintf("%-12s | %20s | %20s | %20s\n",
    "------------", "--------------------", "--------------------", "--------------------"))

fmt_se <- function(se) sprintf("%.4f (%.4f)", mean(se, na.rm = TRUE), median(se, na.rm = TRUE))

cat(sprintf("%-12s | %20s | %20s | %20s\n", "Overall",
    fmt_se(rl_joint$se_delta), fmt_se(rl_noaa$se_delta), fmt_se(rl_adcirc$se_delta)))
cat(sprintf("%-12s | %20s | %20s | %20s\n", "Atlantic",
    fmt_se(rl_joint$se_delta[atlantic]), fmt_se(rl_noaa$se_delta[atlantic]),
    fmt_se(rl_adcirc$se_delta[atlantic])))
cat(sprintf("%-12s | %20s | %20s | %20s\n", "Gulf",
    fmt_se(rl_joint$se_delta[gulf]), fmt_se(rl_noaa$se_delta[gulf]),
    fmt_se(rl_adcirc$se_delta[gulf])))
cat("\n")

# ══════════════════════════════════════════════════════════════
# TABLE 4: SE Ratio Summary
# ══════════════════════════════════════════════════════════════
ratio <- rl_noaa$se_delta / rl_joint$se_delta
ratio <- ratio[is.finite(ratio)]

cat("================================================================================\n")
cat("TABLE 4: SE Ratio Summary -- SE(NOAA-only) / SE(joint)\n")
cat("================================================================================\n")
cat(sprintf("  Min:    %.3f\n", min(ratio)))
cat(sprintf("  Q25:    %.3f\n", quantile(ratio, 0.25)))
cat(sprintf("  Median: %.3f\n", median(ratio)))
cat(sprintf("  Mean:   %.3f\n", mean(ratio)))
cat(sprintf("  Q75:    %.3f\n", quantile(ratio, 0.75)))
cat(sprintf("  Max:    %.3f\n", max(ratio)))
cat(sprintf("  > 1.0:  %.1f%% (%d / %d)\n",
    100 * mean(ratio > 1), sum(ratio > 1), length(ratio)))
cat(sprintf("  > 1.5:  %.1f%% (%d / %d)\n",
    100 * mean(ratio > 1.5), sum(ratio > 1.5), length(ratio)))
cat(sprintf("  > 2.0:  %.1f%% (%d / %d)\n",
    100 * mean(ratio > 2), sum(ratio > 2), length(ratio)))
cat("\n")
