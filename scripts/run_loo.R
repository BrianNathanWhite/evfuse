#!/usr/bin/env Rscript
# LOO-CV comparison: Joint vs NOAA-only at 29 NOAA sites
devtools::load_all()

# ── Load fitted models ──────────────────────────────────────
joint_model <- readRDS("data-raw/model_6dim_best.rds")
noaa_model  <- readRDS("data-raw/model_noaa_only.rds")

# ── Compute LOO-CV ──────────────────────────────────────────
cat("Computing LOO-CV for joint model (387 obs)...\n")
loo_joint <- loo_cv(joint_model)

cat("Computing LOO-CV for NOAA-only model (87 obs)...\n")
loo_noaa <- loo_cv(noaa_model)

# ── Summaries ───────────────────────────────────────────────
sum_joint <- loo_summary(loo_joint, r = 100)
sum_noaa  <- loo_summary(loo_noaa, r = 100)

# ══════════════════════════════════════════════════════════════
# TABLE: LOO-CV Parameter Accuracy
# ══════════════════════════════════════════════════════════════
cat("\n")
cat("════════════════════════════════════════════════════════════\n")
cat("LOO-CV: Parameter Prediction Accuracy at 29 NOAA Sites\n")
cat("════════════════════════════════════════════════════════════\n")
cat(sprintf("%-12s %-10s %10s %10s\n", "Parameter", "Model", "RMSE", "MAD"))
cat(sprintf("%-12s %-10s %10s %10s\n", "---------", "-----", "----", "---"))

for (j in 1:3) {
  cat(sprintf("%-12s %-10s %10.4f %10.4f\n",
      sum_joint$param_stats$parameter[j], "Joint",
      sum_joint$param_stats$rmse[j], sum_joint$param_stats$mad[j]))
  cat(sprintf("%-12s %-10s %10.4f %10.4f\n",
      "", "NOAA-only",
      sum_noaa$param_stats$rmse[j], sum_noaa$param_stats$mad[j]))
}
cat("\n")

# ══════════════════════════════════════════════════════════════
# TABLE: LOO-CV Return Level Accuracy
# ══════════════════════════════════════════════════════════════
cat("════════════════════════════════════════════════════════════\n")
cat("LOO-CV: 100-Year Return Level Accuracy\n")
cat("════════════════════════════════════════════════════════════\n")
cat(sprintf("%-15s %10s %10s %10s\n", "Model", "RMSE", "MAD", "Mean SE"))
cat(sprintf("%-15s %10s %10s %10s\n", "-----", "----", "---", "-------"))
cat(sprintf("%-15s %10.4f %10.4f %10.4f\n",
    "Joint", sum_joint$rl_rmse, sum_joint$rl_mad, sum_joint$rl_mean_se))
cat(sprintf("%-15s %10.4f %10.4f %10.4f\n",
    "NOAA-only", sum_noaa$rl_rmse, sum_noaa$rl_mad, sum_noaa$rl_mean_se))
cat("\n")

# ══════════════════════════════════════════════════════════════
# TABLE: LOO-CV Log Predictive Density
# ══════════════════════════════════════════════════════════════
cat("════════════════════════════════════════════════════════════\n")
cat("LOO-CV: Log Predictive Density (higher = better)\n")
cat("════════════════════════════════════════════════════════════\n")
cat(sprintf("%-15s %12s %12s\n", "Model", "Total LPD", "Mean LPD"))
cat(sprintf("%-15s %12s %12s\n", "-----", "---------", "--------"))
cat(sprintf("%-15s %12.2f %12.4f\n",
    "Joint", sum_joint$total_lpd, sum_joint$mean_lpd))
cat(sprintf("%-15s %12.2f %12.4f\n",
    "NOAA-only", sum_noaa$total_lpd, sum_noaa$mean_lpd))
cat("\n")

# ══════════════════════════════════════════════════════════════
# Per-site comparison
# ══════════════════════════════════════════════════════════════
cat("════════════════════════════════════════════════════════════\n")
cat("Per-Site LOO LPD Comparison\n")
cat("════════════════════════════════════════════════════════════\n")
n_joint_better <- sum(loo_joint$loo_lpd > loo_noaa$loo_lpd)
cat(sprintf("Joint better at %d / %d sites (%.0f%%)\n",
    n_joint_better, 29, 100 * n_joint_better / 29))
cat(sprintf("Mean LPD difference (Joint - NOAA): %.4f\n",
    mean(loo_joint$loo_lpd - loo_noaa$loo_lpd)))
cat("\n")

# ── Save results ────────────────────────────────────────────
loo_results <- list(
  joint = loo_joint,
  noaa  = loo_noaa,
  summary_joint = sum_joint,
  summary_noaa  = sum_noaa
)
saveRDS(loo_results, "data-raw/loo_cv_results.rds")
cat("Saved data-raw/loo_cv_results.rds\n")
