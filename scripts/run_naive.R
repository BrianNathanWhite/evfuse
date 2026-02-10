#!/usr/bin/env Rscript
# Fit naive baseline models (NOAA-only and ADCIRC-only) for comparison
devtools::load_all()

# ── 1. Load data ─────────────────────────────────────────────
df <- read.csv("data-raw/combined_max.csv")
dat <- load_data(df)
D <- compute_distances(dat$sites)

# ── 2. Stage 1 ───────────────────────────────────────────────
stage1 <- fit_gev_all(dat)
cat("Stage 1:", sum(stage1$converged), "/", dat$n_sites, "sites converged\n")

# ── 3. Bootstrap (same seed as main pipeline) ────────────────
bs <- bootstrap_W(dat, B = 500, seed = 42)
cat("Bootstrap failures:", bs$n_failures, "\n")

# ── 4. Fit NOAA-only model ───────────────────────────────────
cat("\n========================================\n")
cat("NOAA-only model (p=3, L=29)\n")
cat("========================================\n")
set.seed(2026)
noaa_model <- fit_naive_model(stage1, dat, bs$W_bs, D,
                               source = "NOAA", lambda = 300)
cat("  rho:", round(noaa_model$rho, 1), "\n")
cat("  beta:", round(noaa_model$beta, 4), "\n")
saveRDS(noaa_model, "data-raw/model_noaa_only.rds")

# ── 5. Fit ADCIRC-only model ─────────────────────────────────
cat("\n========================================\n")
cat("ADCIRC-only model (p=3, L=100)\n")
cat("========================================\n")
set.seed(2026)
adcirc_model <- fit_naive_model(stage1, dat, bs$W_bs, D,
                                 source = "ADCIRC", lambda = 300)
cat("  rho:", round(adcirc_model$rho, 1), "\n")
cat("  beta:", round(adcirc_model$beta, 4), "\n")
saveRDS(adcirc_model, "data-raw/model_adcirc_only.rds")

# ── 6. Krige to prediction grid ──────────────────────────────
grid <- read.csv("data-raw/prediction_grid.csv")

cat("\n=== Kriging: NOAA-only ===\n")
preds_noaa <- predict_krig_naive(noaa_model, grid)
cat("  mu:        ", round(range(preds_noaa$noaa_mean[, "mu"]), 3), "\n")
cat("  log_sigma: ", round(range(preds_noaa$noaa_mean[, "log_sigma"]), 3), "\n")
cat("  xi:        ", round(range(preds_noaa$noaa_mean[, "xi"]), 3), "\n")
saveRDS(preds_noaa, "data-raw/predictions_noaa_only.rds")

cat("\n=== Kriging: ADCIRC-only ===\n")
preds_adcirc <- predict_krig_naive(adcirc_model, grid)
cat("  mu:        ", round(range(preds_adcirc$noaa_mean[, "mu"]), 3), "\n")
cat("  log_sigma: ", round(range(preds_adcirc$noaa_mean[, "log_sigma"]), 3), "\n")
cat("  xi:        ", round(range(preds_adcirc$noaa_mean[, "xi"]), 3), "\n")
saveRDS(preds_adcirc, "data-raw/predictions_adcirc_only.rds")

# ── 7. Return levels ─────────────────────────────────────────
cat("\n=== 100-year return levels ===\n")
rl_noaa <- compute_return_levels(preds_noaa, r = 100, method = "both",
                                  n_sim = 5000, seed = 123)
rl_adcirc <- compute_return_levels(preds_adcirc, r = 100, method = "both",
                                    n_sim = 5000, seed = 123)
saveRDS(rl_noaa, "data-raw/rl_noaa_only.rds")
saveRDS(rl_adcirc, "data-raw/rl_adcirc_only.rds")

cat(sprintf("NOAA-only:   RL range [%.3f, %.3f], mean SE %.4f\n",
            min(rl_noaa$return_level), max(rl_noaa$return_level),
            mean(rl_noaa$se_delta, na.rm = TRUE)))
cat(sprintf("ADCIRC-only: RL range [%.3f, %.3f], mean SE %.4f\n",
            min(rl_adcirc$return_level), max(rl_adcirc$return_level),
            mean(rl_adcirc$se_delta, na.rm = TRUE)))

# ── 8. Load joint model for comparison ────────────────────────
if (file.exists("data-raw/model_6dim_best.rds")) {
  joint <- readRDS("data-raw/model_6dim_best.rds")
  rl_joint <- readRDS("data-raw/return_levels_100yr.rds")
  cat(sprintf("\nJoint 6-dim: RL range [%.3f, %.3f], mean SE %.4f\n",
              min(rl_joint$return_level), max(rl_joint$return_level),
              mean(rl_joint$se_delta, na.rm = TRUE)))
  cat(sprintf("  NLL: %.4f (33 params, 129 sites)\n", joint$optim_result$value))
}

cat(sprintf("\nNOAA-only:   NLL: %.4f (12 params, 29 sites)\n",
            noaa_model$optim_result$value))
cat(sprintf("ADCIRC-only: NLL: %.4f (12 params, 100 sites)\n",
            adcirc_model$optim_result$value))

cat("\nDone. Outputs saved to data-raw/\n")
