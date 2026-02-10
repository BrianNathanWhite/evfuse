#!/usr/bin/env Rscript
# Full evfuse pipeline: Stage 1 → Bootstrap → Stage 2 → Kriging → Return levels

devtools::load_all()

# ── 1. Load data ─────────────────────────────────────────────
df <- read.csv("data-raw/combined_max.csv")
dat <- load_data(df)
D <- compute_distances(dat$sites)

# ── 2. Stage 1: per-site GEV fits ───────────────────────────
stage1 <- fit_gev_all(dat)
cat("Stage 1:", sum(stage1$converged), "/", dat$n_sites, "sites converged\n")

# ── 3. Bootstrap W ──────────────────────────────────────────
bs <- bootstrap_W(dat, B = 500, seed = 42)
cat("Bootstrap failures:", bs$n_failures, "\n")

# Taper with lambda = 500 km
W_tap <- taper_W(bs$W_bs, D, lambda = 500)

# ── 4. Stage 2: spatial model fit ───────────────────────────
fit <- fit_spatial_model(stage1, dat, W_tap, D,
                         control = list(maxit = 2000, trace = 1))

cat("\nStage 2 converged:", fit$optim_result$convergence == 0, "\n")
cat("NLL:", fit$optim_result$value, "\n")
cat("Beta:", round(fit$beta, 4), "\n")
cat("Rho:", round(fit$rho, 1), "\n")

saveRDS(fit, "data-raw/model_6dim.rds")

# ── 5. Kriging on prediction grid ───────────────────────────
grid <- read.csv("data-raw/prediction_grid.csv")
preds <- predict_krig(fit, grid)

cat("\nKriging predictions (NOAA params) at", nrow(grid), "sites:\n")
cat("  mu:        ", round(range(preds$noaa_mean[,"mu"]), 3), "\n")
cat("  log_sigma: ", round(range(preds$noaa_mean[,"log_sigma"]), 3), "\n")
cat("  xi:        ", round(range(preds$noaa_mean[,"xi"]), 3), "\n")

saveRDS(preds, "data-raw/predictions_grid.rds")

# ── 6. Return levels ────────────────────────────────────────
rl_100 <- compute_return_levels(preds, r = 100, method = "both",
                                 n_sim = 5000, seed = 123)

cat("\n100-year return levels:\n")
cat("  Range:", round(range(rl_100$return_level), 3), "\n")
cat("  Mean SE (delta):", round(mean(rl_100$se_delta, na.rm = TRUE), 4), "\n")

saveRDS(rl_100, "data-raw/return_levels_100yr.rds")

cat("\nDone. Outputs saved to data-raw/\n")
