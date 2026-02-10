#!/usr/bin/env Rscript
# Full evfuse pipeline: reproducible from scratch
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
W_tap <- taper_W(bs$W_bs, D, lambda = 300)

# ── 4. Stage 2: spatial model fit + multi-start ──────────────
fit <- fit_spatial_model(stage1, dat, W_tap, D,
                         control = list(maxit = 2000, trace = 1))

# Retry with perturbed rho if convergence fails
max_retries <- 5
attempt <- 1
while (fit$optim_result$convergence != 0 && attempt <= max_retries) {
  cat(sprintf("Retry %d: perturbing rho...\n", attempt))
  set.seed(100 + attempt)
  start_retry <- list(beta = fit$beta, A = fit$A,
                      rho = fit$rho * runif(6, 0.5, 2))
  fit <- fit_spatial_model(stage1, dat, W_tap, D,
                           start = start_retry,
                           control = list(maxit = 2000, trace = 1))
  attempt <- attempt + 1
}

# Multi-start from the converged solution
set.seed(2026)
n_starts <- 20
best_fit <- fit
best_nll <- fit$optim_result$value
cat(sprintf("\nInitial fit NLL: %.4f\n", best_nll))

for (i in seq_len(n_starts)) {
  rho_i <- exp(runif(6, log(50), log(5000)))
  start_i <- list(beta = fit$beta, A = fit$A, rho = rho_i)
  fit_i <- tryCatch(
    suppressWarnings(
      fit_spatial_model(stage1, dat, W_tap, D,
                        start = start_i,
                        control = list(maxit = 2000, trace = 0))
    ),
    error = function(e) NULL
  )
  if (is.null(fit_i)) next
  nll_i <- fit_i$optim_result$value
  cat(sprintf("Start %2d: NLL = %.4f  rho = %s\n",
              i, nll_i, paste(round(fit_i$rho, 1), collapse = " ")))
  if (nll_i < best_nll) {
    best_fit <- fit_i
    best_nll <- nll_i
  }
}
cat(sprintf("\nBest NLL: %.4f\n", best_nll))
cat("Best rho:", round(best_fit$rho, 1), "\n")
saveRDS(best_fit, "data-raw/model_6dim_best.rds")

# ── 5. Kriging ──────────────────────────────────────────────
grid <- read.csv("data-raw/prediction_grid.csv")
preds <- predict_krig(best_fit, grid)
cat("\nKriging at", nrow(grid), "sites:\n")
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

cat("\nDone. All outputs in data-raw/\n")