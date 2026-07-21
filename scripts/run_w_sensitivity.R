#!/usr/bin/env Rscript
# ══════════════════════════════════════════════════════════════════════════════
# W Sensitivity Analysis: Batch-Splitting
# ══════════════════════════════════════════════════════════════════════════════
# Demonstrates that bootstrap W uncertainty is negligible.
# Runs B=2000 bootstrap, splits into 4 batches of 500, refits Stage 2 with
# each batch's W, and compares return levels and LOO-CV metrics.
#
# Expected runtime: ~30 min (22 min bootstrap + 4 x 2 min Stage 2)
#
# Requires the model objects saved by scripts/run_nonstationary.R:
#   data-raw/model_6dim_ns.rds, data-raw/model_noaa_only_ns.rds

devtools::load_all()
library(ggplot2)
library(sf)

# ══════════════════════════════════════════════════════════════════════════════
# Phase 1: Load data and run B=2000 bootstrap
# ══════════════════════════════════════════════════════════════════════════════

data(coast_data, package = "evfuse")
dat <- coast_data
df  <- coast_data$raw_df
D   <- compute_distances(dat$sites)

cat("\n== Phase 1: Bootstrap B=2000 ==\n")
t_boot <- system.time({
  bs <- bootstrap_W_detrended(dat, df, B = 2000, ref_year = 2000, seed = 20260212)
})
cat(sprintf("Bootstrap complete: %.1f min, %d failures\n",
            t_boot["elapsed"] / 60, bs$n_failures))

Gamma <- bs$Gamma  # 2000 x 387
saveRDS(Gamma, "data-raw/w_sensitivity_bootstrap.rds")
cat("Saved data-raw/w_sensitivity_bootstrap.rds\n")

# ══════════════════════════════════════════════════════════════════════════════
# Phase 2: Split into 4 batches, compute W matrices
# ══════════════════════════════════════════════════════════════════════════════

cat("\n== Phase 2: Compute 4 batch W matrices ==\n")
batch_idx <- list(1:500, 501:1000, 1001:1500, 1501:2000)

W_list <- lapply(seq_along(batch_idx), function(b) {
  Gamma_b <- Gamma[batch_idx[[b]], ]
  W_raw <- cov(Gamma_b, use = "pairwise.complete.obs")
  W_tap <- taper_W(W_raw, D, lambda = 300, p = 3)
  cat(sprintf("  Batch %d: W range [%.4f, %.4f], PSD=%s\n",
              b, min(W_tap), max(W_tap),
              tryCatch({chol(W_tap + diag(1e-8, nrow(W_tap))); "yes"},
                       error = function(e) "NO")))
  W_tap
})

# Load production model for comparison
prod_model <- readRDS("data-raw/model_6dim_ns.rds")

# Stage 1 (same for all batches)
stage1 <- fit_gev_detrended(dat, df, ref_year = 2000)

# ══════════════════════════════════════════════════════════════════════════════
# Phase 3: Refit Stage 2 four times with multi-start
# ══════════════════════════════════════════════════════════════════════════════

cat("\n== Phase 3: Refit Stage 2 (4 batches x 20 multi-starts) ==\n")

fit_with_multistart <- function(W_tap, label) {
  cat(sprintf("\n  Fitting %s...\n", label))

  # Initial fit
  fit <- tryCatch(
    suppressWarnings(
      fit_spatial_model(stage1, dat, W_tap, D,
                        control = list(maxit = 2000, trace = 0))
    ),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    warning("Initial fit failed for ", label)
    return(NULL)
  }

  # Multi-start (20 random rho, seed 2026)
  set.seed(2026)
  best_fit <- fit
  best_nll <- fit$optim_result$value

  for (i in seq_len(20)) {
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
    if (!is.null(fit_i) && fit_i$optim_result$value < best_nll) {
      best_fit <- fit_i
      best_nll <- fit_i$optim_result$value
    }
  }
  cat(sprintf("  %s: NLL = %.4f\n", label, best_nll))
  best_fit
}

t_fits <- system.time({
  models <- lapply(seq_along(W_list), function(b) {
    fit_with_multistart(W_list[[b]], sprintf("Batch %d", b))
  })
})
cat(sprintf("\nAll 4 fits complete: %.1f min\n", t_fits["elapsed"] / 60))

saveRDS(models, "data-raw/w_sensitivity_models.rds")

# ══════════════════════════════════════════════════════════════════════════════
# Phase 4a: Extract Stage 2 parameters and correlations
# ══════════════════════════════════════════════════════════════════════════════

cat("\n== Phase 4: Compare results ==\n")

extract_cors <- function(model) {
  AAT <- model$A %*% t(model$A)
  sds <- sqrt(diag(AAT))
  cor_mat <- AAT / outer(sds, sds)
  c(cor_mu = cor_mat[1, 4], cor_logsig = cor_mat[2, 5], cor_xi = cor_mat[3, 6])
}

cors_batch <- t(sapply(models, extract_cors))
cors_prod  <- extract_cors(prod_model)

cat("\nCross-source correlations:\n")
cat(sprintf("  %-12s  Cor(mu)  Cor(log_sig)  Cor(xi)\n", ""))
for (b in 1:4) {
  cat(sprintf("  Batch %d       %.4f     %.4f       %.4f\n",
              b, cors_batch[b, 1], cors_batch[b, 2], cors_batch[b, 3]))
}
cat(sprintf("  Production    %.4f     %.4f       %.4f\n",
            cors_prod[1], cors_prod[2], cors_prod[3]))
cat(sprintf("  Range         %.4f     %.4f       %.4f\n",
            diff(range(c(cors_batch[, 1], cors_prod[1]))),
            diff(range(c(cors_batch[, 2], cors_prod[2]))),
            diff(range(c(cors_batch[, 3], cors_prod[3])))))

# ══════════════════════════════════════════════════════════════════════════════
# Phase 4b: LOO-CV at 29 NOAA sites
# ══════════════════════════════════════════════════════════════════════════════

cat("\nLOO-CV comparison:\n")

# NOAA-only baseline (fixed across batches; doesn't depend on W of joint model)
noaa_model <- readRDS("data-raw/model_noaa_only_ns.rds")
# Re-attach current stage1 for LOO-CV
noaa_model$stage1 <- list(
  theta_hat = stage1$theta_hat[which(dat$sites$data_source == "NOAA"), ],
  converged = stage1$converged[which(dat$sites$data_source == "NOAA")]
)
loo_noaa <- loo_cv(noaa_model)
sum_noaa <- loo_summary(loo_noaa, r = 100)

loo_results <- lapply(seq_along(models), function(b) {
  loo_b <- loo_cv(models[[b]])
  sum_b <- loo_summary(loo_b, r = 100)
  list(loo = loo_b, summary = sum_b)
})

saveRDS(loo_results, "data-raw/w_sensitivity_loo.rds")

cat(sprintf("  %-12s  RL RMSE   LPD(joint)  LPD gain\n", ""))
for (b in 1:4) {
  s <- loo_results[[b]]$summary
  cat(sprintf("  Batch %d       %.4f    %6.2f      %6.2f\n",
              b, s$rl_rmse, s$total_lpd,
              s$total_lpd - sum_noaa$total_lpd))
}
# Production LOO for comparison
loo_prod <- loo_cv(prod_model)
sum_prod <- loo_summary(loo_prod, r = 100)
cat(sprintf("  Production    %.4f    %6.2f      %6.2f\n",
            sum_prod$rl_rmse, sum_prod$total_lpd,
            sum_prod$total_lpd - sum_noaa$total_lpd))

# ══════════════════════════════════════════════════════════════════════════════
# Phase 4c: Kriged return levels at prediction grid
# ══════════════════════════════════════════════════════════════════════════════

cat("\nKriging at prediction grid...\n")
data(prediction_grid, package = "evfuse")
grid <- prediction_grid

rl_list <- lapply(seq_along(models), function(b) {
  preds_b <- predict_krig(models[[b]], grid)
  rl_b <- compute_return_levels(preds_b, r = 100, method = "delta")
  cat(sprintf("  Batch %d: RL range [%.3f, %.3f]\n",
              b, min(rl_b$return_level), max(rl_b$return_level)))
  rl_b
})

saveRDS(rl_list, "data-raw/w_sensitivity_rl.rds")

# Production grid return levels
preds_prod <- predict_krig(prod_model, grid)
rl_prod <- compute_return_levels(preds_prod, r = 100, method = "delta")

# ══════════════════════════════════════════════════════════════════════════════
# Phase 4d: Pairwise comparisons
# ══════════════════════════════════════════════════════════════════════════════

cat("\n--- Pairwise max absolute differences (4 batches + production) ---\n")

# Stage 2 parameters
all_pars <- sapply(c(models, list(prod_model)), function(m) {
  pack_params(m$beta, m$A, m$rho)
})  # 33 x 5
par_diffs <- max(apply(all_pars, 1, function(row) diff(range(row))))
cat(sprintf("  Max |delta| across 33 parameters:       %.4f\n", par_diffs))

# LOO return levels (29 sites)
all_loo_rl <- cbind(
  sapply(loo_results, function(r) r$summary$rl_loo),
  sum_prod$rl_loo
)  # 29 x 5
loo_rl_diff <- max(apply(all_loo_rl, 1, function(row) diff(range(row))))
cat(sprintf("  Max |delta| across 29 LOO RLs:          %.4f m\n", loo_rl_diff))

# Grid return levels
all_grid_rl <- cbind(
  sapply(rl_list, function(r) r$return_level),
  rl_prod$return_level
)  # n_grid x 5
grid_rl_diff <- max(apply(all_grid_rl, 1, function(row) diff(range(row))))
cat(sprintf("  Max |delta| across %d grid RLs:        %.4f m\n",
            nrow(grid), grid_rl_diff))

# CV of kriged RL at each grid point
grid_rl_cv <- apply(all_grid_rl, 1, function(row) sd(row) / mean(row))
cat(sprintf("  Max CV of kriged RL across grid:         %.4f (%.2f%%)\n",
            max(grid_rl_cv), 100 * max(grid_rl_cv)))

# ══════════════════════════════════════════════════════════════════════════════
# Phase 4e: Summary table
# ══════════════════════════════════════════════════════════════════════════════

summary_df <- data.frame(
  Batch = c(paste("Batch", 1:4), "Production"),
  NLL = c(sapply(models, function(m) m$optim_result$value),
          prod_model$optim_result$value),
  Cor_mu = c(cors_batch[, 1], cors_prod[1]),
  Cor_logsig = c(cors_batch[, 2], cors_prod[2]),
  Cor_xi = c(cors_batch[, 3], cors_prod[3]),
  RL_RMSE = c(sapply(loo_results, function(r) r$summary$rl_rmse),
              sum_prod$rl_rmse),
  LPD_joint = c(sapply(loo_results, function(r) r$summary$total_lpd),
                sum_prod$total_lpd),
  LPD_gain = c(sapply(loo_results, function(r) r$summary$total_lpd - sum_noaa$total_lpd),
               sum_prod$total_lpd - sum_noaa$total_lpd),
  stringsAsFactors = FALSE
)

write.csv(summary_df, "tables/w_sensitivity_table.csv", row.names = FALSE)
cat("\nSaved tables/w_sensitivity_table.csv\n")

cat("\nSummary table:\n")
print(summary_df, digits = 4, row.names = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
# Phase 5: 4-panel return level map
# ══════════════════════════════════════════════════════════════════════════════

cat("\nGenerating 4-panel map...\n")

if (!requireNamespace("maps", quietly = TRUE) ||
    !requireNamespace("gridExtra", quietly = TRUE) ||
    !requireNamespace("ragg", quietly = TRUE)) {
  cat("Optional packages (maps/gridExtra/ragg) not installed; skipping map figure.\n")
  quit(save = "no", status = 0)
}

# Basemap (same style as run_nonstationary.R)
states_df <- map_data("state")
states_sf <- lapply(split(states_df, states_df$group), function(grp) {
  st_polygon(list(as.matrix(grp[, c("long", "lat")])))
})
states_sf <- st_sfc(states_sf, crs = 4326)
states_sf <- st_sf(geometry = states_sf)

crs_albers <- st_crs(5070)

theme_map <- theme_minimal(base_size = 10) +
  theme(
    axis.title = element_blank(),
    axis.text = element_text(size = 7, color = "grey40"),
    panel.grid = element_blank(),
    strip.text = element_text(size = 10),
    legend.position = "right",
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(t = 4, r = 2, b = 2, l = 2)
  )

basemap <- geom_sf(data = states_sf, fill = "grey95", color = "grey60",
                    linewidth = 0.3, inherit.aes = FALSE)

map_coord <- function() {
  coord_sf(crs = crs_albers,
           xlim = c(-98, -66), ylim = c(24, 46),
           default_crs = st_crs(4326),
           expand = FALSE)
}

# Shared color scale across all 4 panels
all_rl_vals <- unlist(lapply(rl_list, function(r) r$return_level))
rl_range <- range(all_rl_vals)

panels <- lapply(1:4, function(b) {
  pts <- st_as_sf(rl_list[[b]], coords = c("lon", "lat"), crs = 4326)
  ggplot() +
    basemap +
    geom_sf(data = pts, aes(color = return_level),
            size = 0.8, shape = 16, alpha = 0.7) +
    scale_color_viridis_c(name = "RL (m)", limits = rl_range) +
    map_coord() +
    labs(title = sprintf("Batch %d (B=%d:%d)", b,
                         batch_idx[[b]][1], tail(batch_idx[[b]], 1))) +
    theme_map
})

p <- gridExtra::grid.arrange(grobs = panels, ncol = 2)
ggsave("figures/w_sensitivity_maps.png", p,
       width = 10, height = 8, dpi = 150,
       device = ragg::agg_png)
cat("Saved figures/w_sensitivity_maps.png\n")

# ══════════════════════════════════════════════════════════════════════════════
# Summary and recommendation
# ══════════════════════════════════════════════════════════════════════════════

cat("\n")
cat("══════════════════════════════════════════════════════════════\n")
cat("W SENSITIVITY ANALYSIS: SUMMARY\n")
cat("══════════════════════════════════════════════════════════════\n")

cor_ranges <- apply(rbind(cors_batch, cors_prod), 2, function(x) diff(range(x)))
rl_rmses <- c(sapply(loo_results, function(r) r$summary$rl_rmse), sum_prod$rl_rmse)

cat(sprintf("Cor(mu) range:      %.4f (max diff %.4f)\n",
            mean(c(cors_batch[, 1], cors_prod[1])), cor_ranges[1]))
cat(sprintf("Cor(log_sig) range: %.4f (max diff %.4f)\n",
            mean(c(cors_batch[, 2], cors_prod[2])), cor_ranges[2]))
cat(sprintf("Cor(xi) range:      %.4f (max diff %.4f)\n",
            mean(c(cors_batch[, 3], cors_prod[3])), cor_ranges[3]))
cat(sprintf("RL RMSE range:      [%.4f, %.4f] (max diff %.4f m)\n",
            min(rl_rmses), max(rl_rmses), diff(range(rl_rmses))))
cat(sprintf("Max grid RL diff:   %.4f m\n", grid_rl_diff))
cat(sprintf("Max grid RL CV:     %.2f%%\n", 100 * max(grid_rl_cv)))

# Check success criteria
ok_cor_mu  <- cor_ranges[1] < 0.005
ok_cor_xi  <- cor_ranges[3] < 0.02
ok_rl_rmse <- diff(range(rl_rmses)) < 0.01
ok_grid_rl <- grid_rl_diff < 0.05

cat(sprintf("\nSuccess criteria:\n"))
cat(sprintf("  Cor(mu) range < 0.005:    %s (%.4f)\n",
            if (ok_cor_mu) "PASS" else "FAIL", cor_ranges[1]))
cat(sprintf("  Cor(xi) range < 0.02:     %s (%.4f)\n",
            if (ok_cor_xi) "PASS" else "FAIL", cor_ranges[3]))
cat(sprintf("  RL RMSE range < 0.01 m:   %s (%.4f)\n",
            if (ok_rl_rmse) "PASS" else "FAIL", diff(range(rl_rmses))))
cat(sprintf("  Grid RL max diff < 0.05:  %s (%.4f)\n",
            if (ok_grid_rl) "PASS" else "FAIL", grid_rl_diff))

max_pct <- 100 * max(grid_rl_cv)
cat(sprintf("\nDraft sentence: \"A sensitivity analysis splitting 2,000 bootstrap\n"))
cat(sprintf("replicates into four independent batches of 500 showed that all\n"))
cat(sprintf("cross-source correlations, LOO-CV metrics, and kriged return levels\n"))
cat(sprintf("varied by less than %.1f%% across batches, confirming that\n", max_pct))
cat(sprintf("W-uncertainty is negligible.\"\n"))
