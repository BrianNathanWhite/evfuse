#!/usr/bin/env Rscript
# ══════════════════════════════════════════════════════════════════════════════
# Delta Method vs. Simulation SE Comparison for 100-year Return Levels
# ══════════════════════════════════════════════════════════════════════════════
# Compares return level SEs from the delta method and Monte Carlo simulation
# at (1) 29 NOAA LOO-CV sites and (2) the full coastal prediction grid.

devtools::load_all()
library(ggplot2)

n_sim <- 10000
seed  <- 123

# ══════════════════════════════════════════════════════════════════════════════
# 1. LOO-CV predictions at 29 NOAA sites
# ══════════════════════════════════════════════════════════════════════════════

loo_data <- readRDS("data-raw/loo_cv_ns.rds")
loo <- loo_data$joint  # evfuse_loo object with loo_mean, loo_cov, sites

# Wrap LOO predictions into a predictions-like object for compute_return_levels
noaa_preds <- structure(list(
  noaa_mean = loo$loo_mean,
  noaa_cov  = loo$loo_cov,
  new_sites = loo$sites
), class = "evfuse_predictions")

rl_noaa <- compute_return_levels(noaa_preds, r = 100, method = "both",
                                  n_sim = n_sim, seed = seed)

# ══════════════════════════════════════════════════════════════════════════════
# 2. Grid predictions
# ══════════════════════════════════════════════════════════════════════════════

preds_grid <- readRDS("data-raw/predictions_grid_ns.rds")
rl_grid <- compute_return_levels(preds_grid, r = 100, method = "both",
                                  n_sim = n_sim, seed = seed)

# ══════════════════════════════════════════════════════════════════════════════
# 3. Comparison statistics
# ══════════════════════════════════════════════════════════════════════════════

compare <- function(label, se_delta, se_sim, xi_hat) {
  abs_diff <- abs(se_delta - se_sim)
  rel_diff <- abs_diff / pmax(se_delta, 1e-8)
  ratio    <- se_sim / se_delta

  cat(sprintf("\n--- %s ---\n", label))
  cat(sprintf("  Max |delta - sim|:    %.4f\n", max(abs_diff)))
  cat(sprintf("  Max relative diff:    %.1f%%\n", 100 * max(rel_diff)))
  cat(sprintf("  Median relative diff: %.1f%%\n", 100 * median(rel_diff)))
  cat(sprintf("  Correlation:          %.4f\n", cor(se_delta, se_sim)))
  cat(sprintf("  Median ratio (sim/delta): %.3f\n", median(ratio)))
  cat(sprintf("  Fraction sim > delta: %.1f%%\n", 100 * mean(se_sim > se_delta)))

  # Divergent points (>5% relative difference)
  div_idx <- which(rel_diff > 0.05)
  if (length(div_idx) > 0) {
    cat(sprintf("  Divergent points (>5%%): %d of %d\n", length(div_idx), length(se_delta)))
    cat(sprintf("  xi at divergent points: range [%.3f, %.3f], median %.3f\n",
                min(xi_hat[div_idx]), max(xi_hat[div_idx]), median(xi_hat[div_idx])))
    cat(sprintf("  xi at all points:       range [%.3f, %.3f], median %.3f\n",
                min(xi_hat), max(xi_hat), median(xi_hat)))
  } else {
    cat("  Divergent points (>5%): none\n")
  }

  list(abs_diff = abs_diff, rel_diff = rel_diff, ratio = ratio, div_idx = div_idx)
}

# NOAA sites
xi_noaa <- loo$loo_mean[, "xi"]
noaa_stats <- compare("29 NOAA sites (LOO-CV)", rl_noaa$se_delta, rl_noaa$se_sim, xi_noaa)

# Grid points
xi_grid <- preds_grid$noaa_mean[, "xi"]
grid_stats <- compare("Prediction grid", rl_grid$se_delta, rl_grid$se_sim, xi_grid)

# ══════════════════════════════════════════════════════════════════════════════
# 4. Figure: scatter plots
# ══════════════════════════════════════════════════════════════════════════════

df_noaa <- data.frame(
  delta = rl_noaa$se_delta,
  sim   = rl_noaa$se_sim,
  xi    = xi_noaa,
  panel = "29 NOAA sites (LOO-CV)"
)
df_grid <- data.frame(
  delta = rl_grid$se_delta,
  sim   = rl_grid$se_sim,
  xi    = xi_grid,
  panel = "Prediction grid"
)
df_all <- rbind(df_noaa, df_grid)

p <- ggplot(df_all, aes(x = delta, y = sim, color = xi)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
  geom_point(size = 1.5, alpha = 0.7) +
  scale_color_viridis_c(name = expression(xi)) +
  facet_wrap(~panel, scales = "free") +
  labs(x = "Delta method SE (m)", y = "Simulation SE (m)",
       title = "100-year Return Level SE: Delta Method vs. Simulation") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(size = 12),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA))

ggsave("figures/se_comparison.png", p, width = 10, height = 5, dpi = 300,
       device = ragg::agg_png)
cat("\nSaved figures/se_comparison.png\n")

# ══════════════════════════════════════════════════════════════════════════════
# 5. Save summary + recommendation
# ══════════════════════════════════════════════════════════════════════════════

out <- c(
  "SE COMPARISON: Delta Method vs. Simulation (100-year Return Levels)",
  sprintf("Monte Carlo draws: %d, seed: %d", n_sim, seed),
  "",
  "--- 29 NOAA sites (LOO-CV) ---",
  sprintf("Max |delta - sim|:       %.4f m", max(noaa_stats$abs_diff)),
  sprintf("Max relative diff:       %.1f%%", 100 * max(noaa_stats$rel_diff)),
  sprintf("Median relative diff:    %.1f%%", 100 * median(noaa_stats$rel_diff)),
  sprintf("Correlation:             %.4f", cor(rl_noaa$se_delta, rl_noaa$se_sim)),
  sprintf("Median ratio (sim/delta): %.3f", median(noaa_stats$ratio)),
  "",
  "--- Prediction grid (%d points) ---",
  sprintf("Max |delta - sim|:       %.4f m", max(grid_stats$abs_diff)),
  sprintf("Max relative diff:       %.1f%%", 100 * max(grid_stats$rel_diff)),
  sprintf("Median relative diff:    %.1f%%", 100 * median(grid_stats$rel_diff)),
  sprintf("Correlation:             %.4f", cor(rl_grid$se_delta, rl_grid$se_sim)),
  sprintf("Median ratio (sim/delta): %.3f", median(grid_stats$ratio)),
  sprintf("Fraction sim > delta:    %.1f%%", 100 * mean(rl_grid$se_sim > rl_grid$se_delta)),
  ""
)

# Fix the grid points line
out[11] <- sprintf("--- Prediction grid (%d points) ---", nrow(rl_grid))

# Recommendation
max_rel <- max(max(noaa_stats$rel_diff), max(grid_stats$rel_diff))
if (max_rel < 0.05) {
  rec <- "RECOMMENDATION: Delta method and simulation agree closely (max relative diff < 5%); either is fine."
} else {
  sim_larger <- mean(rl_grid$se_sim > rl_grid$se_delta)
  if (sim_larger > 0.6) {
    rec <- sprintf(
      "RECOMMENDATION: Simulation SEs are systematically larger (%.0f%% of grid points). They capture nonlinearity in xi that the delta method misses. Consider using simulation SEs.",
      100 * sim_larger)
  } else {
    rec <- sprintf(
      "RECOMMENDATION: Methods diverge at some points (max rel diff %.1f%%) but neither is systematically larger. Delta method is adequate for most purposes.",
      100 * max_rel)
  }
}
out <- c(out, rec)

cat("\n", rec, "\n")

writeLines(out, "tables/se_comparison.txt")
cat("Saved tables/se_comparison.txt\n")
