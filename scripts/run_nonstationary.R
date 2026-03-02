#!/usr/bin/env Rscript
# Full nonstationary pipeline: nonstationary NOAA GEV (Stage 1), joint Stage 2,
# baselines, kriging, return levels, LOO-CV, figures, block CV, taper sensitivity,
# xi bootstrap QQ, PIT calibration, and manuscript numbers summary.
#
# Usage: Rscript scripts/run_nonstationary.R
# Expected runtime: ~15 minutes

devtools::load_all()
library(ggplot2)
library(sf)
library(gridExtra)

dir.create("figures", showWarnings = FALSE)
dir.create("tables", showWarnings = FALSE)


# ══════════════════════════════════════════════════════════════════════════════
# Section 1: Stage 1 + Bootstrap (nonstationary)
# ══════════════════════════════════════════════════════════════════════════════
cat("\n== Section 1: Nonstationary Stage 1 + Bootstrap ==\n")

data(coast_data, package = "evfuse")
dat <- coast_data
df  <- coast_data$raw_df
D   <- compute_distances(dat$sites)
L   <- dat$n_sites  # 129

stage1 <- fit_gev_detrended(dat, df, ref_year = 2000)
cat(sprintf("Stage 1 converged: %d / %d\n", sum(stage1$converged), L))

cat("Running bootstrap (B=500, seed=42)...\n")
bs <- bootstrap_W_detrended(dat, df, B = 500, ref_year = 2000, seed = 42)
cat(sprintf("Bootstrap failures: %d\n", bs$n_failures))
W_bs  <- bs$W_bs
W_tap <- taper_W(W_bs, D, lambda = 300)


# ══════════════════════════════════════════════════════════════════════════════
# Section 2: Stage 2 Joint Model + Multi-start
# ══════════════════════════════════════════════════════════════════════════════
cat("\n== Section 2: Joint 6-dim model (nonstationary) ==\n")

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

# Multi-start
set.seed(2026)
best_fit <- fit
best_nll <- fit$optim_result$value
cat(sprintf("\nInitial NLL: %.4f\n", best_nll))

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
  if (is.null(fit_i)) next
  nll_i <- fit_i$optim_result$value
  cat(sprintf("Start %2d: NLL = %.4f  rho = %s\n",
              i, nll_i, paste(round(fit_i$rho, 1), collapse = " ")))
  if (nll_i < best_nll) {
    best_fit <- fit_i
    best_nll <- nll_i
  }
}
cat(sprintf("\nBest NLL (nonstationary): %.4f\n", best_nll))
cat("Best rho:", round(best_fit$rho, 1), "\n")
saveRDS(best_fit, "data-raw/model_6dim_ns.rds")
cat("Saved data-raw/model_6dim_ns.rds\n")


# ══════════════════════════════════════════════════════════════════════════════
# Section 3: NOAA-only + ADCIRC-only Baselines
# ══════════════════════════════════════════════════════════════════════════════
cat("\n== Section 3: Baseline models ==\n")

# NOAA-only (nonstationary stage1 + bootstrap)
cat("Fitting NOAA-only model...\n")
noaa_model <- fit_naive_model(stage1, dat, W_bs, D,
                               source = "NOAA", lambda = 300, n_starts = 20,
                               control = list(maxit = 2000, trace = 0))
cat(sprintf("NOAA-only NLL: %.4f\n", noaa_model$optim_result$value))
saveRDS(noaa_model, "data-raw/model_noaa_only_ns.rds")
cat("Saved data-raw/model_noaa_only_ns.rds\n")

# ADCIRC-only (stationary Stage 1 + bootstrap, same for all pipelines)
cat("Fitting ADCIRC-only model...\n")
adcirc_model <- fit_naive_model(stage1, dat, W_bs, D,
                                 source = "ADCIRC", lambda = 300, n_starts = 20,
                                 control = list(maxit = 2000, trace = 0))
cat(sprintf("ADCIRC-only NLL: %.4f\n", adcirc_model$optim_result$value))
saveRDS(adcirc_model, "data-raw/model_adcirc_only.rds")


# ══════════════════════════════════════════════════════════════════════════════
# Section 4: Kriging + Return Levels
# ══════════════════════════════════════════════════════════════════════════════
cat("\n== Section 4: Kriging + Return Levels ==\n")

data(prediction_grid, package = "evfuse")
grid <- prediction_grid

# Joint model predictions
preds_joint <- predict_krig(best_fit, grid)
cat(sprintf("Joint kriging at %d grid points\n", nrow(grid)))
saveRDS(preds_joint, "data-raw/predictions_grid_ns.rds")

rl_joint <- compute_return_levels(preds_joint, r = 100, method = "both",
                                   n_sim = 5000, seed = 123)
cat(sprintf("Joint 100-yr RL range: [%.3f, %.3f]\n",
            min(rl_joint$return_level), max(rl_joint$return_level)))
saveRDS(rl_joint, "data-raw/return_levels_ns.rds")

# NOAA-only predictions
preds_noaa <- predict_krig_naive(noaa_model, grid)
saveRDS(preds_noaa, "data-raw/predictions_noaa_only_ns.rds")

rl_noaa <- compute_return_levels(preds_noaa, r = 100, method = "both",
                                  n_sim = 5000, seed = 123)
cat(sprintf("NOAA-only 100-yr RL range: [%.3f, %.3f]\n",
            min(rl_noaa$return_level), max(rl_noaa$return_level)))
saveRDS(rl_noaa, "data-raw/rl_noaa_only_ns.rds")

# ADCIRC-only predictions
preds_adcirc <- predict_krig_naive(adcirc_model, grid)
saveRDS(preds_adcirc, "data-raw/predictions_adcirc_only.rds")

rl_adcirc <- compute_return_levels(preds_adcirc, r = 100, method = "both",
                                    n_sim = 5000, seed = 123)
cat(sprintf("ADCIRC-only 100-yr RL range: [%.3f, %.3f]\n",
            min(rl_adcirc$return_level), max(rl_adcirc$return_level)))
saveRDS(rl_adcirc, "data-raw/rl_adcirc_only.rds")


# ══════════════════════════════════════════════════════════════════════════════
# Section 5: LOO-CV
# ══════════════════════════════════════════════════════════════════════════════
cat("\n== Section 5: LOO-CV ==\n")

loo_joint <- loo_cv(best_fit)
sum_joint <- loo_summary(loo_joint, r = 100)

loo_noaa <- loo_cv(noaa_model)
sum_noaa <- loo_summary(loo_noaa, r = 100)

n_wins <- sum(loo_joint$loo_lpd > loo_noaa$loo_lpd)
lpd_gain <- sum_joint$total_lpd - sum_noaa$total_lpd
rl_pct <- 100 * (1 - sum_joint$rl_rmse / sum_noaa$rl_rmse)

cat(sprintf("Joint:     RL RMSE = %.4f, Total LPD = %.2f\n",
            sum_joint$rl_rmse, sum_joint$total_lpd))
cat(sprintf("NOAA-only: RL RMSE = %.4f, Total LPD = %.2f\n",
            sum_noaa$rl_rmse, sum_noaa$total_lpd))
cat(sprintf("Gain: LPD +%.2f, RL RMSE reduction %.1f%%, wins %d/29\n",
            lpd_gain, rl_pct, n_wins))

# Save LOO-CV results
saveRDS(list(joint = loo_joint, noaa = loo_noaa,
             summary_joint = sum_joint, summary_noaa = sum_noaa),
        "data-raw/loo_cv_ns.rds")
cat("Saved data-raw/loo_cv_ns.rds\n")


# ══════════════════════════════════════════════════════════════════════════════
# Section 6: Generate 10 Figures
# ══════════════════════════════════════════════════════════════════════════════
cat("\n== Section 6: Generating figures ==\n")

sites <- best_fit$dat$sites

# ── Basemap ────────────────────────────────────────────────────────────────
states_df <- map_data("state")
states_sf <- lapply(split(states_df, states_df$group), function(grp) {
  st_polygon(list(as.matrix(grp[, c("long", "lat")])))
})
states_sf <- st_sfc(states_sf, crs = 4326)
states_sf <- st_sf(geometry = states_sf)

crs_albers <- st_crs(5070)

theme_map <- theme_minimal(base_size = 14) +
  theme(
    axis.title = element_blank(),
    axis.text = element_text(size = 10, color = "grey40"),
    panel.grid = element_blank(),
    strip.text = element_text(size = 14),
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(t = 8, r = 5, b = 5, l = 5)
  )

map_coord <- function() {
  coord_sf(crs = crs_albers,
           xlim = c(-98, -66), ylim = c(24, 46),
           default_crs = st_crs(4326),
           expand = FALSE)
}

basemap <- geom_sf(data = states_sf, fill = "grey95", color = "grey60",
                    linewidth = 0.3, inherit.aes = FALSE)


# ── Figure 1: Stage 1 QQ plots (2x2, base R) ──────────────────────────────
cat("  Figure 1: stage1_qq_plots.png\n")

noaa_idx   <- which(sites$data_source == "NOAA")
adcirc_idx <- which(sites$data_source == "ADCIRC")

noaa_xi  <- stage1$theta_hat[noaa_idx, "xi"]
noaa_lat <- sites$lat[noaa_idx]

# (a) Gulf Coast high xi
gulf_mask <- noaa_lat < 31
site_a <- noaa_idx[gulf_mask][which.max(noaa_xi[gulf_mask])]
# (b) Mid-Atlantic xi near zero
midatl_mask <- noaa_lat >= 35 & noaa_lat <= 41
site_b <- noaa_idx[midatl_mask][which.min(abs(noaa_xi[midatl_mask]))]
# (c) New England
ne_mask <- noaa_lat > 41
site_c <- noaa_idx[ne_mask][which.max(noaa_xi[ne_mask])]
# (d) ADCIRC median xi
adcirc_xi <- stage1$theta_hat[adcirc_idx, "xi"]
site_d <- adcirc_idx[which.min(abs(adcirc_xi - median(adcirc_xi)))]

sites_qq <- c(site_a, site_b, site_c, site_d)
qq_labels <- c("NOAA Gulf Coast (high xi)", "NOAA Mid-Atlantic (xi ~ 0)",
               "NOAA New England", "ADCIRC")

gev_quantile <- function(p, mu, sigma, xi) {
  yp <- -log(p)
  if (abs(xi) < 1e-8) mu - sigma * log(yp)
  else mu + (sigma / xi) * (yp^(-xi) - 1)
}

png("figures/stage1_qq_plots.png", width = 2400, height = 2400, res = 300,
    bg = "white")
par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3.5, 1))

for (k in seq_along(sites_qq)) {
  i <- sites_qq[k]
  x <- sort(dat$maxima[[i]])
  n_obs <- length(x)
  is_noaa <- sites$data_source[i] == "NOAA"

  if (is_noaa) {
    # Nonstationary fit: extract mu0, scale, shape
    fit_k <- stage1$fits[[i]]
    pars <- fit_k$results$par
    mu <- pars["mu0"]
    sigma <- pars["scale"]
    xi <- pars["shape"]
  } else {
    fit_k <- stage1$fits[[i]]
    pars <- fit_k$results$par
    mu <- pars["location"]
    sigma <- pars["scale"]
    xi <- pars["shape"]
  }

  pp <- (seq_len(n_obs) - 0.5) / n_obs
  theo <- sapply(pp, function(p) gev_quantile(p, mu, sigma, xi))

  rng <- range(c(x, theo))
  pad <- diff(rng) * 0.05

  plot(theo, x, pch = 16, cex = 0.9, col = "steelblue",
       xlab = "GEV theoretical quantiles (m)",
       ylab = "Observed quantiles (m)",
       xlim = rng + c(-pad, pad), ylim = rng + c(-pad, pad),
       main = sprintf("%s\n%s", qq_labels[k], sites$location[i]),
       cex.main = 0.95, cex.lab = 1.0, cex.axis = 0.9)
  abline(0, 1, col = "red", lwd = 1.5)

  mu_label <- if (is_noaa) "mu0" else "mu"
  legend("topleft",
         legend = sprintf("%s=%.2f  sigma=%.3f  xi=%.3f", mu_label, mu, sigma, xi),
         bty = "n", cex = 0.8)
}
dev.off()
cat("  Saved figures/stage1_qq_plots.png\n")


# ── Figure 2: Stage 1 MLE map (3-panel, mu0 label) ────────────────────────
cat("  Figure 2: stage1_mle_map.png\n")

s1_df <- data.frame(
  lon = sites$lon, lat = sites$lat,
  source = sites$data_source,
  mu = stage1$theta_hat[, "mu"],
  log_sigma = stage1$theta_hat[, "log_sigma"],
  xi = stage1$theta_hat[, "xi"]
)

s1_noaa_df   <- s1_df[s1_df$source == "NOAA", ]
s1_adcirc_df <- s1_df[s1_df$source == "ADCIRC", ]

make_s1_panel <- function(param_col, label, fill_scale, color_scale,
                          show_source_legend = FALSE,
                          show_y = TRUE) {
  p <- ggplot() +
    geom_polygon(data = states_df,
                 aes(x = long, y = lat, group = group),
                 fill = "gray90", colour = "white", linewidth = 0.3) +
    geom_point(data = s1_adcirc_df,
               aes(x = lon, y = lat,
                   fill = .data[[param_col]], color = .data[[param_col]]),
               shape = 24, size = 2.5, stroke = 0.4, alpha = 0.85) +
    geom_point(data = s1_noaa_df,
               aes(x = lon, y = lat, fill = .data[[param_col]]),
               shape = 21, size = 4, stroke = 0.9, color = "black", alpha = 0.9) +
    fill_scale + color_scale +
    coord_cartesian(xlim = c(-98, -66), ylim = c(24, 46), expand = FALSE) +
    labs(title = label, x = "Longitude",
         y = if (show_y) "Latitude" else NULL, fill = NULL) +
    guides(colour = "none") +
    theme_minimal(base_size = 12) +
    theme(
      plot.background   = element_rect(fill = "white", colour = NA),
      panel.background  = element_rect(fill = "white", colour = NA),
      panel.grid.major  = element_line(colour = "gray92", linewidth = 0.3),
      panel.grid.minor  = element_blank(),
      legend.key.height = unit(1.2, "cm"),
      plot.title        = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.margin       = margin(t = 5, r = 5, b = 5, l = if (show_y) 5 else 2)
    )

  if (!show_y) {
    p <- p + theme(axis.text.y = element_blank(),
                   axis.ticks.y = element_blank())
  }

  if (show_source_legend) {
    legend_df <- data.frame(
      x = c(-80, -80), y = c(35, 35),
      Source = factor(c("NOAA", "ADCIRC"), levels = c("NOAA", "ADCIRC"))
    )
    p <- p +
      geom_point(data = legend_df, aes(x = x, y = y, shape = Source),
                 size = NA_real_, alpha = 0) +
      scale_shape_manual(values = c("NOAA" = 21, "ADCIRC" = 24), name = NULL) +
      guides(shape = guide_legend(override.aes = list(
        size = c(4, 2.5), alpha = 1,
        fill = "grey50", color = c("black", "grey50"),
        stroke = c(0.9, 0.4)
      )))
  }
  p
}

mu_lims <- range(s1_df$mu)
ls_lims <- range(s1_df$log_sigma)
xi_lims <- range(s1_df$xi)

# Label mu as mu_0 for nonstationary (Unicode for clean rendering)
p1_mu  <- make_s1_panel("mu", "\u03bc\u2080",
           scale_fill_viridis_c(option = "inferno", begin = 0.1, end = 0.95,
                                limits = mu_lims),
           scale_color_viridis_c(option = "inferno", begin = 0.1, end = 0.95,
                                 limits = mu_lims))
p1_sig <- make_s1_panel("log_sigma", "log \u03c3",
           scale_fill_distiller(palette = "YlGnBu", direction = 1,
                                limits = ls_lims),
           scale_color_distiller(palette = "YlGnBu", direction = 1,
                                 limits = ls_lims),
           show_y = FALSE)
xi_sym <- c(-1, 1) * max(abs(xi_lims))
p1_xi  <- make_s1_panel("xi", "\u03be",
           scale_fill_gradient2(low = "#2166AC", mid = "grey95", high = "#B2182B",
                                midpoint = 0, limits = xi_sym),
           scale_color_gradient2(low = "#2166AC", mid = "grey95", high = "#B2182B",
                                 midpoint = 0, limits = xi_sym),
           show_source_legend = TRUE)

# Compose 2-over-1 layout: mu0 and log_sigma on top, xi centered below
top_row <- gridExtra::arrangeGrob(p1_mu, p1_sig, ncol = 2)
bot_row <- gridExtra::arrangeGrob(grid::nullGrob(), p1_xi, grid::nullGrob(),
                                   ncol = 3, widths = c(1, 2, 1))
s1_combined <- gridExtra::arrangeGrob(top_row, bot_row, nrow = 2)
ggsave("figures/stage1_mle_map.png", s1_combined, width = 20, height = 16,
       dpi = 300, bg = "white", device = ragg::agg_png)
cat("  Saved figures/stage1_mle_map.png\n")


# ── Figure 3: Kriged parameter map (3-panel, mu0,N label) ─────────────────
cat("  Figure 3: kriged_params_map.png\n")

make_krig_panel <- function(param_col, label, scale, show_y = TRUE) {
  kdf <- data.frame(lon = grid$lon, lat = grid$lat,
                     value = preds_joint$noaa_mean[, param_col])

  p <- ggplot() +
    geom_polygon(data = states_df,
                 aes(x = long, y = lat, group = group),
                 fill = "gray90", colour = "white", linewidth = 0.3) +
    geom_point(data = kdf, aes(x = lon, y = lat, color = value),
               size = 1.2, shape = 16, alpha = 0.7) +
    scale +
    coord_cartesian(xlim = c(-98, -66), ylim = c(24, 46), expand = FALSE) +
    labs(title = label, x = "Longitude",
         y = if (show_y) "Latitude" else NULL, color = NULL) +
    theme_minimal(base_size = 12) +
    theme(
      plot.background   = element_rect(fill = "white", colour = NA),
      panel.background  = element_rect(fill = "white", colour = NA),
      panel.grid.major  = element_line(colour = "gray92", linewidth = 0.3),
      panel.grid.minor  = element_blank(),
      legend.key.height = unit(1.2, "cm"),
      plot.title        = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.margin       = margin(t = 5, r = 5, b = 5, l = if (show_y) 5 else 2)
    )

  if (!show_y) {
    p <- p + theme(axis.text.y = element_blank(),
                   axis.ticks.y = element_blank())
  }
  p
}

krig_mu_lims <- range(preds_joint$noaa_mean[, "mu"])
krig_ls_lims <- range(preds_joint$noaa_mean[, "log_sigma"])
krig_xi_lims <- range(preds_joint$noaa_mean[, "xi"])

p2_mu  <- make_krig_panel("mu", "\u03bc\u2080",
           scale_color_viridis_c(option = "inferno", begin = 0.1, end = 0.95,
                                 limits = krig_mu_lims))
p2_sig <- make_krig_panel("log_sigma", "log \u03c3",
           scale_color_distiller(palette = "YlGnBu", direction = 1,
                                 limits = krig_ls_lims),
           show_y = FALSE)
krig_xi_sym <- c(-1, 1) * max(abs(krig_xi_lims))
p2_xi  <- make_krig_panel("xi", "\u03be",
           scale_color_gradient2(low = "#2166AC", mid = "grey95", high = "#B2182B",
                                 midpoint = 0,
                                 limits = krig_xi_sym))

# Save individual panels, then compose 2-over-1 triangle via ImageMagick
# Compose 2-over-1 layout: mu0 and log_sigma on top, xi centered below
top_row <- gridExtra::arrangeGrob(p2_mu, p2_sig, ncol = 2)
bot_row <- gridExtra::arrangeGrob(grid::nullGrob(), p2_xi, grid::nullGrob(),
                                   ncol = 3, widths = c(1, 2, 1))
kp_combined <- gridExtra::arrangeGrob(top_row, bot_row, nrow = 2)
ggsave("figures/kriged_params_map.png", kp_combined, width = 20, height = 16,
       dpi = 300, bg = "white", device = ragg::agg_png)
cat("  Saved figures/kriged_params_map.png\n")


# ── Figure 4: Return levels (2-panel: RL + SE, side-by-side) ─────────────
cat("  Figure 4: return_levels_rl.png + return_levels_se.png\n")

# NOAA MLE return levels for overlay
noaa_rl_df <- data.frame(
  lon = sites$lon[noaa_idx],
  lat = sites$lat[noaa_idx]
)
noaa_rl_df$return_level <- vapply(noaa_idx, function(i) {
  gev_return_level(
    stage1$theta_hat[i, "mu"],
    exp(stage1$theta_hat[i, "log_sigma"]),
    stage1$theta_hat[i, "xi"], r = 100)
}, numeric(1))

rl_lims <- range(c(rl_joint$return_level, noaa_rl_df$return_level))

p3 <- ggplot() +
  geom_polygon(data = states_df,
               aes(x = long, y = lat, group = group),
               fill = "gray90", colour = "white", linewidth = 0.3) +
  geom_point(data = rl_joint,
             aes(x = lon, y = lat, colour = return_level,
                 size = "Prediction grid", shape = "Prediction grid"),
             alpha = 0.7) +
  geom_point(data = noaa_rl_df,
             aes(x = lon, y = lat, colour = return_level,
                 size = "NOAA tidal gauge", shape = "NOAA tidal gauge"),
             stroke = 0.6, alpha = 0.85) +
  scale_color_viridis_c(option = "viridis", direction = -1,
                        name = "RL (m)", limits = rl_lims) +
  scale_size_manual(name = NULL,
                    values = c("NOAA tidal gauge" = 3.5,
                               "Prediction grid" = 1.2),
                    breaks = c("NOAA tidal gauge", "Prediction grid")) +
  scale_shape_manual(name = NULL,
                     values = c("NOAA tidal gauge" = 16,
                                "Prediction grid" = 16),
                     breaks = c("NOAA tidal gauge", "Prediction grid")) +
  coord_cartesian(xlim = c(-98, -66), ylim = c(24, 46), expand = FALSE) +
  labs(x = "Longitude", y = "Latitude",
       title = "100-Year Return Level (Year 2000)") +
  theme_minimal(base_size = 12) +
  theme(
    plot.background   = element_rect(fill = "white", colour = NA),
    panel.background  = element_rect(fill = "white", colour = NA),
    panel.grid.major  = element_line(colour = "gray92", linewidth = 0.3),
    panel.grid.minor  = element_blank(),
    legend.key.height = unit(1.2, "cm"),
    plot.title        = element_text(face = "bold", size = 14, hjust = 0.5)
  )

ggsave("figures/return_levels_rl.png", p3, width = 10, height = 8, dpi = 300,
       bg = "white")
cat("  Saved figures/return_levels_rl.png\n")

p4 <- ggplot() +
  geom_polygon(data = states_df,
               aes(x = long, y = lat, group = group),
               fill = "gray90", colour = "white", linewidth = 0.3) +
  geom_point(data = rl_joint,
             aes(x = lon, y = lat, colour = se_sim),
             size = 1.2, shape = 16, alpha = 0.7) +
  scale_color_viridis_c(option = "viridis", direction = -1,
                        name = "SE (m)",
                        limits = range(rl_joint$se_sim)) +
  coord_cartesian(xlim = c(-98, -66), ylim = c(24, 46), expand = FALSE) +
  labs(x = "Longitude", y = "Latitude",
       title = "Standard Error") +
  theme_minimal(base_size = 12) +
  theme(
    plot.background   = element_rect(fill = "white", colour = NA),
    panel.background  = element_rect(fill = "white", colour = NA),
    panel.grid.major  = element_line(colour = "gray92", linewidth = 0.3),
    panel.grid.minor  = element_blank(),
    legend.key.height = unit(1.2, "cm"),
    plot.title        = element_text(face = "bold", size = 14, hjust = 0.5)
  )

ggsave("figures/return_levels_se.png", p4, width = 10, height = 8, dpi = 300,
       bg = "white")
cat("  Saved figures/return_levels_se.png\n")

# Combine into side-by-side
rl_combined <- gridExtra::arrangeGrob(p3, p4, ncol = 2)
ggsave("figures/return_levels.png", rl_combined, width = 20, height = 8,
       dpi = 300, bg = "white", device = ragg::agg_png)
cat("  Saved figures/return_levels.png (combined)\n")


# ── Figure 5: Return level profile along coast ────────────────────────────
cat("  Figure 5: return_level_profile.png\n")

n_grid <- nrow(grid)
gulf_grid_mask <- grid$lon < -82
gulf_grid_idx  <- which(gulf_grid_mask)
atl_grid_idx   <- which(!gulf_grid_mask)

nn_chain <- function(idx, start) {
  n_ch <- length(idx)
  visited <- logical(n_ch)
  chain <- integer(n_ch)
  chain[1] <- start
  visited[match(start, idx)] <- TRUE
  for (k in 2:n_ch) {
    prev <- chain[k - 1]
    d2 <- (grid$lon[idx] - grid$lon[prev])^2 +
           (grid$lat[idx] - grid$lat[prev])^2
    d2[visited] <- Inf
    best <- which.min(d2)
    chain[k] <- idx[best]
    visited[best] <- TRUE
  }
  ## Fix straggler points left behind by the greedy traversal.
  ## When the coast is concave (e.g. Mississippi delta), the chain can
  ## skip nearby points, then jump back to collect them at the end.
  ## Detect the first large jump (> 2 degrees ~ 220 km, well above
  ## any legitimate step at 15 km grid spacing) and reinsert stragglers
  ## at their optimal positions in the main chain.
  step_d <- sqrt(
    (grid$lon[chain[-1]] - grid$lon[chain[-n_ch]])^2 +
    (grid$lat[chain[-1]] - grid$lat[chain[-n_ch]])^2
  )
  jump_pos <- which(step_d > 2.0) + 1L
  if (length(jump_pos) > 0) {
    cut <- jump_pos[1]
    main <- chain[1:(cut - 1)]
    stragglers <- chain[cut:n_ch]
    for (s in stragglers) {
      n_m <- length(main)
      best_cost <- Inf
      best_after <- n_m
      for (e in seq_len(n_m - 1)) {
        d_es <- sqrt((grid$lon[main[e]] - grid$lon[s])^2 +
                     (grid$lat[main[e]] - grid$lat[s])^2)
        d_sn <- sqrt((grid$lon[s] - grid$lon[main[e + 1]])^2 +
                     (grid$lat[s] - grid$lat[main[e + 1]])^2)
        d_en <- sqrt((grid$lon[main[e]] - grid$lon[main[e + 1]])^2 +
                     (grid$lat[main[e]] - grid$lat[main[e + 1]])^2)
        cost <- d_es + d_sn - d_en
        if (cost < best_cost) {
          best_cost  <- cost
          best_after <- e
        }
      }
      main <- append(main, s, after = best_after)
    }
    chain <- main
  }
  chain
}

gulf_start <- gulf_grid_idx[which.min(grid$lon[gulf_grid_idx])]
gulf_chain <- nn_chain(gulf_grid_idx, gulf_start)

gulf_end <- gulf_chain[length(gulf_chain)]
d2_to_gulf_end <- (grid$lon[atl_grid_idx] - grid$lon[gulf_end])^2 +
                  (grid$lat[atl_grid_idx] - grid$lat[gulf_end])^2
atl_start <- atl_grid_idx[which.min(d2_to_gulf_end)]
atl_chain <- nn_chain(atl_grid_idx, atl_start)

order_idx <- c(gulf_chain, atl_chain)

haversine <- function(lon1, lat1, lon2, lat2) {
  R <- 6371
  dlon <- (lon2 - lon1) * pi / 180
  dlat <- (lat2 - lat1) * pi / 180
  a <- sin(dlat / 2)^2 + cos(lat1 * pi / 180) * cos(lat2 * pi / 180) * sin(dlon / 2)^2
  2 * R * asin(sqrt(a))
}

cum_dist <- numeric(n_grid)
for (k in 2:n_grid) {
  i_k <- order_idx[k]
  j_k <- order_idx[k - 1]
  cum_dist[k] <- cum_dist[k - 1] + haversine(grid$lon[j_k], grid$lat[j_k],
                                               grid$lon[i_k], grid$lat[i_k])
}
dist_km <- numeric(n_grid)
dist_km[order_idx] <- cum_dist

cities <- data.frame(
  label = c("Miami", "Tampa", "New Orleans",
            "Charleston", "NYC", "Boston"),
  lon = c(-80.19, -82.46, -90.07,
          -79.93, -74.01, -71.06),
  lat = c(25.76, 27.95, 29.95,
          32.78, 40.71, 42.36),
  stringsAsFactors = FALSE
)
city_dist <- vapply(seq_len(nrow(cities)), function(cc) {
  d2 <- (grid$lon - cities$lon[cc])^2 + (grid$lat - cities$lat[cc])^2
  dist_km[which.min(d2)]
}, numeric(1))
cities$dist_km <- city_dist

profile_df <- rbind(
  data.frame(dist = dist_km, rl = rl_joint$return_level,
             se = rl_joint$se_sim, model = "Joint"),
  data.frame(dist = dist_km, rl = rl_noaa$return_level,
             se = rl_noaa$se_sim, model = "NOAA-only"),
  data.frame(dist = dist_km, rl = rl_adcirc$return_level,
             se = rl_adcirc$se_sim, model = "ADCIRC-only")
)
profile_df$lower <- profile_df$rl - 1.96 * profile_df$se
profile_df$upper <- profile_df$rl + 1.96 * profile_df$se
profile_df <- profile_df[order(profile_df$model, profile_df$dist), ]

cities_profile <- cities[order(cities$dist_km), ]
rownames(cities_profile) <- NULL

p5 <- ggplot(profile_df, aes(x = dist, y = rl, color = model, fill = model)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = c("Joint" = "#E41A1C", "NOAA-only" = "#377EB8",
                                 "ADCIRC-only" = "#4DAF4A"), name = NULL) +
  scale_fill_manual(values = c("Joint" = "#E41A1C", "NOAA-only" = "#377EB8",
                                "ADCIRC-only" = "#4DAF4A"), name = NULL) +
  scale_x_continuous(breaks = NULL) +
  labs(x = "Distance along coast", y = "100-year return level at year 2000 (m)") +
  coord_cartesian(clip = "off",
                  xlim = c(0, max(cities$dist_km) * 1.03)) +
  theme_minimal(base_size = 11) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.position = "top",
    plot.margin = margin(5, 10, 55, 5, unit = "pt"),
    axis.title.x = element_text(margin = margin(t = 40))
  )

y_range <- range(profile_df$upper, profile_df$lower, na.rm = TRUE)
y_base  <- y_range[1] - 0.08 * diff(y_range)
y_hi    <- y_base - 0.10 * diff(y_range)
y_lo    <- y_base - 0.22 * diff(y_range)

for (ci in seq_len(nrow(cities_profile))) {
  cx <- cities_profile$dist_km[ci]
  y_label <- if (ci %% 2 == 1) y_hi else y_lo
  p5 <- p5 +
    annotate("segment", x = cx, xend = cx, y = y_base, yend = y_label + 0.02 * diff(y_range),
             colour = "grey40", linewidth = 0.4) +
    annotate("text", x = cx, y = y_label, label = cities_profile$label[ci],
             size = 10 / .pt, colour = "grey20", hjust = 0.5, vjust = 1)
}

ggsave("figures/return_level_profile.png", p5,
       width = 4000 / 300, height = 2000 / 300, dpi = 300, bg = "white")
cat("  Saved figures/return_level_profile.png\n")


# ── Figure 6: SE ratio map ────────────────────────────────────────────────
cat("  Figure 6: se_ratio_map.png\n")

ratio_df <- data.frame(
  lon = rl_joint$lon, lat = rl_joint$lat,
  ratio = rl_noaa$se_sim / rl_joint$se_sim
)
ratio_df <- ratio_df[is.finite(ratio_df$ratio), ]

p6 <- ggplot() +
  geom_polygon(data = states_df,
               aes(x = long, y = lat, group = group),
               fill = "gray90", colour = "white", linewidth = 0.3) +
  geom_point(data = ratio_df, aes(x = lon, y = lat, color = ratio),
             size = 1.2, shape = 16, alpha = 0.7) +
  scale_color_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                        midpoint = 1, limits = c(0.5, 2.0),
                        oob = scales::squish,
                        name = "SE ratio") +
  coord_cartesian(xlim = c(-98, -66), ylim = c(24, 46), expand = FALSE) +
  labs(title = "SE(NOAA-only) / SE(joint)",
       x = "Longitude", y = "Latitude") +
  theme_minimal(base_size = 12) +
  theme(
    plot.background   = element_rect(fill = "white", colour = NA),
    panel.background  = element_rect(fill = "white", colour = NA),
    panel.grid.major  = element_line(colour = "gray92", linewidth = 0.3),
    panel.grid.minor  = element_blank(),
    legend.key.height = unit(1.2, "cm"),
    plot.title        = element_text(face = "bold", size = 14, hjust = 0.5)
  )

ggsave("figures/se_ratio_map.png", p6,
       width = 10, height = 8, dpi = 300, bg = "white")
cat("  Saved figures/se_ratio_map.png\n")


# ── Figure 7: RL ratio map ────────────────────────────────────────────────
cat("  Figure 7: rl_ratio_map.png\n")

rl_ratio_df <- data.frame(
  lon = rl_joint$lon, lat = rl_joint$lat,
  ratio = rl_joint$return_level / rl_noaa$return_level
)
rl_ratio_df <- rl_ratio_df[is.finite(rl_ratio_df$ratio), ]

p7 <- ggplot() +
  geom_polygon(data = states_df,
               aes(x = long, y = lat, group = group),
               fill = "gray90", colour = "white", linewidth = 0.3) +
  geom_point(data = rl_ratio_df, aes(x = lon, y = lat, color = ratio),
             size = 1.2, shape = 16, alpha = 0.7) +
  scale_color_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                        midpoint = 1, limits = c(0.8, 1.3),
                        oob = scales::squish, name = "RL ratio") +
  coord_cartesian(xlim = c(-98, -66), ylim = c(24, 46), expand = FALSE) +
  labs(title = "RL(joint) / RL(NOAA-only)",
       x = "Longitude", y = "Latitude") +
  theme_minimal(base_size = 12) +
  theme(
    plot.background   = element_rect(fill = "white", colour = NA),
    panel.background  = element_rect(fill = "white", colour = NA),
    panel.grid.major  = element_line(colour = "gray92", linewidth = 0.3),
    panel.grid.minor  = element_blank(),
    legend.key.height = unit(1.2, "cm"),
    plot.title        = element_text(face = "bold", size = 14, hjust = 0.5)
  )

ggsave("figures/rl_ratio_map.png", p7,
       width = 10, height = 8, dpi = 300, bg = "white")
cat("  Saved figures/rl_ratio_map.png\n")


# ── Figure 8: Comparison return levels (3-panel) ──────────────────────────
cat("  Figure 8: comparison_return_levels.png\n")

comp_rl_lims <- range(c(rl_joint$return_level, rl_noaa$return_level,
                         rl_adcirc$return_level))

make_rl_panel <- function(rl_data, title) {
  ggplot() +
    geom_polygon(data = states_df,
                 aes(x = long, y = lat, group = group),
                 fill = "gray90", colour = "white", linewidth = 0.3) +
    geom_point(data = rl_data, aes(x = lon, y = lat, color = return_level),
               size = 1.2, shape = 16, alpha = 0.7) +
    scale_color_viridis_c(option = "viridis", direction = -1,
                          name = "RL (m)", limits = comp_rl_lims) +
    coord_cartesian(xlim = c(-98, -66), ylim = c(24, 46), expand = FALSE) +
    labs(title = title, x = "Longitude", y = "Latitude") +
    theme_minimal(base_size = 12) +
    theme(
      plot.background   = element_rect(fill = "white", colour = NA),
      panel.background  = element_rect(fill = "white", colour = NA),
      panel.grid.major  = element_line(colour = "gray92", linewidth = 0.3),
      panel.grid.minor  = element_blank(),
      legend.key.height = unit(1.2, "cm"),
      plot.title        = element_text(face = "bold", size = 14, hjust = 0.5)
    )
}

p8a <- make_rl_panel(rl_joint, "Joint (6-dim)")
p8b <- make_rl_panel(rl_noaa, "NOAA-only")
p8c <- make_rl_panel(rl_adcirc, "ADCIRC-only")

p8 <- grid.arrange(p8a, p8b, p8c, nrow = 1)
ggsave("figures/comparison_return_levels.png", p8,
       width = 16, height = 5.5, dpi = 300, bg = "white")
cat("  Saved figures/comparison_return_levels.png\n")


# ── Figure 9: NOAA site comparison scatter (3-panel, mu0 label) ───────────
cat("  Figure 9: noaa_site_comparison.png\n")

noaa_sites <- sites[noaa_idx, ]
preds_joint_noaa <- predict_krig(best_fit, noaa_sites)
preds_noaa_noaa  <- predict_krig_naive(noaa_model, noaa_sites)

scatter_df <- data.frame(
  mu_joint = preds_joint_noaa$noaa_mean[, "mu"],
  mu_noaa  = preds_noaa_noaa$noaa_mean[, "mu"],
  ls_joint = preds_joint_noaa$noaa_mean[, "log_sigma"],
  ls_noaa  = preds_noaa_noaa$noaa_mean[, "log_sigma"],
  xi_joint = preds_joint_noaa$noaa_mean[, "xi"],
  xi_noaa  = preds_noaa_noaa$noaa_mean[, "xi"]
)

scatter_theme <- theme_minimal(base_size = 11) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(size = 14, hjust = 0.5),
    plot.margin = margin(t = 12, r = 5, b = 5, l = 5)
  )

p9a <- ggplot(scatter_df, aes(x = mu_noaa, y = mu_joint)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 2, alpha = 0.8) +
  labs(x = "NOAA-only", y = "Joint", title = "\u03bc\u2080") +
  scatter_theme + coord_fixed()

p9b <- ggplot(scatter_df, aes(x = ls_noaa, y = ls_joint)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 2, alpha = 0.8) +
  labs(x = "NOAA-only", y = "Joint", title = "log \u03c3") +
  scatter_theme + coord_fixed()

p9c <- ggplot(scatter_df, aes(x = xi_noaa, y = xi_joint)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 2, alpha = 0.8) +
  labs(x = "NOAA-only", y = "Joint", title = "\u03be") +
  scatter_theme + coord_fixed()

p9 <- grid.arrange(p9a, p9b, p9c, nrow = 1)
ggsave("figures/noaa_site_comparison.png", p9,
       width = 14, height = 5, dpi = 300, bg = "white", device = ragg::agg_png)
cat("  Saved figures/noaa_site_comparison.png\n")


# ── Figure 10: PIT histograms (4-panel) ───────────────────────────────────
cat("  Figure 10: pit_histograms.png\n")

n_noaa <- length(noaa_idx)
params <- c("mu", "log_sigma", "xi")

z_param   <- matrix(NA_real_, n_noaa, 3, dimnames = list(NULL, params))
pit_param <- matrix(NA_real_, n_noaa, 3, dimnames = list(NULL, params))

for (k in seq_len(n_noaa)) {
  for (j in 1:3) {
    se <- sqrt(loo_joint$loo_cov[[k]][j, j])
    z_param[k, j] <- (loo_joint$observed[k, j] - loo_joint$loo_mean[k, j]) / se
    pit_param[k, j] <- pnorm(z_param[k, j])
  }
}

rl_obs_loo <- rl_loo_pred <- se_rl_loo <- numeric(n_noaa)
for (k in seq_len(n_noaa)) {
  rl_obs_loo[k] <- gev_return_level(
    loo_joint$observed[k, "mu"],
    exp(loo_joint$observed[k, "log_sigma"]),
    loo_joint$observed[k, "xi"], 100)
  rl_loo_pred[k] <- gev_return_level(
    loo_joint$loo_mean[k, "mu"],
    exp(loo_joint$loo_mean[k, "log_sigma"]),
    loo_joint$loo_mean[k, "xi"], 100)
  grad_k <- rl_gradient(
    loo_joint$loo_mean[k, "mu"],
    loo_joint$loo_mean[k, "log_sigma"],
    loo_joint$loo_mean[k, "xi"], 100)
  se_rl_loo[k] <- sqrt(max(drop(t(grad_k) %*% loo_joint$loo_cov[[k]] %*% grad_k), 0))
}
z_rl   <- (rl_obs_loo - rl_loo_pred) / se_rl_loo
pit_rl <- pnorm(z_rl)

n_bins <- 5
expected_count <- n_noaa / n_bins

pit_df <- data.frame(
  value = c(pit_param[, 1], pit_param[, 2], pit_param[, 3], pit_rl),
  parameter = factor(
    rep(c("\u03bc\u2080", "log \u03c3", "\u03be", "100-yr RL"), each = n_noaa),
    levels = c("\u03bc\u2080", "log \u03c3", "\u03be", "100-yr RL")
  )
)

p_pit <- ggplot(pit_df, aes(x = value)) +
  geom_histogram(breaks = seq(0, 1, length.out = n_bins + 1),
                 fill = "grey70", color = "black", linewidth = 0.3) +
  geom_hline(yintercept = expected_count, linetype = "dashed",
             color = "red", linewidth = 0.5) +
  facet_wrap(~parameter, nrow = 1) +
  scale_x_continuous(breaks = c(0, 0.5, 1)) +
  scale_y_continuous(breaks = function(x) seq(0, max(x), by = 2)) +
  labs(x = "PIT value", y = "Count") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(size = 14),
    axis.line = element_line(color = "black", linewidth = 0.3),
    plot.margin = margin(t = 8, r = 5, b = 5, l = 5)
  )

ggsave("figures/pit_histograms.png", p_pit,
       width = 10, height = 3.2, dpi = 300, bg = "white", device = ragg::agg_png)
cat("  Saved figures/pit_histograms.png\n")


# ══════════════════════════════════════════════════════════════════════════════
# Section 7: Block CV
# ══════════════════════════════════════════════════════════════════════════════
cat("\n== Section 7: Block CV ==\n")

blocks <- list(
  "Gulf Coast"   = c("8775870", "8761724", "8735180", "8729840"),
  "Florida"      = c("8727520", "8724580", "8725110", "8720030", "8723970"),
  "Southeast"    = c("8670870", "8665530", "8661070", "8656483", "8651370"),
  "Mid-Atlantic" = c("8635750", "8557380", "8536110", "8534720", "8518750", "8531680"),
  "New England"  = c("8467150", "8461490", "8510560", "8443970", "8447930",
                      "8418150", "8449130", "8413320", "8410140")
)

noaa_idx_all <- which(dat$sites$data_source == "NOAA")

subset_dat <- function(dat_in, keep_idx) {
  structure(
    list(
      sites    = dat_in$sites[keep_idx, ],
      maxima   = dat_in$maxima[keep_idx],
      n_noaa   = sum(dat_in$sites$data_source[keep_idx] == "NOAA"),
      n_adcirc = sum(dat_in$sites$data_source[keep_idx] == "ADCIRC"),
      n_sites  = length(keep_idx),
      source_params = dat_in$source_params,
      p = dat_in$p
    ),
    class = "evfuse_data"
  )
}

subset_stage1 <- function(s1, keep_idx) {
  structure(
    list(
      theta_hat = s1$theta_hat[keep_idx, , drop = FALSE],
      vcov_list = s1$vcov_list[keep_idx],
      fits      = s1$fits[keep_idx],
      converged = s1$converged[keep_idx],
      log_scale = s1$log_scale
    ),
    class = "evfuse_stage1"
  )
}

subset_W_bs_idx <- function(W_full, L_full, keep_idx) {
  w_idx <- c(keep_idx, keep_idx + L_full, keep_idx + 2 * L_full)
  W_full[w_idx, w_idx]
}

# Preallocate
bcv_j <- list(
  pred_mean = matrix(NA_real_, n_noaa, 3),
  pred_cov  = vector("list", n_noaa),
  observed  = stage1$theta_hat[noaa_idx_all, ],
  lpd       = numeric(n_noaa),
  block     = character(n_noaa),
  rl_pred   = numeric(n_noaa),
  rl_obs    = numeric(n_noaa),
  rl_se     = numeric(n_noaa)
)
bcv_n <- bcv_j
colnames(bcv_j$pred_mean) <- colnames(bcv_n$pred_mean) <- c("mu", "log_sigma", "xi")
colnames(bcv_j$observed)  <- colnames(bcv_n$observed)  <- c("mu", "log_sigma", "xi")

for (k in seq_len(n_noaa)) {
  bcv_j$rl_obs[k] <- gev_return_level(
    bcv_j$observed[k, 1], exp(bcv_j$observed[k, 2]), bcv_j$observed[k, 3], 100)
  bcv_n$rl_obs[k] <- bcv_j$rl_obs[k]
}

for (b in seq_along(blocks)) {
  block_name <- names(blocks)[b]
  held_out_locs <- blocks[[b]]
  held_out_idx  <- match(held_out_locs, dat$sites$location)
  retained_idx  <- setdiff(seq_len(L), held_out_idx)
  noaa_pos      <- match(held_out_locs, dat$sites$location[noaa_idx_all])

  cat(sprintf("  Fold %d: %s (%d held out)\n", b, block_name, length(held_out_idx)))

  dat_sub    <- subset_dat(dat, retained_idx)
  stage1_sub <- subset_stage1(stage1, retained_idx)
  D_sub      <- D[retained_idx, retained_idx]
  W_bs_sub   <- subset_W_bs_idx(W_bs, L, retained_idx)
  W_tap_sub  <- taper_W(W_bs_sub, D_sub, lambda = 300)
  held_out_sites <- dat$sites[held_out_idx, ]

  # Joint model
  cat("    Joint model...")
  start0 <- list(beta = best_fit$beta, A = best_fit$A, rho = best_fit$rho)
  bcv_fit <- tryCatch(
    suppressWarnings(
      fit_spatial_model(stage1_sub, dat_sub, W_tap_sub, D_sub,
                        start = start0,
                        control = list(maxit = 2000, trace = 0))
    ),
    error = function(e) NULL
  )
  bcv_nll <- if (!is.null(bcv_fit)) bcv_fit$optim_result$value else Inf

  set.seed(2026)
  for (s in seq_len(20)) {
    rho_s <- exp(runif(6, log(50), log(5000)))
    start_s <- list(beta = best_fit$beta, A = best_fit$A, rho = rho_s)
    fit_s <- tryCatch(
      suppressWarnings(
        fit_spatial_model(stage1_sub, dat_sub, W_tap_sub, D_sub,
                          start = start_s,
                          control = list(maxit = 2000, trace = 0))
      ),
      error = function(e) NULL
    )
    if (!is.null(fit_s) && fit_s$optim_result$value < bcv_nll) {
      bcv_nll <- fit_s$optim_result$value
      bcv_fit <- fit_s
    }
  }
  cat(sprintf(" NLL = %.4f\n", bcv_nll))

  preds_bcv_j <- predict_krig(bcv_fit, held_out_sites)

  # NOAA-only model
  cat("    NOAA-only model...")
  noaa_fit <- tryCatch(
    suppressWarnings(
      fit_naive_model(stage1_sub, dat_sub, W_bs_sub, D_sub,
                      source = "NOAA", lambda = 300, n_starts = 20,
                      control = list(maxit = 2000, trace = 0))
    ),
    error = function(e) { cat(" FAILED\n"); NULL }
  )
  if (!is.null(noaa_fit)) {
    cat(sprintf(" NLL = %.4f\n", noaa_fit$optim_result$value))
  }
  preds_bcv_n <- if (!is.null(noaa_fit)) predict_krig_naive(noaa_fit, held_out_sites) else NULL

  # Store results
  for (k in seq_along(held_out_idx)) {
    row <- noaa_pos[k]
    obs_k <- stage1$theta_hat[held_out_idx[k], ]

    mu_j  <- preds_bcv_j$noaa_mean[k, ]
    cov_j <- preds_bcv_j$noaa_cov[[k]]
    bcv_j$pred_mean[row, ] <- mu_j
    bcv_j$pred_cov[[row]]  <- cov_j
    bcv_j$block[row]       <- block_name

    resid_k <- obs_k - mu_j
    cov_reg <- cov_j + diag(1e-8, 3)
    R_k <- chol(cov_reg)
    bcv_j$lpd[row] <- -0.5 * (3 * log(2 * pi) + 2 * sum(log(diag(R_k))) +
      sum(resid_k * backsolve(R_k, backsolve(R_k, resid_k, transpose = TRUE))))

    bcv_j$rl_pred[row] <- gev_return_level(mu_j[1], exp(mu_j[2]), mu_j[3], 100)
    grad_j <- rl_gradient(mu_j[1], mu_j[2], mu_j[3], 100)
    bcv_j$rl_se[row] <- sqrt(max(drop(t(grad_j) %*% cov_j %*% grad_j), 0))

    if (!is.null(preds_bcv_n)) {
      mu_n  <- preds_bcv_n$noaa_mean[k, ]
      cov_n <- preds_bcv_n$noaa_cov[[k]]
      bcv_n$pred_mean[row, ] <- mu_n
      bcv_n$pred_cov[[row]]  <- cov_n
      bcv_n$block[row]       <- block_name

      resid_n <- obs_k - mu_n
      cov_reg_n <- cov_n + diag(1e-8, 3)
      R_n <- chol(cov_reg_n)
      bcv_n$lpd[row] <- -0.5 * (3 * log(2 * pi) + 2 * sum(log(diag(R_n))) +
        sum(resid_n * backsolve(R_n, backsolve(R_n, resid_n, transpose = TRUE))))

      bcv_n$rl_pred[row] <- gev_return_level(mu_n[1], exp(mu_n[2]), mu_n[3], 100)
      grad_n <- rl_gradient(mu_n[1], mu_n[2], mu_n[3], 100)
      bcv_n$rl_se[row] <- sqrt(max(drop(t(grad_n) %*% cov_n %*% grad_n), 0))
    }
  }
}

# Block CV aggregate metrics
bcv_rmse_j <- bcv_rmse_n <- numeric(3)
for (j in 1:3) {
  bcv_rmse_j[j] <- sqrt(mean((bcv_j$pred_mean[, j] - bcv_j$observed[, j])^2))
  bcv_rmse_n[j] <- sqrt(mean((bcv_n$pred_mean[, j] - bcv_n$observed[, j])^2))
}
bcv_rl_rmse_j <- sqrt(mean((bcv_j$rl_pred - bcv_j$rl_obs)^2))
bcv_rl_rmse_n <- sqrt(mean((bcv_n$rl_pred - bcv_n$rl_obs)^2))
bcv_lpd_j     <- sum(bcv_j$lpd)
bcv_lpd_n     <- sum(bcv_n$lpd)
bcv_wins      <- sum(bcv_j$lpd > bcv_n$lpd)
bcv_rl_pct    <- 100 * (1 - bcv_rl_rmse_j / bcv_rl_rmse_n)

# Excl Gulf
gulf_bcv_idx <- which(bcv_j$block == "Gulf Coast")
ex_idx <- setdiff(seq_len(n_noaa), gulf_bcv_idx)
bcv_rl_rmse_j_ex <- sqrt(mean((bcv_j$rl_pred[ex_idx] - bcv_j$rl_obs[ex_idx])^2))
bcv_rl_rmse_n_ex <- sqrt(mean((bcv_n$rl_pred[ex_idx] - bcv_n$rl_obs[ex_idx])^2))
bcv_lpd_j_ex     <- sum(bcv_j$lpd[ex_idx])
bcv_lpd_n_ex     <- sum(bcv_n$lpd[ex_idx])
bcv_wins_ex      <- sum(bcv_j$lpd[ex_idx] > bcv_n$lpd[ex_idx])
bcv_rl_pct_ex    <- 100 * (1 - bcv_rl_rmse_j_ex / bcv_rl_rmse_n_ex)

cat(sprintf("\nBlock CV: Joint RL RMSE = %.4f, NOAA-only RL RMSE = %.4f\n",
            bcv_rl_rmse_j, bcv_rl_rmse_n))
cat(sprintf("Block CV: RL reduction = %.1f%% (all 29), %.1f%% (excl Gulf)\n",
            bcv_rl_pct, bcv_rl_pct_ex))
cat(sprintf("Block CV: Joint wins %d/29, LPD gain (excl Gulf) = %.2f\n",
            bcv_wins, bcv_lpd_j_ex - bcv_lpd_n_ex))

saveRDS(list(joint = bcv_j, noaa = bcv_n, blocks = blocks),
        "data-raw/block_cv_ns.rds")
cat("Saved data-raw/block_cv_ns.rds\n")


# ══════════════════════════════════════════════════════════════════════════════
# Section 8: PIT Calibration Table
# ══════════════════════════════════════════════════════════════════════════════
cat("\n== Section 8: PIT Calibration ==\n")

nominal <- c(0.50, 0.80, 0.95)
cov_labels <- c("mu", "log_sigma", "xi", "100-yr RL")
coverage <- matrix(NA_real_, 4, 3,
                   dimnames = list(cov_labels, paste0(100 * nominal, "%")))

for (j in 1:3) {
  for (a in seq_along(nominal)) {
    z_crit <- qnorm(1 - (1 - nominal[a]) / 2)
    coverage[j, a] <- mean(abs(z_param[, j]) < z_crit)
  }
}
for (a in seq_along(nominal)) {
  z_crit <- qnorm(1 - (1 - nominal[a]) / 2)
  coverage[4, a] <- mean(abs(z_rl) < z_crit)
}

cat("\nCoverage Table (nominal vs. empirical, n = 29)\n")
cat(sprintf("%-12s %8s %8s %8s\n", "Parameter", "50%", "80%", "95%"))
cat(sprintf("%-12s %8s %8s %8s\n", "---------", "---", "---", "---"))
for (i_cov in 1:4) {
  counts <- round(coverage[i_cov, ] * n_noaa)
  cat(sprintf("%-12s %4d/29 %4d/29 %4d/29\n",
      cov_labels[i_cov], counts[1], counts[2], counts[3]))
}
cat(sprintf("%-12s %8s %8s %8s\n", "(expected)", "14.5", "23.2", "27.6"))

cat("\nKS Tests for Uniformity of PIT Values\n")
cat(sprintf("%-12s %8s %10s\n", "Parameter", "D", "p-value"))
cat(sprintf("%-12s %8s %10s\n", "---------", "---", "-------"))
for (j in 1:3) {
  ks <- ks.test(pit_param[, j], "punif")
  cat(sprintf("%-12s %8.3f %10.3f\n", params[j], ks$statistic, ks$p.value))
}
ks_rl <- ks.test(pit_rl, "punif")
cat(sprintf("%-12s %8.3f %10.3f\n", "100-yr RL", ks_rl$statistic, ks_rl$p.value))
ks_all <- ks.test(c(pit_param), "punif")
cat(sprintf("%-12s %8.3f %10.3f\n", "Combined", ks_all$statistic, ks_all$p.value))


# ══════════════════════════════════════════════════════════════════════════════
# Section 9: Taper Sensitivity
# ══════════════════════════════════════════════════════════════════════════════
cat("\n== Section 9: Taper Sensitivity ==\n")

extract_cors <- function(model) {
  AAT <- model$A %*% t(model$A)
  sds <- sqrt(diag(AAT))
  cor_mat <- AAT / outer(sds, sds)
  c(mu = cor_mat[1, 4], log_sigma = cor_mat[2, 5], xi = cor_mat[3, 6])
}

taper_lambdas <- c(150, 300, 500)
taper_results <- data.frame(
  lambda = taper_lambdas,
  nll = NA_real_, rl_rmse_j = NA_real_, rl_rmse_n = NA_real_,
  rl_pct = NA_real_, lpd_gain = NA_real_,
  cor_mu = NA_real_, cor_ls = NA_real_, cor_xi = NA_real_
)

for (tl in seq_along(taper_lambdas)) {
  lam <- taper_lambdas[tl]
  cat(sprintf("  lambda = %d km\n", lam))

  if (lam == 300) {
    # Primary model: reuse existing
    taper_results$nll[tl]      <- best_nll
    taper_results$rl_rmse_j[tl] <- sum_joint$rl_rmse
    taper_results$rl_rmse_n[tl] <- sum_noaa$rl_rmse
    taper_results$rl_pct[tl]   <- rl_pct
    taper_results$lpd_gain[tl] <- lpd_gain
    cors_300 <- extract_cors(best_fit)
    taper_results$cor_mu[tl] <- cors_300["mu"]
    taper_results$cor_ls[tl] <- cors_300["log_sigma"]
    taper_results$cor_xi[tl] <- cors_300["xi"]
    next
  }

  W_tap_l <- taper_W(W_bs, D, lambda = lam)

  # Fit joint model with this lambda
  fit_l <- tryCatch(
    suppressWarnings(
      fit_spatial_model(stage1, dat, W_tap_l, D,
                        start = list(beta = best_fit$beta, A = best_fit$A,
                                     rho = best_fit$rho),
                        control = list(maxit = 2000, trace = 0))
    ),
    error = function(e) NULL
  )
  nll_l <- if (!is.null(fit_l)) fit_l$optim_result$value else Inf

  set.seed(2026)
  for (s in seq_len(20)) {
    rho_s <- exp(runif(6, log(50), log(5000)))
    start_s <- list(beta = best_fit$beta, A = best_fit$A, rho = rho_s)
    fit_s <- tryCatch(
      suppressWarnings(
        fit_spatial_model(stage1, dat, W_tap_l, D,
                          start = start_s,
                          control = list(maxit = 2000, trace = 0))
      ),
      error = function(e) NULL
    )
    if (!is.null(fit_s) && fit_s$optim_result$value < nll_l) {
      nll_l <- fit_s$optim_result$value
      fit_l <- fit_s
    }
  }

  taper_results$nll[tl] <- nll_l

  # LOO-CV
  loo_l <- loo_cv(fit_l)
  sum_l <- loo_summary(loo_l, r = 100)
  taper_results$rl_rmse_j[tl] <- sum_l$rl_rmse
  taper_results$rl_rmse_n[tl] <- sum_noaa$rl_rmse  # same NOAA-only baseline
  taper_results$rl_pct[tl]    <- 100 * (1 - sum_l$rl_rmse / sum_noaa$rl_rmse)
  taper_results$lpd_gain[tl]  <- sum_l$total_lpd - sum_noaa$total_lpd

  cors_l <- extract_cors(fit_l)
  taper_results$cor_mu[tl] <- cors_l["mu"]
  taper_results$cor_ls[tl] <- cors_l["log_sigma"]
  taper_results$cor_xi[tl] <- cors_l["xi"]
}

cat("\nTaper Sensitivity Table\n")
cat(sprintf("%-8s %10s %10s %12s %10s %8s %8s %8s\n",
    "lambda", "NLL", "RL RMSE", "RL Redn %", "LPD gain", "cor_mu", "cor_ls", "cor_xi"))
for (tl in seq_len(nrow(taper_results))) {
  cat(sprintf("%-8d %10.4f %10.4f %11.1f%% %10.2f %8.3f %8.3f %8.3f\n",
      taper_results$lambda[tl], taper_results$nll[tl],
      taper_results$rl_rmse_j[tl], taper_results$rl_pct[tl],
      taper_results$lpd_gain[tl],
      taper_results$cor_mu[tl], taper_results$cor_ls[tl], taper_results$cor_xi[tl]))
}


# ══════════════════════════════════════════════════════════════════════════════
# Section 10: xi Bootstrap QQ Plots
# ══════════════════════════════════════════════════════════════════════════════
cat("\n== Section 10: xi Bootstrap QQ plots ==\n")

# Select 4 NOAA sites by xi value
noaa_xi_vals <- stage1$theta_hat[noaa_idx, "xi"]
noaa_xi_order <- order(noaa_xi_vals)

# Pick: lowest xi, 1st quartile, median, highest xi
qi <- c(1, ceiling(n_noaa / 4), ceiling(n_noaa / 2), n_noaa)
xi_qq_sites <- noaa_idx[noaa_xi_order[qi]]

png("figures/xi_bootstrap_qq.png", width = 2400, height = 2400, res = 300,
    bg = "white")
par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3.5, 1))

for (k in seq_along(xi_qq_sites)) {
  i <- xi_qq_sites[k]

  # Extract xi bootstrap replicates from Gamma: param 3, site i -> column (2*L + i)
  xi_boot <- bs$Gamma[, 2 * L + i]
  xi_boot <- xi_boot[!is.na(xi_boot)]

  xi_hat <- stage1$theta_hat[i, "xi"]

  qqnorm(xi_boot, pch = 16, cex = 0.7, col = "steelblue",
         main = sprintf("%s\n(xi = %.3f, n_boot = %d)",
                        sites$location[i], xi_hat, length(xi_boot)),
         cex.main = 0.9, cex.lab = 1.0, cex.axis = 0.9)
  qqline(xi_boot, col = "red", lwd = 1.5)
}
dev.off()
cat("Saved figures/xi_bootstrap_qq.png\n")


# ══════════════════════════════════════════════════════════════════════════════
# Section 11: Manuscript Numbers Summary
# ══════════════════════════════════════════════════════════════════════════════
cat("\n== Section 11: Manuscript Numbers Summary ==\n")

cors <- extract_cors(best_fit)

# Losing LOO sites
lpd_diff <- loo_joint$loo_lpd - loo_noaa$loo_lpd
losing_idx <- which(lpd_diff < 0)

# mu RMSE reduction
mu_rmse_pct <- 100 * (1 - sum_joint$param_stats$rmse[1] / sum_noaa$param_stats$rmse[1])

# Mean SE reduction
mean_se_j <- mean(rl_joint$se_sim, na.rm = TRUE)
mean_se_n <- mean(rl_noaa$se_sim, na.rm = TRUE)
se_pct <- 100 * (1 - mean_se_j / mean_se_n)

# NN distances
nn_dist <- apply(D[noaa_idx, adcirc_idx], 1, min)

summary_text <- function() {
  cat("\n")
  cat("════════════════════════════════════════════════════════════════════\n")
  cat("NONSTATIONARY PIPELINE: MANUSCRIPT NUMBERS SUMMARY\n")
  cat("════════════════════════════════════════════════════════════════════\n\n")

  cat("── Model Fit ──────────────────────────────────────────────────────\n")
  cat(sprintf("Joint 6-dim NLL:    %.4f (33 params, 129 sites)\n", best_nll))
  cat(sprintf("NOAA-only NLL:      %.4f (12 params, 29 sites)\n",
              noaa_model$optim_result$value))
  cat(sprintf("ADCIRC-only NLL:    %.4f (12 params, 100 sites)\n",
              adcirc_model$optim_result$value))
  cat(sprintf("Joint rho:          %s\n", paste(round(best_fit$rho, 1), collapse = ", ")))
  cat(sprintf("Joint beta:         %s\n",
              paste(round(best_fit$beta, 4), collapse = ", ")))
  cat("\n")

  cat("── Cross-Source Correlations ──────────────────────────────────────\n")
  cat(sprintf("Cor(mu0_N, mu0_A):     %.3f\n", cors["mu"]))
  cat(sprintf("Cor(log_sig_N, log_sig_A): %.3f\n", cors["log_sigma"]))
  cat(sprintf("Cor(xi_N, xi_A):       %.3f\n", cors["xi"]))
  cat("\n")

  cat("── LOO-CV ────────────────────────────────────────────────────────\n")
  cat(sprintf("Joint:     mu RMSE = %.4f, RL RMSE = %.4f, LPD = %.2f\n",
              sum_joint$param_stats$rmse[1], sum_joint$rl_rmse, sum_joint$total_lpd))
  cat(sprintf("NOAA-only: mu RMSE = %.4f, RL RMSE = %.4f, LPD = %.2f\n",
              sum_noaa$param_stats$rmse[1], sum_noaa$rl_rmse, sum_noaa$total_lpd))
  cat(sprintf("Joint wins: %d / 29\n", n_wins))
  cat(sprintf("LPD gain:   %.2f\n", lpd_gain))
  cat(sprintf("RL RMSE reduction: %.1f%%\n", rl_pct))
  cat(sprintf("mu RMSE reduction: %.1f%%\n", mu_rmse_pct))
  cat(sprintf("Mean SE reduction: %.1f%%\n", se_pct))

  if (length(losing_idx) > 0) {
    cat("\nLosing LOO sites:\n")
    for (idx in losing_idx) {
      cat(sprintf("  %s: diff = %+.4f\n",
                  loo_joint$sites$location[idx], lpd_diff[idx]))
    }
  }
  cat("\n")

  cat("── Block CV ──────────────────────────────────────────────────────\n")
  cat(sprintf("Joint RL RMSE:     %.4f\n", bcv_rl_rmse_j))
  cat(sprintf("NOAA-only RL RMSE: %.4f\n", bcv_rl_rmse_n))
  cat(sprintf("RL RMSE reduction: %.1f%% (all 29), %.1f%% (excl Gulf)\n",
              bcv_rl_pct, bcv_rl_pct_ex))
  cat(sprintf("Joint wins: %d/29 (all), %d/25 (excl Gulf)\n",
              bcv_wins, bcv_wins_ex))
  cat(sprintf("LPD gain (excl Gulf): %.2f\n", bcv_lpd_j_ex - bcv_lpd_n_ex))
  cat("\n")

  cat("── Taper Sensitivity ─────────────────────────────────────────────\n")
  for (tl in seq_len(nrow(taper_results))) {
    cat(sprintf("  lambda=%d: NLL=%.4f, RL RMSE=%.4f, reduction=%.1f%%, LPD gain=%.2f\n",
        taper_results$lambda[tl], taper_results$nll[tl],
        taper_results$rl_rmse_j[tl], taper_results$rl_pct[tl],
        taper_results$lpd_gain[tl]))
  }
  cat("\n")

  cat("── PIT Calibration ───────────────────────────────────────────────\n")
  for (i_cov in 1:4) {
    counts <- round(coverage[i_cov, ] * n_noaa)
    cat(sprintf("  %-12s: 50%%=%d/29, 80%%=%d/29, 95%%=%d/29\n",
        cov_labels[i_cov], counts[1], counts[2], counts[3]))
  }
  cat("\n")

  cat("── Nearest-Neighbor Distances ────────────────────────────────────\n")
  cat(sprintf("Median NOAA->ADCIRC: %.1f km\n", median(nn_dist)))
  cat(sprintf("Range: %.1f - %.1f km\n", min(nn_dist), max(nn_dist)))
  cat("\n")

  cat("── Return Levels ─────────────────────────────────────────────────\n")
  cat(sprintf("Joint RL range:     [%.3f, %.3f] m\n",
              min(rl_joint$return_level), max(rl_joint$return_level)))
  cat(sprintf("Joint mean SE:      %.4f m\n", mean_se_j))
  cat(sprintf("NOAA-only mean SE:  %.4f m\n", mean_se_n))
  cat("\n")
}

summary_text()

# Save to file
sink("tables/nonstationary_summary.txt")
summary_text()
sink()
cat("Saved tables/nonstationary_summary.txt\n")


cat("\n══════════════════════════════════════════════════════════════════════════════\n")
cat("DONE. All nonstationary pipeline outputs generated.\n")
cat("══════════════════════════════════════════════════════════════════════════════\n")
