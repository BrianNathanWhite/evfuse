#!/usr/bin/env Rscript
# Generate publication-quality result plots for evfuse analysis

library(ggplot2)
library(sf)
library(gridExtra)
devtools::load_all()

dir.create("figures", showWarnings = FALSE)

# ── Load data ────────────────────────────────────────────────
model <- readRDS("data-raw/model_6dim_best.rds")
preds <- readRDS("data-raw/predictions_grid.rds")
rl    <- readRDS("data-raw/return_levels_100yr.rds")

sites   <- model$dat$sites
stage1  <- model$stage1
grid    <- read.csv("data-raw/prediction_grid.csv")

# ── Basemap: US state boundaries as sf ───────────────────────
states_df <- map_data("state")
states_sf <- lapply(split(states_df, states_df$group), function(grp) {
  st_polygon(list(as.matrix(grp[, c("long", "lat")])))
})
states_sf <- st_sfc(states_sf, crs = 4326)
states_sf <- st_sf(geometry = states_sf)

crs_albers <- st_crs(5070)

# ── Common theme ─────────────────────────────────────────────
theme_map <- theme_minimal(base_size = 11) +
  theme(
    axis.title = element_blank(),
    axis.text = element_text(size = 8, color = "grey40"),
    panel.grid = element_blank(),
    strip.text = element_text(size = 12),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

map_coord <- function() {
  coord_sf(crs = crs_albers,
           xlim = c(-98, -66), ylim = c(24, 46),
           default_crs = st_crs(4326),
           expand = FALSE)
}

basemap <- geom_sf(data = states_sf, fill = "grey95", color = "grey60",
                    linewidth = 0.3, inherit.aes = FALSE)

# ── 1. Stage 1 MLE map ──────────────────────────────────────
s1_df <- data.frame(
  lon = sites$lon,
  lat = sites$lat,
  source = sites$data_source,
  mu = stage1$theta_hat[, "mu"],
  log_sigma = stage1$theta_hat[, "log_sigma"],
  xi = stage1$theta_hat[, "xi"]
)

# Split into NOAA and ADCIRC sf objects
s1_noaa <- st_as_sf(s1_df[s1_df$source == "NOAA", ],
                     coords = c("lon", "lat"), crs = 4326)
s1_adcirc <- st_as_sf(s1_df[s1_df$source == "ADCIRC", ],
                       coords = c("lon", "lat"), crs = 4326)

make_s1_panel <- function(param_col, label, fill_scale, color_scale,
                          show_source_legend = FALSE) {
  p <- ggplot() +
    basemap +
    # ADCIRC: small filled triangles, thin matching outline
    geom_sf(data = s1_adcirc, aes(fill = .data[[param_col]],
            color = .data[[param_col]]),
            shape = 24, size = 1.8, stroke = 0.3, alpha = 0.85) +
    # NOAA: larger, black outline, foreground
    geom_sf(data = s1_noaa, aes(fill = .data[[param_col]]),
            shape = 21, size = 3, stroke = 0.8, color = "black", alpha = 0.9) +
    fill_scale +
    color_scale +
    map_coord() +
    labs(title = label, fill = NULL) +
    guides(colour = "none") +
    theme_map +
    theme(legend.key.height = unit(0.8, "cm"),
          plot.title = element_text(size = 12, hjust = 0.5))

  if (show_source_legend) {
    # Invisible geom_point (non-sf) to generate shape legend without sf_grob issues
    legend_df <- data.frame(
      x = c(-80, -80), y = c(35, 35),
      Source = factor(c("NOAA", "ADCIRC"), levels = c("NOAA", "ADCIRC"))
    )
    p <- p +
      geom_point(data = legend_df, aes(x = x, y = y, shape = Source),
                 size = NA_real_, alpha = 0) +
      scale_shape_manual(values = c("NOAA" = 21, "ADCIRC" = 24),
                         name = NULL) +
      guides(shape = guide_legend(override.aes = list(
        size = c(3, 1.8), alpha = 1,
        fill = "grey50", color = c("black", "grey50"),
        stroke = c(0.8, 0.3)
      )))
  }
  p
}

mu_lims <- range(s1_df$mu)
ls_lims <- range(s1_df$log_sigma)
xi_lims <- range(s1_df$xi)

p1_mu  <- make_s1_panel("mu", expression(mu),
           scale_fill_viridis_c(option = "inferno", begin = 0.1, end = 0.95,
                                limits = mu_lims),
           scale_color_viridis_c(option = "inferno", begin = 0.1, end = 0.95,
                                 limits = mu_lims))
p1_sig <- make_s1_panel("log_sigma", expression(log~sigma),
           scale_fill_distiller(palette = "YlGnBu", direction = 1,
                                limits = ls_lims),
           scale_color_distiller(palette = "YlGnBu", direction = 1,
                                 limits = ls_lims))
p1_xi  <- make_s1_panel("xi", expression(xi),
           scale_fill_gradient2(low = "#2166AC", mid = "grey95", high = "#B2182B",
                                midpoint = mean(xi_lims), limits = xi_lims),
           scale_color_gradient2(low = "#2166AC", mid = "grey95", high = "#B2182B",
                                 midpoint = mean(xi_lims), limits = xi_lims),
           show_source_legend = TRUE)

p1 <- grid.arrange(p1_mu, p1_sig, p1_xi, nrow = 1)
ggsave("figures/stage1_mle_map.png", p1, width = 14, height = 6, dpi = 300,
       bg = "white")
cat("Saved figures/stage1_mle_map.png\n")

# ── 2. Kriged parameter map ─────────────────────────────────
make_krig_panel <- function(param_col, label, scale) {
  df <- data.frame(lon = grid$lon, lat = grid$lat,
                    value = preds$noaa_mean[, param_col])
  pts <- st_as_sf(df, coords = c("lon", "lat"), crs = 4326)

  ggplot() +
    basemap +
    geom_sf(data = pts, aes(color = value), size = 0.9, shape = 16, alpha = 0.7) +
    scale +
    map_coord() +
    labs(title = label, color = NULL) +
    theme_map +
    theme(legend.key.height = unit(0.8, "cm"),
          plot.title = element_text(size = 12, hjust = 0.5))
}

krig_mu_lims <- range(preds$noaa_mean[, "mu"])
krig_ls_lims <- range(preds$noaa_mean[, "log_sigma"])
krig_xi_lims <- range(preds$noaa_mean[, "xi"])

p2_mu  <- make_krig_panel("mu", expression(mu),
           scale_color_viridis_c(option = "inferno", begin = 0.1, end = 0.95,
                                 limits = krig_mu_lims))
p2_sig <- make_krig_panel("log_sigma", expression(log~sigma),
           scale_color_distiller(palette = "YlGnBu", direction = 1,
                                 limits = krig_ls_lims))
p2_xi  <- make_krig_panel("xi", expression(xi),
           scale_color_gradient2(low = "#2166AC", mid = "grey95", high = "#B2182B",
                                 midpoint = mean(krig_xi_lims),
                                 limits = krig_xi_lims))

p2 <- grid.arrange(p2_mu, p2_sig, p2_xi, nrow = 1)
ggsave("figures/kriged_params_map.png", p2, width = 14, height = 6, dpi = 300,
       bg = "white")
cat("Saved figures/kriged_params_map.png\n")

# ── 3. 100-year return level map ────────────────────────────
noaa_idx <- which(sites$data_source == "NOAA")
noaa_rl_df <- data.frame(
  lon = sites$lon[noaa_idx],
  lat = sites$lat[noaa_idx]
)
noaa_rl_df$return_level <- vapply(noaa_idx, function(i) {
  gev_return_level(
    stage1$theta_hat[i, "mu"],
    exp(stage1$theta_hat[i, "log_sigma"]),
    stage1$theta_hat[i, "xi"],
    r = 100
  )
}, numeric(1))
noaa_pts <- st_as_sf(noaa_rl_df, coords = c("lon", "lat"), crs = 4326)
rl_pts <- st_as_sf(rl, coords = c("lon", "lat"), crs = 4326)

rl_lims <- range(c(rl$return_level, noaa_rl_df$return_level))

p3 <- ggplot() +
  basemap +
  geom_sf(data = rl_pts, aes(color = return_level),
          size = 0.9, shape = 16, alpha = 0.7) +
  geom_sf(data = noaa_pts, aes(fill = return_level),
          shape = 21, size = 2.5, stroke = 0.5, color = "black", alpha = 0.85) +
  scale_color_viridis_c(option = "inferno", begin = 0.1, end = 0.95,
                        aesthetics = c("colour", "fill"),
                        name = "RL (m)", limits = rl_lims) +
  map_coord() +
  theme_map +
  theme(legend.key.height = unit(1.2, "cm"))

# ── 4. Return level uncertainty map ─────────────────────────
p4 <- ggplot() +
  basemap +
  geom_sf(data = rl_pts, aes(color = se_delta),
          size = 0.9, shape = 16, alpha = 0.7) +
  scale_color_viridis_c(option = "inferno", name = "SE (m)",
                        limits = range(rl$se_delta)) +
  map_coord() +
  theme_map +
  theme(legend.key.height = unit(1.2, "cm"))

p3_4 <- grid.arrange(p3, p4, nrow = 1)
ggsave("figures/return_levels.png", p3_4, width = 14, height = 6, dpi = 300,
       bg = "white")
cat("Saved figures/return_levels.png\n")

cat("\nAll plots saved to figures/\n")
