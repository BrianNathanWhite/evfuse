#!/usr/bin/env Rscript
# Comparison plots and tables: Joint vs NOAA-only vs ADCIRC-only
library(ggplot2)
library(sf)
library(gridExtra)
devtools::load_all()

dir.create("figures", showWarnings = FALSE)

# ── Load all model outputs ────────────────────────────────────
joint_model  <- readRDS("data-raw/model_6dim_best.rds")
noaa_model   <- readRDS("data-raw/model_noaa_only.rds")
adcirc_model <- readRDS("data-raw/model_adcirc_only.rds")

preds_joint  <- readRDS("data-raw/predictions_grid.rds")
preds_noaa   <- readRDS("data-raw/predictions_noaa_only.rds")
preds_adcirc <- readRDS("data-raw/predictions_adcirc_only.rds")

rl_joint  <- readRDS("data-raw/return_levels_100yr.rds")
rl_noaa   <- readRDS("data-raw/rl_noaa_only.rds")
rl_adcirc <- readRDS("data-raw/rl_adcirc_only.rds")

grid   <- read.csv("data-raw/prediction_grid.csv")
sites  <- joint_model$dat$sites
stage1 <- joint_model$stage1

# ── Basemap ───────────────────────────────────────────────────
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
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(size = 14, hjust = 0.5)
  )

map_coord <- function() {
  coord_sf(crs = crs_albers,
           xlim = c(-98, -66), ylim = c(24, 46),
           default_crs = st_crs(4326),
           expand = FALSE)
}

basemap <- geom_sf(data = states_sf, fill = "grey95", color = "grey60",
                    linewidth = 0.3, inherit.aes = FALSE)

# ══════════════════════════════════════════════════════════════
# TABLE 1: Model fit summary
# ══════════════════════════════════════════════════════════════
cat("\n")
cat("════════════════════════════════════════════════════════════\n")
cat("TABLE 1: Model Fit Summary\n")
cat("════════════════════════════════════════════════════════════\n")
cat(sprintf("%-15s %6s %6s %8s %10s\n", "Model", "p", "Sites", "Params", "NLL"))
cat(sprintf("%-15s %6s %6s %8s %10s\n", "-----", "--", "-----", "------", "---"))
cat(sprintf("%-15s %6d %6d %8d %10.4f\n", "Joint (6-dim)",
            6, 129, 33, joint_model$optim_result$value))
cat(sprintf("%-15s %6d %6d %8d %10.4f\n", "NOAA-only",
            3, 29, 12, noaa_model$optim_result$value))
cat(sprintf("%-15s %6d %6d %8d %10.4f\n", "ADCIRC-only",
            3, 100, 12, adcirc_model$optim_result$value))
cat("\n")

# ══════════════════════════════════════════════════════════════
# TABLE 2: NOAA site prediction accuracy (joint vs NOAA-only)
# ══════════════════════════════════════════════════════════════
# Krige both models at the 29 NOAA site locations
noaa_idx <- which(sites$data_source == "NOAA")
noaa_sites <- sites[noaa_idx, ]

# Joint model predictions at NOAA sites
preds_joint_noaa <- predict_krig(joint_model, noaa_sites)
# NOAA-only model predictions at NOAA sites
preds_noaa_noaa <- predict_krig_naive(noaa_model, noaa_sites)

# Stage 1 MLEs (truth for comparison)
theta_noaa <- stage1$theta_hat[noaa_idx, ]

cat("════════════════════════════════════════════════════════════\n")
cat("TABLE 2: Prediction Accuracy at 29 NOAA Sites\n")
cat("         (Kriging prediction vs Stage 1 MLE)\n")
cat("════════════════════════════════════════════════════════════\n")
cat(sprintf("%-12s %-10s %10s %10s\n", "Parameter", "Model", "RMSE", "MAD"))
cat(sprintf("%-12s %-10s %10s %10s\n", "---------", "-----", "----", "---"))

for (param in c("mu", "log_sigma", "xi")) {
  j <- match(param, c("mu", "log_sigma", "xi"))
  obs <- theta_noaa[, j]

  pred_j <- preds_joint_noaa$noaa_mean[, param]
  pred_n <- preds_noaa_noaa$noaa_mean[, param]

  rmse_j <- sqrt(mean((pred_j - obs)^2))
  rmse_n <- sqrt(mean((pred_n - obs)^2))
  mad_j <- mean(abs(pred_j - obs))
  mad_n <- mean(abs(pred_n - obs))

  cat(sprintf("%-12s %-10s %10.4f %10.4f\n", param, "Joint", rmse_j, mad_j))
  cat(sprintf("%-12s %-10s %10.4f %10.4f\n", "", "NOAA-only", rmse_n, mad_n))
}
cat("\n")

# ══════════════════════════════════════════════════════════════
# TABLE 3: Mean SE of 100-year return levels by region
# ══════════════════════════════════════════════════════════════
atlantic <- grid$lat > 35
gulf <- grid$lat <= 35

cat("════════════════════════════════════════════════════════════\n")
cat("TABLE 3: Mean SE of 100-Year Return Levels\n")
cat("════════════════════════════════════════════════════════════\n")
cat(sprintf("%-15s %10s %10s %10s\n", "Model", "Overall", "Atlantic", "Gulf"))
cat(sprintf("%-15s %10s %10s %10s\n", "-----", "-------", "--------", "----"))
cat(sprintf("%-15s %10.4f %10.4f %10.4f\n", "Joint",
            mean(rl_joint$se_delta, na.rm = TRUE),
            mean(rl_joint$se_delta[atlantic], na.rm = TRUE),
            mean(rl_joint$se_delta[gulf], na.rm = TRUE)))
cat(sprintf("%-15s %10.4f %10.4f %10.4f\n", "NOAA-only",
            mean(rl_noaa$se_delta, na.rm = TRUE),
            mean(rl_noaa$se_delta[atlantic], na.rm = TRUE),
            mean(rl_noaa$se_delta[gulf], na.rm = TRUE)))
cat(sprintf("%-15s %10.4f %10.4f %10.4f\n", "ADCIRC-only",
            mean(rl_adcirc$se_delta, na.rm = TRUE),
            mean(rl_adcirc$se_delta[atlantic], na.rm = TRUE),
            mean(rl_adcirc$se_delta[gulf], na.rm = TRUE)))
cat("\n")

# ══════════════════════════════════════════════════════════════
# PLOT 1: Comparison return level map (3-panel)
# ══════════════════════════════════════════════════════════════
rl_lims <- range(c(rl_joint$return_level, rl_noaa$return_level,
                    rl_adcirc$return_level))

make_rl_panel <- function(rl_df, title) {
  pts <- st_as_sf(rl_df, coords = c("lon", "lat"), crs = 4326)
  ggplot() +
    basemap +
    geom_sf(data = pts, aes(color = return_level),
            size = 1.2, shape = 16, alpha = 0.7) +
    scale_color_viridis_c(option = "inferno", begin = 0.1, end = 0.95,
                          name = "RL (m)", limits = rl_lims) +
    map_coord() +
    labs(title = title) +
    theme_map +
    theme(legend.key.height = unit(1.2, "cm"))
}

p1a <- make_rl_panel(rl_joint, "Joint (6-dim)")
p1b <- make_rl_panel(rl_noaa, "NOAA-only")
p1c <- make_rl_panel(rl_adcirc, "ADCIRC-only")

p1 <- grid.arrange(p1a, p1b, p1c, nrow = 1)
ggsave("figures/comparison_return_levels.png", p1,
       width = 16, height = 5.5, dpi = 300, bg = "white")
cat("Saved figures/comparison_return_levels.png\n")

# ══════════════════════════════════════════════════════════════
# PLOT 2: Comparison SE map (3-panel)
# ══════════════════════════════════════════════════════════════
se_lims <- range(c(rl_joint$se_delta, rl_noaa$se_delta,
                    rl_adcirc$se_delta), na.rm = TRUE)

make_se_panel <- function(rl_df, title) {
  pts <- st_as_sf(rl_df, coords = c("lon", "lat"), crs = 4326)
  ggplot() +
    basemap +
    geom_sf(data = pts, aes(color = se_delta),
            size = 1.2, shape = 16, alpha = 0.7) +
    scale_color_viridis_c(option = "inferno", name = "SE (m)",
                          limits = se_lims) +
    map_coord() +
    labs(title = title) +
    theme_map +
    theme(legend.key.height = unit(1.2, "cm"))
}

p2a <- make_se_panel(rl_joint, "Joint (6-dim)")
p2b <- make_se_panel(rl_noaa, "NOAA-only")
p2c <- make_se_panel(rl_adcirc, "ADCIRC-only")

p2 <- grid.arrange(p2a, p2b, p2c, nrow = 1)
ggsave("figures/comparison_se.png", p2,
       width = 16, height = 5.5, dpi = 300, bg = "white")
cat("Saved figures/comparison_se.png\n")

# ══════════════════════════════════════════════════════════════
# PLOT 3: SE ratio map (NOAA-only SE / Joint SE)
# ══════════════════════════════════════════════════════════════
ratio_df <- data.frame(
  lon = rl_joint$lon,
  lat = rl_joint$lat,
  ratio = rl_noaa$se_delta / rl_joint$se_delta
)
ratio_df <- ratio_df[is.finite(ratio_df$ratio), ]
ratio_pts <- st_as_sf(ratio_df, coords = c("lon", "lat"), crs = 4326)

p3 <- ggplot() +
  basemap +
  geom_sf(data = ratio_pts, aes(color = ratio),
          size = 1.2, shape = 16, alpha = 0.7) +
  scale_color_viridis_c(option = "plasma", begin = 0.05, end = 0.95,
                        name = "SE ratio",
                        limits = range(ratio_df$ratio)) +
  map_coord() +
  labs(title = "SE(NOAA-only) / SE(joint)") +
  theme_map +
  theme(legend.key.height = unit(1.2, "cm"))

ggsave("figures/se_ratio_map.png", p3,
       width = 10, height = 8, dpi = 300, bg = "white")
cat("Saved figures/se_ratio_map.png\n")

# ══════════════════════════════════════════════════════════════
# PLOT 3b: RL ratio map (Joint RL / NOAA-only RL)
# ══════════════════════════════════════════════════════════════
rl_ratio_df <- data.frame(
  lon = rl_joint$lon,
  lat = rl_joint$lat,
  ratio = rl_joint$return_level / rl_noaa$return_level
)
rl_ratio_df <- rl_ratio_df[is.finite(rl_ratio_df$ratio), ]
rl_ratio_pts <- st_as_sf(rl_ratio_df, coords = c("lon", "lat"), crs = 4326)

p3b <- ggplot() +
  basemap +
  geom_sf(data = rl_ratio_pts, aes(color = ratio),
          size = 1.2, shape = 16, alpha = 0.7) +
  scale_color_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                        midpoint = 1, limits = c(0.8, 1.3),
                        oob = scales::squish,
                        name = "RL ratio") +
  map_coord() +
  labs(title = "RL(joint) / RL(NOAA-only)") +
  theme_map +
  theme(legend.key.height = unit(1.2, "cm"))

ggsave("figures/rl_ratio_map.png", p3b,
       width = 10, height = 8, dpi = 300, bg = "white")
cat("Saved figures/rl_ratio_map.png\n")

# ══════════════════════════════════════════════════════════════
# PLOT 4: Scatter at 29 NOAA sites (joint vs NOAA-only)
# ══════════════════════════════════════════════════════════════
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
    plot.title = element_text(size = 12, hjust = 0.5)
  )

p4a <- ggplot(scatter_df, aes(x = mu_noaa, y = mu_joint)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 2, alpha = 0.8) +
  labs(x = "NOAA-only", y = "Joint", title = expression(mu)) +
  scatter_theme + coord_fixed()

p4b <- ggplot(scatter_df, aes(x = ls_noaa, y = ls_joint)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 2, alpha = 0.8) +
  labs(x = "NOAA-only", y = "Joint", title = expression(log~sigma)) +
  scatter_theme + coord_fixed()

p4c <- ggplot(scatter_df, aes(x = xi_noaa, y = xi_joint)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 2, alpha = 0.8) +
  labs(x = "NOAA-only", y = "Joint", title = expression(xi)) +
  scatter_theme + coord_fixed()

p4 <- grid.arrange(p4a, p4b, p4c, nrow = 1)
ggsave("figures/noaa_site_comparison.png", p4,
       width = 14, height = 5, dpi = 300, bg = "white")
cat("Saved figures/noaa_site_comparison.png\n")

# ══════════════════════════════════════════════════════════════
# PLOT 5: Return level profile by distance along coast
# ══════════════════════════════════════════════════════════════

# Order grid points along the coast using nearest-neighbor chains
# within two geographic regions (Gulf, then Atlantic) to ensure
# Houston → New Orleans → Tampa → Key West → Miami → ... → Maine.
n_grid <- nrow(grid)
gulf_mask <- grid$lon < -82
gulf_idx <- which(gulf_mask)
atl_idx  <- which(!gulf_mask)

nn_chain <- function(idx, start) {
  n <- length(idx)
  visited <- logical(n)
  chain <- integer(n)
  chain[1] <- start
  visited[match(start, idx)] <- TRUE
  for (k in 2:n) {
    prev <- chain[k - 1]
    d2 <- (grid$lon[idx] - grid$lon[prev])^2 +
           (grid$lat[idx] - grid$lat[prev])^2
    d2[visited] <- Inf
    best <- which.min(d2)
    chain[k] <- idx[best]
    visited[best] <- TRUE
  }
  chain
}

gulf_start <- gulf_idx[which.min(grid$lon[gulf_idx])]
gulf_chain <- nn_chain(gulf_idx, gulf_start)

# Start Atlantic chain from the point nearest the end of Gulf chain
gulf_end <- gulf_chain[length(gulf_chain)]
d2_to_gulf_end <- (grid$lon[atl_idx] - grid$lon[gulf_end])^2 +
                  (grid$lat[atl_idx] - grid$lat[gulf_end])^2
atl_start <- atl_idx[which.min(d2_to_gulf_end)]
atl_chain <- nn_chain(atl_idx, atl_start)

order_idx <- c(gulf_chain, atl_chain)

# Compute cumulative geodesic distance along the ordered chain
haversine <- function(lon1, lat1, lon2, lat2) {
  R <- 6371
  dlon <- (lon2 - lon1) * pi / 180
  dlat <- (lat2 - lat1) * pi / 180
  a <- sin(dlat / 2)^2 + cos(lat1 * pi / 180) * cos(lat2 * pi / 180) * sin(dlon / 2)^2
  2 * R * asin(sqrt(a))
}

cum_dist <- numeric(n_grid)
for (k in 2:n_grid) {
  i <- order_idx[k]
  j <- order_idx[k - 1]
  cum_dist[k] <- cum_dist[k - 1] + haversine(grid$lon[j], grid$lat[j],
                                               grid$lon[i], grid$lat[i])
}

# Map: original grid index -> cumulative distance
dist_km <- numeric(n_grid)
dist_km[order_idx] <- cum_dist

# Reference cities (lon, lat, label)
cities <- data.frame(
  label = c("Key West", "Tampa", "New Orleans", "Houston",
            "Miami", "Charleston", "Norfolk", "NYC", "Boston", "Portland"),
  lon = c(-81.78, -82.46, -90.07, -95.36,
          -80.19, -79.93, -76.29, -74.01, -71.06, -70.26),
  lat = c(24.56, 27.95, 29.95, 29.76,
          25.76, 32.78, 36.85, 40.71, 42.36, 43.66),
  stringsAsFactors = FALSE
)
# Find nearest grid point to each city
city_dist <- vapply(seq_len(nrow(cities)), function(c) {
  d2 <- (grid$lon - cities$lon[c])^2 + (grid$lat - cities$lat[c])^2
  dist_km[which.min(d2)]
}, numeric(1))
cities$dist_km <- city_dist

# Build profile data frame
profile_df <- rbind(
  data.frame(dist = dist_km, rl = rl_joint$return_level,
             se = rl_joint$se_delta, model = "Joint"),
  data.frame(dist = dist_km, rl = rl_noaa$return_level,
             se = rl_noaa$se_delta, model = "NOAA-only"),
  data.frame(dist = dist_km, rl = rl_adcirc$return_level,
             se = rl_adcirc$se_delta, model = "ADCIRC-only")
)
profile_df$lower <- profile_df$rl - 1.96 * profile_df$se
profile_df$upper <- profile_df$rl + 1.96 * profile_df$se
profile_df <- profile_df[order(profile_df$model, profile_df$dist), ]

p5 <- ggplot(profile_df, aes(x = dist, y = rl, color = model, fill = model)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = c("Joint" = "#E41A1C", "NOAA-only" = "#377EB8",
                                 "ADCIRC-only" = "#4DAF4A"),
                      name = NULL) +
  scale_fill_manual(values = c("Joint" = "#E41A1C", "NOAA-only" = "#377EB8",
                                "ADCIRC-only" = "#4DAF4A"),
                     name = NULL) +
  scale_x_continuous(
    breaks = cities$dist_km,
    labels = cities$label
  ) +
  labs(x = "Distance along coast (Houston to Maine)",
       y = "100-year return level (m)") +
  theme_minimal(base_size = 11) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9)
  )

ggsave("figures/return_level_profile.png", p5,
       width = 12, height = 6, dpi = 300, bg = "white")
cat("Saved figures/return_level_profile.png\n")

cat("\nAll comparison plots saved to figures/\n")
