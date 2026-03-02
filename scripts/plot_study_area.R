#!/usr/bin/env Rscript
# plot_study_area.R: Study area map for evfuse package
# Shows NOAA tidal gauge and ADCIRC simulation site locations along the
# US East and Gulf coasts.

library(ggplot2)

devtools::load_all()
data(coast_data)

# ---------- site data ---------------------------------------------------------
sites <- coast_data$sites
noaa   <- sites[sites$data_source == "NOAA",  ]
adcirc <- sites[sites$data_source == "ADCIRC",]

# ---------- basemap -----------------------------------------------------------
states <- ggplot2::map_data("state")

# ---------- city labels -------------------------------------------------------
cities <- data.frame(
  city = c("Houston", "New Orleans", "Tampa", "Key West", "Miami",
           "Charleston", "Norfolk", "New York", "Boston", "Portland ME"),
  lon  = c(-95.37, -90.07, -82.46, -81.78, -80.19,
           -79.93, -76.29, -74.01, -71.06, -70.26),
  lat  = c( 29.76,  29.95,  27.95,  24.56,  25.76,
            32.78,  36.85,  40.71,  42.36,  43.66)
)

# Manual nudge offsets (degrees) so labels don't overlap points.
#                Houston  NewOrl  Tampa  KeyW   Miami  Charl  Norf   NY     Bos    Port
cities$nudge_x <- c(-0.8,  0.0,  1.0, -1.5,  1.2,   1.2,   1.2,   1.5,   1.5,   1.5)
cities$nudge_y <- c( 0.6,  0.7,  0.5,  0.3,  0.5,  -0.5,  -0.5,  -0.5,   0.5,   0.5)

# ---------- plot --------------------------------------------------------------
p <- ggplot() +
  # State polygons (basemap)
  geom_polygon(data = states,
               aes(x = long, y = lat, group = group),
               fill = "gray90", colour = "white", linewidth = 0.3) +

  # ADCIRC sites (triangles, plotted first so NOAA overlays)
  geom_point(data = adcirc,
             aes(x = lon, y = lat, colour = "ADCIRC simulation",
                 shape = "ADCIRC simulation"),
             size = 2) +

  # NOAA sites (circles, on top)
  geom_point(data = noaa,
             aes(x = lon, y = lat, colour = "NOAA tidal gauge",
                 shape = "NOAA tidal gauge"),
             size = 3) +

  # City markers (star shape for clear distinction from data points)
  geom_point(data = cities,
             aes(x = lon, y = lat),
             shape = 8, size = 2, stroke = 0.6, colour = "black") +

  # City labels (manual nudge since ggrepel not available)
  geom_text(data = cities,
            aes(x = lon + nudge_x, y = lat + nudge_y, label = city),
            size = 3, colour = "gray20") +

  # Scales
  scale_colour_manual(
    name = NULL,
    values = c("NOAA tidal gauge"   = "#1f4e79",
               "ADCIRC simulation"  = "#d95f02"),
    breaks = c("NOAA tidal gauge", "ADCIRC simulation")
  ) +
  scale_shape_manual(
    name = NULL,
    values = c("NOAA tidal gauge"  = 16,   # filled circle
               "ADCIRC simulation" = 17),   # filled triangle
    breaks = c("NOAA tidal gauge", "ADCIRC simulation")
  ) +

  # Coordinate limits
  coord_cartesian(xlim = c(-98, -66), ylim = c(24, 46), expand = FALSE) +

  # Labels
  labs(x = "Longitude", y = "Latitude") +

  # Theme
  theme_minimal(base_size = 12) +
  theme(
    plot.background   = element_rect(fill = "white", colour = NA),
    panel.background  = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(colour = "gray92", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    legend.position  = "inside",
    legend.position.inside = c(0.82, 0.18),
    legend.background = element_rect(fill = alpha("white", 0.85),
                                      colour = NA),
    legend.key = element_rect(fill = NA, colour = NA),
    plot.title = element_text(face = "bold", size = 14)
  )

# ---------- save --------------------------------------------------------------
out_dir <- file.path(getwd(), "figures")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

ggsave(file.path(out_dir, "study_area_map.png"), plot = p,
       width = 10, height = 8, dpi = 300, bg = "white")

message("Saved to figures/study_area_map.png")
