#!/usr/bin/env Rscript
# Trend diagnostics at 29 NOAA sites
library(ggplot2)
library(sf)
devtools::load_all()

dir.create("figures", showWarnings = FALSE)

# ── Load data ───────────────────────────────────────────────
df <- read.csv("data-raw/combined_max.csv")
dat <- load_data(df)
noaa_df <- df[df$data_source == "NOAA", ]
noaa_sites <- dat$sites[dat$sites$data_source == "NOAA", ]
locations <- unique(noaa_df$location)

# ── Run diagnostics at each site ────────────────────────────
results <- vector("list", length(locations))

for (i in seq_along(locations)) {
  loc <- locations[i]
  site_df <- noaa_df[noaa_df$location == loc, ]
  site_df <- site_df[order(site_df$year), ]
  site_info <- noaa_sites[noaa_sites$location == loc, ]

  mk <- mann_kendall_test(site_df$max_sea_level, site_df$year)
  gev <- tryCatch(
    gev_trend_test(site_df$max_sea_level, site_df$year),
    error = function(e) {
      warning("GEV trend test failed at site ", loc, ": ", e$message)
      list(mu0 = NA, mu1 = NA, lrt_stat = NA, lrt_pvalue = NA,
           aic_stat = NA, aic_nonstat = NA)
    }
  )

  # Sen's slope CI via normal approximation on S
  n <- nrow(site_df)
  var_S <- n * (n - 1) * (2 * n + 5) / 18
  se_slope <- sqrt(var_S) / (n * (n - 1) / 2) *
    sd(site_df$max_sea_level) / sd(site_df$year)

  results[[i]] <- data.frame(
    location = loc,
    lon = site_info$lon[1],
    lat = site_info$lat[1],
    n_years = nrow(site_df),
    mk_tau = mk$tau,
    mk_pvalue = mk$p_value,
    sens_slope = mk$sens_slope,
    sens_slope_se = se_slope,
    mu1_hat = gev$mu1,
    lrt_pvalue = gev$lrt_pvalue,
    aic_stat = gev$aic_stat,
    aic_nonstat = gev$aic_nonstat,
    stringsAsFactors = FALSE
  )
}

tab <- do.call(rbind, results)
tab <- tab[order(tab$lat), ]
rownames(tab) <- NULL

# ══════════════════════════════════════════════════════════════
# Summary table
# ══════════════════════════════════════════════════════════════
cat("\n")
cat("════════════════════════════════════════════════════════════════════════════════════════════\n")
cat("Trend Diagnostics at 29 NOAA Sites\n")
cat("════════════════════════════════════════════════════════════════════════════════════════════\n")
cat(sprintf("%-10s %6s %6s %7s %10s %9s %10s %10s %10s\n",
    "Site", "Lat", "N", "MK_tau", "MK_p", "Sen_slope", "mu1_hat", "LRT_p", "dAIC"))
cat(sprintf("%-10s %6s %6s %7s %10s %9s %10s %10s %10s\n",
    "----", "---", "--", "------", "----", "---------", "-------", "-----", "----"))

for (i in seq_len(nrow(tab))) {
  daic <- tab$aic_nonstat[i] - tab$aic_stat[i]
  mk_flag <- ifelse(tab$mk_pvalue[i] < 0.05, "*", " ")
  lrt_flag <- ifelse(!is.na(tab$lrt_pvalue[i]) && tab$lrt_pvalue[i] < 0.05, "*", " ")
  cat(sprintf("%-10s %6.1f %6d %7.3f %9.4f%s %9.4f %10.5f %9.4f%s %10.2f\n",
      tab$location[i], tab$lat[i], tab$n_years[i],
      tab$mk_tau[i], tab$mk_pvalue[i], mk_flag,
      tab$sens_slope[i], tab$mu1_hat[i],
      tab$lrt_pvalue[i], lrt_flag, daic))
}

n_mk_sig <- sum(tab$mk_pvalue < 0.05)
n_lrt_sig <- sum(tab$lrt_pvalue < 0.05, na.rm = TRUE)
n_aic_nonstat <- sum(tab$aic_nonstat < tab$aic_stat, na.rm = TRUE)

cat("\n")
cat(sprintf("Mann-Kendall significant (p < 0.05): %d / %d sites\n", n_mk_sig, nrow(tab)))
cat(sprintf("GEV LRT significant (p < 0.05):      %d / %d sites\n", n_lrt_sig, nrow(tab)))
cat(sprintf("AIC favors nonstationary:             %d / %d sites\n", n_aic_nonstat, nrow(tab)))
cat(sprintf("Median Sen's slope: %.4f m/yr\n", median(tab$sens_slope)))
cat(sprintf("Mean |Sen's slope|: %.4f m/yr\n", mean(abs(tab$sens_slope))))
cat("\n")

# ── Save results ────────────────────────────────────────────
saveRDS(tab, "data-raw/trend_diagnostics.rds")

sink("data-raw/trend_diagnostics.txt")
cat("Trend Diagnostics at 29 NOAA Sites\n")
cat("Site, Lat, N, MK_tau, MK_pvalue, Sens_slope, mu1_hat, LRT_pvalue, AIC_stat, AIC_nonstat\n")
for (i in seq_len(nrow(tab))) {
  cat(sprintf("%s, %.4f, %d, %.4f, %.4f, %.5f, %.5f, %.4f, %.2f, %.2f\n",
      tab$location[i], tab$lat[i], tab$n_years[i],
      tab$mk_tau[i], tab$mk_pvalue[i], tab$sens_slope[i],
      tab$mu1_hat[i], tab$lrt_pvalue[i],
      tab$aic_stat[i], tab$aic_nonstat[i]))
}
sink()
cat("Saved data-raw/trend_diagnostics.rds\n")
cat("Saved data-raw/trend_diagnostics.txt\n")

# ══════════════════════════════════════════════════════════════
# PLOT 1: Sen's slope dot plot ordered by latitude
# ══════════════════════════════════════════════════════════════
tab$sig <- tab$mk_pvalue < 0.05
tab$location_f <- factor(tab$location, levels = tab$location)
tab$ci_lower <- tab$sens_slope - 1.96 * tab$sens_slope_se
tab$ci_upper <- tab$sens_slope + 1.96 * tab$sens_slope_se

p1 <- ggplot(tab, aes(x = sens_slope, y = location_f, color = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper),
                 height = 0.3, linewidth = 0.4) +
  geom_point(size = 2.5) +
  scale_color_manual(values = c("FALSE" = "grey50", "TRUE" = "#E41A1C"),
                     labels = c("p >= 0.05", "p < 0.05"),
                     name = "Mann-Kendall") +
  labs(x = "Sen's slope (m/year)", y = "NOAA site (ordered by latitude)") +
  theme_minimal(base_size = 11) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.position = "top"
  )

ggsave("figures/trend_summary.png", p1,
       width = 8, height = 8, dpi = 300, bg = "white")
cat("Saved figures/trend_summary.png\n")

# ══════════════════════════════════════════════════════════════
# PLOT 2: Map of NOAA sites colored by Sen's slope
# ══════════════════════════════════════════════════════════════
states_df <- map_data("state")
states_sf <- lapply(split(states_df, states_df$group), function(grp) {
  st_polygon(list(as.matrix(grp[, c("long", "lat")])))
})
states_sf <- st_sfc(states_sf, crs = 4326)
states_sf <- st_sf(geometry = states_sf)
crs_albers <- st_crs(5070)

slope_lim <- max(abs(tab$sens_slope)) * c(-1, 1)

p2 <- ggplot() +
  geom_sf(data = states_sf, fill = "grey95", color = "grey60",
          linewidth = 0.3, inherit.aes = FALSE) +
  geom_point(data = tab, aes(x = lon, y = lat, color = sens_slope,
                              shape = sig), size = 3.5) +
  scale_color_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                        midpoint = 0, limits = slope_lim,
                        name = "Sen's slope\n(m/year)") +
  scale_shape_manual(values = c("FALSE" = 1, "TRUE" = 16),
                     labels = c("p >= 0.05", "p < 0.05"),
                     name = "Mann-Kendall") +
  coord_sf(crs = crs_albers,
           xlim = c(-98, -66), ylim = c(24, 46),
           default_crs = st_crs(4326), expand = FALSE) +
  theme_minimal(base_size = 11) +
  theme(
    axis.title = element_blank(),
    axis.text = element_text(size = 8, color = "grey40"),
    panel.grid = element_blank(),
    legend.position = "right",
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

ggsave("figures/trend_map.png", p2,
       width = 10, height = 7, dpi = 300, bg = "white")
cat("Saved figures/trend_map.png\n")
