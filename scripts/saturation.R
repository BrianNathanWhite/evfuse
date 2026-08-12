#!/usr/bin/env Rscript
# ADCIRC subsampling saturation curve (Figure 7, Section 4.5).
#
# Quantifies how fusion benefits depend on ADCIRC network density under the
# NONSTATIONARY Stage 1 pipeline (the model reported in the paper). Subsamples
# the 100 ADCIRC sites to n in {5, 10, 15, 25, 50, 75, 100} by maximin design,
# refits the joint model at each density (20 multi-starts + warm start), and
# computes LOO-CV metrics at the 29 NOAA sites against the NOAA-only baseline.
#
# Requires the model objects saved by scripts/run_nonstationary.R:
#   data-raw/model_6dim_ns.rds, data-raw/model_noaa_only_ns.rds
#
# The theoretical saturation scale uses L ~ 6,000 km: the NOAA general
# coastline length for the study domain (Atlantic 3,330 km + Gulf 2,625 km
# = 5,955 km; NOAA shoreline statistics), not the fractal shoreline length,
# with rho_min from the fitted joint model.
#
# Usage: Rscript scripts/saturation.R   (~30-60 min)

devtools::load_all()
library(ggplot2)
source("scripts/fig_theme.R")

# ── 1. Load pre-fitted nonstationary models ─────────────────
joint_full <- readRDS("data-raw/model_6dim_ns.rds")
dat    <- joint_full$dat
stage1 <- joint_full$stage1
D      <- joint_full$D
L      <- dat$n_sites  # 129

# Rebuild the raw observed-space bootstrap W exactly as the pipeline does
# (deterministic with seed 42). NOTE: joint_full$W_tap is the EMBEDDED
# 774-dim (6L) version and must not be subset as if it were 387-dim.
data(coast_data, package = "evfuse")
df <- coast_data$raw_df
cat("Running bootstrap (B=500, seed=42)...\n")
bs   <- bootstrap_W_detrended(dat, df, B = 500, ref_year = 2000, seed = 42)
W_bs <- bs$W_bs
stopifnot(nrow(W_bs) == 3 * L)

noaa_full <- readRDS("data-raw/model_noaa_only_ns.rds")
loo_noaa  <- loo_cv(noaa_full)
sum_noaa  <- loo_summary(loo_noaa, r = 100)

cat(sprintf("NOAA-only baseline (nonstationary): RL RMSE = %.4f, Total LPD = %.2f\n\n",
    sum_noaa$rl_rmse, sum_noaa$total_lpd))

# ── 2. Maximin subsampling ──────────────────────────────────
# Greedy deletion: repeatedly remove the ADCIRC site whose removal
# maximizes the minimum nearest-neighbor distance among remaining sites.
adcirc_idx <- which(dat$sites$data_source == "ADCIRC")
noaa_idx   <- which(dat$sites$data_source == "NOAA")
n_adcirc_full <- length(adcirc_idx)

D_adcirc <- D[adcirc_idx, adcirc_idx]

maximin_subsample <- function(D_full, n_target) {
  n <- nrow(D_full)
  if (n_target >= n) return(seq_len(n))
  active <- seq_len(n)
  while (length(active) > n_target) {
    best_remove <- NA
    best_min_nn <- -Inf
    for (i in seq_along(active)) {
      remain <- active[-i]
      D_r <- D_full[remain, remain]
      diag(D_r) <- Inf
      min_nn <- min(apply(D_r, 1, min))
      if (min_nn > best_min_nn) {
        best_min_nn <- min_nn
        best_remove <- i
      }
    }
    active <- active[-best_remove]
  }
  active
}

# n = 100 runs first: it must reproduce the production model's LOO-CV
# exactly (identity check), else something is wrong with the subsetting.
n_values <- c(100, 75, 50, 25, 15, 10, 5)
cat("Computing maximin subsamples...\n")
selections <- list()
selections[["100"]] <- seq_len(n_adcirc_full)
for (nt in rev(n_values[n_values < 100])) {
  cat(sprintf("  n=%d...", nt))
  selections[[as.character(nt)]] <- maximin_subsample(D_adcirc, nt)
  cat(" done\n")
}

# ── 3. Subsetting helpers ──────────────────────────────────
subset_dat <- function(dat, keep_idx) {
  out <- dat
  out$sites    <- dat$sites[keep_idx, ]
  out$maxima   <- dat$maxima[keep_idx]
  out$n_noaa   <- sum(dat$sites$data_source[keep_idx] == "NOAA")
  out$n_adcirc <- sum(dat$sites$data_source[keep_idx] == "ADCIRC")
  out$n_sites  <- length(keep_idx)
  out
}

subset_stage1 <- function(stage1, keep_idx) {
  out <- stage1
  n <- nrow(stage1$theta_hat)
  for (nm in names(out)) {
    x <- out[[nm]]
    if (is.matrix(x) && nrow(x) == n)           out[[nm]] <- x[keep_idx, , drop = FALSE]
    else if (is.list(x) && length(x) == n)      out[[nm]] <- x[keep_idx]
    else if (is.atomic(x) && length(x) == n)    out[[nm]] <- x[keep_idx]
  }
  out
}

subset_W_idx <- function(W, L, keep_idx) {
  w_idx <- c(keep_idx, keep_idx + L, keep_idx + 2 * L)
  W[w_idx, w_idx]
}

# ── 4. Main loop over ADCIRC densities ─────────────────────
results <- data.frame(
  n_adcirc           = n_values,
  mean_spacing_km    = NA_real_,
  mu_rmse_joint      = NA_real_,
  logsig_rmse_joint  = NA_real_,
  xi_rmse_joint      = NA_real_,
  rl_rmse_joint      = NA_real_,
  rl_rmse_noaa       = sum_noaa$rl_rmse,
  rl_rmse_reduction_pct = NA_real_,
  total_lpd_joint    = NA_real_,
  total_lpd_noaa     = sum_noaa$total_lpd,
  lpd_gain           = NA_real_,
  joint_wins         = NA_integer_
)
site_selections <- list()

for (r in seq_along(n_values)) {
  nt <- n_values[r]
  sel_adcirc <- selections[[as.character(nt)]]
  keep_idx <- sort(c(noaa_idx, adcirc_idx[sel_adcirc]))

  D_sel <- D_adcirc[sel_adcirc, sel_adcirc]
  diag(D_sel) <- Inf
  results$mean_spacing_km[r] <- mean(apply(D_sel, 1, min))

  site_selections[[r]] <- data.frame(
    n_adcirc   = nt,
    station_id = dat$sites$location[adcirc_idx[sel_adcirc]],
    lon        = dat$sites$lon[adcirc_idx[sel_adcirc]],
    lat        = dat$sites$lat[adcirc_idx[sel_adcirc]]
  )

  cat(sprintf("== n_ADCIRC = %d (%d total sites, mean spacing = %.0f km) ==\n",
      nt, length(keep_idx), results$mean_spacing_km[r]))

  dat_sub    <- subset_dat(dat, keep_idx)
  stage1_sub <- subset_stage1(stage1, keep_idx)
  D_sub      <- D[keep_idx, keep_idx]
  W_tap_sub  <- taper_W(subset_W_idx(W_bs, L, keep_idx), D_sub, lambda = 300)

  cat("  Fitting joint model...")
  start0 <- list(beta = joint_full$beta, A = joint_full$A, rho = joint_full$rho)
  best_fit <- tryCatch(
    suppressWarnings(
      fit_spatial_model(stage1_sub, dat_sub, W_tap_sub, D_sub,
                        start = start0,
                        control = list(maxit = 2000, trace = 0))
    ),
    error = function(e) NULL
  )
  best_nll <- if (!is.null(best_fit)) best_fit$optim_result$value else Inf

  set.seed(2026)
  for (s in seq_len(20)) {
    rho_s <- exp(runif(6, log(50), log(5000)))
    start_s <- list(beta = joint_full$beta, A = joint_full$A, rho = rho_s)
    fit_s <- tryCatch(
      suppressWarnings(
        fit_spatial_model(stage1_sub, dat_sub, W_tap_sub, D_sub,
                          start = start_s,
                          control = list(maxit = 2000, trace = 0))
      ),
      error = function(e) NULL
    )
    if (!is.null(fit_s) && fit_s$optim_result$value < best_nll) {
      best_nll <- fit_s$optim_result$value
      best_fit <- fit_s
    }
  }
  cat(sprintf(" NLL = %.4f\n", best_nll))

  cat("  Computing LOO-CV...")
  loo_j <- tryCatch(loo_cv(best_fit), error = function(e) NULL)
  if (!is.null(loo_j)) {
    sum_j <- loo_summary(loo_j, r = 100)
    results$mu_rmse_joint[r]     <- sum_j$param_stats$rmse[1]
    results$logsig_rmse_joint[r] <- sum_j$param_stats$rmse[2]
    results$xi_rmse_joint[r]     <- sum_j$param_stats$rmse[3]
    results$rl_rmse_joint[r]     <- sum_j$rl_rmse
    results$total_lpd_joint[r]   <- sum_j$total_lpd
    results$lpd_gain[r]          <- sum_j$total_lpd - sum_noaa$total_lpd
    results$joint_wins[r]        <- sum(loo_j$loo_lpd > loo_noaa$loo_lpd)
    results$rl_rmse_reduction_pct[r] <- 100 * (1 - sum_j$rl_rmse / sum_noaa$rl_rmse)
    cat(sprintf(" RL RMSE = %.4f, reduction = %.1f%%, LPD gain = %.2f\n",
        sum_j$rl_rmse, results$rl_rmse_reduction_pct[r], results$lpd_gain[r]))
  } else {
    cat(" LOO-CV FAILED\n")
  }

  # Identity check: the full-network "subset" must reproduce the production
  # model's LOO-CV (paper Table 2). Abort loudly rather than produce a curve
  # inconsistent with the paper.
  if (nt == n_adcirc_full) {
    loo_prod <- loo_summary(loo_cv(joint_full), r = 100)
    if (abs(results$rl_rmse_joint[r] - loo_prod$rl_rmse) > 1e-4) {
      stop(sprintf(
        "n=100 identity check FAILED: refit RL RMSE %.6f vs production %.6f",
        results$rl_rmse_joint[r], loo_prod$rl_rmse))
    }
    cat(sprintf("  Identity check passed: n=100 RL RMSE matches production (%.4f)\n",
        loo_prod$rl_rmse))
  }
}

results <- results[order(results$n_adcirc), ]

# ── 5. Summary output ──────────────────────────────────────
rho_min <- min(joint_full$rho)
output_text <- function() {
  cat("\n")
  cat("======================================================================\n")
  cat("  ADCIRC Subsampling Saturation Curve (nonstationary pipeline)\n")
  cat("======================================================================\n")
  cat(sprintf("  NOAA-only baseline: RL RMSE = %.4f, Total LPD = %.2f\n\n",
      sum_noaa$rl_rmse, sum_noaa$total_lpd))

  cat(sprintf("  %-8s %12s %9s %9s %12s %10s %6s\n",
      "n_ADCIRC", "Mean Spacing", "mu RMSE", "RL RMSE", "RL Reduction", "LPD Gain", "Wins"))
  cat(sprintf("  %-8s %12s %9s %9s %12s %10s %6s\n",
      "--------", "------------", "-------", "-------", "------------", "--------", "----"))
  for (r in seq_along(n_values)) {
    cat(sprintf("  %8d %10.0f km %9.4f %9.4f %11.1f%% %10.2f %4d/29\n",
        results$n_adcirc[r], results$mean_spacing_km[r],
        results$mu_rmse_joint[r], results$rl_rmse_joint[r],
        results$rl_rmse_reduction_pct[r], results$lpd_gain[r],
        results$joint_wins[r]))
  }

  cat(sprintf("\n  Theoretical saturation scale: n_eff = L/rho_min = 6000/%.0f = %.0f\n",
      rho_min, 6000 / rho_min))
  cat("  (L ~ 6,000 km: NOAA general coastline length for the domain,\n")
  cat("   Atlantic 3,330 km + Gulf 2,625 km = 5,955 km)\n")
  cat(sprintf("  Smallest rho in full nonstationary model: %.1f km\n", rho_min))

  cat("\n  Detail by parameter RMSE:\n")
  cat(sprintf("  %-8s %9s %9s %9s\n", "n_ADCIRC", "mu", "log_sig", "xi"))
  cat(sprintf("  %-8s %9s %9s %9s\n", "--------", "-------", "-------", "-------"))
  for (r in seq_along(n_values)) {
    cat(sprintf("  %8d %9.4f %9.4f %9.4f\n",
        results$n_adcirc[r], results$mu_rmse_joint[r],
        results$logsig_rmse_joint[r], results$xi_rmse_joint[r]))
  }
  cat("\n")
}

output_text()

sink("tables/saturation_curve_summary.txt")
output_text()
sink()
cat("Saved tables/saturation_curve_summary.txt\n")

write.csv(results, "tables/saturation_curve.csv", row.names = FALSE)
cat("Saved tables/saturation_curve.csv\n")

site_sel_df <- do.call(rbind, site_selections)
write.csv(site_sel_df, "tables/saturation_site_selections.csv", row.names = FALSE)
cat("Saved tables/saturation_site_selections.csv\n")

saveRDS(list(results = results, selections = selections,
             noaa_summary = sum_noaa),
        "data-raw/saturation_results_ns.rds")
cat("Saved data-raw/saturation_results_ns.rds\n")

# ── 6. Figure 7 ────────────────────────────────────────────
n_eff <- round(6000 / rho_min)  # theoretical: coastline length / smallest range
n_emp <- 75                     # empirical plateau
ymax  <- 40

p <- ggplot(results, aes(x = n_adcirc, y = rl_rmse_reduction_pct)) +
  geom_vline(xintercept = n_eff, linetype = "dashed", colour = "grey50",
             linewidth = 0.5) +
  geom_vline(xintercept = n_emp, linetype = "dashed", colour = "grey50",
             linewidth = 0.5) +
  annotate("text", x = n_eff + 1.5, y = ymax * 0.18,
           label = paste0("italic(n)[eff] %~~% ", n_eff),
           parse = TRUE, hjust = 0, size = 3.2, colour = "grey40") +
  annotate("text", x = n_eff + 1.5, y = ymax * 0.10,
           label = "(heuristic)", hjust = 0, size = 3.0, colour = "grey50") +
  annotate("text", x = n_emp + 1.5, y = ymax * 0.18,
           label = paste0("italic(n)[plateau] %~~% ", n_emp),
           parse = TRUE, hjust = 0, size = 3.2, colour = "grey40") +
  annotate("text", x = n_emp + 1.5, y = ymax * 0.10,
           label = "(empirical)", hjust = 0, size = 3.0, colour = "grey50") +
  geom_line(linewidth = 0.8, colour = "#2166AC") +
  geom_point(size = 2.5, colour = "#2166AC") +
  scale_x_continuous(breaks = n_values) +
  scale_y_continuous(limits = c(0, ymax),
                     labels = function(x) paste0(x, "%")) +
  labs(x = expression(n[ADCIRC]),
       y = "100-yr RL RMSE reduction vs. NOAA-only") +
  theme_bw_nogrid(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.margin = margin(10, 15, 10, 10))

ggsave("figures/saturation_curve.png", p,
       width = 6.5, height = 4, dpi = 300, bg = "white")
cat("Saved figures/saturation_curve.png\n")
