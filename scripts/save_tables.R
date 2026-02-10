#!/usr/bin/env Rscript
# Save comparison tables as formatted text and RDS
devtools::load_all()

# ── Load all model outputs ────────────────────────────────────
joint    <- readRDS("data-raw/model_6dim_best.rds")
noaa_m   <- readRDS("data-raw/model_noaa_only.rds")
adcirc_m <- readRDS("data-raw/model_adcirc_only.rds")

rl_joint  <- readRDS("data-raw/return_levels_100yr.rds")
rl_noaa   <- readRDS("data-raw/rl_noaa_only.rds")
rl_adcirc <- readRDS("data-raw/rl_adcirc_only.rds")

grid   <- read.csv("data-raw/prediction_grid.csv")
sites  <- joint$dat$sites
stage1 <- joint$stage1

# ── Krige at NOAA sites for Table 2 ──────────────────────────
noaa_idx <- which(sites$data_source == "NOAA")
noaa_sites <- sites[noaa_idx, ]

cat("Kriging at 29 NOAA sites...\n")
preds_joint_noaa <- predict_krig(joint, noaa_sites)
preds_noaa_noaa  <- predict_krig_naive(noaa_m, noaa_sites)
theta_noaa <- stage1$theta_hat[noaa_idx, ]

rl_joint_noaa <- vapply(seq_len(nrow(noaa_sites)), function(k) {
  gev_return_level(preds_joint_noaa$noaa_mean[k, "mu"],
                   exp(preds_joint_noaa$noaa_mean[k, "log_sigma"]),
                   preds_joint_noaa$noaa_mean[k, "xi"], r = 100)
}, numeric(1))
rl_noaa_noaa <- vapply(seq_len(nrow(noaa_sites)), function(k) {
  gev_return_level(preds_noaa_noaa$noaa_mean[k, "mu"],
                   exp(preds_noaa_noaa$noaa_mean[k, "log_sigma"]),
                   preds_noaa_noaa$noaa_mean[k, "xi"], r = 100)
}, numeric(1))
rl_obs <- vapply(seq_along(noaa_idx), function(k) {
  gev_return_level(theta_noaa[k, "mu"], exp(theta_noaa[k, "log_sigma"]),
                   theta_noaa[k, "xi"], r = 100)
}, numeric(1))

# ══════════════════════════════════════════════════════════════
# TABLE 1: Model Fit Summary
# ══════════════════════════════════════════════════════════════
tab1 <- data.frame(
  Model = c("Joint", "NOAA_only", "ADCIRC_only"),
  Sites = c(129, 29, 100),
  Num_Parameters = c(33, 12, 12),
  NLL = c(joint$optim_result$value, noaa_m$optim_result$value, adcirc_m$optim_result$value),
  rho = c(paste(round(joint$rho, 1), collapse = ", "),
          paste(round(noaa_m$rho, 1), collapse = ", "),
          paste(round(adcirc_m$rho, 1), collapse = ", ")),
  stringsAsFactors = FALSE
)

# ══════════════════════════════════════════════════════════════
# TABLE 2: NOAA Site Comparison
# ══════════════════════════════════════════════════════════════
params <- c("mu", "log_sigma", "xi")
tab2_rows <- list()
for (param in params) {
  j <- match(param, params)
  obs <- theta_noaa[, j]
  pred_j <- preds_joint_noaa$noaa_mean[, param]
  pred_n <- preds_noaa_noaa$noaa_mean[, param]
  tab2_rows[[param]] <- data.frame(
    Parameter = param,
    Joint_RMSE = sqrt(mean((pred_j - obs)^2)),
    NOAA_only_RMSE = sqrt(mean((pred_n - obs)^2)),
    Joint_MAD = mean(abs(pred_j - obs)),
    NOAA_only_MAD = mean(abs(pred_n - obs)),
    stringsAsFactors = FALSE
  )
}
tab2_rows[["100yr_RL"]] <- data.frame(
  Parameter = "100yr_RL",
  Joint_RMSE = sqrt(mean((rl_joint_noaa - rl_obs)^2)),
  NOAA_only_RMSE = sqrt(mean((rl_noaa_noaa - rl_obs)^2)),
  Joint_MAD = mean(abs(rl_joint_noaa - rl_obs)),
  NOAA_only_MAD = mean(abs(rl_noaa_noaa - rl_obs)),
  stringsAsFactors = FALSE
)
tab2 <- do.call(rbind, tab2_rows)
rownames(tab2) <- NULL

# ══════════════════════════════════════════════════════════════
# TABLE 3: Mean SE by Region
# ══════════════════════════════════════════════════════════════
atlantic <- grid$lon > -82
gulf     <- grid$lon <= -82

tab3 <- data.frame(
  Region = c("Overall", "Atlantic", "Gulf"),
  Joint_mean_SE = c(mean(rl_joint$se_delta, na.rm = TRUE),
                    mean(rl_joint$se_delta[atlantic], na.rm = TRUE),
                    mean(rl_joint$se_delta[gulf], na.rm = TRUE)),
  NOAA_only_mean_SE = c(mean(rl_noaa$se_delta, na.rm = TRUE),
                         mean(rl_noaa$se_delta[atlantic], na.rm = TRUE),
                         mean(rl_noaa$se_delta[gulf], na.rm = TRUE)),
  ADCIRC_only_mean_SE = c(mean(rl_adcirc$se_delta, na.rm = TRUE),
                           mean(rl_adcirc$se_delta[atlantic], na.rm = TRUE),
                           mean(rl_adcirc$se_delta[gulf], na.rm = TRUE)),
  stringsAsFactors = FALSE
)

# ══════════════════════════════════════════════════════════════
# TABLE 4: SE Ratio Summary
# ══════════════════════════════════════════════════════════════
ratio <- rl_noaa$se_delta / rl_joint$se_delta
ratio <- ratio[is.finite(ratio)]

tab4 <- data.frame(
  Statistic = c("min", "Q25", "median", "mean", "Q75", "max", "pct_gt_1.5"),
  Value = c(min(ratio), quantile(ratio, 0.25), median(ratio), mean(ratio),
            quantile(ratio, 0.75), max(ratio), 100 * mean(ratio > 1.5)),
  stringsAsFactors = FALSE
)
rownames(tab4) <- NULL

# ── Print to console ──────────────────────────────────────────
cat("\n")
cat("Table 1: Model Fit Summary\n")
print(tab1, row.names = FALSE, right = FALSE)

cat("\nTable 2: Prediction Accuracy at 29 NOAA Sites\n")
print(tab2, row.names = FALSE, right = FALSE)

cat("\nTable 3: Mean SE of 100-Year Return Levels by Region\n")
print(tab3, row.names = FALSE, right = FALSE)

cat("\nTable 4: SE Ratio Summary -- SE(NOAA-only) / SE(joint)\n")
print(tab4, row.names = FALSE, right = FALSE)

# ── Save as formatted text ────────────────────────────────────
sink("data-raw/comparison_tables.txt")
cat("Table 1: Model Fit Summary\n")
cat("Model, Sites, Num_Parameters, NLL, rho\n")
for (i in seq_len(nrow(tab1))) {
  cat(sprintf("%s, %d, %d, %.4f, %s\n",
      tab1$Model[i], tab1$Sites[i], tab1$Num_Parameters[i], tab1$NLL[i], tab1$rho[i]))
}
cat("\nTable 2: Prediction Accuracy at 29 NOAA Sites (kriging vs Stage 1 MLE)\n")
cat("Parameter, Joint_RMSE, NOAA_only_RMSE, Joint_MAD, NOAA_only_MAD\n")
for (i in seq_len(nrow(tab2))) {
  cat(sprintf("%s, %.4f, %.4f, %.4f, %.4f\n",
      tab2$Parameter[i], tab2$Joint_RMSE[i], tab2$NOAA_only_RMSE[i],
      tab2$Joint_MAD[i], tab2$NOAA_only_MAD[i]))
}
cat("\nTable 3: Mean SE of 100-Year Return Levels by Region\n")
cat("Region, Joint_mean_SE, NOAA_only_mean_SE, ADCIRC_only_mean_SE\n")
for (i in seq_len(nrow(tab3))) {
  cat(sprintf("%s, %.4f, %.4f, %.4f\n",
      tab3$Region[i], tab3$Joint_mean_SE[i], tab3$NOAA_only_mean_SE[i],
      tab3$ADCIRC_only_mean_SE[i]))
}
cat("\nTable 4: SE Ratio Summary -- SE(NOAA-only) / SE(joint)\n")
cat("Statistic, Value\n")
for (i in seq_len(nrow(tab4))) {
  cat(sprintf("%s, %.4f\n", tab4$Statistic[i], tab4$Value[i]))
}
sink()

# ── Save as RDS ───────────────────────────────────────────────
tables <- list(
  model_fit = tab1,
  noaa_site_accuracy = tab2,
  se_by_region = tab3,
  se_ratio = tab4
)
saveRDS(tables, "data-raw/comparison_tables.rds")

cat("\nSaved data-raw/comparison_tables.txt\n")
cat("Saved data-raw/comparison_tables.rds\n")
