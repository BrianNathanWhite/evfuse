#!/usr/bin/env Rscript
# Detrended pipeline: nonstationary NOAA GEV, then full Stage 2 + LOO-CV
# Compares stationary vs detrended results side-by-side
devtools::load_all()

# ── 1. Load data ─────────────────────────────────────────────
df <- read.csv("data-raw/combined_max.csv")
dat <- load_data(df)
D <- compute_distances(dat$sites)

# ── 2. Detrended Stage 1 ────────────────────────────────────
cat("\n== Stage 1: Detrended GEV (NOAA nonstationary, ADCIRC stationary) ==\n")
stage1_dt <- fit_gev_detrended(dat, df, ref_year = 2000)
cat(sprintf("Converged: %d / %d\n", sum(stage1_dt$converged), dat$n_sites))

# ── 3. Detrended Bootstrap W ────────────────────────────────
cat("\n== Bootstrap W (detrended, B=500) ==\n")
bs_dt <- bootstrap_W_detrended(dat, df, B = 500, ref_year = 2000, seed = 42)
cat(sprintf("Bootstrap failures: %d\n", bs_dt$n_failures))
W_tap_dt <- taper_W(bs_dt$W_bs, D, lambda = 300)

# ── 4. Stage 2: Joint model + multi-start ───────────────────
cat("\n== Stage 2: Joint 6-dim model (detrended) ==\n")
fit_dt <- fit_spatial_model(stage1_dt, dat, W_tap_dt, D,
                             control = list(maxit = 2000, trace = 1))

# Retry if needed
max_retries <- 5
attempt <- 1
while (fit_dt$optim_result$convergence != 0 && attempt <= max_retries) {
  cat(sprintf("Retry %d: perturbing rho...\n", attempt))
  set.seed(100 + attempt)
  start_retry <- list(beta = fit_dt$beta, A = fit_dt$A,
                      rho = fit_dt$rho * runif(6, 0.5, 2))
  fit_dt <- fit_spatial_model(stage1_dt, dat, W_tap_dt, D,
                               start = start_retry,
                               control = list(maxit = 2000, trace = 1))
  attempt <- attempt + 1
}

# Multi-start
set.seed(2026)
best_fit <- fit_dt
best_nll <- fit_dt$optim_result$value
cat(sprintf("\nInitial NLL: %.4f\n", best_nll))

for (i in seq_len(20)) {
  rho_i <- exp(runif(6, log(50), log(5000)))
  start_i <- list(beta = fit_dt$beta, A = fit_dt$A, rho = rho_i)
  fit_i <- tryCatch(
    suppressWarnings(
      fit_spatial_model(stage1_dt, dat, W_tap_dt, D,
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
cat(sprintf("\nBest NLL (detrended): %.4f\n", best_nll))
saveRDS(best_fit, "data-raw/model_6dim_detrended.rds")

# ── 5. LOO-CV (detrended) ───────────────────────────────────
cat("\n== LOO-CV (detrended joint model) ==\n")
loo_dt <- loo_cv(best_fit)
sum_dt <- loo_summary(loo_dt, r = 100)

# Also LOO-CV for detrended NOAA-only model
cat("Fitting NOAA-only detrended model for comparison...\n")
noaa_model_dt <- fit_naive_model(stage1_dt, dat, bs_dt$W_bs, D,
                                  source = "NOAA", lambda = 300, n_starts = 20,
                                  control = list(maxit = 2000, trace = 0))
saveRDS(noaa_model_dt, "data-raw/model_noaa_only_detrended.rds")
loo_noaa_dt <- loo_cv(noaa_model_dt)
sum_noaa_dt <- loo_summary(loo_noaa_dt, r = 100)

# ══════════════════════════════════════════════════════════════
# Load stationary results for comparison
# ══════════════════════════════════════════════════════════════
stat_fit <- readRDS("data-raw/model_6dim_best.rds")
loo_stat <- readRDS("data-raw/loo_cv_results.rds")
sum_stat_joint <- loo_stat$summary_joint
sum_stat_noaa  <- loo_stat$summary_noaa

# ── Cross-source correlations ────────────────────────────────
extract_cors <- function(model) {
  AAT <- model$A %*% t(model$A)
  sds <- sqrt(diag(AAT))
  cor_mat <- AAT / outer(sds, sds)
  c(mu    = cor_mat[1, 4],
    log_sigma = cor_mat[2, 5],
    xi    = cor_mat[3, 6])
}

cors_stat <- extract_cors(stat_fit)
cors_dt   <- extract_cors(best_fit)

# ══════════════════════════════════════════════════════════════
# Print comparison table
# ══════════════════════════════════════════════════════════════
cat("\n")
cat("════════════════════════════════════════════════════════════════════\n")
cat("Stationary vs Detrended Pipeline: Side-by-Side Comparison\n")
cat("════════════════════════════════════════════════════════════════════\n\n")

cat("── Cross-source correlations (from A matrix) ─────────────────────\n")
cat(sprintf("%-20s %12s %12s\n", "Parameter pair", "Stationary", "Detrended"))
cat(sprintf("%-20s %12s %12s\n", "--------------", "----------", "---------"))
cat(sprintf("%-20s %12.3f %12.3f\n", "Cor(mu_N, mu_A)", cors_stat["mu"], cors_dt["mu"]))
cat(sprintf("%-20s %12.3f %12.3f\n", "Cor(lsig_N, lsig_A)", cors_stat["log_sigma"], cors_dt["log_sigma"]))
cat(sprintf("%-20s %12.3f %12.3f\n", "Cor(xi_N, xi_A)", cors_stat["xi"], cors_dt["xi"]))
cat("\n")

cat("── Model fit ──────────────────────────────────────────────────────\n")
cat(sprintf("%-20s %12s %12s\n", "Metric", "Stationary", "Detrended"))
cat(sprintf("%-20s %12s %12s\n", "------", "----------", "---------"))
cat(sprintf("%-20s %12.4f %12.4f\n", "Joint NLL", stat_fit$optim_result$value, best_nll))
cat("\n")

cat("── LOO-CV at 29 NOAA sites (joint model) ─────────────────────────\n")
cat(sprintf("%-20s %12s %12s\n", "Metric", "Stationary", "Detrended"))
cat(sprintf("%-20s %12s %12s\n", "------", "----------", "---------"))
cat(sprintf("%-20s %12.4f %12.4f\n", "mu RMSE",
    sum_stat_joint$param_stats$rmse[1], sum_dt$param_stats$rmse[1]))
cat(sprintf("%-20s %12.4f %12.4f\n", "log_sigma RMSE",
    sum_stat_joint$param_stats$rmse[2], sum_dt$param_stats$rmse[2]))
cat(sprintf("%-20s %12.4f %12.4f\n", "xi RMSE",
    sum_stat_joint$param_stats$rmse[3], sum_dt$param_stats$rmse[3]))
cat(sprintf("%-20s %12.4f %12.4f\n", "100yr RL RMSE",
    sum_stat_joint$rl_rmse, sum_dt$rl_rmse))
cat(sprintf("%-20s %12.2f %12.2f\n", "Total LPD",
    sum_stat_joint$total_lpd, sum_dt$total_lpd))
cat(sprintf("%-20s %12.4f %12.4f\n", "Mean LPD",
    sum_stat_joint$mean_lpd, sum_dt$mean_lpd))
cat("\n")

cat("── LOO-CV: Joint vs NOAA-only (per pipeline) ─────────────────────\n")
cat(sprintf("%-20s %12s %12s\n", "Metric", "Stationary", "Detrended"))
cat(sprintf("%-20s %12s %12s\n", "------", "----------", "---------"))

n_stat_wins <- sum(loo_stat$joint$loo_lpd > loo_stat$noaa$loo_lpd)
n_dt_wins   <- sum(loo_dt$loo_lpd > loo_noaa_dt$loo_lpd)
cat(sprintf("%-20s %12s %12s\n", "Joint wins (LPD)",
    sprintf("%d / 29", n_stat_wins), sprintf("%d / 29", n_dt_wins)))

lpd_diff_stat <- sum_stat_joint$total_lpd - sum_stat_noaa$total_lpd
lpd_diff_dt   <- sum_dt$total_lpd - sum_noaa_dt$total_lpd
cat(sprintf("%-20s %12.2f %12.2f\n", "LPD gain (J - N)", lpd_diff_stat, lpd_diff_dt))

rl_pct_stat <- 100 * (1 - sum_stat_joint$rl_rmse / sum_stat_noaa$rl_rmse)
rl_pct_dt   <- 100 * (1 - sum_dt$rl_rmse / sum_noaa_dt$rl_rmse)
cat(sprintf("%-20s %11.1f%% %11.1f%%\n", "RL RMSE reduction", rl_pct_stat, rl_pct_dt))
cat("\n")

# ── Save all detrended results ──────────────────────────────
detrended_results <- list(
  stage1 = stage1_dt,
  model = best_fit,
  loo_joint = loo_dt,
  loo_noaa = loo_noaa_dt,
  summary_joint = sum_dt,
  summary_noaa = sum_noaa_dt,
  cors = cors_dt
)
saveRDS(detrended_results, "data-raw/detrended_results.rds")
cat("Saved data-raw/model_6dim_detrended.rds\n")
cat("Saved data-raw/model_noaa_only_detrended.rds\n")
cat("Saved data-raw/detrended_results.rds\n")
