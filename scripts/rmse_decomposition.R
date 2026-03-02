#!/usr/bin/env Rscript
# rmse_decomposition.R
#
# Decompose 100-yr return level LOO-CV improvement by parameter.
# For each GEV parameter (mu, log_sigma, xi), substitute only that
# parameter from the joint model into the NOAA-only prediction and
# measure the RMSE reduction. Also breaks down by geographic region.
#
# Outputs: printed summary (Table 3 in manuscript).

library(devtools)
load_all()

# ── Load fitted models ────────────────────────────────────────────────────────
cat("Loading models...\n")
dat <- evfuse::coast_data
D <- compute_distances(dat$sites)

model_joint <- readRDS("data-raw/model_6dim_ns.rds")
model_noaa  <- readRDS("data-raw/model_noaa_only_ns.rds")

# ── LOO-CV ────────────────────────────────────────────────────────────────────
loo_j <- loo_cv(model_joint)
loo_n <- loo_cv(model_noaa)

n_noaa <- nrow(loo_j$loo_mean)

# Compute return levels from observed, joint LOO, and NOAA-only LOO
rl_obs <- rl_joint <- rl_noaa <- numeric(n_noaa)
# Hybrid: replace only one param from joint, keep rest at NOAA-only
rl_hybrid_mu <- rl_hybrid_logsig <- rl_hybrid_xi <- numeric(n_noaa)

for (k in seq_len(n_noaa)) {
  mu_obs <- loo_j$observed[k, 1]
  ls_obs <- loo_j$observed[k, 2]
  xi_obs <- loo_j$observed[k, 3]

  mu_j <- loo_j$loo_mean[k, 1]
  ls_j <- loo_j$loo_mean[k, 2]
  xi_j <- loo_j$loo_mean[k, 3]

  mu_n <- loo_n$loo_mean[k, 1]
  ls_n <- loo_n$loo_mean[k, 2]
  xi_n <- loo_n$loo_mean[k, 3]

  rl_obs[k]   <- gev_return_level(mu_obs, exp(ls_obs), xi_obs, 100)
  rl_joint[k] <- gev_return_level(mu_j, exp(ls_j), xi_j, 100)
  rl_noaa[k]  <- gev_return_level(mu_n, exp(ls_n), xi_n, 100)

  # Hybrids: swap one param from joint into NOAA-only prediction
  rl_hybrid_mu[k]     <- gev_return_level(mu_j, exp(ls_n), xi_n, 100)
  rl_hybrid_logsig[k] <- gev_return_level(mu_n, exp(ls_j), xi_n, 100)
  rl_hybrid_xi[k]     <- gev_return_level(mu_n, exp(ls_n), xi_j, 100)
}

rmse <- function(x) sqrt(mean(x^2))

rmse_noaa  <- rmse(rl_noaa - rl_obs)
rmse_joint <- rmse(rl_joint - rl_obs)
rmse_hyb_mu <- rmse(rl_hybrid_mu - rl_obs)
rmse_hyb_ls <- rmse(rl_hybrid_logsig - rl_obs)
rmse_hyb_xi <- rmse(rl_hybrid_xi - rl_obs)

cat(sprintf("100-yr RL LOO-CV RMSE:\n"))
cat(sprintf("  NOAA-only:             %.3f m\n", rmse_noaa))
cat(sprintf("  Joint (all params):    %.3f m (%.1f%% reduction)\n",
            rmse_joint, 100 * (1 - rmse_joint / rmse_noaa)))
cat(sprintf("  Hybrid (mu from J):    %.3f m (%.1f%% reduction)\n",
            rmse_hyb_mu, 100 * (1 - rmse_hyb_mu / rmse_noaa)))
cat(sprintf("  Hybrid (logsig from J):%.3f m (%.1f%% reduction)\n",
            rmse_hyb_ls, 100 * (1 - rmse_hyb_ls / rmse_noaa)))
cat(sprintf("  Hybrid (xi from J):    %.3f m (%.1f%% reduction)\n",
            rmse_hyb_xi, 100 * (1 - rmse_hyb_xi / rmse_noaa)))
cat("\n")

# Regional decomposition
sites_loo <- loo_j$sites
regions <- ifelse(sites_loo$lon > -82 & sites_loo$lat < 31, "Gulf",
           ifelse(sites_loo$lat < 35.5, "SE_Atlantic",
           ifelse(sites_loo$lat < 41, "Mid_Atlantic", "New_England")))

cat("Regional 100-yr RL RMSE:\n")
cat(sprintf("%-15s %5s %8s %8s %8s\n", "Region", "n", "NOAA", "Joint", "Reduction"))
for (reg in c("Gulf", "SE_Atlantic", "Mid_Atlantic", "New_England")) {
  idx <- which(regions == reg)
  if (length(idx) == 0) next
  r_n <- rmse(rl_noaa[idx] - rl_obs[idx])
  r_j <- rmse(rl_joint[idx] - rl_obs[idx])
  cat(sprintf("%-15s %5d %8.3f %8.3f %7.1f%%\n",
              reg, length(idx), r_n, r_j, 100 * (1 - r_j / r_n)))
}
cat("\n")
