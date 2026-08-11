# M4 check: bootstrap W with calendar years resampled JOINTLY across sources
# (headline analysis resamples years independently across sources, zeroing the
# cross-source block of W in expectation). If shared-storm Stage-1 errors are
# being absorbed into the signal covariance, the cross-source correlations
# should drop under this scheme.
devtools::load_all(quiet = TRUE)

data(coast_data)
dat <- coast_data
df  <- coast_data$raw_df
D   <- compute_distances(dat$sites)
L   <- dat$n_sites
ref_year <- 2000
B <- 500

model_joint <- readRDS("data-raw/model_6dim_ns.rds")

cat("Stage 1 (deterministic refit)...\n")
stage1 <- fit_gev_detrended(dat, df, ref_year = ref_year)

# ── Shared-years bootstrap (mirrors bootstrap_W_detrended except the draw) ──
set.seed(42)
noaa_idx <- which(dat$sites$data_source == "NOAA")
site_years <- vector("list", L); site_maxima <- vector("list", L)
for (i in seq_len(L)) {
  loc <- dat$sites$location[i]
  site_df <- df[as.character(df$location) == loc, ]
  site_df <- site_df[order(site_df$year), ]
  site_years[[i]] <- site_df$year
  site_maxima[[i]] <- site_df$max_sea_level
}
years_all <- 1979:2021
Lp <- L * 3
Gamma <- matrix(NA_real_, nrow = B, ncol = Lp)
n_failures <- 0

cat(sprintf("Shared-years bootstrap: %d replicates...\n", B))
for (b in seq_len(B)) {
  yrs_b <- sample(years_all, length(years_all), replace = TRUE)
  theta_b <- rep(NA_real_, Lp)
  for (i in seq_len(L)) {
    m <- match(yrs_b, site_years[[i]]); m <- m[!is.na(m)]
    if (length(m) < 20) { n_failures <- n_failures + 1; next }
    x_boot <- site_maxima[[i]][m]
    if (dat$sites$data_source[i] == "NOAA") {
      d_boot <- data.frame(x = x_boot, year_c = site_years[[i]][m] - ref_year)
      fit <- tryCatch(extRemes::fevd(d_boot$x, data = d_boot, type = "GEV",
                                     location.fun = ~year_c),
                      error = function(e) NULL)
      if (is.null(fit) || fit$results$convergence != 0) { n_failures <- n_failures + 1; next }
      pars <- fit$results$par
      theta_i <- c(pars["mu0"], log(pars["scale"]), pars["shape"])
    } else {
      fit <- tryCatch(extRemes::fevd(x_boot, type = "GEV", method = "MLE"),
                      error = function(e) NULL)
      if (is.null(fit) || fit$results$convergence != 0) { n_failures <- n_failures + 1; next }
      pars <- fit$results$par
      theta_i <- c(pars["location"], log(pars["scale"]), pars["shape"])
    }
    for (j in 1:3) theta_b[(j - 1) * L + i] <- theta_i[j]
  }
  Gamma[b, ] <- theta_b
  if (b %% 100 == 0) cat(sprintf("  %d / %d\n", b, B))
}
W_bs_shared <- cov(Gamma, use = "pairwise.complete.obs")
cat(sprintf("Failures: %d / %d fits\n", n_failures, B * L))

W_tap_shared <- taper_W(W_bs_shared, D, lambda = 300)

# ── Cross-source block magnitude: shared vs headline ──
src_of <- dat$sites$data_source
idx_meta <- expand.grid(site = seq_len(L), param = 1:3)
col_src <- src_of[idx_meta$site]
cross_mask <- outer(col_src, col_src, function(a, b) a != b)
blk <- function(W) { v <- W[cross_mask]; c(mean_abs = mean(abs(v)), max_abs = max(abs(v)), frob = sqrt(sum(v^2))) }
cat("\nCross-source block of tapered W (shared-years scheme):\n"); print(round(blk(W_tap_shared), 5))
# headline: embedded 774 W_tap from the fitted model, restricted to observed comps
obs <- evfuse:::build_observation_structure(stage1, dat)
W_head_obs <- model_joint$W_tap[obs$obs_idx, obs$obs_idx]
op <- ((obs$obs_idx - 1) %/% L) + 1; os <- ((obs$obs_idx - 1) %% L) + 1
cross_mask_h <- outer(src_of[os], src_of[os], function(a, b) a != b)
v <- W_head_obs[cross_mask_h]
cat("Cross-source block of tapered W (headline, independent scheme):\n")
print(round(c(mean_abs = mean(abs(v)), max_abs = max(abs(v)), frob = sqrt(sum(v^2))), 5))

# ── Stage 2 refit under shared-years W: headline optimum start + 7 random ──
fits <- list()
fits[[1]] <- tryCatch(suppressWarnings(fit_spatial_model(
  stage1, dat, W_tap_shared, D,
  start = list(beta = model_joint$beta, A = model_joint$A, rho = model_joint$rho),
  control = list(maxit = 2000, trace = 0))), error = function(e) NULL)
set.seed(2026)
for (s in 2:8) {
  st <- list(beta = model_joint$beta, A = model_joint$A,
             rho = exp(runif(6, log(50), log(5000))))
  fits[[s]] <- tryCatch(suppressWarnings(fit_spatial_model(
    stage1, dat, W_tap_shared, D, start = st,
    control = list(maxit = 2000, trace = 0))), error = function(e) NULL)
}
nlls <- vapply(fits, function(f) if (is.null(f)) Inf else f$optim_result$value, 0)
cat("\nMulti-start NLLs (shared-years W):", round(nlls, 3), "\n")
best <- fits[[which.min(nlls)]]

cors <- function(A) { S <- A %*% t(A)
  c(mu = S[1,4]/sqrt(S[1,1]*S[4,4]), logsig = S[2,5]/sqrt(S[2,2]*S[5,5]),
    xi = S[3,6]/sqrt(S[3,3]*S[6,6])) }
cat("\nCross-source correlations, shared-years W: ", round(cors(best$A), 3), "\n")
cat("Cross-source correlations, headline:        ", round(cors(model_joint$A), 3), "\n")

# ── LOO comparison under shared-years W ──
loo_j <- loo_cv(best); sum_j <- loo_summary(loo_j, r = 100)
noaa_shared <- fit_naive_model(stage1, dat, W_bs_shared, D, source = "NOAA",
                               lambda = 300, n_starts = 10,
                               control = list(maxit = 2000, trace = 0))
loo_n <- loo_cv(noaa_shared); sum_n <- loo_summary(loo_n, r = 100)
cat(sprintf("\nLOO RL RMSE joint=%.4f noaa=%.4f reduction=%.1f%% (headline: 0.4670 / 0.7245 / 35.5%%)\n",
            sum_j$rl_rmse, sum_n$rl_rmse, 100 * (1 - sum_j$rl_rmse / sum_n$rl_rmse)))
saveRDS(list(W_tap_shared = W_tap_shared, best = best, cors_shared = cors(best$A),
             nlls = nlls, sum_j = sum_j, sum_n = sum_n, n_failures = n_failures),
        "data-raw/shared_years_W_results.rds")
cat("DONE m4\n")
