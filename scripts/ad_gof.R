#!/usr/bin/env Rscript
# Anderson-Darling goodness-of-fit for the Stage 1 GEV fits (Section S2).
#
# The A^2 statistic is applied to the PIT of each site's annual maxima under
# the fitted GEV, matching the Stage 1 model: NONSTATIONARY at NOAA sites
# (mu(t) = mu0 + mu1*(year - ref_year)), stationary at ADCIRC sites. Because
# the GEV parameters are estimated, the standard "parameters known" null does
# not apply; the p-value is obtained by parametric bootstrap (simulate from the
# fitted GEV, refit, recompute A^2), which is calibrated for estimated
# parameters (Stephens 1977; Stute et al. 1993).
#
# Usage: Rscript scripts/ad_gof.R   (~3 min, uses bundled coast_data)

devtools::load_all()

B        <- 2999
SEED0    <- 20260601      # site i uses SEED0 + i (deterministic, core-independent)
N_CORES  <- max(1L, parallel::detectCores() - 2L)
ALPHA    <- 0.05
REF_YEAR <- 2000

data(coast_data, package = "evfuse")
dat    <- coast_data
df     <- coast_data$raw_df
stage1 <- fit_gev_detrended(dat, df, ref_year = REF_YEAR)

gev_pit <- function(x, mu, sigma, xi) {        # mu may be a vector (per year)
  if (abs(xi) < 1e-8) {
    u <- exp(-exp(-(x - mu) / sigma))
  } else {
    t_val <- pmax(1 + xi * (x - mu) / sigma, 1e-15)
    u <- exp(-t_val^(-1 / xi))
  }
  pmax(pmin(u, 1 - 1e-15), 1e-15)
}
ad_A2 <- function(u) {
  u <- sort(pmax(pmin(u, 1 - 1e-15), 1e-15)); n <- length(u); i <- seq_len(n)
  -n - (1 / n) * sum((2 * i - 1) * (log(u) + log(1 - u[n + 1 - i])))
}
gev_sim <- function(mu, sigma, xi) {           # one draw per element of mu
  u <- runif(length(mu))
  if (abs(xi) < 1e-8) mu - sigma * log(-log(u))
  else                mu + sigma * ((-log(u))^(-xi) - 1) / xi
}

one_site <- function(i) {
  if (is.null(stage1$fits[[i]]) || !stage1$converged[i]) return(NULL)
  loc     <- dat$sites$location[i]
  is_noaa <- dat$sites$data_source[i] == "NOAA"
  par     <- stage1$fits[[i]]$results$par
  set.seed(SEED0 + i)

  if (is_noaa) {
    sdf <- df[as.character(df$location) == loc, ]
    sdf <- sdf[order(sdf$year), ]
    yc  <- sdf$year - REF_YEAR
    x   <- sdf$max_sea_level
    mu_t <- par["mu0"] + par["mu1"] * yc
    A2_obs <- ad_A2(gev_pit(x, mu_t, par["scale"], par["shape"]))
    A2_boot <- numeric(B); ok <- logical(B)
    for (b in seq_len(B)) {
      xb  <- gev_sim(mu_t, par["scale"], par["shape"])
      bdf <- data.frame(max_sea_level = xb, year_c = yc)
      fit <- tryCatch(suppressWarnings(extRemes::fevd(bdf$max_sea_level, data = bdf,
                       type = "GEV", location.fun = ~year_c)), error = function(e) NULL)
      if (is.null(fit) || fit$results$convergence != 0) { ok[b] <- FALSE; next }
      pb   <- fit$results$par
      mu_b <- pb["mu0"] + pb["mu1"] * yc
      A2_boot[b] <- ad_A2(gev_pit(xb, mu_b, pb["scale"], pb["shape"])); ok[b] <- TRUE
    }
  } else {
    x  <- dat$maxima[[i]]
    A2_obs <- ad_A2(gev_pit(x, par["location"], par["scale"], par["shape"]))
    A2_boot <- numeric(B); ok <- logical(B)
    for (b in seq_len(B)) {
      xb  <- gev_sim(rep(par["location"], length(x)), par["scale"], par["shape"])
      fit <- tryCatch(suppressWarnings(extRemes::fevd(xb, type = "GEV", method = "MLE")),
                      error = function(e) NULL)
      if (is.null(fit) || fit$results$convergence != 0) { ok[b] <- FALSE; next }
      pb <- fit$results$par
      A2_boot[b] <- ad_A2(gev_pit(xb, pb["location"], pb["scale"], pb["shape"])); ok[b] <- TRUE
    }
  }
  nb <- sum(ok)
  data.frame(site = loc, source = dat$sites$data_source[i], n = length(x),
             xi = round(as.numeric(par["shape"]), 3), A2 = round(A2_obs, 4),
             p_boot = round((1 + sum(A2_boot[ok] >= A2_obs)) / (nb + 1), 4),
             stringsAsFactors = FALSE)
}

cat(sprintf("Parametric-bootstrap AD goodness-of-fit: %d sites, B=%d, cores=%d\n",
            dat$n_sites, B, N_CORES))
res  <- parallel::mclapply(seq_len(dat$n_sites), one_site, mc.cores = N_CORES)
res  <- do.call(rbind, res[!vapply(res, is.null, logical(1))])
res  <- res[order(res$p_boot), ]
bonf <- ALPHA / nrow(res)

cat(sprintf("\nReject GEV null @ a=%.2f: %d / %d (expected by chance %.1f)\n",
            ALPHA, sum(res$p_boot < ALPHA), nrow(res), ALPHA * nrow(res)))
cat(sprintf("Reject after Bonferroni (%.2g): %d / %d; smallest p = %.4f\n",
            bonf, sum(res$p_boot < bonf), nrow(res), min(res$p_boot)))
cat(sprintf("By source @ a=%.2f: NOAA %d/%d | ADCIRC %d/%d\n", ALPHA,
            sum(res$p_boot[res$source == "NOAA"] < ALPHA), sum(res$source == "NOAA"),
            sum(res$p_boot[res$source == "ADCIRC"] < ALPHA), sum(res$source == "ADCIRC")))
