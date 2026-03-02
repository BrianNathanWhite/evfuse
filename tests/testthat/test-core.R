test_that("gev_return_level matches known values", {
  # For Gumbel (xi = 0): RL_r = mu - sigma * log(-log(1 - 1/r))
  mu <- 10
  sigma <- 2
  xi <- 0
  r <- 100

  rl <- gev_return_level(mu, sigma, xi, r)
  expected <- mu - sigma * log(-log(1 - 1/r))
  expect_equal(rl, expected, tolerance = 1e-6)
})

test_that("gev_return_level handles positive xi", {
  mu <- 10
  sigma <- 2
  xi <- 0.2
  r <- 100

  rl <- gev_return_level(mu, sigma, xi, r)
  # Manual: mu - (sigma/xi)(1 - (-log(1-1/r))^{-xi})
  p <- 1 / r
  yp <- -log(1 - p)
  expected <- mu - (sigma / xi) * (1 - yp^(-xi))
  expect_equal(rl, expected, tolerance = 1e-6)
})

test_that("gev_exceedance_prob is inverse of return level", {
  mu <- 10
  sigma <- 2
  xi <- 0.1
  r <- 50

  rl <- gev_return_level(mu, sigma, xi, r)
  prob <- gev_exceedance_prob(rl, mu, sigma, xi)
  expect_equal(prob, 1 / r, tolerance = 1e-6)
})

test_that("wendland2 has compact support", {
  expect_equal(wendland2(0), 1)
  expect_equal(wendland2(1), 0)
  expect_equal(wendland2(1.5), 0)
  expect_true(wendland2(0.5) > 0)
})

test_that("pack/unpack params are inverses", {
  p <- 6
  beta <- rnorm(p)
  A <- matrix(0, p, p)
  A[lower.tri(A, diag = TRUE)] <- rnorm(p * (p + 1) / 2)
  rho <- abs(rnorm(p)) + 1

  par <- pack_params(beta, A, rho, p)
  out <- unpack_params(par, p)

  expect_equal(out$beta, beta)
  expect_equal(out$A, A)
  expect_equal(out$rho, rho)
})

test_that("build_sigma is symmetric positive semi-definite", {
  p <- 6
  L <- 5
  D <- as.matrix(dist(matrix(rnorm(L * 2), L, 2)))
  A <- diag(0.5, p)
  A[lower.tri(A)] <- rnorm(p * (p - 1) / 2) * 0.1
  rho <- rep(1, p)

  Sigma <- build_sigma(A, rho, D, p)

  # Symmetric
  expect_equal(Sigma, t(Sigma), tolerance = 1e-10)

  # PSD: all eigenvalues >= 0
  eigs <- eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eigs >= -1e-10))
})

test_that("exp_cov_matrix produces valid correlation matrix", {
  D <- matrix(c(0, 1, 2, 1, 0, 1.5, 2, 1.5, 0), 3, 3)
  rho <- 1

  C <- exp_cov_matrix(D, rho)
  expect_equal(diag(C), rep(1, 3))
  expect_equal(C, t(C))
  expect_true(all(C >= 0 & C <= 1))
})

test_that("fit_gev_all returns valid output structure", {
  skip_if_not_installed("extRemes")

  # Create minimal test data (2 sites, 30 years each)
  set.seed(42)
  df <- data.frame(
    lon = rep(c(-80, -81), each = 30),
    lat = rep(c(25, 26), each = 30),
    location = rep(c("site1", "site2"), each = 30),
    year = rep(1990:2019, 2),
    max_sea_level = c(
      extRemes::revd(30, loc = 1, scale = 0.2, shape = 0.1, type = "GEV"),
      extRemes::revd(30, loc = 2, scale = 0.3, shape = 0.05, type = "GEV")
    ),
    data_source = rep(c("NOAA", "ADCIRC"), each = 30)
  )

  dat <- load_data(df)
  stage1 <- fit_gev_all(dat)

  # Structure checks

expect_s3_class(stage1, "evfuse_stage1")
  expect_equal(nrow(stage1$theta_hat), 2)
  expect_equal(ncol(stage1$theta_hat), 3)
  expect_equal(length(stage1$vcov_list), 2)
  expect_equal(length(stage1$fits), 2)
  expect_equal(length(stage1$converged), 2)
  expect_true(stage1$log_scale)

  # All should converge
  expect_true(all(stage1$converged))

  # Vcov should be 3x3 and positive definite
  for (V in stage1$vcov_list) {
    expect_equal(dim(V), c(3, 3))
    eigs <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
    expect_true(all(eigs > 0))
  }
})

test_that("fit_gev_all delta method transforms correctly", {
  skip_if_not_installed("extRemes")

  # Create test data
  set.seed(123)
  df <- data.frame(
    lon = -80, lat = 25, location = "test",
    year = 1990:2019,
    max_sea_level = extRemes::revd(30, loc = 1, scale = 0.5, shape = 0.1, type = "GEV"),
    data_source = "NOAA"
  )

  dat <- load_data(df)

  # Fit with and without log scale
  stage1_log <- fit_gev_all(dat, log_scale = TRUE)
  stage1_raw <- fit_gev_all(dat, log_scale = FALSE)

  # Check that log(sigma) transformation is correct
  expect_equal(
    unname(exp(stage1_log$theta_hat[1, "log_sigma"])),
    unname(stage1_raw$theta_hat[1, "sigma"]),
    tolerance = 1e-10
  )

  # mu and xi should be identical
  expect_equal(stage1_log$theta_hat[1, "mu"], stage1_raw$theta_hat[1, "mu"])
  expect_equal(stage1_log$theta_hat[1, "xi"], stage1_raw$theta_hat[1, "xi"])
})

test_that("stage 2 analytic gradient matches numerical finite differences", {
  skip_if_not_installed("extRemes")

  set.seed(42)
  n_noaa <- 3
  n_adcirc <- 4
  n_sites <- n_noaa + n_adcirc
  n_years <- 30

  df <- data.frame(
    lon = rep(seq(-80, -78, length.out = n_sites), each = n_years),
    lat = rep(seq(25, 27, length.out = n_sites), each = n_years),
    location = rep(paste0("site", 1:n_sites), each = n_years),
    year = rep(1990:(1990 + n_years - 1), n_sites),
    max_sea_level = c(
      unlist(lapply(1:n_noaa, function(i)
        extRemes::revd(n_years, loc = 1 + 0.1 * i, scale = 0.2,
                       shape = 0.1, type = "GEV")
      )),
      unlist(lapply(1:n_adcirc, function(i)
        extRemes::revd(n_years, loc = 2 + 0.1 * i, scale = 0.3,
                       shape = 0.05, type = "GEV")
      ))
    ),
    data_source = rep(c(rep("NOAA", n_noaa), rep("ADCIRC", n_adcirc)),
                      each = n_years)
  )

  dat <- load_data(df)
  stage1 <- fit_gev_all(dat)
  D <- compute_distances(dat$sites)

  # Simple W_tap: small diagonal in observed space, then embed
  W_tap <- embed_W(diag(0.01, n_sites * 3), dat)

  result <- suppressWarnings(
    fit_spatial_model(stage1, dat, W_tap, D,
                      check_gradient = TRUE,
                      control = list(maxit = 1, trace = 0))
  )

  expect_true(!is.null(result$grad_check))
  expect_lt(result$grad_check$max_rel_err, 1e-4)
})

test_that("rl_gradient finite difference check", {
  mu <- 5
  log_sigma <- log(2)
  xi <- 0.15
  r <- 100
  eps <- 1e-5

  grad <- rl_gradient(mu, log_sigma, xi, r)

  # Numerical gradient
  f <- function(par) gev_return_level(par[1], exp(par[2]), par[3], r)
  par0 <- c(mu, log_sigma, xi)
  num_grad <- vapply(1:3, function(j) {
    par_up <- par0; par_up[j] <- par_up[j] + eps
    par_dn <- par0; par_dn[j] <- par_dn[j] - eps
    (f(par_up) - f(par_dn)) / (2 * eps)
  }, numeric(1))

  expect_equal(grad, num_grad, tolerance = 1e-4)
})

# ============================================================================
# Nonstationary covariate functions
# ============================================================================

# Helper: create test data with a covariate
make_ns_test_data <- function(n_noaa = 3, n_adcirc = 4, n_years = 30,
                               seed = 42) {
  set.seed(seed)
  n_sites <- n_noaa + n_adcirc

  # SST-like covariate: shared across NOAA sites, with a trend
  sst_vals <- 25 + 0.05 * (1:n_years) + rnorm(n_years, 0, 0.3)

  rows <- list()
  for (i in seq_len(n_sites)) {
    is_noaa <- i <= n_noaa
    src <- if (is_noaa) "NOAA" else "ADCIRC"
    loc_name <- paste0("site", i)
    lon <- -80 + (i - 1) * 0.5
    lat <- 25 + (i - 1) * 0.3

    if (is_noaa) {
      # Nonstationary: mu depends on SST
      mu0 <- 1 + 0.1 * i
      mu1 <- 0.15  # location sensitivity to SST
      sigma <- 0.2
      xi <- 0.1
      mu_t <- mu0 + mu1 * (sst_vals - mean(sst_vals))
      maxima <- mu_t + sigma * extRemes::revd(n_years, loc = 0, scale = 1,
                                                shape = xi, type = "GEV")
    } else {
      # ADCIRC: stationary
      mu0 <- 2 + 0.1 * (i - n_noaa)
      sigma <- 0.3
      xi <- 0.05
      maxima <- extRemes::revd(n_years, loc = mu0, scale = sigma,
                                shape = xi, type = "GEV")
    }

    rows[[i]] <- data.frame(
      lon = lon, lat = lat, location = loc_name,
      year = 1990:(1990 + n_years - 1),
      max_sea_level = maxima,
      data_source = src,
      sst_warm = sst_vals,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

test_that("fit_gev_ns returns valid output structure", {
  skip_if_not_installed("extRemes")

  df <- make_ns_test_data()
  dat <- load_data(df, source_params = list(NOAA = 1:4, ADCIRC = 5:7))
  stage1 <- fit_gev_ns(dat, df, covariate = "sst_warm")

  # Class and dimensions
  expect_s3_class(stage1, "evfuse_stage1")
  expect_equal(nrow(stage1$theta_hat), 7)
  expect_equal(ncol(stage1$theta_hat), 4)
  expect_equal(colnames(stage1$theta_hat), c("mu0", "mu1", "log_sigma", "xi"))

  # Metadata
  expect_true(stage1$nonstationary)
  expect_equal(stage1$covariate, "sst_warm")
  expect_true(is.numeric(stage1$ref_value))

  # All should converge
  expect_true(all(stage1$converged))
})

test_that("fit_gev_ns: NOAA has 4x4 vcov, ADCIRC has 3x3 vcov", {
  skip_if_not_installed("extRemes")

  df <- make_ns_test_data()
  dat <- load_data(df, source_params = list(NOAA = 1:4, ADCIRC = 5:7))
  stage1 <- fit_gev_ns(dat, df, covariate = "sst_warm")

  noaa_idx <- which(dat$sites$data_source == "NOAA")
  adcirc_idx <- which(dat$sites$data_source == "ADCIRC")

  # NOAA: 4 params, 4x4 vcov
  for (i in noaa_idx) {
    expect_equal(dim(stage1$vcov_list[[i]]), c(4, 4))
    expect_false(any(is.na(stage1$theta_hat[i, ])))
    eigs <- eigen(stage1$vcov_list[[i]], symmetric = TRUE, only.values = TRUE)$values
    expect_true(all(eigs > 0))
  }

  # ADCIRC: 3 params in cols 1-3, col 4 = NA, 3x3 vcov
  for (i in adcirc_idx) {
    expect_equal(dim(stage1$vcov_list[[i]]), c(3, 3))
    expect_false(any(is.na(stage1$theta_hat[i, 1:3])))
    expect_true(is.na(stage1$theta_hat[i, 4]))
  }
})

test_that("fit_gev_ns: custom ref_value is respected", {
  skip_if_not_installed("extRemes")

  df <- make_ns_test_data()
  dat <- load_data(df, source_params = list(NOAA = 1:4, ADCIRC = 5:7))

  stage1_a <- fit_gev_ns(dat, df, covariate = "sst_warm", ref_value = 25.0)
  stage1_b <- fit_gev_ns(dat, df, covariate = "sst_warm", ref_value = 26.0)

  expect_equal(stage1_a$ref_value, 25.0)
  expect_equal(stage1_b$ref_value, 26.0)

  # Different ref_value -> different mu0 but same mu1
  noaa_1 <- which(dat$sites$data_source == "NOAA")[1]
  expect_false(isTRUE(all.equal(
    stage1_a$theta_hat[noaa_1, "mu0"],
    stage1_b$theta_hat[noaa_1, "mu0"]
  )))
  expect_equal(
    stage1_a$theta_hat[noaa_1, "mu1"],
    stage1_b$theta_hat[noaa_1, "mu1"],
    tolerance = 1e-4
  )
})

test_that("build_observation_structure handles asymmetric theta_hat", {
  skip_if_not_installed("extRemes")

  df <- make_ns_test_data()
  dat <- load_data(df, source_params = list(NOAA = 1:4, ADCIRC = 5:7))
  stage1 <- fit_gev_ns(dat, df, covariate = "sst_warm")

  obs <- build_observation_structure(stage1, dat)

  noaa_idx <- which(dat$sites$data_source == "NOAA")
  adcirc_idx <- which(dat$sites$data_source == "ADCIRC")

  # Expected observation count: NOAA * 4 + ADCIRC * 3
  expected_n <- length(noaa_idx) * 4 + length(adcirc_idx) * 3
  expect_equal(length(obs$obs_idx), expected_n)
  expect_equal(length(obs$theta_obs), expected_n)

  # No NAs in theta_obs
  expect_false(any(is.na(obs$theta_obs)))
})

test_that("embed_W handles heterogeneous observation dimensions", {
  set.seed(1)
  L <- 5
  p_bs <- 4
  p_full <- 7

  # Simulate a 4L x 4L W_tap with NaN for "ADCIRC" mu1 (position 2)
  W_tap <- crossprod(matrix(rnorm(p_bs * L * p_bs * L), p_bs * L))

  # Sites: first 2 are "NOAA", last 3 are "ADCIRC"
  sites_df <- data.frame(
    location = paste0("s", 1:L),
    data_source = c(rep("NOAA", 2), rep("ADCIRC", 3)),
    lon = seq(-80, -78, length.out = L),
    lat = seq(25, 27, length.out = L)
  )
  dat <- list(
    n_sites = L,
    sites = sites_df,
    source_params = list(NOAA = 1:4, ADCIRC = 5:7),
    p = 7
  )

  # Set ADCIRC mu1 entries to NA (as cov() produces for all-NA columns)
  adcirc_sites <- 3:5
  mu1_pos <- L + adcirc_sites
  W_tap[mu1_pos, ] <- NA
  W_tap[, mu1_pos] <- NA

  W_full <- embed_W(W_tap, dat)

  # Correct dimensions
  expect_equal(dim(W_full), c(L * p_full, L * p_full))

  # NOAA site 1: params 1-4 map to full positions 1-4
  noaa_1 <- 1
  expect_equal(W_full[noaa_1, noaa_1], W_tap[noaa_1, noaa_1])
  expect_equal(W_full[L + noaa_1, L + noaa_1], W_tap[L + noaa_1, L + noaa_1])

  # ADCIRC site 3: mu0 (W_tap pos 3) -> full param 5 (pos 4*L+3)
  adcirc_1 <- 3
  expect_equal(W_full[4 * L + adcirc_1, 4 * L + adcirc_1],
               W_tap[adcirc_1, adcirc_1])
  # ADCIRC log_sigma (W_tap pos 2*L+3) -> full param 6 (pos 5*L+3)
  expect_equal(W_full[5 * L + adcirc_1, 5 * L + adcirc_1],
               W_tap[2 * L + adcirc_1, 2 * L + adcirc_1])
  # ADCIRC xi (W_tap pos 3*L+3) -> full param 7 (pos 6*L+3)
  expect_equal(W_full[6 * L + adcirc_1, 6 * L + adcirc_1],
               W_tap[3 * L + adcirc_1, 3 * L + adcirc_1])

  # ADCIRC mu1 in full space (pos L+3) should be 0 (not mapped)
  expect_equal(W_full[L + adcirc_1, L + adcirc_1], 0)

  # Cross-covariance: NOAA mu0 x ADCIRC mu0
  expect_equal(W_full[noaa_1, 4 * L + adcirc_1],
               W_tap[noaa_1, adcirc_1])
})

test_that("embed_W backward compatible with uniform p_obs", {
  set.seed(2)
  L <- 4

  sites_df <- data.frame(
    location = paste0("s", 1:L),
    data_source = c(rep("NOAA", 2), rep("ADCIRC", 2)),
    lon = seq(-80, -78, length.out = L),
    lat = seq(25, 27, length.out = L)
  )
  dat <- list(
    n_sites = L,
    sites = sites_df,
    source_params = list(NOAA = 1:3, ADCIRC = 4:6),
    p = 6
  )

  W_tap <- crossprod(matrix(rnorm(3 * L * 3 * L), 3 * L))
  W_full <- embed_W(W_tap, dat)

  expect_equal(dim(W_full), c(6 * L, 6 * L))

  # NOAA site 1, param 1 -> full position 1
  expect_equal(W_full[1, 1], W_tap[1, 1])
  # ADCIRC site 3, param 1 -> full position 3*L+3
  expect_equal(W_full[3 * L + 3, 3 * L + 3], W_tap[3, 3])
})

test_that("bootstrap_W_ns produces correct structure", {
  skip_if_not_installed("extRemes")

  df <- make_ns_test_data(n_noaa = 2, n_adcirc = 3, n_years = 30)
  dat <- load_data(df, source_params = list(NOAA = 1:4, ADCIRC = 5:7))
  L <- dat$n_sites

  bs <- suppressWarnings(
    bootstrap_W_ns(dat, df, covariate = "sst_warm", B = 20, seed = 42)
  )

  # Gamma dimensions
  expect_equal(dim(bs$Gamma), c(20, 4 * L))

  # W_bs dimensions
  expect_equal(dim(bs$W_bs), c(4 * L, 4 * L))

  # ADCIRC mu1 columns in Gamma should be all NA
  adcirc_idx <- which(dat$sites$data_source == "ADCIRC")
  mu1_cols <- L + adcirc_idx
  expect_true(all(is.na(bs$Gamma[, mu1_cols])))

  # ADCIRC mu1 diagonal in W_bs should be NA
  expect_true(all(is.na(diag(bs$W_bs)[mu1_cols])))

  # Non-NA entries should be finite
  non_na <- !is.na(bs$W_bs)
  expect_true(all(is.finite(bs$W_bs[non_na])))
})

test_that("compute_return_levels_ns delta method gradient is correct", {
  # Test 4-dim gradient against numerical finite differences
  mu0 <- 1.5
  mu1 <- 0.15
  log_sigma <- log(0.3)
  xi <- 0.1
  sst <- 2.0
  r <- 100
  eps <- 1e-5

  # Manual 4-dim gradient
  mu <- mu0 + mu1 * sst
  grad3 <- rl_gradient(mu, log_sigma, xi, r)
  grad4 <- c(grad3[1], grad3[1] * sst, grad3[2], grad3[3])

  # Numerical gradient w.r.t. (mu0, mu1, log_sigma, xi)
  f <- function(par) {
    m <- par[1] + par[2] * sst
    gev_return_level(m, exp(par[3]), par[4], r)
  }
  par0 <- c(mu0, mu1, log_sigma, xi)
  num_grad <- vapply(1:4, function(j) {
    p_up <- par0; p_up[j] <- p_up[j] + eps
    p_dn <- par0; p_dn[j] <- p_dn[j] - eps
    (f(p_up) - f(p_dn)) / (2 * eps)
  }, numeric(1))

  expect_equal(grad4, num_grad, tolerance = 1e-4)
})

test_that("compute_return_levels_ns produces valid output", {
  # Create a mock predictions object
  n_new <- 3
  predictions <- list(
    noaa_mean = matrix(c(
      1.0, 0.10, log(0.3), 0.10,
      1.5, 0.15, log(0.4), 0.05,
      2.0, 0.20, log(0.2), 0.15
    ), nrow = n_new, byrow = TRUE),
    noaa_cov = lapply(1:n_new, function(k) {
      V <- diag(c(0.01, 0.001, 0.01, 0.005))
      V[1, 2] <- V[2, 1] <- 0.001
      V
    }),
    new_sites = data.frame(lon = c(-80, -81, -82), lat = c(25, 26, 27))
  )

  rl <- compute_return_levels_ns(predictions, covariate_value = 2.0,
                                   r = 100, seed = 123)

  # Output structure
  expect_equal(nrow(rl), n_new)
  expect_true(all(c("return_level", "se_delta", "ci_lower_delta",
                     "ci_upper_delta") %in% names(rl)))

  # Return levels should be finite and positive (sea levels)
  expect_true(all(is.finite(rl$return_level)))
  expect_true(all(rl$se_delta > 0))

  # CI should bracket the point estimate
  expect_true(all(rl$ci_lower_delta < rl$return_level))
  expect_true(all(rl$ci_upper_delta > rl$return_level))
})

test_that("compute_return_levels_ns: scalar vs vector covariate", {
  predictions <- list(
    noaa_mean = matrix(c(1.0, 0.1, log(0.3), 0.1,
                          1.5, 0.2, log(0.4), 0.05), nrow = 2, byrow = TRUE),
    noaa_cov = lapply(1:2, function(k) diag(0.01, 4)),
    new_sites = data.frame(lon = c(-80, -81), lat = c(25, 26))
  )

  # Scalar: same SST for both sites
  rl_scalar <- compute_return_levels_ns(predictions, covariate_value = 1.5,
                                          r = 100, method = "delta")
  # Vector: different SST per site
  rl_vector <- compute_return_levels_ns(predictions, covariate_value = c(1.5, 1.5),
                                          r = 100, method = "delta")

  expect_equal(rl_scalar$return_level, rl_vector$return_level)
  expect_equal(rl_scalar$se_delta, rl_vector$se_delta)
})

test_that("full nonstationary pipeline: Stage 1 through spatial model", {
  skip_if_not_installed("extRemes")

  df <- make_ns_test_data(n_noaa = 3, n_adcirc = 4, n_years = 30)
  dat <- load_data(df, source_params = list(NOAA = 1:4, ADCIRC = 5:7))
  D <- compute_distances(dat$sites)

  # Stage 1
  stage1 <- fit_gev_ns(dat, df, covariate = "sst_warm")
  expect_true(all(stage1$converged))

  # Build W_tap (use vcov-based diagonal for speed instead of bootstrap)
  L <- dat$n_sites
  W_diag <- diag(0.01, 4 * L)
  # Set ADCIRC mu1 to NA to simulate bootstrap output
  adcirc_idx <- which(dat$sites$data_source == "ADCIRC")
  mu1_pos <- L + adcirc_idx
  W_diag[mu1_pos, ] <- NA
  W_diag[, mu1_pos] <- NA

  # Embed and fit
  W_full <- embed_W(W_diag, dat)
  expect_equal(dim(W_full), c(7 * L, 7 * L))

  # Stage 2: should run with p=7
  result <- suppressWarnings(
    fit_spatial_model(stage1, dat, W_full, D,
                      control = list(maxit = 5, trace = 0))
  )

  expect_s3_class(result, "evfuse_model")
  expect_equal(result$p, 7)
  expect_equal(length(result$beta), 7)
  expect_equal(dim(result$A), c(7, 7))
  expect_equal(length(result$rho), 7)
  expect_equal(length(result$optim_result$par), 42)  # 7 + 28 + 7
})
