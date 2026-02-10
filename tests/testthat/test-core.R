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
