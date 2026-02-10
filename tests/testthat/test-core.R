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
