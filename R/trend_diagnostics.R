#' Mann-Kendall trend test with Sen's slope
#'
#' Tests for a monotonic trend in a time series using Kendall's tau
#' statistic. Implements the test directly from rank correlations
#' (no external package required). Sen's slope is the median of all
#' pairwise slopes.
#'
#' @param annual_maxima Numeric vector of annual maxima.
#' @param years Numeric vector of corresponding years (same length).
#' @return A list with components:
#'   \describe{
#'     \item{tau}{Kendall's tau statistic.}
#'     \item{S}{Mann-Kendall S statistic.}
#'     \item{p_value}{Two-sided p-value (normal approximation).}
#'     \item{sens_slope}{Sen's slope estimate (units per year).}
#'     \item{sens_intercept}{Sen's intercept estimate.}
#'   }
#' @export
mann_kendall_test <- function(annual_maxima, years) {
  # Remove NAs pairwise
  ok <- !is.na(annual_maxima) & !is.na(years)
  x <- annual_maxima[ok]
  t <- years[ok]
  n <- length(x)

  # Mann-Kendall S statistic: sum of sign(x[j] - x[i]) for all j > i
  S <- 0
  for (i in seq_len(n - 1)) {
    for (j in (i + 1):n) {
      S <- S + sign(x[j] - x[i])
    }
  }

  # Variance of S (no ties correction for simplicity)
  var_S <- n * (n - 1) * (2 * n + 5) / 18

  # Normal approximation with continuity correction
  if (S > 0) {
    Z <- (S - 1) / sqrt(var_S)
  } else if (S < 0) {
    Z <- (S + 1) / sqrt(var_S)
  } else {
    Z <- 0
  }

  p_value <- 2 * pnorm(-abs(Z))
  tau <- S / (n * (n - 1) / 2)

  # Sen's slope: median of all pairwise slopes
  slopes <- numeric(0)
  for (i in seq_len(n - 1)) {
    for (j in (i + 1):n) {
      if (t[j] != t[i]) {
        slopes <- c(slopes, (x[j] - x[i]) / (t[j] - t[i]))
      }
    }
  }
  sens_slope <- median(slopes)
  sens_intercept <- median(x - sens_slope * t)

  list(
    tau = tau,
    S = S,
    p_value = p_value,
    sens_slope = sens_slope,
    sens_intercept = sens_intercept
  )
}

#' GEV trend test via likelihood ratio
#'
#' Fits stationary GEV(mu, sigma, xi) and nonstationary
#' GEV(mu0 + mu1*year, sigma, xi) via \code{extRemes::fevd}, then
#' compares with a likelihood ratio test (chi-squared, df=1).
#'
#' @param annual_maxima Numeric vector of annual maxima.
#' @param years Numeric vector of corresponding years (same length).
#' @return A list with components:
#'   \describe{
#'     \item{mu0}{Intercept of nonstationary location.}
#'     \item{mu1}{Slope of nonstationary location (units per year).}
#'     \item{lrt_stat}{Likelihood ratio test statistic.}
#'     \item{lrt_pvalue}{P-value from chi-squared(df=1).}
#'     \item{aic_stat}{AIC of stationary model.}
#'     \item{aic_nonstat}{AIC of nonstationary model.}
#'     \item{nllh_stat}{Negative log-likelihood of stationary model.}
#'     \item{nllh_nonstat}{Negative log-likelihood of nonstationary model.}
#'   }
#' @export
gev_trend_test <- function(annual_maxima, years) {
  ok <- !is.na(annual_maxima) & !is.na(years)
  x <- annual_maxima[ok]
  yr <- years[ok]
  d <- data.frame(x = x, year = yr)

  fit0 <- extRemes::fevd(x, data = d, type = "GEV")
  fit1 <- extRemes::fevd(x, data = d, type = "GEV", location.fun = ~year)

  nllh0 <- fit0$results$value
  nllh1 <- fit1$results$value
  n_par0 <- length(fit0$results$par)
  n_par1 <- length(fit1$results$par)

  # LRT: deviance = -2 * (loglik_stat - loglik_nonstat)
  #     = -2 * (-nllh0 - (-nllh1)) = 2 * (nllh0 - nllh1)
  deviance <- 2 * (nllh0 - nllh1)
  lrt_pvalue <- pchisq(deviance, df = n_par1 - n_par0, lower.tail = FALSE)

  aic0 <- 2 * n_par0 + 2 * nllh0
  aic1 <- 2 * n_par1 + 2 * nllh1

  list(
    mu0 = unname(fit1$results$par["mu0"]),
    mu1 = unname(fit1$results$par["mu1"]),
    lrt_stat = deviance,
    lrt_pvalue = lrt_pvalue,
    aic_stat = aic0,
    aic_nonstat = aic1,
    nllh_stat = nllh0,
    nllh_nonstat = nllh1
  )
}
