#!/usr/bin/env Rscript
# ══════════════════════════════════════════════════════════════════════════════
# Gradient Benchmark: Analytic vs Numerical (Central Finite Differences)
# ══════════════════════════════════════════════════════════════════════════════
# Controlled comparison of optim() with analytic gradient vs numerical gradient
# on the Stage 2 likelihood. Both methods start from the same perturbed point
# near the production optimum.
#
# Output: tables/gradient_benchmark.txt

devtools::load_all()

# ══════════════════════════════════════════════════════════════════════════════
# 1. Load production model and reconstruct objective function
# ══════════════════════════════════════════════════════════════════════════════

model <- readRDS("data-raw/model_6dim_ns.rds")
prod_nll <- model$optim_result$value
cat(sprintf("Production model NLL: %.4f\n", prod_nll))

L <- model$dat$n_sites
p <- model$p
D <- model$D
W_tap <- model$W_tap  # already in full 774x774 space

# Observation structure
obs <- build_observation_structure(model$stage1, model$dat)
n_obs <- length(obs$obs_idx)
obs_param <- ((obs$obs_idx - 1) %/% L) + 1
obs_site  <- ((obs$obs_idx - 1) %% L) + 1
param_groups <- lapply(1:p, function(j) which(obs_param == j))
param_sites  <- lapply(1:p, function(j) obs_site[param_groups[[j]]])

W_tap_obs <- W_tap[obs$obs_idx, obs$obs_idx]

cat(sprintf("Problem size: L=%d sites, p=%d params, n_obs=%d, n_par=33\n",
            L, p, n_obs))

# ══════════════════════════════════════════════════════════════════════════════
# 2. Construct nll / nll_grad closures (replicated from fit_spatial_model)
# ══════════════════════════════════════════════════════════════════════════════

cache <- new.env(parent = emptyenv())
cache$par <- NULL

compute_core <- function(par) {
  if (!is.null(cache$par) && identical(par, cache$par)) {
    return(cache$result)
  }

  params <- unpack_params(par, p = p)
  A <- params$A
  rho <- params$rho
  beta <- params$beta

  Omega <- vector("list", p)
  Omega_dot <- vector("list", p)
  for (i in seq_len(p)) {
    Omega[[i]] <- exp_cov_matrix(D, rho[i])
    Omega_dot[[i]] <- D / rho[i]^2 * Omega[[i]]
  }

  Sigma_obs <- matrix(0, n_obs, n_obs)
  for (i in seq_len(p)) {
    a_vec <- A[obs_param, i]
    Sigma_obs <- Sigma_obs + tcrossprod(a_vec) * Omega[[i]][obs_site, obs_site]
  }

  V_obs <- Sigma_obs + W_tap_obs + diag(1e-6, n_obs)
  mu_obs <- rep(beta, each = L)[obs$obs_idx]
  resid <- obs$theta_obs - mu_obs

  R_chol <- tryCatch(chol(V_obs), error = function(e) NULL)
  if (is.null(R_chol)) {
    cache$par <- par; cache$result <- NULL; return(NULL)
  }

  alpha <- backsolve(R_chol, backsolve(R_chol, resid, transpose = TRUE))
  log_det <- 2 * sum(log(diag(R_chol)))
  nll_val <- 0.5 * (n_obs * log(2 * pi) + log_det + sum(resid * alpha))

  result <- list(A = A, rho = rho, alpha = alpha, R_chol = R_chol,
                 nll_val = nll_val, Omega = Omega, Omega_dot = Omega_dot)
  cache$par <- par
  cache$result <- result
  result
}

nll <- function(par) {
  core <- compute_core(par)
  if (is.null(core)) return(1e20)
  core$nll_val
}

nll_grad <- function(par) {
  core <- compute_core(par)
  if (is.null(core)) return(rep(0, length(par)))

  A <- core$A
  alpha <- core$alpha
  R_chol <- core$R_chol
  Omega <- core$Omega
  Omega_dot <- core$Omega_dot

  V_inv <- chol2inv(R_chol)
  Q <- V_inv - tcrossprod(alpha)

  compute_H <- function(Omega_mat) {
    H <- matrix(0, p, p)
    for (j1 in seq_len(p)) {
      idx1 <- param_groups[[j1]]
      if (length(idx1) == 0) next
      s1 <- param_sites[[j1]]
      for (j2 in j1:p) {
        idx2 <- param_groups[[j2]]
        if (length(idx2) == 0) next
        s2 <- param_sites[[j2]]
        val <- sum(Q[idx1, idx2, drop = FALSE] *
                     Omega_mat[s1, s2, drop = FALSE])
        H[j1, j2] <- val
        if (j1 != j2) H[j2, j1] <- val
      }
    }
    H
  }

  H <- lapply(Omega, compute_H)
  H_dot <- lapply(Omega_dot, compute_H)

  grad_beta <- numeric(p)
  for (k in seq_len(p)) {
    idx_k <- param_groups[[k]]
    if (length(idx_k) > 0) grad_beta[k] <- -sum(alpha[idx_k])
  }

  grad_A <- matrix(0, p, p)
  for (a in seq_len(p)) {
    for (b in seq_len(a)) {
      grad_A[a, b] <- sum(H[[b]][a, ] * A[, b])
    }
  }
  grad_A_vec <- grad_A[lower.tri(grad_A, diag = TRUE)]

  rho <- core$rho
  grad_rho <- numeric(p)
  for (k in seq_len(p)) {
    Ak <- A[, k]
    grad_rho[k] <- 0.5 * drop(crossprod(Ak, H_dot[[k]] %*% Ak)) * rho[k]
  }

  c(grad_beta, grad_A_vec, grad_rho)
}

# ══════════════════════════════════════════════════════════════════════════════
# 3. Sanity checks
# ══════════════════════════════════════════════════════════════════════════════

prod_par <- pack_params(model$beta, model$A, model$rho, p = p)
nll_at_opt <- nll(prod_par)
cat(sprintf("NLL at production optimum: %.4f (saved: %.4f, diff: %.2e)\n",
            nll_at_opt, prod_nll, abs(nll_at_opt - prod_nll)))

cache$par <- NULL  # clear cache before gradient check
grad_check <- verify_gradient(prod_par, nll, nll_grad, eps = 1e-5)
cat(sprintf("Gradient check at optimum: max rel error = %.2e\n",
            grad_check$max_rel_err))

# ══════════════════════════════════════════════════════════════════════════════
# 4. Perturbed starting point
# ══════════════════════════════════════════════════════════════════════════════

set.seed(99)
start_par <- prod_par * (1 + rnorm(33, 0, 0.10))
nll_at_start <- nll(start_par)
cat(sprintf("NLL at perturbed start: %.4f\n", nll_at_start))

# ══════════════════════════════════════════════════════════════════════════════
# 5. Run A: Analytic gradient
# ══════════════════════════════════════════════════════════════════════════════

cat("\n== Run A: Analytic gradient ==\n")
cache$par <- NULL
time_A <- system.time({
  result_A <- optim(start_par, nll, gr = nll_grad, method = "L-BFGS-B",
                    lower = rep(-Inf, 33), upper = rep(Inf, 33),
                    control = list(trace = 1, REPORT = 1, maxit = 1000))
})
cat(sprintf("  Time: %.1f s | NLL: %.4f | Conv: %d | Fn: %d | Gr: %d\n",
            time_A["elapsed"], result_A$value, result_A$convergence,
            result_A$counts["function"], result_A$counts["gradient"]))

# ══════════════════════════════════════════════════════════════════════════════
# 6. Run B: Numerical gradient (maxit = 1000)
# ══════════════════════════════════════════════════════════════════════════════

cat("\n== Run B: Numerical gradient (maxit=1000) ==\n")
cache$par <- NULL
time_B <- system.time({
  result_B <- optim(start_par, nll, gr = NULL, method = "L-BFGS-B",
                    lower = rep(-Inf, 33), upper = rep(Inf, 33),
                    control = list(trace = 1, REPORT = 1, maxit = 1000))
})
cat(sprintf("  Time: %.1f s | NLL: %.4f | Conv: %d | Fn: %d | Gr: %d\n",
            time_B["elapsed"], result_B$value, result_B$convergence,
            result_B$counts["function"], result_B$counts["gradient"]))

# ══════════════════════════════════════════════════════════════════════════════
# 7. Run C: Numerical gradient (maxit = 5000, convergence test)
# ══════════════════════════════════════════════════════════════════════════════

cat("\n== Run C: Numerical gradient (maxit=5000) ==\n")
cache$par <- NULL
time_C <- system.time({
  result_C <- optim(start_par, nll, gr = NULL, method = "L-BFGS-B",
                    lower = rep(-Inf, 33), upper = rep(Inf, 33),
                    control = list(trace = 1, REPORT = 10, maxit = 5000))
})
cat(sprintf("  Time: %.1f s | NLL: %.4f | Conv: %d | Fn: %d | Gr: %d\n",
            time_C["elapsed"], result_C$value, result_C$convergence,
            result_C$counts["function"], result_C$counts["gradient"]))

# ══════════════════════════════════════════════════════════════════════════════
# 8. Comparison summary
# ══════════════════════════════════════════════════════════════════════════════

nll_diff_AB <- abs(result_A$value - result_B$value)
nll_diff_AC <- abs(result_A$value - result_C$value)
par_diff_AB <- max(abs(result_A$par - result_B$par))
par_diff_AC <- max(abs(result_A$par - result_C$par))

# Per-iteration time: wall / gradient evaluations (each L-BFGS-B iteration = 1 grad eval)
iter_time_A <- time_A["elapsed"] / max(result_A$counts["gradient"], 1)
iter_time_B <- time_B["elapsed"] / max(result_B$counts["gradient"], 1)
iter_time_C <- time_C["elapsed"] / max(result_C$counts["gradient"], 1)

speedup_B <- time_B["elapsed"] / max(time_A["elapsed"], 0.01)
speedup_C <- time_C["elapsed"] / max(time_A["elapsed"], 0.01)

summary_lines <- c(
  "══════════════════════════════════════════════════════════════",
  "GRADIENT BENCHMARK: ANALYTIC vs NUMERICAL",
  "══════════════════════════════════════════════════════════════",
  "",
  sprintf("Problem: %d observations, 33 parameters", n_obs),
  sprintf("Starting NLL: %.4f (production optimum: %.4f)", nll_at_start, prod_nll),
  "",
  "--- Single-start comparison ---",
  "",
  sprintf("%-28s %12s %12s %12s",
          "", "Analytic", "Numerical", "Num (5000)"),
  sprintf("%-28s %12s %12s %12s",
          "", "(Run A)", "(Run B)", "(Run C)"),
  sprintf("%-28s %12s %12s %12s",
          "---", "---", "---", "---"),
  sprintf("%-28s %12.1f %12.1f %12.1f",
          "Wall time (s)", time_A["elapsed"], time_B["elapsed"], time_C["elapsed"]),
  sprintf("%-28s %12d %12d %12d",
          "Function evaluations", result_A$counts["function"],
          result_B$counts["function"], result_C$counts["function"]),
  sprintf("%-28s %12d %12d %12d",
          "Gradient evaluations", result_A$counts["gradient"],
          result_B$counts["gradient"], result_C$counts["gradient"]),
  sprintf("%-28s %12.4f %12.4f %12.4f",
          "Final NLL", result_A$value, result_B$value, result_C$value),
  sprintf("%-28s %12d %12d %12d",
          "Convergence (0=success)", result_A$convergence,
          result_B$convergence, result_C$convergence),
  sprintf("%-28s %12.3f %12.3f %12.3f",
          "Per-iteration time (s)", iter_time_A, iter_time_B, iter_time_C),
  "",
  sprintf("%-28s %12s %12.2e %12.2e",
          "|NLL_A - NLL_x|", "--", nll_diff_AB, nll_diff_AC),
  sprintf("%-28s %12s %12.2e %12.2e",
          "max |par_A - par_x|", "--", par_diff_AB, par_diff_AC),
  sprintf("%-28s %12s %12.1fx %12.1fx",
          "Total speedup (vs A)", "--", speedup_B, speedup_C),
  sprintf("%-28s %12s %12.1fx %12.1fx",
          "Per-iteration cost ratio", "--",
          iter_time_B / max(iter_time_A, 1e-6),
          iter_time_C / max(iter_time_A, 1e-6)),
  ""
)

# Check if Run B converged
if (result_B$convergence == 1) {
  summary_lines <- c(summary_lines,
    sprintf("Note: Run B did not converge (maxit=1000 reached)."),
    sprintf("  NLL gap from analytic optimum: %.4f", nll_diff_AB),
    sprintf("  Still decreasing: check Run C results.")
  )
}
if (result_C$convergence == 1) {
  summary_lines <- c(summary_lines,
    sprintf("Note: Run C did not converge (maxit=5000 reached)."),
    sprintf("  NLL gap from analytic optimum: %.4f", nll_diff_AC)
  )
} else if (result_C$convergence == 0) {
  summary_lines <- c(summary_lines,
    sprintf("Run C converged after %d gradient evaluations (%.1f s).",
            result_C$counts["gradient"], time_C["elapsed"]),
    sprintf("  Same optimum as analytic: %s (NLL diff = %.2e)",
            if (nll_diff_AC < 0.01) "YES" else "NO", nll_diff_AC)
  )
}

# Print to console
cat("\n")
for (line in summary_lines) cat(line, "\n")

# Save to file
writeLines(summary_lines, "tables/gradient_benchmark.txt")
cat("\nSaved to tables/gradient_benchmark.txt\n")

# ══════════════════════════════════════════════════════════════════════════════
# 9. Multi-start comparison (3 starts), if single-start numerical < 10 min
# ══════════════════════════════════════════════════════════════════════════════

if (time_B["elapsed"] < 600) {
  cat("\n== Multi-start comparison (3 starts) ==\n")
  n_starts <- 3
  set.seed(2026)

  # Generate 3 random rho vectors (log-uniform[50, 5000] km)
  rho_starts <- lapply(seq_len(n_starts), function(i) {
    exp(runif(p, log(50), log(5000)))
  })

  # Use production beta and A as base, vary only rho
  base_beta <- model$beta
  base_A <- model$A

  ms_time_A <- 0
  ms_time_B <- 0
  ms_best_A <- Inf
  ms_best_B <- Inf

  for (i in seq_len(n_starts)) {
    par_i <- pack_params(base_beta, base_A, rho_starts[[i]], p = p)

    # Analytic
    cache$par <- NULL
    t_a <- system.time({
      r_a <- tryCatch(
        optim(par_i, nll, gr = nll_grad, method = "L-BFGS-B",
              lower = rep(-Inf, 33), upper = rep(Inf, 33),
              control = list(maxit = 2000, trace = 0)),
        error = function(e) NULL
      )
    })
    if (!is.null(r_a)) {
      ms_time_A <- ms_time_A + t_a["elapsed"]
      if (r_a$value < ms_best_A) ms_best_A <- r_a$value
      cat(sprintf("  Start %d analytic:  %.1f s, NLL=%.4f, conv=%d\n",
                  i, t_a["elapsed"], r_a$value, r_a$convergence))
    }

    # Numerical
    cache$par <- NULL
    t_b <- system.time({
      r_b <- tryCatch(
        optim(par_i, nll, gr = NULL, method = "L-BFGS-B",
              lower = rep(-Inf, 33), upper = rep(Inf, 33),
              control = list(maxit = 2000, trace = 0)),
        error = function(e) NULL
      )
    })
    if (!is.null(r_b)) {
      ms_time_B <- ms_time_B + t_b["elapsed"]
      if (r_b$value < ms_best_B) ms_best_B <- r_b$value
      cat(sprintf("  Start %d numerical: %.1f s, NLL=%.4f, conv=%d\n",
                  i, t_b["elapsed"], r_b$value, r_b$convergence))
    }
  }

  ms_lines <- c(
    "",
    "--- Multi-start comparison (3 starts) ---",
    "",
    sprintf("%-28s %12s %12s", "", "Analytic", "Numerical"),
    sprintf("%-28s %12s %12s", "---", "---", "---"),
    sprintf("%-28s %12.1f %12.1f", "Total time (s)", ms_time_A, ms_time_B),
    sprintf("%-28s %12.4f %12.4f", "Best NLL", ms_best_A, ms_best_B),
    sprintf("%-28s %12.1fx", "Multi-start speedup", ms_time_B / max(ms_time_A, 0.01)),
    sprintf("%-28s %12.1f %12.1f",
            "Extrapolated 20-start (s)", ms_time_A / 3 * 20, ms_time_B / 3 * 20),
    sprintf("%-28s %12.1f %12.1f",
            "Extrapolated 20-start (min)", ms_time_A / 3 * 20 / 60,
            ms_time_B / 3 * 20 / 60)
  )

  for (line in ms_lines) cat(line, "\n")

  # Append to file
  all_lines <- c(readLines("tables/gradient_benchmark.txt"), ms_lines)
  writeLines(all_lines, "tables/gradient_benchmark.txt")
  cat("\nUpdated tables/gradient_benchmark.txt\n")

} else {
  cat(sprintf("\nSkipping multi-start: single numerical run took %.0f s (> 600 s limit)\n",
              time_B["elapsed"]))

  skip_lines <- c(
    "",
    "--- Multi-start comparison ---",
    sprintf("Skipped: single numerical run took %.0f s (> 10 min limit)", time_B["elapsed"]),
    sprintf("Extrapolated 20-start numerical time: ~%.0f min",
            time_B["elapsed"] / 60 * 20)
  )
  all_lines <- c(readLines("tables/gradient_benchmark.txt"), skip_lines)
  writeLines(all_lines, "tables/gradient_benchmark.txt")
}
