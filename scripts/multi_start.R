#!/usr/bin/env Rscript
# Multi-start optimization for stage 2 spatial model
#
# Generates random starting rho values (log-uniform on [50, 5000] km),
# keeps beta and A from the current best fit, and runs fit_spatial_model
# from each. Prints a summary table and saves the best model.

devtools::load_all()

# --- Load current best model ---
model5 <- readRDS("data-raw/model_6dim.rds")

cat("Current model NLL:", model5$optim_result$value, "\n")
cat("Current rho:", round(model5$rho, 1), "\n\n")

# --- Extract fixed components ---
stage1 <- model5$stage1
dat    <- model5$dat
W_tap  <- model5$W_tap
D      <- model5$D
beta0  <- model5$beta
A0     <- model5$A

# --- Multi-start settings ---
n_starts <- 20
set.seed(2026)

# Log-uniform on [50, 5000]: exp(runif(6, log(50), log(5000)))
rho_starts <- lapply(seq_len(n_starts), function(i) {
  exp(runif(6, log(50), log(5000)))
})

# --- Run fits ---
results <- vector("list", n_starts)
nlls    <- numeric(n_starts)

cat(sprintf("%-5s  %-12s  %s\n", "Start", "NLL", "rho"))
cat(strrep("-", 80), "\n")

for (i in seq_len(n_starts)) {
  start_i <- list(beta = beta0, A = A0, rho = rho_starts[[i]])

  fit <- tryCatch(
    suppressWarnings(
      fit_spatial_model(stage1, dat, W_tap, D,
                        start = start_i,
                        control = list(maxit = 2000, trace = 0))
    ),
    error = function(e) {
      message(sprintf("  Start %2d: ERROR — %s", i, e$message))
      NULL
    }
  )

  if (is.null(fit)) {
    nlls[i] <- Inf
    next
  }

  results[[i]] <- fit
  nlls[i] <- fit$optim_result$value

  rho_str <- paste(sprintf("%7.1f", fit$rho), collapse = "  ")
  cat(sprintf("%-5d  %-12.4f  %s\n", i, nlls[i], rho_str))
}

# --- Include the original model for comparison ---
cat(strrep("-", 80), "\n")
rho_str <- paste(sprintf("%7.1f", model5$rho), collapse = "  ")
cat(sprintf("%-5s  %-12.4f  %s\n", "orig", model5$optim_result$value, rho_str))

# --- Pick the best ---
best_idx <- which.min(nlls)
best_nll <- nlls[best_idx]
orig_nll <- model5$optim_result$value

cat("\n")
if (best_nll < orig_nll) {
  cat(sprintf("Best: start %d (NLL = %.4f, improvement = %.4f)\n",
              best_idx, best_nll, orig_nll - best_nll))
  cat("Rho:", round(results[[best_idx]]$rho, 1), "\n")

  saveRDS(results[[best_idx]], "data-raw/model_6dim_best.rds")
  cat("Saved to data-raw/model_6dim_best.rds\n")
} else {
  cat(sprintf("Original model is still best (NLL = %.4f)\n", orig_nll))
}

# --- Save all results for inspection ---
summary_df <- data.frame(
  start = seq_len(n_starts),
  nll = nlls,
  converged = vapply(results, function(r) {
    if (is.null(r)) NA else r$optim_result$convergence == 0
  }, logical(1)),
  do.call(rbind, lapply(seq_len(n_starts), function(i) {
    if (is.null(results[[i]])) rep(NA, 6)
    else results[[i]]$rho
  }))
)
names(summary_df)[4:9] <- paste0("rho", 1:6)
saveRDS(summary_df, "data-raw/multi_start_summary.rds")
cat("Summary saved to data-raw/multi_start_summary.rds\n")
