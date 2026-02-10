# evfuse - Spatial Extremes Data Fusion

R package for dissertation on spatial extremes data fusion. Implements a two-stage frequentist approach to jointly model NOAA tidal gauge and ADCIRC simulation sea level data via a 6-dimensional Gaussian process with coregionalization structure, following Russell et al. (2019).

## Pipeline

1. **Stage 1**: Fit GEV(mu, log_sigma, xi) independently at each site via `extRemes::fevd`
2. **Stage 2**: Fit a 6-dim GP (3 NOAA + 3 ADCIRC GEV params) with coregionalization (beta, A, rho) via `optim`
3. **Kriging**: Predict with partial observations (NOAA sites observe dims 1-3, ADCIRC sites observe dims 4-6)
4. **Return levels**: Delta method and simulation CIs

## Data

- Source: `data-raw/combined_max.csv`
- 129 sites: 29 NOAA, 100 ADCIRC
- 34-43 years of annual maxima per site
- Load via:
  ```r
  df <- read.csv("data-raw/combined_max.csv")
  dat <- load_data(df)
  ```

## Mathematical Details

### Parameter Structure
- 6 parameters per site: (mu_NOAA, log_sigma_NOAA, xi_NOAA, mu_ADCIRC, log_sigma_ADCIRC, xi_ADCIRC)
- Stacked across L=129 sites giving 774-dimensional theta vector
- Each site observes only 3 of 6 dimensions
- All parameters in log_sigma space (unconstrained optimization)

### Stage 2 Likelihood
Multivariate normal: Theta_hat ~ N(beta ⊗ 1_L, Sigma_{A,rho} + W_tap)

Only observed components enter likelihood (selection matrix per site).

### Covariance Structure
```
Sigma_{A,rho} = (I_L ⊗ A) [sum_i (e_i e_i^T ⊗ Omega_i(rho))] (I_L ⊗ A)^T
```
- A is lower triangular 6x6
- Omega_i is exponential covariance with range rho_i
- p=6 exponential correlation matrices in block-diagonal structure

### Tapering Matrix (W_tap)
- Uses Wendland2 taper: `T_tap = 1_p 1_p^T ⊗ C_W2(D/lambda)`
- Applied via Hadamard product to bootstrap covariance from stage 1
- Adds bootstrap uncertainty while maintaining sparsity

### Kriging
- Standard GP conditioning with partial observations
- Extract upper-left 3x3 (NOAA params) from predicted 6x6 covariance at ADCIRC locations

## Development Practices

- Use `devtools::load_all()` and `devtools::test()` to iterate
- Run `test()` after every change
- Don't install packages without asking first
- If a fix involves changing mathematical structure (not just R syntax), explain what's changing and why before implementing

## Package Dependencies

- extRemes
- Matrix
- stats
- testthat (suggests)

## Current Status

### Validated
- `load_data()`: 129 sites, 29 NOAA, 100 ADCIRC, 34-43 years
- Stage 1 GEV fitting: 129/129 sites converge, all vcov matrices PSD
- Stage 2 `fit_spatial_model()` with analytic gradient and log-rho reparameterization
- `embed_W()`: maps 387x387 bootstrap covariance into 774x774 space
- `verify_gradient()`: analytic vs numerical finite-difference check (max rel error ~9e-6)
- 33 tests passing

### Key Implementation Notes
- `extRemes::fevd` returns observed information matrix directly (positive definite), not Hessian of negative log-likelihood
- Use `solve(fit$results$hessian)` not `solve(-fit$results$hessian)` for vcov
- Delta method applied for log_sigma transformation: J = diag(1, 1/sigma, 1)
- **Rho is optimized in log space**: `pack_params` stores `log(rho)`, `unpack_params` exponentiates back. Gradient uses chain rule (`grad_log_rho = grad_rho * rho`). No box constraints needed.
- **Analytic gradient**: contracts Q = V_inv - alpha*alpha^T with Omega_i into 6x6 H matrices. dSigma/dA_{ab} and dSigma/drho_k are rank-1 Kronecker products, so each gradient component is a dot product on H. See `data-raw/evfuse_gradient.pdf`.
- **Sigma_obs built directly** in 387x387 observed space, never forming the full 774x774 matrix.
- **nll/nll_grad share a cache** keyed on the parameter vector to avoid redundant Cholesky when L-BFGS-B calls fn then gr at the same point.

### Saved Models
- `data-raw/model_6dim.rds`: current best fit (NLL = -77.52, convergence 0)
  - rho = (3052, 260, 893, 362, 870, 1000) — rho[6]=1000 looks stuck

## Next Steps

- Run `scripts/multi_start.R` to search for better local optima (20 random rho starts, log-uniform [50, 5000] km)
- Bootstrap for W_tap
- Kriging predictions and return level maps
