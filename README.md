# evfuse

Fusing sparse observations and dense simulations for spatial extreme value
analysis. Implements the two-stage frequentist framework from:

> White, B. N., Blanton, B., Luettich, R., & Smith, R. L. Fusing Sparse
> Observations and Dense Simulations for Spatial Extreme Value Analysis:
> Application to U.S. Coastal Sea Levels. *In preparation.*

The package was developed for fusing NOAA tidal gauge observations with
ADCIRC hydrodynamic simulations of coastal sea levels, but the framework is
general: any application with annual maxima from multiple spatial data
sources (observed, modeled, satellite, reanalysis, etc.) can use `evfuse`
by specifying the source-to-parameter mapping via `source_params`.

## Installation

```r
# install.packages("devtools")
devtools::install_github("BrianNathanWhite/evfuse")
```

## Quick Start

```r
library(evfuse)

# Bundled dataset: 29 NOAA + 100 ADCIRC sites along U.S. Gulf/Atlantic coasts
data(gulf_data)
D <- compute_distances(gulf_data$sites)

# Stage 1: Fit GEV at each site
stage1 <- fit_gev_all(gulf_data)

# Bootstrap measurement uncertainty
bs <- bootstrap_W(gulf_data, B = 500, seed = 42)
W_tap <- taper_W(bs$W_bs, D, lambda = 300)

# Stage 2: Fit 6-dim joint GP model
model <- fit_spatial_model(stage1, gulf_data, W_tap, D)

# Predict at new locations
preds <- predict_krig(model, new_sites)
rl <- compute_return_levels(preds, r = 100)
```

## Pipeline Overview

1. **Stage 1** -- Fit GEV(mu, log sigma, xi) independently at each site via
   maximum likelihood (`fit_gev_all()`).

2. **Stage 2** -- Estimate a p-dimensional linear model of coregionalization:
   mean vector beta, lower-triangular cross-covariance factor A, and range
   parameters rho, all via `optim` with analytic gradients
   (`fit_spatial_model()`). A selection matrix accommodates partial
   observations, where each site provides estimates from only one source.

3. **Kriging** -- Predict at unobserved locations, leveraging information
   from all data sources even though each site observes only a subset of
   the p parameters (`predict_krig()`).

4. **Return levels** -- Delta-method and simulation-based confidence
   intervals for r-year return levels (`compute_return_levels()`).

## Applying to Your Own Data

Your input data frame needs columns: `lon`, `lat`, `location`, `year`,
`max_sea_level`, and `data_source`. The `source_params` argument controls
the mapping from data sources to model parameters:

```r
# Default: 2 sources, 6-dim model
dat <- load_data(df)  # source_params = list(NOAA = 1:3, ADCIRC = 4:6)

# Custom: 3 sources, 9-dim model
dat <- load_data(df, source_params = list(
  gauge      = 1:3,
  satellite  = 4:6,
  reanalysis = 7:9
))
```

Each source contributes 3 GEV parameters (mu, log sigma, xi) per site. The
framework jointly models all sources, learning cross-source correlations
that allow dense but potentially biased simulations to improve estimation
at sparse observation locations.

## Key Results (U.S. Coastal Application)

- Cross-source correlations: r = 0.99 (location), r = 0.95 (shape)
- 68% reduction in location parameter LOO-CV prediction error
- 35% reduction in 100-year return level RMSE vs. NOAA-only model
- Fusion gains robust to detrending for sea level rise

## References

White, B. N., Blanton, B., Luettich, R., & Smith, R. L. Fusing Sparse
Observations and Dense Simulations for Spatial Extreme Value Analysis:
Application to U.S. Coastal Sea Levels. *In preparation.*

Russell, B. T., Cooley, D. S., Porter, W. C., Reich, B. J., & Heald, C. L.
(2019). Data fusion for environmental and human health applications.
*International Statistical Review*, 87, S261--S287.
