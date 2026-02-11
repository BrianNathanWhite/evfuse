# evfuse

Spatial data fusion for extreme value analysis. Implements a two-stage
frequentist framework for jointly modeling multiple spatial data sources with
GEV (Generalized Extreme Value) margins via a multivariate Gaussian process
with coregionalization structure, following Russell et al. (2019).

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

2. **Stage 2** -- Estimate a p-dimensional GP with coregionalization: mean
   vector beta, lower-triangular cross-covariance factor A, and range
   parameters rho, all via `optim` with analytic gradients
   (`fit_spatial_model()`).

3. **Kriging** -- Predict at unobserved locations using partial observations:
   each data source observes a subset of the p parameters
   (`predict_krig()`).

4. **Return levels** -- Delta-method and simulation-based confidence
   intervals for r-year return levels (`compute_return_levels()`).

## Customization

The `source_params` argument controls the mapping from data sources to model
parameters. The default `list(NOAA = 1:3, ADCIRC = 4:6)` defines a 6-dim
model where NOAA sites observe parameters 1-3 and ADCIRC sites observe 4-6.
This can be adapted for other source configurations:

```r
# Example: 3 sources, 9-dim model
dat <- load_data(df, source_params = list(
  gauge = 1:3,
  satellite = 4:6,
  model = 7:9
))
```

## References

Russell, B. T., Cooley, D. S., Porter, W. C., Reich, B. J., & Heald, C. L.
(2019). Data fusion for environmental and human health applications.
*International Statistical Review*, 87, S261--S287.
