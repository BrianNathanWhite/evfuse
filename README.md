# evfuse

Fusing sparse observations and dense simulations for spatial extreme value
analysis. Implements the two-stage frequentist framework from:

> White, B. N., Blanton, B., Luettich, R., & Smith, R. L. Fusing Sparse
> Observations and Dense Simulations for Spatial Extreme Value Analysis:
> Application to U.S. Coastal Sea Levels. *arXiv preprint*, 2026.
> [arXiv:2603.03247](https://arxiv.org/abs/2603.03247)

Developed for fusing NOAA tide gauge observations with ADCIRC hydrodynamic
simulations, but the framework is general: any application with annual
maxima from multiple spatial data sources can use `evfuse` by specifying
the source-to-parameter mapping via `source_params`. See the
[tutorial vignette](vignettes/evfuse-tutorial.Rmd) for details.

## Installation

```r
# install.packages("devtools")
devtools::install_github("BrianNathanWhite/evfuse")
```

Requires R >= 3.5. Dependencies (`extRemes`, `Matrix`) are installed
automatically.

## Quick Start

```r
library(evfuse)

data(coast_data)
D <- compute_distances(coast_data$sites)

# Stage 1: site-wise GEV fits
stage1 <- fit_gev_all(coast_data)

# Bootstrap measurement uncertainty
bs <- bootstrap_W(coast_data, B = 500, seed = 42)
W_tap <- taper_W(bs$W_bs, D, lambda = 300)

# Stage 2: joint GP model
model <- fit_spatial_model(stage1, coast_data, W_tap, D)

# Predict at new locations
new_sites <- data.frame(lon = c(-90.0, -81.5), lat = c(30.0, 31.5))
preds <- predict_krig(model, new_sites)
rl <- compute_return_levels(preds, r = 100)
rl$return_level
rl$se_sim
```

## Reproducing the Paper

```bash
git clone https://github.com/BrianNathanWhite/evfuse.git
cd evfuse
Rscript scripts/run_nonstationary.R
```

This runs the full pipeline end-to-end (~15 min): Stage 1 fitting with
linear trend at NOAA sites, bootstrap, Stage 2 coregionalization, kriging,
return level maps, LOO-CV, block CV, and all manuscript figures. Output
goes to `figures/` and `tables/`.

Additional standalone scripts:

```bash
Rscript scripts/run_trends.R           # Trend diagnostics
Rscript scripts/plot_study_area.R      # Study area map (Figure 1)
Rscript scripts/simulation_study.R     # Parameter recovery (§4.6.4)
Rscript scripts/rmse_decomposition.R   # RMSE by parameter/region (Table 3)
Rscript scripts/baseline_comparisons.R # Bias correction baselines (§5.1)
Rscript scripts/gradient_benchmark.R   # Analytic vs numerical gradient
```

## References

White, B. N., Blanton, B., Luettich, R., & Smith, R. L. Fusing Sparse
Observations and Dense Simulations for Spatial Extreme Value Analysis:
Application to U.S. Coastal Sea Levels. *arXiv preprint*, 2026.
[arXiv:2603.03247](https://arxiv.org/abs/2603.03247)

Russell, B. T., Risser, M. D., Smith, R. L., & Kunkel, K. E. (2020).
Investigating the association between late spring Gulf of Mexico sea
surface temperatures and U.S. Gulf Coast precipitation extremes with
focus on Hurricane Harvey. *Environmetrics*, 31(2), e2595.
