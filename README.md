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
Rscript -e 'install.packages(c("ggplot2", "sf", "gridExtra", "maps", "ragg"))'
Rscript scripts/run_nonstationary.R
```

The `install.packages()` line covers the figure dependencies, which are
Suggests and therefore not installed automatically with the package.

This runs the full pipeline end-to-end (~15 min): Stage 1 fitting with
linear trend at NOAA sites, bootstrap, Stage 2 coregionalization, kriging,
return level maps, LOO-CV, block CV, and most manuscript figures. Output
goes to `figures/` and `tables/`.

Additional scripts produce the remaining figures and analyses (most read
the fitted models written to `data-raw/` by the main pipeline above, so
run it first):

```bash
Rscript scripts/ad_gof.R               # GEV goodness-of-fit, bootstrap (Section S2)
Rscript scripts/run_trends.R           # Trend diagnostics
Rscript scripts/plot_study_area.R      # Study area map (Figure 1)
Rscript scripts/simulation_study.R     # Parameter recovery (§4.6.4)
Rscript scripts/rmse_decomposition.R   # RMSE by parameter/region (Table 3)
Rscript scripts/baseline_comparisons.R # Bias correction baselines (§5.1)
Rscript scripts/gradient_benchmark.R   # Analytic vs numerical gradient
Rscript scripts/saturation.R           # ADCIRC subsampling curve (Figure 7)
Rscript scripts/se_comparison.R        # Delta vs simulation SEs (Figure S6)
Rscript scripts/combine_ratio_maps.R   # Assemble Figure 6 from its two panels
Rscript scripts/run_w_sensitivity.R    # Bootstrap batch stability (§4.6.3)
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
