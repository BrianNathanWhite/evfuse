# U.S. coastal sea level data

A pre-loaded `evfuse_data` object containing annual maximum sea levels
at 129 sites along the U.S. East and Gulf Coasts: 29 NOAA tidal gauge
stations (34–43 years of observations) and 100 ADCIRC numerical
simulation points (43 years each).

## Usage

``` r
coast_data
```

## Format

An `evfuse_data` list with components:

- sites:

  Data frame (129 rows) with lon, lat, location, data_source.

- maxima:

  Named list of 129 numeric vectors of annual maxima (meters).

- n_noaa:

  29

- n_adcirc:

  100

- n_sites:

  129

- source_params:

  `list(NOAA = 1:3, ADCIRC = 4:6)`

- p:

  6

- raw_df:

  Data frame (5503 rows) with lon, lat, location, year, max_sea_level,
  data_source. The long-format source data, needed by
  [`fit_gev_detrended`](https://briannathanwhite.github.io/evfuse/reference/fit_gev_detrended.md)
  and
  [`bootstrap_W_detrended`](https://briannathanwhite.github.io/evfuse/reference/bootstrap_W_detrended.md)
  for nonstationary GEV fitting.

## Source

NOAA CO-OPS tidal gauges and ADCIRC storm surge simulations. See White
et al. for details.
