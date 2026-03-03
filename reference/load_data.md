# Load and validate sea level data

Reads a data frame with columns: lon, lat, location, year,
max_sea_level, data_source. Validates required columns and types, then
returns a structured list ready for stage-1 GEV fitting.

## Usage

``` r
load_data(df, source_params = list(NOAA = 1:3, ADCIRC = 4:6))
```

## Arguments

- df:

  A data frame with the required columns.

- source_params:

  Named list mapping data source labels to parameter indices in the
  joint model. Default: `list(NOAA = 1:3, ADCIRC = 4:6)`.

## Value

An `evfuse_data` list with components:

- sites:

  Data frame of unique sites with lon, lat, location, data_source.

- maxima:

  Named list of numeric vectors of annual maxima, keyed by location.

- n_noaa:

  Number of NOAA sites.

- n_adcirc:

  Number of ADCIRC sites.

- n_sites:

  Total number of sites.

- source_params:

  The source-to-parameter mapping.

- p:

  Total number of parameters per site in the joint model.
