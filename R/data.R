#' Gulf Coast sea level data
#'
#' A pre-loaded \code{evfuse_data} object containing annual maximum sea levels
#' at 129 sites along the U.S. Gulf and Atlantic coasts: 29 NOAA tidal gauge
#' stations (34--43 years of observations) and 100 ADCIRC numerical simulation
#' points (43 years each).
#'
#' @format An \code{evfuse_data} list with components:
#' \describe{
#'   \item{sites}{Data frame (129 rows) with lon, lat, location, data_source.}
#'   \item{maxima}{Named list of 129 numeric vectors of annual maxima (meters).}
#'   \item{n_noaa}{29}
#'   \item{n_adcirc}{100}
#'   \item{n_sites}{129}
#'   \item{source_params}{\code{list(NOAA = 1:3, ADCIRC = 4:6)}}
#'   \item{p}{6}
#' }
#' @source NOAA CO-OPS tidal gauges and ADCIRC storm surge simulations.
"gulf_data"

#' Load and validate sea level data
#'
#' Reads a data frame with columns: lon, lat, location, year, max_sea_level,
#' data_source. Validates required columns and types, then returns a structured
#' list ready for stage-1 GEV fitting.
#'
#' @param df A data frame with the required columns.
#' @param source_params Named list mapping data source labels to parameter
#'   indices in the joint model. Default: \code{list(NOAA = 1:3, ADCIRC = 4:6)}.
#' @return An \code{evfuse_data} list with components:
#'   \describe{
#'     \item{sites}{Data frame of unique sites with lon, lat, location, data_source.}
#'     \item{maxima}{Named list of numeric vectors of annual maxima, keyed by location.}
#'     \item{n_noaa}{Number of NOAA sites.}
#'     \item{n_adcirc}{Number of ADCIRC sites.}
#'     \item{n_sites}{Total number of sites.}
#'     \item{source_params}{The source-to-parameter mapping.}
#'     \item{p}{Total number of parameters per site in the joint model.}
#'   }
#' @export
load_data <- function(df, source_params = list(NOAA = 1:3, ADCIRC = 4:6)) {
  required_cols <- c("lon", "lat", "location", "year", "max_sea_level", "data_source")
  missing <- setdiff(required_cols, names(df))
  if (length(missing) > 0) {
    stop("Missing required columns: ", paste(missing, collapse = ", "))
  }

  # Validate data_source values against source_params
  valid_sources <- names(source_params)
  bad_sources <- setdiff(unique(df$data_source), valid_sources)
  if (length(bad_sources) > 0) {
    stop("Invalid data_source values: ", paste(bad_sources, collapse = ", "),
         ". Must be one of: ", paste(valid_sources, collapse = ", "), ".")
  }

  # Derive total parameter dimension
  p <- max(unlist(source_params))

  # Ensure location is character (numeric IDs cause indexing issues)
  df$location <- as.character(df$location)

  # Check for NA in max_sea_level
  n_na <- sum(is.na(df$max_sea_level))
  if (n_na > 0) {
    message("Dropping ", n_na, " rows with NA max_sea_level.")
    df <- df[!is.na(df$max_sea_level), ]
  }

  # Extract unique sites (take first row per location for coordinates)
  sites <- df[!duplicated(df$location), c("lon", "lat", "location", "data_source")]
  rownames(sites) <- NULL

  # Order: sources in order of source_params, then by location
  sites <- sites[order(factor(sites$data_source, levels = valid_sources),
                        sites$location), ]
  rownames(sites) <- NULL

  # Extract annual maxima as a named list
  maxima <- split(df$max_sea_level, df$location)
  # Reorder to match sites
  maxima <- maxima[sites$location]

  # Per-source counts
  source_counts <- vapply(valid_sources, function(s) sum(sites$data_source == s), integer(1))
  n_noaa <- if ("NOAA" %in% valid_sources) source_counts[["NOAA"]] else 0L
  n_adcirc <- if ("ADCIRC" %in% valid_sources) source_counts[["ADCIRC"]] else 0L

  count_str <- paste(sprintf("%d %s", source_counts, names(source_counts)), collapse = ", ")
  message(sprintf("Loaded %d sites: %s", nrow(sites), count_str))
  message(sprintf("Years of data per site: %d to %d",
                  min(lengths(maxima)), max(lengths(maxima))))

  structure(
    list(
      sites = sites,
      maxima = maxima,
      n_noaa = n_noaa,
      n_adcirc = n_adcirc,
      n_sites = nrow(sites),
      source_params = source_params,
      p = p
    ),
    class = "evfuse_data"
  )
}

#' Compute pairwise great-circle distances between sites
#'
#' Uses the Haversine formula. Returns distances in kilometers.
#'
#' @param sites Data frame with lon and lat columns (in degrees).
#' @return Symmetric distance matrix (km).
#' @export
compute_distances <- function(sites) {
  n <- nrow(sites)
  lon_rad <- sites$lon * pi / 180
  lat_rad <- sites$lat * pi / 180
  R <- 6371  # Earth radius in km

  D <- matrix(0, n, n)
  for (i in seq_len(n - 1)) {
    for (j in (i + 1):n) {
      dlat <- lat_rad[j] - lat_rad[i]
      dlon <- lon_rad[j] - lon_rad[i]
      a <- sin(dlat / 2)^2 + cos(lat_rad[i]) * cos(lat_rad[j]) * sin(dlon / 2)^2
      D[i, j] <- 2 * R * asin(sqrt(a))
      D[j, i] <- D[i, j]
    }
  }
  D
}
