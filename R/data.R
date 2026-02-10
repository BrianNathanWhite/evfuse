#' Load and validate sea level data
#'
#' Reads a data frame with columns: lon, lat, location, year, max_sea_level,
#' data_source (ADCIRC or NOAA). Validates required columns and types, then
#' returns a structured list ready for stage-1 GEV fitting.
#'
#' @param df A data frame with the required columns.
#' @return A list with components:
#'   \describe{
#'     \item{sites}{Data frame of unique sites with lon, lat, location, data_source.}
#'     \item{maxima}{Named list of numeric vectors of annual maxima, keyed by location.}
#'     \item{n_noaa}{Number of NOAA sites.}
#'     \item{n_adcirc}{Number of ADCIRC sites.}
#'     \item{n_sites}{Total number of sites.}
#'   }
#' @export
load_data <- function(df) {
  required_cols <- c("lon", "lat", "location", "year", "max_sea_level", "data_source")
  missing <- setdiff(required_cols, names(df))
  if (length(missing) > 0) {
    stop("Missing required columns: ", paste(missing, collapse = ", "))
  }

  # Validate data_source values
  valid_sources <- c("NOAA", "ADCIRC")
  bad_sources <- setdiff(unique(df$data_source), valid_sources)
  if (length(bad_sources) > 0) {
    stop("Invalid data_source values: ", paste(bad_sources, collapse = ", "),
         ". Must be 'NOAA' or 'ADCIRC'.")
  }

  # Check for NA in max_sea_level
  n_na <- sum(is.na(df$max_sea_level))
  if (n_na > 0) {
    message("Dropping ", n_na, " rows with NA max_sea_level.")
    df <- df[!is.na(df$max_sea_level), ]
  }

  # Extract unique sites (take first row per location for coordinates)
  sites <- df[!duplicated(df$location), c("lon", "lat", "location", "data_source")]
  rownames(sites) <- NULL

  # Order: NOAA first, then ADCIRC
  sites <- sites[order(factor(sites$data_source, levels = c("NOAA", "ADCIRC")),
                        sites$location), ]
  rownames(sites) <- NULL

  # Extract annual maxima as a named list
  maxima <- split(df$max_sea_level, df$location)
  # Reorder to match sites
  maxima <- maxima[sites$location]

  n_noaa <- sum(sites$data_source == "NOAA")
  n_adcirc <- sum(sites$data_source == "ADCIRC")

  message(sprintf("Loaded %d sites: %d NOAA, %d ADCIRC", nrow(sites), n_noaa, n_adcirc))
  message(sprintf("Years of data per site: %d to %d",
                  min(lengths(maxima)), max(lengths(maxima))))

  structure(
    list(
      sites = sites,
      maxima = maxima,
      n_noaa = n_noaa,
      n_adcirc = n_adcirc,
      n_sites = nrow(sites)
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
