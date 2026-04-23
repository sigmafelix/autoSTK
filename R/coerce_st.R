# coerce_st.R — internal helpers for ST class coercion
# All functions are internal (not exported); used throughout the package
# to avoid repeating the conversion chain across files.
#
# Supported input classes:
#   spacetime: STFDF, STSDF, STIDF
#   sftime:    sftime
#   cubble:    cubble_df (nested or temporal face, requires the cubble package)


# ---------------------------------------------------------------------------
# cubble → sftime conversion
# ---------------------------------------------------------------------------

# Check whether x is a cubble (nested or temporal face).
.is_cubble <- function(x) {
  inherits(x, "cubble_df")
}

#' Convert a cubble object to an sftime object
#'
#' Handles both nested and temporal faces of a \pkg{cubble} object.  The
#' resulting \code{sftime} has one row per (location × time) observation with
#' the spatial geometry attached to every row.
#'
#' @param x A \code{cubble_df} object (from the \pkg{cubble} package).
#' @param time_col Character. Name of the time column in the temporal face.
#'   If \code{NULL} (default), the function tries to detect it automatically
#'   by looking for POSIXct/Date columns.
#' @param key_col Character. Name of the site identifier column shared by the
#'   spatial and temporal faces.  Detected automatically from
#'   \code{cubble::key_vars()} when \code{NULL}.
#' @return An \code{sftime} object with a \code{geometry} column (sfc_POINT).
#' @importFrom cubble is_cubble_spatial face_temporal
#' @export
cubble_to_sftime <- function(x, time_col = NULL, key_col = NULL) {
  if (!requireNamespace("cubble", quietly = TRUE))
    stop("Package 'cubble' is required to convert cubble objects. ",
         "Install it with: install.packages('cubble')")

  # ---- Ensure we work with the temporal (long) face -----------------------
  if (cubble::is_cubble_spatial(x)) {
    x_long <- cubble::face_temporal(x)
  } else {
    x_long <- x
  }

  # ---- Detect key column --------------------------------------------------
  # cubble stores the key as attr(x, "key"), a one-row-per-site tibble whose
  # first non-.rows column is the site identifier.
  if (is.null(key_col)) {
    key_attr <- attr(x_long, "key")
    key_candidates <- setdiff(names(key_attr), ".rows")
    if (length(key_candidates) > 0L) {
      key_col <- key_candidates[[1L]]
    } else {
      stop("Cannot detect key column from cubble object. ",
           "Pass key_col explicitly.")
    }
  }

  # ---- Detect time column -------------------------------------------------
  # cubble stores the time index as attr(x, "index") — a plain character string.
  # Fall back to scanning for POSIXct/Date columns when the attribute is absent.
  df_long <- as.data.frame(x_long)

  if (is.null(time_col)) {
    idx_attr <- attr(x_long, "index")
    # Use the cubble index attribute only if the column is actually POSIXct/Date
    if (!is.null(idx_attr) && nchar(idx_attr) > 0L &&
        idx_attr %in% names(df_long) &&
        inherits(df_long[[idx_attr]], c("POSIXct", "POSIXlt", "Date"))) {
      time_col <- idx_attr
      message("cubble_to_sftime: using '", time_col, "' as the time column.")
    } else {
      # Fallback: scan for POSIXct/Date columns not used as the key
      candidate_cols <- names(df_long)[
        vapply(df_long, function(col)
          inherits(col, c("POSIXct", "POSIXlt", "Date")), logical(1L))
      ]
      candidate_cols <- setdiff(candidate_cols, key_col)
      if (length(candidate_cols) == 0L)
        stop("Cannot auto-detect a time column (POSIXct/Date) in the cubble ",
             "temporal face. Pass time_col explicitly.")
      time_col <- candidate_cols[[1L]]
      message("cubble_to_sftime: using '", time_col, "' as the time column.")
    }
  }

  # ---- Extract spatial geometry from the spatial (nested) face ------------
  if (cubble::is_cubble_spatial(x)) {
    x_nested <- x
  } else {
    x_nested <- cubble::face_spatial(x)
  }

  # Retrieve geometry: cubble keeps it as sf geometry on the nested face.
  nested_df <- sf::st_as_sf(as.data.frame(x_nested))
  geom_col  <- attr(nested_df, "sf_column")

  if (is.null(geom_col) || !geom_col %in% names(nested_df))
    stop("The cubble nested face does not carry an sf geometry column. ",
         "Make sure the cubble was created from an sf object.")

  # Keep only key + geometry for the join
  spatial_lookup <- nested_df[, c(key_col, geom_col), drop = FALSE]

  # ---- Join geometry onto temporal face -----------------------------------
  df_long[[key_col]] <- as.character(df_long[[key_col]])
  spatial_lookup[[key_col]] <- as.character(spatial_lookup[[key_col]])

  merged <- merge(df_long, spatial_lookup, by = key_col, all.x = TRUE)

  if (anyNA(sf::st_geometry(sf::st_as_sf(merged,
                                          sf_column_name = geom_col))))
    warning("Some rows in the temporal face have no matching geometry after ",
            "the join. Check that key_col ('", key_col, "') values match ",
            "between spatial and temporal faces.")

  # ---- Build sftime -------------------------------------------------------
  merged_sf <- sf::st_as_sf(merged, sf_column_name = geom_col)

  sftime::st_as_sftime(
    merged_sf,
    time_column_name = time_col,
    sf_column_name   = geom_col
  )
}


# ---------------------------------------------------------------------------
# Unified .to_stsdf / .to_stfdf  (now cubble-aware)
# ---------------------------------------------------------------------------

# Convert any ST*DF, sftime, or cubble object to STSDF.
.to_stsdf <- function(x) {
  if (inherits(x, "STSDF")) return(x)
  if (inherits(x, "STFDF")) return(as(x, "STSDF"))
  if (inherits(x, "STIDF")) return(as(as(x, "STFDF"), "STSDF"))
  if (inherits(x, "sftime")) return(as(as(as(x, "STIDF"), "STFDF"), "STSDF"))
  if (.is_cubble(x))         return(.to_stsdf(cubble_to_sftime(x)))
  stop("Cannot coerce object of class '", class(x)[1L], "' to STSDF. ",
       "Supported classes: STFDF, STSDF, STIDF, sftime, cubble_df.")
}

# Convert any ST*DF, sftime, or cubble object to STFDF.
.to_stfdf <- function(x) {
  if (inherits(x, "STFDF")) return(x)
  if (inherits(x, "STSDF")) return(as(x, "STFDF"))
  if (inherits(x, "STIDF")) return(as(x, "STFDF"))
  if (inherits(x, "sftime")) return(as(as(x, "STIDF"), "STFDF"))
  if (.is_cubble(x))         return(.to_stfdf(cubble_to_sftime(x)))
  stop("Cannot coerce object of class '", class(x)[1L], "' to STFDF. ",
       "Supported classes: STFDF, STSDF, STIDF, sftime, cubble_df.")
}

# Extract the numeric time values (seconds since epoch or similar) from
# a spacetime ST* object's @time slot.
.extract_time_numeric <- function(x) {
  as.numeric(zoo::index(x@time))
}

# Extract all pairwise (spacelag, timelag) for the observations in an STSDF.
# Returns a data.frame with columns spacelag and timelag (n_obs^2 rows).
.obs_lag_pairs <- function(stsdf) {
  sp_idx <- stsdf@index[, 1]
  t_idx  <- stsdf@index[, 2]

  all_sp_coords <- sp::coordinates(stsdf@sp)
  all_t_vals    <- .extract_time_numeric(stsdf)

  sp_dist_full <- as.matrix(stats::dist(all_sp_coords))
  t_dist_full  <- abs(outer(all_t_vals, all_t_vals, "-"))

  sp_dists_obs <- sp_dist_full[sp_idx, sp_idx]
  t_dists_obs  <- t_dist_full[t_idx,  t_idx]

  data.frame(
    spacelag = as.vector(sp_dists_obs),
    timelag  = as.vector(t_dists_obs)
  )
}
