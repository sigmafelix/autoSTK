# coerce_st.R — internal helpers for ST class coercion
# All functions are internal (not exported); used throughout the package
# to avoid repeating the conversion chain across files.

# Convert any ST*DF or sftime object to STSDF.
.to_stsdf <- function(x) {
  if (inherits(x, "STSDF")) return(x)
  if (inherits(x, "STFDF")) return(as(x, "STSDF"))
  if (inherits(x, "STIDF")) return(as(as(x, "STFDF"), "STSDF"))
  if (inherits(x, "sftime")) return(as(as(as(x, "STIDF"), "STFDF"), "STSDF"))
  stop("Cannot coerce object of class '", class(x)[1], "' to STSDF")
}

# Convert any ST*DF or sftime object to STFDF.
.to_stfdf <- function(x) {
  if (inherits(x, "STFDF")) return(x)
  if (inherits(x, "STSDF")) return(as(x, "STFDF"))
  if (inherits(x, "STIDF")) return(as(x, "STFDF"))
  if (inherits(x, "sftime")) return(as(as(x, "STIDF"), "STFDF"))
  stop("Cannot coerce object of class '", class(x)[1], "' to STFDF")
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

  sp_dist_full <- as.matrix(dist(all_sp_coords))
  t_dist_full  <- abs(outer(all_t_vals, all_t_vals, "-"))

  sp_dists_obs <- sp_dist_full[sp_idx, sp_idx]
  t_dists_obs  <- t_dist_full[t_idx,  t_idx]

  data.frame(
    spacelag = as.vector(sp_dists_obs),
    timelag  = as.vector(t_dists_obs)
  )
}
