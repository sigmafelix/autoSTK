# param_space.R — semantically-aware parameter bounds for vgmST models.
#
# The key improvement over the existing 0.25x / 1.5x heuristic:
# parameter names returned by gstat::extractPar() are used to classify each
# parameter as a sill/psill, nugget, range, anisotropy ratio, or k-factor,
# and each class gets bounds derived from the empirical variogram statistics.

#' Compute parameter bounds for vgmST model optimisation
#'
#' @param model_template A \code{vgmST} object (initial guess).
#' @param stva_emp Empirical ST variogram (\code{StVariogram}/data.frame) from
#'   \code{setSTI} or \code{gstat::variogramST}.
#' @param sill_scale Numeric (2). Upper bound multiplier applied to
#'   \code{max(stva_emp$gamma)} for sill/psill/nugget parameters.
#' @param range_scale Numeric (3). Upper bound multiplier applied to the
#'   maximum spatial lag for range parameters.
#' @param ani_scale Numeric (20). Controls the breadth of the stAni search
#'   interval relative to the data extent.
#' @return A named list with elements \code{lower} and \code{upper} (both
#'   named numeric vectors matching \code{gstat::extractPar(model_template)}).
#' @export
st_param_bounds <- function(model_template, stva_emp,
                             sill_scale  = 2.0,
                             range_scale = 3.0,
                             ani_scale   = 20.0) {
  init_par  <- extractPar(model_template)
  par_names <- names(init_par)

  maxgamma <- max(stva_emp$gamma, na.rm = TRUE)
  maxspl   <- max(stva_emp$spacelag, na.rm = TRUE)
  maxtl    <- as.numeric(max(stva_emp$timelag, na.rm = TRUE))

  lower <- setNames(numeric(length(init_par)), par_names)
  upper <- setNames(numeric(length(init_par)), par_names)

  for (i in seq_along(init_par)) {
    nm <- if (!is.null(par_names)) par_names[i] else ""
    v  <- init_par[i]

    if (grepl("sill|psill|nugget", nm, ignore.case = TRUE)) {
      lower[i] <- 0
      upper[i] <- maxgamma * sill_scale

    } else if (grepl("^range\\.", nm, ignore.case = TRUE) ||
               identical(nm, "range")) {
      # Spatial range
      lower[i] <- 1e-6
      upper[i] <- maxspl * range_scale

    } else if (grepl("^range\\.t", nm, ignore.case = TRUE)) {
      # Temporal range
      lower[i] <- 1e-6
      upper[i] <- maxtl * range_scale

    } else if (grepl("^(stAni|anis)", nm, ignore.case = TRUE) ||
               identical(nm, "anis")) {
      # Anisotropy ratio: space units per time unit
      # Lower: at least a small positive ratio
      # Upper: allow a wide range up to ani_scale * max_space / max_time
      lower[i] <- max(maxspl / (maxtl * ani_scale), 1e-6)
      upper[i] <- maxspl * ani_scale / max(maxtl, 1e-6)

    } else {
      # k (productSum) or unknown: relative bounds from initial value
      if (v > 0) {
        lower[i] <- v * 0.01
        upper[i] <- v * 5.0
      } else {
        lower[i] <- 1e-6
        upper[i] <- maxgamma * sill_scale
      }
    }
  }

  # Safety: ensure lower >= 0 and upper > lower
  lower <- pmax(lower, 0)
  upper <- pmax(upper, lower + 1e-6)

  list(lower = lower, upper = upper)
}
