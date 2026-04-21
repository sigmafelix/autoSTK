# Simulated Annealing Optimiser for ST Variogram Fitting

Uses `stats::optim(method = "SANN")` to globally search the WLS
objective surface for a `vgmST` model, then refines the best solution
with a short L-BFGS-B run. Parameters are transformed to an
unconstrained space via a scaled logistic map so that SANN can explore
freely without explicit boundary enforcement.

## Usage

``` r
.optimize_stv_sa(
  stva_emp,
  model_template,
  bounds,
  maxit = 5000L,
  temp = 10,
  tmax = 10L,
  lbfgsb_refine = TRUE,
  lbfgsb_maxit = 500L
)
```

## Arguments

- stva_emp:

  Empirical ST variogram (`StVariogram`).

- model_template:

  Initial `vgmST` model (defines type/structure).

- bounds:

  List with `lower` and `upper` named numeric vectors (from
  [`st_param_bounds`](https://sigmafelix.github.io/autoSTK/reference/st_param_bounds.md)).

- maxit:

  Integer. Maximum number of SANN function evaluations. Default 5000.

- temp:

  Numeric. Initial temperature for the SANN cooling schedule. Default
  10.

- tmax:

  Integer. Number of function evaluations at each temperature level
  before cooling. Default 10.

- lbfgsb_refine:

  Logical. If `TRUE` (default), run a final short L-BFGS-B step from the
  SANN solution.

- lbfgsb_maxit:

  Integer. Maximum L-BFGS-B iterations for the refinement step. Default
  500.

## Value

Best-fitting `vgmST` object found.

## See also

[`autofitVariogramST`](https://sigmafelix.github.io/autoSTK/reference/autofitVariogramST.md),
[`st_param_bounds`](https://sigmafelix.github.io/autoSTK/reference/st_param_bounds.md)
