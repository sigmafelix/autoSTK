# Genetic Algorithm Optimiser for ST Variogram Fitting

Uses the GA package's real-valued genetic algorithm to globally search
the WLS objective surface, then refines the solution with a short
L-BFGS-B run. Falls back to a built-in Differential Evolution
(`DE/rand/1/bin`) implementation when GA is not installed.

## Usage

``` r
.optimize_stv_ga(
  stva_emp,
  model_template,
  bounds,
  popSize = 50L,
  maxiter = 200L,
  run = 30L,
  seed = 42L,
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

- popSize:

  Integer. Population size. Default 50.

- maxiter:

  Integer. Maximum number of generations. Default 200.

- run:

  Integer. Stop after this many generations without improvement. Default
  30.

- seed:

  Integer. Random seed for reproducibility. Default 42.

- lbfgsb_refine:

  Logical. Refine the GA solution with L-BFGS-B. Default `TRUE`.

- lbfgsb_maxit:

  Integer. Maximum L-BFGS-B iterations for refinement. Default 500.

## Value

Best-fitting `vgmST` object found.

## See also

[`autofitVariogramST`](https://sigmafelix.github.io/autoSTK/reference/autofitVariogramST.md),
[`st_param_bounds`](https://sigmafelix.github.io/autoSTK/reference/st_param_bounds.md)
