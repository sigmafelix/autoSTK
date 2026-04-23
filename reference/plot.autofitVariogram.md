# Plot the automatically fitted variogram

Plot the automatically fitted variogram

## Usage

``` r
# S3 method for class 'autofitVariogram'
plot(
  x,
  plotit = TRUE,
  title = "Experimental variogram and fitted variogram model",
  ...
)
```

## Arguments

- x:

  A result object of autofitVariogram.

- plotit:

  boolean. Print graph or not.

- title:

  character. the title of the plot.

- ...:

  passed to xyplot

## Value

A lattice::xyplot object.
