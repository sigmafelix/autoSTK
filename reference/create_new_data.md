# Create new spatial data for interpolation

Create new spatial data for interpolation

## Usage

``` r
create_new_data(obj, gen_mode = "chull", npoints = 10000, return_class = "sf")
```

## Arguments

- obj:

  A Spatial\*DataFrame.

- gen_mode:

  character. One of 'rect' (rectangular) and 'chull' (convex hull).

- npoints:

  integer. the number of points that will be generated

- return_class:

  character(1). One of 'sf' and 'sp'

## Value

A sf (return_class == 'sf') or SpatialPointsDataFrame (return_class ==
'sp').
