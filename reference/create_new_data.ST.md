# Generate a new spatiotemporal points for the spatiotemporal prediction and interpolation

Generate a new spatiotemporal points for the spatiotemporal prediction
and interpolation

## Usage

``` r
create_new_data.ST(obj, form, gen_mode = "chull", npoints = 10000, forward = 6)
```

## Arguments

- obj:

  a sftime object.

- form:

  formula.

- gen_mode:

  character. The type of shape by the surface is generated. One of
  'rect' (rectangular) and 'chull' (convex hull).

- npoints:

  integer. the number of points that will be generated

- forward:

  integer. the time length of the data will generate ahead of the last
  time point of the input data. If NULL is passed, the spatiotemporal
  interpolation mode in obj will be conducted.

## Value

A sftime object.
