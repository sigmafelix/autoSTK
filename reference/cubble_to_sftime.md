# Convert a cubble object to an sftime object

Handles both nested and temporal faces of a cubble object. The resulting
`sftime` has one row per (location × time) observation with the spatial
geometry attached to every row.

## Usage

``` r
cubble_to_sftime(x, time_col = NULL, key_col = NULL)
```

## Arguments

- x:

  A `cubble_df` object (from the cubble package).

- time_col:

  Character. Name of the time column in the temporal face. If `NULL`
  (default), the function tries to detect it automatically by looking
  for POSIXct/Date columns.

- key_col:

  Character. Name of the site identifier column shared by the spatial
  and temporal faces. Detected automatically from `cubble::key_vars()`
  when `NULL`.

## Value

An `sftime` object with a `geometry` column (sfc_POINT).
