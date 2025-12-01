# Extract X, Y or Z from a stanreg object

Extract X, Y or Z from a stanreg object

## Usage

``` r
get_y(object, ...)

get_x(object, ...)

get_z(object, ...)
```

## Arguments

- object:

  A fitted model object returned by one of the rstanarm modeling
  functions. See
  [`stanreg-objects`](https://mc-stan.org/rstanarm/reference/stanreg-objects.md).

- ...:

  Other arguments passed to methods. For a `stanmvreg` object this can
  be an integer `m` specifying the submodel.

## Value

For `get_x` and `get_z`, a matrix. For `get_y`, either a vector or a
matrix, depending on how the response variable was specified.
