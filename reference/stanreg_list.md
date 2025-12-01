# Create lists of fitted model objects, combine them, or append new models to existing lists of models.

Create lists of fitted model objects, combine them, or append new models
to existing lists of models.

## Usage

``` r
stanreg_list(..., model_names = NULL)

stanmvreg_list(..., model_names = NULL)

stanjm_list(..., model_names = NULL)

# S3 method for class 'stanreg_list'
print(x, ...)
```

## Arguments

- ...:

  Objects to combine into a `"stanreg_list"`, `"stanmvreg_list"`, or
  `"stanjm_list"`. Can be fitted model objects, existing `"stan*_list"`
  objects to combine, or one existing `"stan*_list"` object followed by
  fitted model objects to append to the list.

- model_names:

  Optionally, a character vector of model names. If not specified then
  the names are inferred from the name of the objects passed in via
  `...`. These model names are used, for example, when printing the
  results of the `loo_compare.stanreg_list` and
  `loo_model_weights.stanreg_list` methods.

- x:

  The object to print.

## Value

A list of class `"stanreg_list"`, `"stanmvreg_list"`, or
`"stanjm_list"`, containing the fitted model objects and some metadata
stored as attributes.

## See also

[`loo_model_weights`](https://mc-stan.org/rstanarm/reference/loo.stanreg.md)
for usage of `stanreg_list`.
