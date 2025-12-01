# Predict method for stanreg objects

This method is primarily intended to be used only for models fit using
optimization. For models fit using MCMC or one of the variational
approximations, see
[`posterior_predict`](https://mc-stan.org/rstanarm/reference/posterior_predict.stanreg.md).

## Usage

``` r
# S3 method for class 'stanreg'
predict(
  object,
  ...,
  newdata = NULL,
  type = c("link", "response"),
  se.fit = FALSE
)
```

## Arguments

- object:

  A fitted model object returned by one of the rstanarm modeling
  functions. See
  [`stanreg-objects`](https://mc-stan.org/rstanarm/reference/stanreg-objects.md).

- ...:

  Ignored.

- newdata:

  Optionally, a data frame in which to look for variables with which to
  predict. If omitted, the model matrix is used.

- type:

  The type of prediction. The default `'link'` is on the scale of the
  linear predictors; the alternative `'response'` is on the scale of the
  response variable.

- se.fit:

  A logical scalar indicating if standard errors should be returned. The
  default is `FALSE`.

## Value

A vector if `se.fit` is `FALSE` and a list if `se.fit` is `TRUE`.

## See also

[`posterior_predict`](https://mc-stan.org/rstanarm/reference/posterior_predict.stanreg.md)
