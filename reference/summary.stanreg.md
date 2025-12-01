# Summary method for stanreg objects

Summaries of parameter estimates and MCMC convergence diagnostics (Monte
Carlo error, effective sample size, Rhat).

## Usage

``` r
# S3 method for class 'stanreg'
summary(
  object,
  pars = NULL,
  regex_pars = NULL,
  probs = c(0.1, 0.5, 0.9),
  ...,
  digits = 1
)

# S3 method for class 'summary.stanreg'
print(x, digits = max(1, attr(x, "print.digits")), ...)

# S3 method for class 'summary.stanreg'
as.data.frame(x, ...)

# S3 method for class 'stanmvreg'
summary(object, pars = NULL, regex_pars = NULL, probs = NULL, ..., digits = 3)

# S3 method for class 'summary.stanmvreg'
print(x, digits = max(1, attr(x, "print.digits")), ...)
```

## Arguments

- object:

  A fitted model object returned by one of the rstanarm modeling
  functions. See
  [`stanreg-objects`](https://mc-stan.org/rstanarm/reference/stanreg-objects.md).

- pars:

  An optional character vector specifying a subset of parameters to
  display. Parameters can be specified by name or several shortcuts can
  be used. Using `pars="beta"` will restrict the displayed parameters to
  only the regression coefficients (without the intercept). `"alpha"`
  can also be used as a shortcut for `"(Intercept)"`. If the model has
  varying intercepts and/or slopes they can be selected using
  `pars = "varying"`.

  In addition, for `stanmvreg` objects there are some additional
  shortcuts available. Using `pars = "long"` will display the parameter
  estimates for the longitudinal submodels only (excluding
  group-specific pparameters, but including auxiliary parameters). Using
  `pars = "event"` will display the parameter estimates for the event
  submodel only, including any association parameters. Using
  `pars = "assoc"` will display only the association parameters. Using
  `pars = "fixef"` will display all fixed effects, but not the random
  effects or the auxiliary parameters. `pars` and `regex_pars` are set
  to `NULL` then all fixed effect regression coefficients are selected,
  as well as any auxiliary parameters and the log posterior.

  If `pars` is `NULL` all parameters are selected for a `stanreg`
  object, while for a `stanmvreg` object all fixed effect regression
  coefficients are selected as well as any auxiliary parameters and the
  log posterior. See **Examples**.

- regex_pars:

  An optional character vector of [regular
  expressions](https://rdrr.io/r/base/grep.html) to use for parameter
  selection. `regex_pars` can be used in place of `pars` or in addition
  to `pars`. Currently, all functions that accept a `regex_pars`
  argument ignore it for models fit using optimization.

- probs:

  For models fit using MCMC or one of the variational algorithms, an
  optional numeric vector of probabilities passed to
  [`quantile`](https://rdrr.io/r/stats/quantile.html).

- ...:

  Currently ignored.

- digits:

  Number of digits to use for formatting numbers when printing. When
  calling `summary`, the value of digits is stored as the
  `"print.digits"` attribute of the returned object.

- x:

  An object of class `"summary.stanreg"`.

## Value

The `summary` method returns an object of class `"summary.stanreg"` (or
`"summary.stanmvreg"`, inheriting `"summary.stanreg"`), which is a
matrix of summary statistics and diagnostics, with attributes storing
information for use by the `print` method. The `print` method for
`summary.stanreg` or `summary.stanmvreg` objects is called for its side
effect and just returns its input. The `as.data.frame` method for
`summary.stanreg` objects converts the matrix to a data.frame,
preserving row and column names but dropping the `print`-related
attributes.

## Details

### mean_PPD diagnostic

Summary statistics are also reported for `mean_PPD`, the sample average
posterior predictive distribution of the outcome. This is useful as a
quick diagnostic. A useful heuristic is to check if `mean_PPD` is
plausible when compared to `mean(y)`. If it is plausible then this does
*not* mean that the model is good in general (only that it can reproduce
the sample mean), however if `mean_PPD` is implausible then it is a sign
that something is wrong (severe model misspecification, problems with
the data, computational issues, etc.).

## See also

[`prior_summary`](https://mc-stan.org/rstanarm/reference/prior_summary.stanreg.md)
to extract or print a summary of the priors used for a particular model.

## Examples

``` r
if (.Platform$OS.type != "windows" || .Platform$r_arch != "i386") {
if (!exists("example_model")) example(example_model) 
summary(example_model, probs = c(0.1, 0.9))

# These produce the same output for this example, 
# but the second method can be used for any model
summary(example_model, pars = c("(Intercept)", "size", 
                                paste0("period", 2:4)))
summary(example_model, pars = c("alpha", "beta"))

# Only show parameters varying by group
summary(example_model, pars = "varying")
as.data.frame(summary(example_model, pars = "varying"))
}
#>                               mean       mcse        sd         10%         50%
#> b[(Intercept) herd:1]   0.62722341 0.01608742 0.4415459  0.07595455  0.62715031
#> b[(Intercept) herd:2]  -0.36427769 0.01525883 0.4424714 -0.94483006 -0.33901142
#> b[(Intercept) herd:3]   0.38045628 0.01263177 0.3766349 -0.07852149  0.37502550
#> b[(Intercept) herd:4]   0.03841076 0.01803683 0.4905975 -0.58000948  0.04103359
#> b[(Intercept) herd:5]  -0.26015881 0.01415271 0.4160125 -0.80240461 -0.23959978
#> b[(Intercept) herd:6]  -0.45988461 0.01559916 0.4388333 -1.02158661 -0.43434658
#> b[(Intercept) herd:7]   0.92618148 0.01426487 0.4360411  0.40131569  0.89970975
#> b[(Intercept) herd:8]   0.51032101 0.01925821 0.5221492 -0.13943996  0.52135539
#> b[(Intercept) herd:9]  -0.25893282 0.01666986 0.5406258 -0.90647896 -0.24139284
#> b[(Intercept) herd:10] -0.62914073 0.01481049 0.4468560 -1.22576496 -0.60609808
#> b[(Intercept) herd:11] -0.14637179 0.01431972 0.4129541 -0.71249718 -0.11485555
#> b[(Intercept) herd:12] -0.04363117 0.01846334 0.5122429 -0.69627174 -0.02758330
#> b[(Intercept) herd:13] -0.78627355 0.01656960 0.4573965 -1.36763552 -0.75501245
#> b[(Intercept) herd:14]  1.00700584 0.01707785 0.4423126  0.43654877  0.98415992
#> b[(Intercept) herd:15] -0.60198255 0.01381385 0.4725656 -1.22122313 -0.57879172
#>                                90% n_eff      Rhat
#> b[(Intercept) herd:1]   1.18994683   753 1.0005647
#> b[(Intercept) herd:2]   0.18168238   841 0.9992947
#> b[(Intercept) herd:3]   0.86572288   889 0.9996979
#> b[(Intercept) herd:4]   0.64919383   740 0.9988455
#> b[(Intercept) herd:5]   0.24399929   864 0.9985012
#> b[(Intercept) herd:6]   0.06249890   791 1.0008678
#> b[(Intercept) herd:7]   1.52063980   934 0.9988557
#> b[(Intercept) herd:8]   1.15757378   735 1.0005324
#> b[(Intercept) herd:9]   0.41947464  1052 0.9983954
#> b[(Intercept) herd:10] -0.09920854   910 0.9994387
#> b[(Intercept) herd:11]  0.35927443   832 0.9998923
#> b[(Intercept) herd:12]  0.59755834   770 1.0019334
#> b[(Intercept) herd:13] -0.25753682   762 1.0027890
#> b[(Intercept) herd:14]  1.59355623   671 1.0003160
#> b[(Intercept) herd:15] -0.01508412  1170 0.9985240
```
