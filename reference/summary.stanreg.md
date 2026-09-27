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
#> b[(Intercept) herd:1]   0.63476660 0.02111340 0.4430841  0.08082475  0.62485778
#> b[(Intercept) herd:2]  -0.36108534 0.01944247 0.4543615 -0.93856571 -0.33217068
#> b[(Intercept) herd:3]   0.39225026 0.01660559 0.3709991 -0.09050650  0.41577756
#> b[(Intercept) herd:4]   0.05745181 0.02079834 0.4996894 -0.59872409  0.07318090
#> b[(Intercept) herd:5]  -0.24667459 0.01714514 0.4263930 -0.80419571 -0.22214783
#> b[(Intercept) herd:6]  -0.47564865 0.01650612 0.4401350 -1.02829154 -0.47863528
#> b[(Intercept) herd:7]   0.93534077 0.02486005 0.4530910  0.37907438  0.90319516
#> b[(Intercept) herd:8]   0.52173407 0.02212551 0.5389708 -0.17825159  0.52102931
#> b[(Intercept) herd:9]  -0.27611669 0.02164574 0.5320603 -0.95877054 -0.26097761
#> b[(Intercept) herd:10] -0.63113193 0.01714923 0.4464421 -1.21605263 -0.61235169
#> b[(Intercept) herd:11] -0.13922484 0.01876186 0.4139470 -0.65302783 -0.14146057
#> b[(Intercept) herd:12] -0.06108168 0.02016889 0.5503159 -0.76784490 -0.06076016
#> b[(Intercept) herd:13] -0.83508068 0.02172721 0.5044752 -1.49593445 -0.76716353
#> b[(Intercept) herd:14]  1.01052351 0.02523687 0.4708093  0.41451925  0.98582412
#> b[(Intercept) herd:15] -0.61524326 0.01422651 0.4498592 -1.20843121 -0.57588912
#>                                90% n_eff      Rhat
#> b[(Intercept) herd:1]   1.22749000   440 1.0109947
#> b[(Intercept) herd:2]   0.21543789   546 1.0061973
#> b[(Intercept) herd:3]   0.83565863   499 1.0015294
#> b[(Intercept) herd:4]   0.70046845   577 1.0038669
#> b[(Intercept) herd:5]   0.28415157   618 0.9989187
#> b[(Intercept) herd:6]   0.07958805   711 0.9990246
#> b[(Intercept) herd:7]   1.51805591   332 1.0109432
#> b[(Intercept) herd:8]   1.18438214   593 1.0010381
#> b[(Intercept) herd:9]   0.36176788   604 1.0051933
#> b[(Intercept) herd:10] -0.07264965   678 0.9995690
#> b[(Intercept) herd:11]  0.36728593   487 1.0008207
#> b[(Intercept) herd:12]  0.61987868   744 0.9990974
#> b[(Intercept) herd:13] -0.25235638   539 1.0002926
#> b[(Intercept) herd:14]  1.62110902   348 1.0048214
#> b[(Intercept) herd:15] -0.07666799  1000 0.9988275
```
