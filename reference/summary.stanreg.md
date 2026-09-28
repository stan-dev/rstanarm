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
#> b[(Intercept) herd:1]   0.62300438 0.01620054 0.4588130  0.03796704  0.62611637
#> b[(Intercept) herd:2]  -0.36481538 0.01489670 0.4327220 -0.90822480 -0.35023927
#> b[(Intercept) herd:3]   0.37725127 0.01568539 0.3696453 -0.08760737  0.37827066
#> b[(Intercept) herd:4]   0.03678426 0.01566984 0.4916136 -0.58336494  0.02393670
#> b[(Intercept) herd:5]  -0.26173939 0.01476787 0.4332863 -0.83041946 -0.22689161
#> b[(Intercept) herd:6]  -0.48400034 0.01559027 0.4606396 -1.09863009 -0.45470874
#> b[(Intercept) herd:7]   0.93232301 0.01603749 0.4173846  0.41634352  0.93039901
#> b[(Intercept) herd:8]   0.53408380 0.02178227 0.5098291 -0.09138152  0.52287568
#> b[(Intercept) herd:9]  -0.25936885 0.01859386 0.5542501 -0.93166615 -0.23831961
#> b[(Intercept) herd:10] -0.62692607 0.01603693 0.4352790 -1.20094641 -0.60839733
#> b[(Intercept) herd:11] -0.17262225 0.01698866 0.4154737 -0.72043939 -0.16538408
#> b[(Intercept) herd:12] -0.08228036 0.01522644 0.5096514 -0.73990528 -0.06503475
#> b[(Intercept) herd:13] -0.82538180 0.01737403 0.4745142 -1.47618453 -0.80554091
#> b[(Intercept) herd:14]  1.02188545 0.01670857 0.4402157  0.45192065  1.02659147
#> b[(Intercept) herd:15] -0.63270314 0.01345009 0.4438263 -1.22398627 -0.60756536
#>                                90% n_eff      Rhat
#> b[(Intercept) herd:1]   1.20556247   802 0.9982976
#> b[(Intercept) herd:2]   0.18783005   844 0.9993093
#> b[(Intercept) herd:3]   0.85919710   555 0.9990585
#> b[(Intercept) herd:4]   0.65275493   984 0.9986838
#> b[(Intercept) herd:5]   0.25863332   861 0.9993085
#> b[(Intercept) herd:6]   0.08783862   873 1.0012477
#> b[(Intercept) herd:7]   1.45930283   677 0.9981217
#> b[(Intercept) herd:8]   1.18441094   548 1.0004232
#> b[(Intercept) herd:9]   0.42257927   889 1.0004938
#> b[(Intercept) herd:10] -0.11416995   737 0.9995577
#> b[(Intercept) herd:11]  0.35144067   598 1.0020647
#> b[(Intercept) herd:12]  0.51117888  1120 0.9988075
#> b[(Intercept) herd:13] -0.25850026   746 1.0002000
#> b[(Intercept) herd:14]  1.61145793   694 0.9987592
#> b[(Intercept) herd:15] -0.07991711  1089 0.9989376
```
