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
#>                               mean       mcse        sd         10%        50%
#> b[(Intercept) herd:1]   0.61366803 0.01684104 0.4358830  0.06018384  0.5997267
#> b[(Intercept) herd:2]  -0.39076127 0.01559322 0.4326705 -0.95744794 -0.3628813
#> b[(Intercept) herd:3]   0.37405608 0.01444953 0.3745823 -0.11806346  0.3760999
#> b[(Intercept) herd:4]   0.02385185 0.01824713 0.4855213 -0.56488796  0.0372756
#> b[(Intercept) herd:5]  -0.26762793 0.01495113 0.4249044 -0.80254997 -0.2595050
#> b[(Intercept) herd:6]  -0.47048746 0.01630673 0.4623683 -1.03788228 -0.4603273
#> b[(Intercept) herd:7]   0.92120225 0.01671739 0.4424558  0.37613467  0.8979430
#> b[(Intercept) herd:8]   0.51262524 0.01892408 0.5139155 -0.11905568  0.5028575
#> b[(Intercept) herd:9]  -0.28377136 0.02003379 0.5725331 -0.98793800 -0.2832377
#> b[(Intercept) herd:10] -0.66338033 0.01346737 0.4170828 -1.21194037 -0.6452184
#> b[(Intercept) herd:11] -0.17299306 0.01554727 0.4223052 -0.70529282 -0.1589002
#> b[(Intercept) herd:12] -0.08972901 0.01820619 0.5175226 -0.72362116 -0.0974024
#> b[(Intercept) herd:13] -0.83293321 0.01577726 0.4797535 -1.48537375 -0.8033615
#> b[(Intercept) herd:14]  1.00932780 0.01561415 0.4575590  0.43742382  0.9785853
#> b[(Intercept) herd:15] -0.64913880 0.01602881 0.4970131 -1.29509249 -0.6365910
#>                                90% n_eff      Rhat
#> b[(Intercept) herd:1]   1.18636950   670 1.0016340
#> b[(Intercept) herd:2]   0.11671173   770 0.9990603
#> b[(Intercept) herd:3]   0.87124036   672 1.0000762
#> b[(Intercept) herd:4]   0.64587534   708 0.9982303
#> b[(Intercept) herd:5]   0.30776694   808 0.9993399
#> b[(Intercept) herd:6]   0.09856199   804 1.0017466
#> b[(Intercept) herd:7]   1.51077314   700 0.9990832
#> b[(Intercept) herd:8]   1.19046545   737 1.0001959
#> b[(Intercept) herd:9]   0.39624101   817 0.9994001
#> b[(Intercept) herd:10] -0.16268401   959 0.9982412
#> b[(Intercept) herd:11]  0.37308955   738 0.9988509
#> b[(Intercept) herd:12]  0.61418839   808 1.0007142
#> b[(Intercept) herd:13] -0.23597422   925 0.9982607
#> b[(Intercept) herd:14]  1.58811402   859 0.9992650
#> b[(Intercept) herd:15] -0.04035512   961 1.0003634
```
