# Print method for stanreg objects

The `print` method for stanreg objects displays a compact summary of the
fitted model. See the **Details** section below for descriptions of the
different components of the printed output. For additional summary
statistics and diagnostics use the
[`summary`](https://mc-stan.org/rstanarm/reference/summary.stanreg.md)
method.

## Usage

``` r
# S3 method for class 'stanreg'
print(x, digits = 1, detail = TRUE, ...)

# S3 method for class 'stanmvreg'
print(x, digits = 3, ...)
```

## Arguments

- x:

  A fitted model object returned by one of the rstanarm modeling
  functions. See
  [`stanreg-objects`](https://mc-stan.org/rstanarm/reference/stanreg-objects.md).

- digits:

  Number of digits to use for formatting numbers.

- detail:

  Logical, defaulting to `TRUE`. If `FALSE` a more minimal summary is
  printed consisting only of the parameter estimates.

- ...:

  Ignored.

## Value

Returns `x`, invisibly.

## Details

### Point estimates

Regardless of the estimation algorithm, point estimates are medians
computed from simulations. For models fit using MCMC (`"sampling"`) the
posterior sample is used. For optimization (`"optimizing"`), the
simulations are generated from the asymptotic Gaussian sampling
distribution of the parameters. For the `"meanfield"` and `"fullrank"`
variational approximations, draws from the variational approximation to
the posterior are used. In all cases, the point estimates reported are
the same as the values returned by
[`coef`](https://mc-stan.org/rstanarm/reference/stanreg-methods.md).

### Uncertainty estimates (MAD_SD)

The standard deviations reported (labeled `MAD_SD` in the print output)
are computed from the same set of draws described above and are
proportional to the median absolute deviation
([`mad`](https://rdrr.io/r/stats/mad.html)) from the median. Compared to
the raw posterior standard deviation, the MAD_SD will be more robust for
long-tailed distributions. These are the same as the values returned by
[`se`](https://mc-stan.org/rstanarm/reference/stanreg-methods.md).

### Additional output

- For GLMs with group-specific terms (see
  [`stan_glmer`](https://mc-stan.org/rstanarm/reference/stan_glmer.md))
  the printed output also shows point estimates of the standard
  deviations of the group effects (and correlations if there are both
  intercept and slopes that vary by group).

- For analysis of variance models (see
  [`stan_aov`](https://mc-stan.org/rstanarm/reference/stan_lm.md))
  models, an ANOVA-like table is also displayed.

- For joint longitudinal and time-to-event (see
  [`stan_jm`](https://mc-stan.org/rstanarm/reference/stan_jm.md)) models
  the estimates are presented separately for each of the distinct
  submodels.

## See also

[`summary.stanreg`](https://mc-stan.org/rstanarm/reference/summary.stanreg.md),
[`stanreg-methods`](https://mc-stan.org/rstanarm/reference/stanreg-methods.md)
