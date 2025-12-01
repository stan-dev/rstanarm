# Methods for stanreg objects

The methods documented on this page are actually some of the least
important methods defined for
[stanreg](https://mc-stan.org/rstanarm/reference/stanreg-objects.md)
objects. The most important methods are documented separately, each with
its own page. Links to those pages are provided in the **See Also**
section, below.

## Usage

``` r
# S3 method for class 'stanmvreg'
nobs(object, ...)

# S3 method for class 'stanreg'
coef(object, ...)

# S3 method for class 'stanreg'
confint(object, parm, level = 0.95, ...)

# S3 method for class 'stanreg'
fitted(object, ...)

# S3 method for class 'stanreg'
nobs(object, ...)

# S3 method for class 'stanreg'
residuals(object, ...)

# S3 method for class 'stanreg'
se(object, ...)

# S3 method for class 'stanreg'
update(object, formula., ..., evaluate = TRUE)

# S3 method for class 'stanreg'
vcov(object, correlation = FALSE, ...)

# S3 method for class 'stanreg'
fixef(object, ...)

# S3 method for class 'stanreg'
ngrps(object, ...)

# S3 method for class 'stanreg'
nsamples(object, ...)

# S3 method for class 'stanreg'
ranef(object, ...)

# S3 method for class 'stanreg'
sigma(object, ...)

# S3 method for class 'stanreg'
VarCorr(x, sigma = 1, ...)
```

## Arguments

- object, x:

  A fitted model object returned by one of the rstanarm modeling
  functions. See
  [`stanreg-objects`](https://mc-stan.org/rstanarm/reference/stanreg-objects.md).

- ...:

  Ignored, except by the `update` method. See
  [`update`](https://rdrr.io/r/stats/update.html).

- parm:

  For `confint`, an optional character vector of parameter names.

- level:

  For `confint`, a scalar between \\0\\ and \\1\\ indicating the
  confidence level to use.

- formula., evaluate:

  See [`update`](https://rdrr.io/r/stats/update.html).

- correlation:

  For `vcov`, if `FALSE` (the default) the covariance matrix is
  returned. If `TRUE`, the correlation matrix is returned instead.

- sigma:

  Ignored (included for compatibility with
  [`VarCorr`](https://rdrr.io/pkg/nlme/man/VarCorr.html)).

## Details

The methods documented on this page are similar to the methods defined
for objects of class 'lm', 'glm', 'glmer', etc. However there are a few
key differences:

- `residuals`:

  Residuals are *always* of type `"response"` (not `"deviance"`
  residuals or any other type). However, in the case of
  [`stan_polr`](https://mc-stan.org/rstanarm/reference/stan_polr.md)
  with more than two response categories, the residuals are the
  difference between the latent utility and its linear predictor.

- `coef`:

  Medians are used for point estimates. See the *Point estimates*
  section in
  [`print.stanreg`](https://mc-stan.org/rstanarm/reference/print.stanreg.md)
  for more details.

- `se`:

  The `se` function returns standard errors based on
  [`mad`](https://rdrr.io/r/stats/mad.html). See the *Uncertainty
  estimates* section in
  [`print.stanreg`](https://mc-stan.org/rstanarm/reference/print.stanreg.md)
  for more details.

- `confint`:

  For models fit using optimization, confidence intervals are returned
  via a call to
  [`confint.default`](https://rdrr.io/r/stats/confint.html). If
  `algorithm` is `"sampling"`, `"meanfield"`, or `"fullrank"`, the
  `confint` will throw an error because the
  [`posterior_interval`](https://mc-stan.org/rstanarm/reference/posterior_interval.stanreg.md)
  function should be used to compute Bayesian uncertainty intervals.

- `nsamples`:

  The number of draws from the posterior distribution obtained

## See also

- The
  [`print`](https://mc-stan.org/rstanarm/reference/print.stanreg.md),
  [`summary`](https://mc-stan.org/rstanarm/reference/summary.stanreg.md),
  and
  [`prior_summary`](https://mc-stan.org/rstanarm/reference/prior_summary.stanreg.md)
  methods for stanreg objects for information on the fitted model.

- [`launch_shinystan`](https://mc-stan.org/rstanarm/reference/launch_shinystan.stanreg.md)
  to use the ShinyStan GUI to explore a fitted rstanarm model.

- The [`plot`](https://mc-stan.org/rstanarm/reference/plot.stanreg.md)
  method to plot estimates and diagnostics.

- The
  [`pp_check`](https://mc-stan.org/rstanarm/reference/pp_check.stanreg.md)
  method for graphical posterior predictive checking.

- The
  [`posterior_predict`](https://mc-stan.org/rstanarm/reference/posterior_predict.stanreg.md)
  and
  [`predictive_error`](https://mc-stan.org/rstanarm/reference/predictive_error.stanreg.md)
  methods for predictions and predictive errors.

- The
  [`posterior_interval`](https://mc-stan.org/rstanarm/reference/posterior_interval.stanreg.md)
  and
  [`predictive_interval`](https://mc-stan.org/rstanarm/reference/predictive_interval.stanreg.md)
  methods for uncertainty intervals for model parameters and
  predictions.

- The [`loo`](https://mc-stan.org/rstanarm/reference/loo.stanreg.md),
  [`kfold`](https://mc-stan.org/rstanarm/reference/kfold.stanreg.md),
  and
  [`log_lik`](https://mc-stan.org/rstanarm/reference/log_lik.stanreg.md)
  methods for leave-one-out or K-fold cross-validation, model
  comparison, and computing the log-likelihood of (possibly new) data.

- The
  [`as.matrix`](https://mc-stan.org/rstanarm/reference/as.matrix.stanreg.md),
  `as.data.frame`, and `as.array` methods to access posterior draws.
