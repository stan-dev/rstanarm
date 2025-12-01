# Methods for stanmvreg objects

S3 methods for
[stanmvreg](https://mc-stan.org/rstanarm/reference/stanreg-objects.md)
objects. There are also several methods (listed in **See Also**, below)
with their own individual help pages. The main difference between these
methods and the
[stanreg](https://mc-stan.org/rstanarm/reference/stanreg-methods.md)
methods is that the methods described here generally include an
additional argument `m` which allows the user to specify which submodel
they wish to return the result for. If the argument `m` is set to `NULL`
then the result will generally be a named list with each element of the
list containing the result for one of the submodels.

## Usage

``` r
# S3 method for class 'stanmvreg'
coef(object, m = NULL, ...)

# S3 method for class 'stanmvreg'
fitted(object, m = NULL, ...)

# S3 method for class 'stanmvreg'
residuals(object, m = NULL, ...)

# S3 method for class 'stanmvreg'
se(object, m = NULL, ...)

# S3 method for class 'stanmvreg'
formula(x, fixed.only = FALSE, random.only = FALSE, m = NULL, ...)

# S3 method for class 'stanmvreg'
update(object, formula., ..., evaluate = TRUE)

# S3 method for class 'stanjm'
update(object, formulaLong., formulaEvent., ..., evaluate = TRUE)

# S3 method for class 'stanmvreg'
fixef(object, m = NULL, remove_stub = TRUE, ...)

# S3 method for class 'stanmvreg'
ngrps(object, ...)

# S3 method for class 'stanmvreg'
ranef(object, m = NULL, ...)

# S3 method for class 'stanmvreg'
sigma(object, m = NULL, ...)
```

## Arguments

- object, x:

  A fitted model object returned by one of the multivariate rstanarm
  modelling functions. See
  [`stanreg-objects`](https://mc-stan.org/rstanarm/reference/stanreg-objects.md).

- m:

  Integer specifying the number or name of the submodel

- ...:

  Ignored, except by the `update` method. See
  [`update`](https://rdrr.io/r/stats/update.html).

- fixed.only:

  A logical specifying whether to only retain the fixed effect part of
  the longitudinal submodel formulas

- random.only:

  A logical specifying whether to only retain the random effect part of
  the longitudinal submodel formulas

- formula.:

  An updated formula for the model. For a multivariate model `formula.`
  should be a list of formulas, as described for the `formula` argument
  in
  [`stan_mvmer`](https://mc-stan.org/rstanarm/reference/stan_mvmer.md).

- evaluate:

  See [`update`](https://rdrr.io/r/stats/update.html).

- formulaLong., formulaEvent.:

  An updated formula for the longitudinal or event submodel, when
  `object` was estimated using
  [`stan_jm`](https://mc-stan.org/rstanarm/reference/stan_jm.md). For a
  multivariate joint model `formulaLong.` should be a list of formulas,
  as described for the `formulaLong` argument in
  [`stan_jm`](https://mc-stan.org/rstanarm/reference/stan_jm.md).

- remove_stub:

  Logical specifying whether to remove the string identifying the
  submodel (e.g. `y1|`, `y2|`, `Long1|`, `Long2|`, `Event|`) from each
  of the parameter names.

## Details

Most of these methods are similar to the methods defined for objects of
class 'lm', 'glm', 'glmer', etc. However there are a few exceptions:

- `coef`:

  Medians are used for point estimates. See the *Point estimates*
  section in
  [`print.stanmvreg`](https://mc-stan.org/rstanarm/reference/print.stanreg.md)
  for more details. `coef` returns a list equal to the length of the
  number of submodels. The first elements of the list are the
  coefficients from each of the fitted longitudinal submodels and are
  the same layout as those returned by `coef` method of the lme4
  package, that is, the sum of the random and fixed effects coefficients
  for each explanatory variable for each level of each grouping factor.
  The final element of the returned list is a vector of fixed effect
  coefficients from the event submodel.

- `se`:

  The `se` function returns standard errors based on
  [`mad`](https://rdrr.io/r/stats/mad.html). See the *Uncertainty
  estimates* section in
  [`print.stanmvreg`](https://mc-stan.org/rstanarm/reference/print.stanreg.md)
  for more details.

- `confint`:

  Not supplied, since the
  [`posterior_interval`](https://mc-stan.org/rstanarm/reference/posterior_interval.stanreg.md)
  function should be used instead to compute Bayesian uncertainty
  intervals.

- `residuals`:

  Residuals are *always* of type `"response"` (not `"deviance"`
  residuals or any other type).

## See also

- The
  [`print`](https://mc-stan.org/rstanarm/reference/print.stanreg.md),
  [`summary`](https://mc-stan.org/rstanarm/reference/summary.stanreg.md),
  and
  [`prior_summary`](https://mc-stan.org/rstanarm/reference/prior_summary.stanreg.md)
  methods for `stanmvreg` objects for information on the fitted model.

- The [`plot`](https://mc-stan.org/rstanarm/reference/plot.stanreg.md)
  method to plot estimates and diagnostics.

- The
  [`pp_check`](https://mc-stan.org/rstanarm/reference/pp_check.stanreg.md)
  method for graphical posterior predictive checking of the longitudinal
  or glmer submodels.

- The [`ps_check`](https://mc-stan.org/rstanarm/reference/ps_check.md)
  method for graphical posterior predictive checking of the event
  submodel.

- The
  [`posterior_traj`](https://mc-stan.org/rstanarm/reference/posterior_traj.md)
  for predictions for the longitudinal submodel (for models estimated
  using [`stan_jm`](https://mc-stan.org/rstanarm/reference/stan_jm.md)),
  as well as it's associated
  [`plot`](https://mc-stan.org/rstanarm/reference/plot.predict.stanjm.md)
  method.

- The
  [`posterior_survfit`](https://mc-stan.org/rstanarm/reference/posterior_survfit.md)
  for predictions for the event submodel, including so-called "dynamic"
  predictions (for models estimated using
  [`stan_jm`](https://mc-stan.org/rstanarm/reference/stan_jm.md)), as
  well as it's associated
  [`plot`](https://mc-stan.org/rstanarm/reference/plot.survfit.stanjm.md)
  method.

- The
  [`posterior_predict`](https://mc-stan.org/rstanarm/reference/posterior_predict.stanreg.md)
  for predictions for the glmer submodel (for models estimated using
  [`stan_mvmer`](https://mc-stan.org/rstanarm/reference/stan_mvmer.md)).

- The
  [`posterior_interval`](https://mc-stan.org/rstanarm/reference/posterior_interval.stanreg.md)
  for uncertainty intervals for model parameters.

- The [`loo`](https://mc-stan.org/rstanarm/reference/loo.stanreg.md),
  and
  [`log_lik`](https://mc-stan.org/rstanarm/reference/log_lik.stanreg.md)
  methods for leave-one-out model comparison, and computing the
  log-likelihood of (possibly new) data.

- The
  [`as.matrix`](https://mc-stan.org/rstanarm/reference/as.matrix.stanreg.md),
  `as.data.frame`, and `as.array` methods to access posterior draws.

Other S3 methods for stanmvreg objects, which have separate
documentation, including
[`print.stanmvreg`](https://mc-stan.org/rstanarm/reference/print.stanreg.md),
and
[`summary.stanmvreg`](https://mc-stan.org/rstanarm/reference/summary.stanreg.md).

Also
[`posterior_interval`](https://mc-stan.org/rstanarm/reference/posterior_interval.stanreg.md)
for an alternative to `confint`, and `posterior_predict`,
`posterior_traj` and `posterior_survfit` for predictions based on the
fitted joint model.
