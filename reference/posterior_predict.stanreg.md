# Draw from posterior predictive distribution

The posterior predictive distribution is the distribution of the outcome
implied by the model after using the observed data to update our beliefs
about the unknown parameters in the model. Simulating data from the
posterior predictive distribution using the observed predictors is
useful for checking the fit of the model. Drawing from the posterior
predictive distribution at interesting values of the predictors also
lets us visualize how a manipulation of a predictor affects (a function
of) the outcome(s). With new observations of predictor variables we can
use the posterior predictive distribution to generate predicted
outcomes.

## Usage

``` r
# S3 method for class 'stanreg'
posterior_predict(
  object,
  newdata = NULL,
  draws = NULL,
  re.form = NULL,
  fun = NULL,
  seed = NULL,
  offset = NULL,
  ...
)

# S3 method for class 'stanmvreg'
posterior_predict(
  object,
  m = 1,
  newdata = NULL,
  draws = NULL,
  re.form = NULL,
  fun = NULL,
  seed = NULL,
  offset = NULL,
  ...
)
```

## Arguments

- object:

  A fitted model object returned by one of the rstanarm modeling
  functions. See
  [`stanreg-objects`](https://mc-stan.org/rstanarm/reference/stanreg-objects.md).

- newdata:

  Optionally, a data frame in which to look for variables with which to
  predict. If omitted, the model matrix is used. If `newdata` is
  provided and any variables were transformed (e.g. rescaled) in the
  data used to fit the model, then these variables must also be
  transformed in `newdata`. This only applies if variables were
  transformed before passing the data to one of the modeling functions
  and *not* if transformations were specified inside the model formula.
  Also see the Note section below for a note about using the `newdata`
  argument with with binomial models.

- draws:

  An integer indicating the number of draws to return. The default and
  maximum number of draws is the size of the posterior sample.

- re.form:

  If `object` contains
  [`group-level`](https://mc-stan.org/rstanarm/reference/stan_glmer.md)
  parameters, a formula indicating which group-level parameters to
  condition on when making predictions. `re.form` is specified in the
  same form as for
  [`predict.merMod`](https://rdrr.io/pkg/lme4/man/predict.merMod.html).
  The default, `NULL`, indicates that all estimated group-level
  parameters are conditioned on. To refrain from conditioning on any
  group-level parameters, specify `NA` or `~0`. The `newdata` argument
  may include new *levels* of the grouping factors that were specified
  when the model was estimated, in which case the resulting posterior
  predictions marginalize over the relevant variables.

- fun:

  An optional function to apply to the results. `fun` is found by a call
  to [`match.fun`](https://rdrr.io/r/base/match.fun.html) and so can be
  specified as a function object, a string naming a function, etc.

- seed:

  An optional [`seed`](https://rdrr.io/r/base/Random.html) to use.

- offset:

  A vector of offsets. Only required if `newdata` is specified and an
  `offset` argument was specified when fitting the model.

- ...:

  For `stanmvreg` objects, argument `m` can be specified indicating the
  submodel for which you wish to obtain predictions.

- m:

  Integer specifying the number or name of the submodel

## Value

A `draws` by `nrow(newdata)` matrix of simulations from the posterior
predictive distribution. Each row of the matrix is a vector of
predictions generated using a single draw of the model parameters from
the posterior distribution.

## Note

For binomial models with a number of trials greater than one (i.e., not
Bernoulli models), if `newdata` is specified then it must include all
variables needed for computing the number of binomial trials to use for
the predictions. For example if the left-hand side of the model formula
is `cbind(successes, failures)` then both `successes` and `failures`
must be in `newdata`. The particular values of `successes` and
`failures` in `newdata` do not matter so long as their sum is the
desired number of trials. If the left-hand side of the model formula
were `cbind(successes, trials - successes)` then both `trials` and
`successes` would need to be in `newdata`, probably with `successes` set
to `0` and `trials` specifying the number of trials. See the Examples
section below and the *How to Use the rstanarm Package* for examples.

For models estimated with
[`stan_clogit`](https://mc-stan.org/rstanarm/reference/stan_clogit.md),
the number of successes per stratum is ostensibly fixed by the research
design. Thus, when doing posterior prediction with new data, the
`data.frame` passed to the `newdata` argument must contain an outcome
variable and a stratifying factor, both with the same name as in the
original `data.frame`. Then, the posterior predictions will condition on
this outcome in the new data.

## See also

[`pp_check`](https://mc-stan.org/rstanarm/reference/pp_check.stanreg.md)
for graphical posterior predictive checks. Examples of posterior
predictive checking can also be found in the rstanarm vignettes and
demos.

[`predictive_error`](https://mc-stan.org/rstanarm/reference/predictive_error.stanreg.md)
and
[`predictive_interval`](https://mc-stan.org/rstanarm/reference/predictive_interval.stanreg.md).

## Examples

``` r
if (.Platform$OS.type != "windows" || .Platform$r_arch != "i386") {
if (!exists("example_model")) example(example_model)
yrep <- posterior_predict(example_model)
table(yrep)

# \donttest{
# Using newdata
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)
dat <- data.frame(counts, treatment, outcome)
fit3 <- stan_glm(
  counts ~ outcome + treatment, 
  data = dat,
  family = poisson(link="log"),
  prior = normal(0, 1, autoscale = FALSE), 
  prior_intercept = normal(0, 5, autoscale = FALSE),
  refresh = 0
)
nd <- data.frame(treatment = factor(rep(1,3)), outcome = factor(1:3))
ytilde <- posterior_predict(fit3, nd, draws = 500)
print(dim(ytilde))  # 500 by 3 matrix (draws by nrow(nd))

ytilde <- data.frame(
  count = c(ytilde),
  outcome = rep(nd$outcome, each = 500)
)
ggplot2::ggplot(ytilde, ggplot2::aes(x=outcome, y=count)) +
  ggplot2::geom_boxplot() +
  ggplot2::ylab("predicted count")


# Using newdata with a binomial model.
# example_model is binomial so we need to set
# the number of trials to use for prediction.
# This could be a different number for each
# row of newdata or the same for all rows.
# Here we'll use the same value for all.
nd <- lme4::cbpp
print(formula(example_model))  # cbind(incidence, size - incidence) ~ ...
nd$size <- max(nd$size) + 1L   # number of trials
nd$incidence <- 0  # set to 0 so size - incidence = number of trials
ytilde <- posterior_predict(example_model, newdata = nd)


# Using fun argument to transform predictions
mtcars2 <- mtcars
mtcars2$log_mpg <- log(mtcars2$mpg)
fit <- stan_glm(log_mpg ~ wt, data = mtcars2, refresh = 0)
ytilde <- posterior_predict(fit, fun = exp)
# }
}
#> [1] 500   3
#> cbind(incidence, size - incidence) ~ size + period + (1 | herd)
```
