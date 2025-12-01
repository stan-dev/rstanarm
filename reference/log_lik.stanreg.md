# Pointwise log-likelihood matrix

For models fit using MCMC only, the `log_lik` method returns the \\S\\
by \\N\\ pointwise log-likelihood matrix, where \\S\\ is the size of the
posterior sample and \\N\\ is the number of data points, or in the case
of the `stanmvreg` method (when called on
[`stan_jm`](https://mc-stan.org/rstanarm/reference/stan_jm.md) model
objects) an \\S\\ by \\Npat\\ matrix where \\Npat\\ is the number of
individuals.

## Usage

``` r
# S3 method for class 'stanreg'
log_lik(object, newdata = NULL, offset = NULL, ...)

# S3 method for class 'stanmvreg'
log_lik(object, m = 1, newdata = NULL, ...)

# S3 method for class 'stanjm'
log_lik(object, newdataLong = NULL, newdataEvent = NULL, ...)
```

## Arguments

- object:

  A fitted model object returned by one of the rstanarm modeling
  functions. See
  [`stanreg-objects`](https://mc-stan.org/rstanarm/reference/stanreg-objects.md).

- newdata:

  An optional data frame of new data (e.g. holdout data) to use when
  evaluating the log-likelihood. See the description of `newdata` for
  [`posterior_predict`](https://mc-stan.org/rstanarm/reference/posterior_predict.stanreg.md).

- offset:

  A vector of offsets. Only required if `newdata` is specified and an
  `offset` was specified when fitting the model.

- ...:

  Currently ignored.

- m:

  Integer specifying the number or name of the submodel

- newdataLong, newdataEvent:

  Optional data frames containing new data (e.g. holdout data) to use
  when evaluating the log-likelihood for a model estimated using
  [`stan_jm`](https://mc-stan.org/rstanarm/reference/stan_jm.md). If the
  fitted model was a multivariate joint model (i.e. more than one
  longitudinal outcome), then `newdataLong` is allowed to be a list of
  data frames. If supplying new data, then `newdataEvent` should also
  include variables corresponding to the event time and event indicator
  as these are required for evaluating the log likelihood for the event
  submodel. For more details, see the description of `newdataLong` and
  `newdataEvent` for
  [`posterior_survfit`](https://mc-stan.org/rstanarm/reference/posterior_survfit.md).

## Value

For the `stanreg` and `stanmvreg` methods an \\S\\ by \\N\\ matrix,
where \\S\\ is the size of the posterior sample and \\N\\ is the number
of data points. For the `stanjm` method an \\S\\ by \\Npat\\ matrix
where \\Npat\\ is the number of individuals.

## Examples

``` r
if (.Platform$OS.type != "windows" || .Platform$r_arch != "i386") {
# \donttest{
 roaches$roach100 <- roaches$roach1 / 100
 fit <- stan_glm(
    y ~ roach100 + treatment + senior,
    offset = log(exposure2),
    data = roaches,
    family = poisson(link = "log"),
    prior = normal(0, 2.5),
    prior_intercept = normal(0, 10),
    iter = 500, # just to speed up example,
    refresh = 0
 )
 ll <- log_lik(fit)
 dim(ll)
 all.equal(ncol(ll), nobs(fit))

 # using newdata argument
 nd <- roaches[1:2, ]
 nd$treatment[1:2] <- c(0, 1)
 ll2 <- log_lik(fit, newdata = nd, offset = c(0, 0))
 head(ll2)
 dim(ll2)
 all.equal(ncol(ll2), nrow(nd))
# }
}
#> [1] TRUE
```
