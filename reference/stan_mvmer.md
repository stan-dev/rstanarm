# Bayesian multivariate generalized linear models with correlated group-specific terms via Stan

Bayesian inference for multivariate GLMs with group-specific
coefficients that are assumed to be correlated across the GLM submodels.

## Usage

``` r
stan_mvmer(
  formula,
  data,
  family = gaussian,
  weights,
  prior = normal(autoscale = TRUE),
  prior_intercept = normal(autoscale = TRUE),
  prior_aux = cauchy(0, 5, autoscale = TRUE),
  prior_covariance = lkj(autoscale = TRUE),
  prior_PD = FALSE,
  algorithm = c("sampling", "meanfield", "fullrank"),
  adapt_delta = NULL,
  max_treedepth = 10L,
  init = "random",
  QR = FALSE,
  sparse = FALSE,
  ...
)
```

## Arguments

- formula:

  A two-sided linear formula object describing both the fixed-effects
  and random-effects parts of the longitudinal submodel similar in vein
  to formula specification in the **lme4** package (see
  [`glmer`](https://rdrr.io/pkg/lme4/man/glmer.html) or the **lme4**
  vignette for details). Note however that the double bar (`||`)
  notation is not allowed when specifying the random-effects parts of
  the formula, and neither are nested grouping factors (e.g.
  `(1 | g1/g2))` or `(1 | g1:g2)`, where `g1`, `g2` are grouping
  factors. For a multivariate GLM this should be a list of such formula
  objects, with each element of the list providing the formula for one
  of the GLM submodels.

- data:

  A data frame containing the variables specified in `formula`. For a
  multivariate GLM, this can be either a single data frame which
  contains the data for all GLM submodels, or it can be a list of data
  frames where each element of the list provides the data for one of the
  GLM submodels.

- family:

  The family (and possibly also the link function) for the GLM
  submodel(s). See [`glmer`](https://rdrr.io/pkg/lme4/man/glmer.html)
  for details. If fitting a multivariate GLM, then this can optionally
  be a list of families, in which case each element of the list
  specifies the family for one of the GLM submodels. In other words, a
  different family can be specified for each GLM submodel.

- weights:

  Same as in [`glm`](https://rdrr.io/r/stats/glm.html), except that when
  fitting a multivariate GLM and a list of data frames is provided in
  `data` then a corresponding list of weights must be provided. If
  weights are provided for one of the GLM submodels, then they must be
  provided for all GLM submodels.

- prior, prior_intercept, prior_aux:

  Same as in
  [`stan_glmer`](https://mc-stan.org/rstanarm/reference/stan_glmer.md)
  except that for a multivariate GLM a list of priors can be provided
  for any of `prior`, `prior_intercept` or `prior_aux` arguments. That
  is, different priors can optionally be specified for each of the GLM
  submodels. If a list is not provided, then the same prior
  distributions are used for each GLM submodel. Note that the
  `"product_normal"` prior is not allowed for `stan_mvmer`.

- prior_covariance:

  Cannot be `NULL`; see
  [`priors`](https://mc-stan.org/rstanarm/reference/priors.md) for more
  information about the prior distributions on covariance matrices. Note
  however that the default prior for covariance matrices in `stan_mvmer`
  is slightly different to that in
  [`stan_glmer`](https://mc-stan.org/rstanarm/reference/stan_glmer.md)
  (the details of which are described on the
  [`priors`](https://mc-stan.org/rstanarm/reference/priors.md) page).

- prior_PD:

  A logical scalar (defaulting to `FALSE`) indicating whether to draw
  from the prior predictive distribution instead of conditioning on the
  outcome.

- algorithm:

  A string (possibly abbreviated) indicating the estimation approach to
  use. Can be `"sampling"` for MCMC (the default), `"optimizing"` for
  optimization, `"meanfield"` for variational inference with independent
  normal distributions, or `"fullrank"` for variational inference with a
  multivariate normal distribution. See
  [`rstanarm-package`](https://mc-stan.org/rstanarm/reference/rstanarm-package.md)
  for more details on the estimation algorithms. NOTE: not all fitting
  functions support all four algorithms.

- adapt_delta:

  Only relevant if `algorithm="sampling"`. See the
  [adapt_delta](https://mc-stan.org/rstanarm/reference/adapt_delta.md)
  help page for details.

- max_treedepth:

  A positive integer specifying the maximum treedepth for the non-U-turn
  sampler. See the `control` argument in
  [`stan`](https://mc-stan.org/rstan/reference/stan.html).

- init:

  The method for generating initial values. See
  [`stan`](https://mc-stan.org/rstan/reference/stan.html).

- QR:

  A logical scalar defaulting to `FALSE`, but if `TRUE` applies a scaled
  [`qr`](https://rdrr.io/r/base/qr.html) decomposition to the design
  matrix. The transformation does not change the likelihood of the data
  but is recommended for computational reasons when there are multiple
  predictors. See the
  [QR-argument](https://mc-stan.org/rstanarm/reference/QR-argument.md)
  documentation page for details on how rstanarm does the transformation
  and important information about how to interpret the prior
  distributions of the model parameters when using `QR=TRUE`.

- sparse:

  A logical scalar (defaulting to `FALSE`) indicating whether to use a
  sparse representation of the design (X) matrix. If `TRUE`, the the
  design matrix is not centered (since that would destroy the sparsity)
  and likewise it is not possible to specify both `QR = TRUE` and
  `sparse = TRUE`. Depending on how many zeros there are in the design
  matrix, setting `sparse = TRUE` may make the code run faster and can
  consume much less RAM.

- ...:

  Further arguments passed to the function in the rstan package
  (`sampling`, `vb`, or `optimizing`), corresponding to the estimation
  method named by `algorithm`. For example, if `algorithm` is
  `"sampling"` it is possible to specify `iter`, `chains`, `cores`, and
  other MCMC controls.

  Another useful argument that can be passed to rstan via `...` is
  `refresh`, which specifies how often to print updates when sampling
  (i.e., show the progress every `refresh` iterations). `refresh=0`
  turns off the iteration updates.

## Value

A [stanmvreg](https://mc-stan.org/rstanarm/reference/stanreg-objects.md)
object is returned.

## Details

The `stan_mvmer` function can be used to fit a multivariate generalized
linear model (GLM) with group-specific terms. The model consists of
distinct GLM submodels, each which contains group-specific terms; within
a grouping factor (for example, patient ID) the grouping-specific terms
are assumed to be correlated across the different GLM submodels. It is
possible to specify a different outcome type (for example a different
family and/or link function) for each of the GLM submodels.  
  
Bayesian estimation of the model is performed via MCMC, in the same way
as for
[`stan_glmer`](https://mc-stan.org/rstanarm/reference/stan_glmer.md).
Also, similar to
[`stan_glmer`](https://mc-stan.org/rstanarm/reference/stan_glmer.md), an
unstructured covariance matrix is used for the group-specific terms
within a given grouping factor, with priors on the terms of a
decomposition of the covariance matrix.See
[`priors`](https://mc-stan.org/rstanarm/reference/priors.md) for more
information about the priors distributions that are available for the
covariance matrices, the regression coefficients and the intercept and
auxiliary parameters.

## See also

[`stan_glmer`](https://mc-stan.org/rstanarm/reference/stan_glmer.md),
[`stan_jm`](https://mc-stan.org/rstanarm/reference/stan_jm.md),
[`stanreg-objects`](https://mc-stan.org/rstanarm/reference/stanreg-objects.md),
[`stanmvreg-methods`](https://mc-stan.org/rstanarm/reference/stanmvreg-methods.md),
[`print.stanmvreg`](https://mc-stan.org/rstanarm/reference/print.stanreg.md),
[`summary.stanmvreg`](https://mc-stan.org/rstanarm/reference/summary.stanreg.md),
[`posterior_predict`](https://mc-stan.org/rstanarm/reference/posterior_predict.stanreg.md),
[`posterior_interval`](https://mc-stan.org/rstanarm/reference/posterior_interval.stanreg.md).

## Examples

``` r
if (.Platform$OS.type != "windows" || .Platform$r_arch !="i386") {
# \donttest{
#####
# A multivariate GLM with two submodels. For the grouping factor 'id', the 
# group-specific intercept from the first submodel (logBili) is assumed to
# be correlated with the group-specific intercept and linear slope in the 
# second submodel (albumin)
f1 <- stan_mvmer(
        formula = list(
          logBili ~ year + (1 | id), 
          albumin ~ sex + year + (year | id)),
        data = pbcLong, 
        # this next line is only to keep the example small in size!
        chains = 1, cores = 1, seed = 12345, iter = 1000)
summary(f1) 

#####
# A multivariate GLM with one bernoulli outcome and one
# gaussian outcome. We will artificially create the bernoulli
# outcome by dichotomising log serum bilirubin
pbcLong$ybern <- as.integer(pbcLong$logBili >= mean(pbcLong$logBili))
f2 <- stan_mvmer(
        formula = list(
          ybern ~ year + (1 | id), 
          albumin ~ sex + year + (year | id)),
        data = pbcLong,
        family = list(binomial, gaussian),
        chains = 1, cores = 1, seed = 12345, iter = 1000)
# }
}
#> Fitting a multivariate glmer model.
#> 
#> Please note the warmup may be much slower than later iterations!
#> 
#> SAMPLING FOR MODEL 'mvmer' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 9.5e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.95 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: Iteration:   1 / 1000 [  0%]  (Warmup)
#> Chain 1: Iteration: 100 / 1000 [ 10%]  (Warmup)
#> Chain 1: Iteration: 200 / 1000 [ 20%]  (Warmup)
#> Chain 1: Iteration: 300 / 1000 [ 30%]  (Warmup)
#> Chain 1: Iteration: 400 / 1000 [ 40%]  (Warmup)
#> Chain 1: Iteration: 500 / 1000 [ 50%]  (Warmup)
#> Chain 1: Iteration: 501 / 1000 [ 50%]  (Sampling)
#> Chain 1: Iteration: 600 / 1000 [ 60%]  (Sampling)
#> Chain 1: Iteration: 700 / 1000 [ 70%]  (Sampling)
#> Chain 1: Iteration: 800 / 1000 [ 80%]  (Sampling)
#> Chain 1: Iteration: 900 / 1000 [ 90%]  (Sampling)
#> Chain 1: Iteration: 1000 / 1000 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 3.67 seconds (Warm-up)
#> Chain 1:                2.033 seconds (Sampling)
#> Chain 1:                5.703 seconds (Total)
#> Chain 1: 
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Fitting a multivariate glmer model.
#> 
#> Please note the warmup may be much slower than later iterations!
#> 
#> SAMPLING FOR MODEL 'mvmer' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 9.8e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.98 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: Iteration:   1 / 1000 [  0%]  (Warmup)
#> Chain 1: Iteration: 100 / 1000 [ 10%]  (Warmup)
#> Chain 1: Iteration: 200 / 1000 [ 20%]  (Warmup)
#> Chain 1: Iteration: 300 / 1000 [ 30%]  (Warmup)
#> Chain 1: Iteration: 400 / 1000 [ 40%]  (Warmup)
#> Chain 1: Iteration: 500 / 1000 [ 50%]  (Warmup)
#> Chain 1: Iteration: 501 / 1000 [ 50%]  (Sampling)
#> Chain 1: Iteration: 600 / 1000 [ 60%]  (Sampling)
#> Chain 1: Iteration: 700 / 1000 [ 70%]  (Sampling)
#> Chain 1: Iteration: 800 / 1000 [ 80%]  (Sampling)
#> Chain 1: Iteration: 900 / 1000 [ 90%]  (Sampling)
#> Chain 1: Iteration: 1000 / 1000 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 4.029 seconds (Warm-up)
#> Chain 1:                2.111 seconds (Sampling)
#> Chain 1:                6.14 seconds (Total)
#> Chain 1: 
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
```
