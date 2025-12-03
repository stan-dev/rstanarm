# Conditional logistic (clogit) regression models via Stan

A model for case-control studies with optional prior distributions for
the coefficients, intercept, and auxiliary parameters.

## Usage

``` r
stan_clogit(
  formula,
  data,
  subset,
  na.action = NULL,
  contrasts = NULL,
  ...,
  strata,
  prior = normal(autoscale = TRUE),
  prior_covariance = decov(),
  prior_PD = FALSE,
  algorithm = c("sampling", "optimizing", "meanfield", "fullrank"),
  adapt_delta = NULL,
  QR = FALSE,
  sparse = FALSE
)
```

## Arguments

- formula, data, subset, na.action, contrasts:

  Same as for [`glmer`](https://rdrr.io/pkg/lme4/man/glmer.html), except
  that any global intercept included in the formula will be dropped. *We
  strongly advise against omitting the `data` argument*. Unless `data`
  is specified (and is a data frame) many post-estimation functions
  (including `update`, `loo`, `kfold`) are not guaranteed to work
  properly.

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

- strata:

  A factor indicating the groups in the data where the number of
  successes (possibly one) is fixed by the research design. It may be
  useful to use [`interaction`](https://rdrr.io/r/base/interaction.html)
  or [`strata`](https://rdrr.io/pkg/survival/man/strata.html) to create
  this factor. However, the `strata` argument must not rely on any
  object besides the `data`
  [`data.frame`](https://rdrr.io/r/base/data.frame.html).

- prior:

  The prior distribution for the (non-hierarchical) regression
  coefficients.

  The default priors are described in the vignette [*Prior Distributions
  for rstanarm
  Models*](https://mc-stan.org/rstanarm/articles/priors.html). If not
  using the default, `prior` should be a call to one of the various
  functions provided by rstanarm for specifying priors. The subset of
  these functions that can be used for the prior on the coefficients can
  be grouped into several "families":

  |                                 |                                 |
  |---------------------------------|---------------------------------|
  | **Family**                      | **Functions**                   |
  | *Student t family*              | `normal`, `student_t`, `cauchy` |
  | *Hierarchical shrinkage family* | `hs`, `hs_plus`                 |
  | *Laplace family*                | `laplace`, `lasso`              |
  | *Product normal family*         | `product_normal`                |

  See the [priors help
  page](https://mc-stan.org/rstanarm/reference/priors.md) for details on
  the families and how to specify the arguments for all of the functions
  in the table above. To omit a prior —i.e., to use a flat (improper)
  uniform prior— `prior` can be set to `NULL`, although this is rarely a
  good idea.

  **Note:** Unless `QR=TRUE`, if `prior` is from the Student t family or
  Laplace family, and if the `autoscale` argument to the function used
  to specify the prior (e.g.
  [`normal`](https://mc-stan.org/rstanarm/reference/priors.md)) is left
  at its default and recommended value of `TRUE`, then the default or
  user-specified prior scale(s) may be adjusted internally based on the
  scales of the predictors. See the [priors help
  page](https://mc-stan.org/rstanarm/reference/priors.md) and the *Prior
  Distributions* vignette for details on the rescaling and the
  [`prior_summary`](https://mc-stan.org/rstanarm/reference/prior_summary.stanreg.md)
  function for a summary of the priors used for a particular model.

- prior_covariance:

  Cannot be `NULL` when lme4-style group-specific terms are included in
  the `formula`. See
  [`decov`](https://mc-stan.org/rstanarm/reference/priors.md) for more
  information about the default arguments. Ignored when there are no
  group-specific terms.

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

## Value

A [stanreg](https://mc-stan.org/rstanarm/reference/stanreg-objects.md)
object is returned for `stan_clogit`.

## Details

The `stan_clogit` function is mostly similar in syntax to
[`clogit`](https://rdrr.io/pkg/survival/man/clogit.html) but rather than
performing maximum likelihood estimation of generalized linear models,
full Bayesian estimation is performed (if `algorithm` is `"sampling"`)
via MCMC. The Bayesian model adds priors (independent by default) on the
coefficients of the GLM.

The `data.frame` passed to the `data` argument must be sorted by the
variable passed to the `strata` argument.

The `formula` may have group-specific terms like in
[`stan_glmer`](https://mc-stan.org/rstanarm/reference/stan_glmer.md) but
should not allow the intercept to vary by the stratifying variable,
since there is no information in the data with which to estimate such
deviations in the intercept.

## See also

[`stanreg-methods`](https://mc-stan.org/rstanarm/reference/stanreg-methods.md)
and [`clogit`](https://rdrr.io/pkg/survival/man/clogit.html).

The vignette for Bernoulli and binomial models, which has more details
on using `stan_clogit`. <https://mc-stan.org/rstanarm/articles/>

## Examples

``` r
if (.Platform$OS.type != "windows" || .Platform$r_arch != "i386") {
dat <- infert[order(infert$stratum), ] # order by strata
post <- stan_clogit(case ~ spontaneous + induced + (1 | education), 
                    strata = stratum,
                    data = dat,
                    subset = parity <= 2,
                    QR = TRUE,
                    chains = 2, iter = 500) # for speed only

nd <- dat[dat$parity > 2, c("case", "spontaneous", "induced", "education", "stratum")]
# next line would fail without case and stratum variables                                 
pr <- posterior_epred(post, newdata = nd) # get predicted probabilities

# not a random variable b/c probabilities add to 1 within strata
all.equal(rep(sum(nd$case), nrow(pr)), rowSums(pr)) 
}
#> 
#> SAMPLING FOR MODEL 'bernoulli' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 4.8e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.48 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: Iteration:   1 / 500 [  0%]  (Warmup)
#> Chain 1: Iteration:  50 / 500 [ 10%]  (Warmup)
#> Chain 1: Iteration: 100 / 500 [ 20%]  (Warmup)
#> Chain 1: Iteration: 150 / 500 [ 30%]  (Warmup)
#> Chain 1: Iteration: 200 / 500 [ 40%]  (Warmup)
#> Chain 1: Iteration: 250 / 500 [ 50%]  (Warmup)
#> Chain 1: Iteration: 251 / 500 [ 50%]  (Sampling)
#> Chain 1: Iteration: 300 / 500 [ 60%]  (Sampling)
#> Chain 1: Iteration: 350 / 500 [ 70%]  (Sampling)
#> Chain 1: Iteration: 400 / 500 [ 80%]  (Sampling)
#> Chain 1: Iteration: 450 / 500 [ 90%]  (Sampling)
#> Chain 1: Iteration: 500 / 500 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 0.218 seconds (Warm-up)
#> Chain 1:                0.075 seconds (Sampling)
#> Chain 1:                0.293 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'bernoulli' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 3.4e-05 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.34 seconds.
#> Chain 2: Adjust your expectations accordingly!
#> Chain 2: 
#> Chain 2: 
#> Chain 2: Iteration:   1 / 500 [  0%]  (Warmup)
#> Chain 2: Iteration:  50 / 500 [ 10%]  (Warmup)
#> Chain 2: Iteration: 100 / 500 [ 20%]  (Warmup)
#> Chain 2: Iteration: 150 / 500 [ 30%]  (Warmup)
#> Chain 2: Iteration: 200 / 500 [ 40%]  (Warmup)
#> Chain 2: Iteration: 250 / 500 [ 50%]  (Warmup)
#> Chain 2: Iteration: 251 / 500 [ 50%]  (Sampling)
#> Chain 2: Iteration: 300 / 500 [ 60%]  (Sampling)
#> Chain 2: Iteration: 350 / 500 [ 70%]  (Sampling)
#> Chain 2: Iteration: 400 / 500 [ 80%]  (Sampling)
#> Chain 2: Iteration: 450 / 500 [ 90%]  (Sampling)
#> Chain 2: Iteration: 500 / 500 [100%]  (Sampling)
#> Chain 2: 
#> Chain 2:  Elapsed Time: 0.22 seconds (Warm-up)
#> Chain 2:                0.082 seconds (Sampling)
#> Chain 2:                0.302 seconds (Total)
#> Chain 2: 
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> [1] TRUE
```
