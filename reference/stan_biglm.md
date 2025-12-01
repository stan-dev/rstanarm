# Bayesian regularized linear but big models via Stan

This is the same model as with
[`stan_lm`](https://mc-stan.org/rstanarm/reference/stan_lm.md) but it
utilizes the output from
[`biglm`](https://rdrr.io/pkg/biglm/man/biglm.html) in the biglm package
in order to proceed when the data is too large to fit in memory.

## Usage

``` r
stan_biglm(
  biglm,
  xbar,
  ybar,
  s_y,
  ...,
  prior = R2(stop("'location' must be specified")),
  prior_intercept = NULL,
  prior_PD = FALSE,
  algorithm = c("sampling", "meanfield", "fullrank"),
  adapt_delta = NULL
)

stan_biglm.fit(
  b,
  R,
  SSR,
  N,
  xbar,
  ybar,
  s_y,
  has_intercept = TRUE,
  ...,
  prior = R2(stop("'location' must be specified")),
  prior_intercept = NULL,
  prior_PD = FALSE,
  algorithm = c("sampling", "meanfield", "fullrank", "optimizing"),
  adapt_delta = NULL,
  importance_resampling = TRUE,
  keep_every = 1
)
```

## Arguments

- biglm:

  The list output by [`biglm`](https://rdrr.io/pkg/biglm/man/biglm.html)
  in the biglm package.

- xbar:

  A numeric vector of column means in the implicit design matrix
  excluding the intercept for the observations included in the model.

- ybar:

  A numeric scalar indicating the mean of the outcome for the
  observations included in the model.

- s_y:

  A numeric scalar indicating the unbiased sample standard deviation of
  the outcome for the observations included in the model.

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

- prior:

  Must be a call to
  [`R2`](https://mc-stan.org/rstanarm/reference/priors.md) with its
  `location` argument specified or `NULL`, which would indicate a
  standard uniform prior for the \\R^2\\.

- prior_intercept:

  Either `NULL` (the default) or a call to
  [`normal`](https://mc-stan.org/rstanarm/reference/priors.md). If a
  [`normal`](https://mc-stan.org/rstanarm/reference/priors.md) prior is
  specified without a `scale`, then the standard deviation is taken to
  be the marginal standard deviation of the outcome divided by the
  square root of the sample size, which is legitimate because the
  marginal standard deviation of the outcome is a primitive parameter
  being estimated.

  **Note:** If using a dense representation of the design matrix —i.e.,
  if the `sparse` argument is left at its default value of `FALSE`— then
  the prior distribution for the intercept is set so it applies to the
  value *when all predictors are centered*. If you prefer to specify a
  prior on the intercept without the predictors being auto-centered,
  then you have to omit the intercept from the
  [`formula`](https://rdrr.io/r/stats/formula.html) and include a column
  of ones as a predictor, in which case some element of `prior`
  specifies the prior on it, rather than `prior_intercept`. Regardless
  of how `prior_intercept` is specified, the reported *estimates* of the
  intercept always correspond to a parameterization without centered
  predictors (i.e., same as in `glm`).

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

- b:

  A numeric vector of OLS coefficients, excluding the intercept

- R:

  A square upper-triangular matrix from the QR decomposition of the
  design matrix, excluding the intercept

- SSR:

  A numeric scalar indicating the sum-of-squared residuals for OLS

- N:

  A integer scalar indicating the number of included observations

- has_intercept:

  A logical scalar indicating whether to add an intercept to the model
  when estimating it.

- importance_resampling:

  Logical scalar indicating whether to use importance resampling when
  approximating the posterior distribution with a multivariate normal
  around the posterior mode, which only applies when `algorithm` is
  `"optimizing"` but defaults to `TRUE` in that case

- keep_every:

  Positive integer, which defaults to 1, but can be higher in order to
  thin the importance sampling realizations and also only apples when
  `algorithm` is `"optimizing"` but defaults to `TRUE` in that case

## Value

The output of both `stan_biglm` and `stan_biglm.fit` is an object of
[`stanfit-class`](https://mc-stan.org/rstan/reference/stanfit-class.html)
rather than
[`stanreg-objects`](https://mc-stan.org/rstanarm/reference/stanreg-objects.md),
which is more limited and less convenient but necessitated by the fact
that `stan_biglm` does not bring the full design matrix into memory.
Without the full design matrix,some of the elements of a
[`stanreg-objects`](https://mc-stan.org/rstanarm/reference/stanreg-objects.md)
object cannot be calculated, such as residuals. Thus, the functions in
the rstanarm package that input
[`stanreg-objects`](https://mc-stan.org/rstanarm/reference/stanreg-objects.md),
such as
[`posterior_predict`](https://mc-stan.org/rstanarm/reference/posterior_predict.stanreg.md)
cannot be used.

## Details

The `stan_biglm` function is intended to be used in the same
circumstances as the [`biglm`](https://rdrr.io/pkg/biglm/man/biglm.html)
function in the biglm package but with an informative prior on the
\\R^2\\ of the regression. Like
[`biglm`](https://rdrr.io/pkg/biglm/man/biglm.html), the memory required
to estimate the model depends largely on the number of predictors rather
than the number of observations. However, `stan_biglm` and
`stan_biglm.fit` have additional required arguments that are not
necessary in [`biglm`](https://rdrr.io/pkg/biglm/man/biglm.html), namely
`xbar`, `ybar`, and `s_y`. If any observations have any missing values
on any of the predictors or the outcome, such observations do not
contribute to these statistics.

## Examples

``` r
if (.Platform$OS.type != "windows" || .Platform$r_arch != "i386") {
# create inputs
ols <- lm(mpg ~ wt + qsec + am, data = mtcars, # all row are complete so ...
          na.action = na.exclude)              # not necessary in this case
b <- coef(ols)[-1]
R <- qr.R(ols$qr)[-1,-1]
SSR <- crossprod(ols$residuals)[1]
not_NA <- !is.na(fitted(ols))
N <- sum(not_NA)
xbar <- colMeans(mtcars[not_NA,c("wt", "qsec", "am")])
y <- mtcars$mpg[not_NA]
ybar <- mean(y)
s_y <- sd(y)
post <- stan_biglm.fit(b, R, SSR, N, xbar, ybar, s_y, prior = R2(.75),
                       # the next line is only to make the example go fast
                       chains = 1, iter = 500, seed = 12345)
cbind(lm = b, stan_lm = rstan::get_posterior_mean(post)[13:15,]) # shrunk
}
#> 
#> SAMPLING FOR MODEL 'lm' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 1.7e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.17 seconds.
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
#> Chain 1:  Elapsed Time: 0.512 seconds (Warm-up)
#> Chain 1:                0.259 seconds (Sampling)
#> Chain 1:                0.771 seconds (Total)
#> Chain 1: 
#> Warning: There were 1 divergent transitions after warmup. See
#> https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#> to find out why this is a problem and how to eliminate them.
#> Warning: Examine the pairs() plot to diagnose sampling problems
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
#>             lm   stan_lm
#> wt   -3.916504 -3.694667
#> qsec  1.225886  1.190556
#> am    2.935837  2.984365
```
