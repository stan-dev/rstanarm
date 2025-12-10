# Summarize the priors used for an rstanarm model

The `prior_summary` method provides a summary of the prior distributions
used for the parameters in a given model. In some cases the
user-specified prior does not correspond exactly to the prior used
internally by rstanarm (see the sections below). Especially in these
cases, but also in general, it can be much more useful to visualize the
priors. Visualizing the priors can be done using the
[`posterior_vs_prior`](https://mc-stan.org/rstanarm/reference/posterior_vs_prior.md)
function, or alternatively by fitting the model with the `prior_PD`
argument set to `TRUE` (to draw from the prior predictive distribution
instead of conditioning on the outcome) and then plotting the
parameters.

## Usage

``` r
# S3 method for class 'stanreg'
prior_summary(object, digits = 2, ...)
```

## Arguments

- object:

  A fitted model object returned by one of the rstanarm modeling
  functions. See
  [`stanreg-objects`](https://mc-stan.org/rstanarm/reference/stanreg-objects.md).

- digits:

  Number of digits to use for rounding.

- ...:

  Currently ignored by the method for stanreg objects.

## Value

A list of class "prior_summary.stanreg", which has its own print method.

## Intercept (after predictors centered)

For rstanarm modeling functions that accept a `prior_intercept`
argument, the specified prior for the intercept term applies to the
intercept after rstanarm internally centers the predictors so they each
have mean zero. The estimate of the intercept returned to the user
correspond to the intercept with the predictors as specified by the user
(unmodified by rstanarm), but when *specifying* the prior the intercept
can be thought of as the expected outcome when the predictors are set to
their means. The only exception to this is for models fit with the
`sparse` argument set to `TRUE` (which is only possible with a subset of
the modeling functions and never the default).

## Adjusted scales

For some models you may see "`adjusted scale`" in the printed output and
adjusted scales included in the object returned by `prior_summary`.
These adjusted scale values are the prior scales actually used by
rstanarm and are computed by adjusting the prior scales specified by the
user to account for the scales of the predictors (as described in the
documentation for the
[`autoscale`](https://mc-stan.org/rstanarm/reference/priors.md)
argument). To disable internal prior scale adjustments set the
`autoscale` argument to `FALSE` when setting a prior using one of the
distributions that accepts an `autoscale` argument. For example,
`normal(0, 5, autoscale=FALSE)` instead of just `normal(0, 5)`.

## Coefficients in Q-space

For the models fit with an rstanarm modeling function that supports the
`QR` argument (see e.g,
[`stan_glm`](https://mc-stan.org/rstanarm/reference/stan_glm.md)), if
`QR` is set to `TRUE` then the prior distributions for the regression
coefficients specified using the `prior` argument are not relative to
the original predictor variables \\X\\ but rather to the variables in
the matrix \\Q\\ obtained from the \\QR\\ decomposition of \\X\\.

In particular, if `prior = normal(location,scale)`, then this prior on
the coefficients in \\Q\\-space can be easily translated into a joint
multivariate normal (MVN) prior on the coefficients on the original
predictors in \\X\\. Letting \\\theta\\ denote the coefficients on \\Q\\
and \\\beta\\ the coefficients on \\X\\ then if \\\theta \sim N(\mu,
\sigma)\\ the corresponding prior on \\\beta\\ is \\\beta \sim MVN(R\mu,
R'R\sigma^2)\\, where \\\mu\\ and \\\sigma\\ are vectors of the
appropriate length. Technically, rstanarm uses a scaled \\QR\\
decomposition to ensure that the columns of the predictor matrix used to
fit the model all have unit scale, when the `autoscale` argument to the
function passed to the `prior` argument is `TRUE` (the default), in
which case the matrices actually used are \\Q^\ast = Q \sqrt{n-1}\\ and
\\R^\ast = \frac{1}{\sqrt{n-1}} R\\. If `autoscale = FALSE` we instead
scale such that the lower-right element of \\R^\ast\\ is \\1\\, which is
useful if you want to specify a prior on the coefficient of the last
predictor in its original units (see the documentation for the
[`QR`](https://mc-stan.org/rstanarm/reference/stan_glm.md) argument).

If you are interested in the prior on \\\beta\\ implied by the prior on
\\\theta\\, we strongly recommend visualizing it as described above in
the **Description** section, which is simpler than working it out
analytically.

## See also

The [priors help page](https://mc-stan.org/rstanarm/reference/priors.md)
and the *Prior Distributions* vignette.

## Examples

``` r
if (.Platform$OS.type != "windows" || .Platform$r_arch != "i386") {
if (!exists("example_model")) example(example_model) 
prior_summary(example_model)

priors <- prior_summary(example_model)
names(priors)
priors$prior$scale
priors$prior$adjusted_scale

# for a glm with adjusted scales (see Details, above), compare 
# the default (rstanarm adjusting the scales) to setting 
# autoscale=FALSE for prior on coefficients
fit <- stan_glm(mpg ~ wt + am, data = mtcars, 
                prior = normal(0, c(2.5, 4)), 
                prior_intercept = normal(0, 5), 
                iter = 10, chains = 1) # only for demonstration 
prior_summary(fit)

fit2 <- update(fit, prior = normal(0, c(2.5, 4), autoscale=FALSE), 
               prior_intercept = normal(0, 5, autoscale=FALSE))
prior_summary(fit2)
}
#> 
#> SAMPLING FOR MODEL 'continuous' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 2e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.2 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: WARNING: No variance estimation is
#> Chain 1:          performed for num_warmup < 20
#> Chain 1: 
#> Chain 1: Iteration: 1 / 10 [ 10%]  (Warmup)
#> Chain 1: Iteration: 2 / 10 [ 20%]  (Warmup)
#> Chain 1: Iteration: 3 / 10 [ 30%]  (Warmup)
#> Chain 1: Iteration: 4 / 10 [ 40%]  (Warmup)
#> Chain 1: Iteration: 5 / 10 [ 50%]  (Warmup)
#> Chain 1: Iteration: 6 / 10 [ 60%]  (Sampling)
#> Chain 1: Iteration: 7 / 10 [ 70%]  (Sampling)
#> Chain 1: Iteration: 8 / 10 [ 80%]  (Sampling)
#> Chain 1: Iteration: 9 / 10 [ 90%]  (Sampling)
#> Chain 1: Iteration: 10 / 10 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 0 seconds (Warm-up)
#> Chain 1:                0 seconds (Sampling)
#> Chain 1:                0 seconds (Total)
#> Chain 1: 
#> Warning: There were 1 chains where the estimated Bayesian Fraction of Missing Information was low. See
#> https://mc-stan.org/misc/warnings.html#bfmi-low
#> Warning: Examine the pairs() plot to diagnose sampling problems
#> Warning: The largest R-hat is 1.9, indicating chains have not mixed.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#r-hat
#> Warning: Markov chains did not converge! Do not analyze results!
#> 
#> SAMPLING FOR MODEL 'continuous' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 2.1e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.21 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: WARNING: No variance estimation is
#> Chain 1:          performed for num_warmup < 20
#> Chain 1: 
#> Chain 1: Iteration: 1 / 10 [ 10%]  (Warmup)
#> Chain 1: Iteration: 2 / 10 [ 20%]  (Warmup)
#> Chain 1: Iteration: 3 / 10 [ 30%]  (Warmup)
#> Chain 1: Iteration: 4 / 10 [ 40%]  (Warmup)
#> Chain 1: Iteration: 5 / 10 [ 50%]  (Warmup)
#> Chain 1: Iteration: 6 / 10 [ 60%]  (Sampling)
#> Chain 1: Iteration: 7 / 10 [ 70%]  (Sampling)
#> Chain 1: Iteration: 8 / 10 [ 80%]  (Sampling)
#> Chain 1: Iteration: 9 / 10 [ 90%]  (Sampling)
#> Chain 1: Iteration: 10 / 10 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 0 seconds (Warm-up)
#> Chain 1:                0 seconds (Sampling)
#> Chain 1:                0 seconds (Total)
#> Chain 1: 
#> Warning: The largest R-hat is 1.9, indicating chains have not mixed.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#r-hat
#> Warning: Markov chains did not converge! Do not analyze results!
#> Priors for model 'fit2' 
#> ------
#> Intercept (after predictors centered)
#>  ~ normal(location = 0, scale = 5)
#> 
#> Coefficients
#>  ~ normal(location = [0,0], scale = [2.5,4.0])
#> 
#> Auxiliary (sigma)
#>   Specified prior:
#>     ~ exponential(rate = 1)
#>   Adjusted prior:
#>     ~ exponential(rate = 0.17)
#> ------
#> See help('prior_summary.stanreg') for more details
```
