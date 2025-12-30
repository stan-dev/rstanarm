# Predictive intervals

For models fit using MCMC (`algorithm="sampling"`) or one of the
variational approximations (`"meanfield"` or `"fullrank"`), the
`predictive_interval` function computes Bayesian predictive intervals.
The method for stanreg objects calls
[`posterior_predict`](https://mc-stan.org/rstanarm/reference/posterior_predict.stanreg.md)
internally, whereas the method for matrices accepts the matrix returned
by `posterior_predict` as input and can be used to avoid multiple calls
to `posterior_predict`.

## Usage

``` r
# S3 method for class 'stanreg'
predictive_interval(
  object,
  prob = 0.9,
  newdata = NULL,
  draws = NULL,
  re.form = NULL,
  fun = NULL,
  seed = NULL,
  offset = NULL,
  ...
)

# S3 method for class 'matrix'
predictive_interval(object, prob = 0.9, ...)

# S3 method for class 'ppd'
predictive_interval(object, prob = 0.9, ...)
```

## Arguments

- object:

  Either a fitted model object returned by one of the rstanarm modeling
  functions (a [stanreg
  object](https://mc-stan.org/rstanarm/reference/stanreg-objects.md))
  or, for the matrix method, a matrix of draws from the posterior
  predictive distribution returned by
  [`posterior_predict`](https://mc-stan.org/rstanarm/reference/posterior_predict.stanreg.md).

- prob:

  A number \\p \in (0,1)\\ indicating the desired probability mass to
  include in the intervals. The default is to report \\90\\% intervals
  (`prob=0.9`) rather than the traditionally used \\95\\% (see Details).

- newdata, draws, fun, offset, re.form, seed:

  Passed to
  [`posterior_predict`](https://mc-stan.org/rstanarm/reference/posterior_predict.stanreg.md).

- ...:

  Currently ignored.

## Value

A matrix with two columns and as many rows as are in `newdata`. If
`newdata` is not provided then the matrix will have as many rows as the
data used to fit the model. For a given value of `prob`, \\p\\, the
columns correspond to the lower and upper \\100p\\% central interval
limits and have the names \\100\alpha/2\\% and \\100(1 - \alpha/2)\\%,
where \\\alpha = 1-p\\. For example, if `prob=0.9` is specified (a
\\90\\% interval), then the column names will be `"5%"` and `"95%"`,
respectively.

## See also

[`predictive_error`](https://mc-stan.org/rstanarm/reference/predictive_error.stanreg.md),
[`posterior_predict`](https://mc-stan.org/rstanarm/reference/posterior_predict.stanreg.md),
[`posterior_interval`](https://mc-stan.org/rstanarm/reference/posterior_interval.stanreg.md)

## Examples

``` r
if (.Platform$OS.type != "windows" || .Platform$r_arch != "i386") {
fit <- stan_glm(mpg ~ wt, data = mtcars, iter = 300)
predictive_interval(fit)
predictive_interval(fit, newdata = data.frame(wt = range(mtcars$wt)), 
                    prob = 0.5)

# stanreg vs matrix methods
preds <- posterior_predict(fit, seed = 123)
all.equal(
  predictive_interval(fit, seed = 123),
  predictive_interval(preds)
)
}
#> 
#> SAMPLING FOR MODEL 'continuous' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 2.2e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.22 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: Iteration:   1 / 300 [  0%]  (Warmup)
#> Chain 1: Iteration:  30 / 300 [ 10%]  (Warmup)
#> Chain 1: Iteration:  60 / 300 [ 20%]  (Warmup)
#> Chain 1: Iteration:  90 / 300 [ 30%]  (Warmup)
#> Chain 1: Iteration: 120 / 300 [ 40%]  (Warmup)
#> Chain 1: Iteration: 150 / 300 [ 50%]  (Warmup)
#> Chain 1: Iteration: 151 / 300 [ 50%]  (Sampling)
#> Chain 1: Iteration: 180 / 300 [ 60%]  (Sampling)
#> Chain 1: Iteration: 210 / 300 [ 70%]  (Sampling)
#> Chain 1: Iteration: 240 / 300 [ 80%]  (Sampling)
#> Chain 1: Iteration: 270 / 300 [ 90%]  (Sampling)
#> Chain 1: Iteration: 300 / 300 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 0.009 seconds (Warm-up)
#> Chain 1:                0.005 seconds (Sampling)
#> Chain 1:                0.014 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'continuous' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 1e-05 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.1 seconds.
#> Chain 2: Adjust your expectations accordingly!
#> Chain 2: 
#> Chain 2: 
#> Chain 2: Iteration:   1 / 300 [  0%]  (Warmup)
#> Chain 2: Iteration:  30 / 300 [ 10%]  (Warmup)
#> Chain 2: Iteration:  60 / 300 [ 20%]  (Warmup)
#> Chain 2: Iteration:  90 / 300 [ 30%]  (Warmup)
#> Chain 2: Iteration: 120 / 300 [ 40%]  (Warmup)
#> Chain 2: Iteration: 150 / 300 [ 50%]  (Warmup)
#> Chain 2: Iteration: 151 / 300 [ 50%]  (Sampling)
#> Chain 2: Iteration: 180 / 300 [ 60%]  (Sampling)
#> Chain 2: Iteration: 210 / 300 [ 70%]  (Sampling)
#> Chain 2: Iteration: 240 / 300 [ 80%]  (Sampling)
#> Chain 2: Iteration: 270 / 300 [ 90%]  (Sampling)
#> Chain 2: Iteration: 300 / 300 [100%]  (Sampling)
#> Chain 2: 
#> Chain 2:  Elapsed Time: 0.009 seconds (Warm-up)
#> Chain 2:                0.004 seconds (Sampling)
#> Chain 2:                0.013 seconds (Total)
#> Chain 2: 
#> 
#> SAMPLING FOR MODEL 'continuous' NOW (CHAIN 3).
#> Chain 3: 
#> Chain 3: Gradient evaluation took 9e-06 seconds
#> Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 0.09 seconds.
#> Chain 3: Adjust your expectations accordingly!
#> Chain 3: 
#> Chain 3: 
#> Chain 3: Iteration:   1 / 300 [  0%]  (Warmup)
#> Chain 3: Iteration:  30 / 300 [ 10%]  (Warmup)
#> Chain 3: Iteration:  60 / 300 [ 20%]  (Warmup)
#> Chain 3: Iteration:  90 / 300 [ 30%]  (Warmup)
#> Chain 3: Iteration: 120 / 300 [ 40%]  (Warmup)
#> Chain 3: Iteration: 150 / 300 [ 50%]  (Warmup)
#> Chain 3: Iteration: 151 / 300 [ 50%]  (Sampling)
#> Chain 3: Iteration: 180 / 300 [ 60%]  (Sampling)
#> Chain 3: Iteration: 210 / 300 [ 70%]  (Sampling)
#> Chain 3: Iteration: 240 / 300 [ 80%]  (Sampling)
#> Chain 3: Iteration: 270 / 300 [ 90%]  (Sampling)
#> Chain 3: Iteration: 300 / 300 [100%]  (Sampling)
#> Chain 3: 
#> Chain 3:  Elapsed Time: 0.008 seconds (Warm-up)
#> Chain 3:                0.004 seconds (Sampling)
#> Chain 3:                0.012 seconds (Total)
#> Chain 3: 
#> 
#> SAMPLING FOR MODEL 'continuous' NOW (CHAIN 4).
#> Chain 4: 
#> Chain 4: Gradient evaluation took 9e-06 seconds
#> Chain 4: 1000 transitions using 10 leapfrog steps per transition would take 0.09 seconds.
#> Chain 4: Adjust your expectations accordingly!
#> Chain 4: 
#> Chain 4: 
#> Chain 4: Iteration:   1 / 300 [  0%]  (Warmup)
#> Chain 4: Iteration:  30 / 300 [ 10%]  (Warmup)
#> Chain 4: Iteration:  60 / 300 [ 20%]  (Warmup)
#> Chain 4: Iteration:  90 / 300 [ 30%]  (Warmup)
#> Chain 4: Iteration: 120 / 300 [ 40%]  (Warmup)
#> Chain 4: Iteration: 150 / 300 [ 50%]  (Warmup)
#> Chain 4: Iteration: 151 / 300 [ 50%]  (Sampling)
#> Chain 4: Iteration: 180 / 300 [ 60%]  (Sampling)
#> Chain 4: Iteration: 210 / 300 [ 70%]  (Sampling)
#> Chain 4: Iteration: 240 / 300 [ 80%]  (Sampling)
#> Chain 4: Iteration: 270 / 300 [ 90%]  (Sampling)
#> Chain 4: Iteration: 300 / 300 [100%]  (Sampling)
#> Chain 4: 
#> Chain 4:  Elapsed Time: 0.007 seconds (Warm-up)
#> Chain 4:                0.004 seconds (Sampling)
#> Chain 4:                0.011 seconds (Total)
#> Chain 4: 
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> [1] TRUE
```
