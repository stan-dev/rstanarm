# Create a `draws` object from a `stanreg` object

Convert a `stanreg` object to a format supported by the
[posterior](https://mc-stan.org/posterior/reference/posterior-package.html)
package.

## Usage

``` r
# S3 method for class 'stanreg'
as_draws(x, ...)

# S3 method for class 'stanreg'
as_draws_matrix(x, ...)

# S3 method for class 'stanreg'
as_draws_array(x, ...)

# S3 method for class 'stanreg'
as_draws_df(x, ...)

# S3 method for class 'stanreg'
as_draws_list(x, ...)

# S3 method for class 'stanreg'
as_draws_rvars(x, ...)
```

## Arguments

- x:

  A `stanreg` object returned by one of the rstanarm modeling functions.

- ...:

  Arguments (e.g., `pars`, `regex_pars`) passed internally to
  [`as.matrix.stanreg`](https://mc-stan.org/rstanarm/reference/as.matrix.stanreg.md)
  or `as.array.stanreg`.

## Value

A `draws` object from the
[posterior](https://mc-stan.org/posterior/reference/posterior-package.html)
package. See the posterior package documentation and vignettes for
details on working with these objects.

## Details

To subset iterations, chains, or draws, use
[`subset_draws`](https://mc-stan.org/posterior/reference/subset_draws.html)
after making the `draws` object. To subset variables use `...` to pass
the `pars` and/or `regex_pars` arguments to `as.matrix.stanreg` or
`as.array.stanreg` (these are called internally by `as_draws.stanreg`),
or use
[`subset_draws`](https://mc-stan.org/posterior/reference/subset_draws.html)
after making the `draws` object.

## Examples

``` r
fit <- stan_glm(mpg ~ wt + as.factor(cyl), data = mtcars)
#> 
#> SAMPLING FOR MODEL 'continuous' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 2.1e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.21 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 1: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 1: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 1: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 1: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 1: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 1: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 1: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 1: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 1: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 1: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 1: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 0.049 seconds (Warm-up)
#> Chain 1:                0.045 seconds (Sampling)
#> Chain 1:                0.094 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'continuous' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 9e-06 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.09 seconds.
#> Chain 2: Adjust your expectations accordingly!
#> Chain 2: 
#> Chain 2: 
#> Chain 2: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 2: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 2: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 2: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 2: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 2: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 2: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 2: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 2: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 2: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 2: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 2: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 2: 
#> Chain 2:  Elapsed Time: 0.047 seconds (Warm-up)
#> Chain 2:                0.044 seconds (Sampling)
#> Chain 2:                0.091 seconds (Total)
#> Chain 2: 
#> 
#> SAMPLING FOR MODEL 'continuous' NOW (CHAIN 3).
#> Chain 3: 
#> Chain 3: Gradient evaluation took 8e-06 seconds
#> Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 0.08 seconds.
#> Chain 3: Adjust your expectations accordingly!
#> Chain 3: 
#> Chain 3: 
#> Chain 3: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 3: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 3: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 3: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 3: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 3: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 3: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 3: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 3: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 3: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 3: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 3: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 3: 
#> Chain 3:  Elapsed Time: 0.048 seconds (Warm-up)
#> Chain 3:                0.048 seconds (Sampling)
#> Chain 3:                0.096 seconds (Total)
#> Chain 3: 
#> 
#> SAMPLING FOR MODEL 'continuous' NOW (CHAIN 4).
#> Chain 4: 
#> Chain 4: Gradient evaluation took 9e-06 seconds
#> Chain 4: 1000 transitions using 10 leapfrog steps per transition would take 0.09 seconds.
#> Chain 4: Adjust your expectations accordingly!
#> Chain 4: 
#> Chain 4: 
#> Chain 4: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 4: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 4: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 4: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 4: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 4: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 4: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 4: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 4: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 4: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 4: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 4: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 4: 
#> Chain 4:  Elapsed Time: 0.047 seconds (Warm-up)
#> Chain 4:                0.04 seconds (Sampling)
#> Chain 4:                0.087 seconds (Total)
#> Chain 4: 
as_draws_matrix(fit) # matrix format combines all chains 
#> # A draws_matrix: 4000 iterations, 1 chains, and 5 variables
#>     variable
#> draw (Intercept)   wt as.factor(cyl)6 as.factor(cyl)8 sigma
#>   1           35 -3.1            -7.4            -7.2   2.5
#>   2           32 -2.5            -4.4            -6.8   2.4
#>   3           33 -2.8            -3.4            -6.5   2.6
#>   4           36 -3.6            -5.1            -5.9   2.8
#>   5           37 -3.8            -5.2            -6.5   2.9
#>   6           31 -2.6            -4.0            -5.7   2.3
#>   7           37 -4.5            -3.3            -4.5   2.8
#>   8           33 -2.8            -3.5            -4.9   2.8
#>   9           34 -3.7            -1.6            -5.3   2.3
#>   10          36 -3.5            -4.8            -6.7   3.2
#> # ... with 3990 more draws
as_draws_df(fit, regex_pars = "cyl")
#> # A draws_df: 1000 iterations, 4 chains, and 2 variables
#>    as.factor(cyl)6 as.factor(cyl)8
#> 1             -7.4            -7.2
#> 2             -4.4            -6.8
#> 3             -3.4            -6.5
#> 4             -5.1            -5.9
#> 5             -5.2            -6.5
#> 6             -4.0            -5.7
#> 7             -3.3            -4.5
#> 8             -3.5            -4.9
#> 9             -1.6            -5.3
#> 10            -4.8            -6.7
#> # ... with 3990 more draws
#> # ... hidden reserved variables {'.chain', '.iteration', '.draw'}
posterior::summarize_draws(as_draws_array(fit))
#> # A tibble: 5 × 10
#>   variable         mean median    sd   mad    q5   q95  rhat ess_bulk ess_tail
#>   <chr>           <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>    <dbl>    <dbl>
#> 1 (Intercept)     34.0   33.9  1.96  1.93  30.8  37.2   1.00    2819.    3016.
#> 2 wt              -3.21  -3.20 0.776 0.749 -4.51 -1.96  1.00    2018.    2385.
#> 3 as.factor(cyl)6 -4.24  -4.25 1.48  1.46  -6.65 -1.82  1.00    1924.    2547.
#> 4 as.factor(cyl)8 -6.04  -6.05 1.72  1.68  -8.88 -3.18  1.00    1767.    2050.
#> 5 sigma            2.65   2.61 0.374 0.360  2.12  3.34  1.00    2744.    2655.
```
