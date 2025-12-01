# Example joint longitudinal and time-to-event model

A model for use in the rstanarm examples related to
[`stan_jm`](https://mc-stan.org/rstanarm/reference/stan_jm.md).

## Format

Calling `example("example_jm")` will run the model in the Examples
section, below, and the resulting stanmvreg object will then be
available in the global environment. The `chains` and `iter` arguments
are specified to make this example be small in size. In practice, we
recommend that they be left unspecified in order to use the default
values or increased if there are convergence problems. The `cores`
argument is optional and on a multicore system, the user may well want
to set that equal to the number of chains being executed.

## Examples

``` r
  # set.seed(123)
  if (.Platform$OS.type != "windows" || .Platform$r_arch !="i386")
  example_jm <- 
     stan_jm(formulaLong = logBili ~ year + (1 | id), 
             dataLong = pbcLong[1:101,],
             formulaEvent = survival::Surv(futimeYears, death) ~ sex + trt, 
             dataEvent = pbcSurv[1:15,],
             time_var = "year",
             # this next line is only to keep the example small in size!
             chains = 1, seed = 12345, iter = 100, refresh = 0)
#> Loading required namespace: data.table
#> Fitting a univariate joint model.
#> 
#> Please note the warmup may be much slower than later iterations!
#> Warning: The largest R-hat is 1.07, indicating chains have not mixed.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#r-hat
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess

```
