# Example model

A model for use in rstanarm examples.

## Format

Calling `example("example_model")` will run the model in the Examples
section, below, and the resulting stanreg object will then be available
in the global environment. The `chains` and `iter` arguments are
specified to make this example be small in size. In practice, we
recommend that they be left unspecified in order to use the default
values (4 and 2000 respectively) or increased if there are convergence
problems. The `cores` argument is optional and on a multicore system,
the user may well want to set that equal to the number of chains being
executed.

## See also

[`cbpp`](https://rdrr.io/pkg/lme4/man/cbpp.html) for a description of
the data.

## Examples

``` r
if (.Platform$OS.type != "windows" || .Platform$r_arch != "i386") {
example_model <- 
  stan_glmer(cbind(incidence, size - incidence) ~ size + period + (1|herd),
             data = lme4::cbpp, family = binomial, QR = TRUE,
             # this next line is only to keep the example small in size!
             chains = 2, cores = 1, seed = 12345, iter = 1000, refresh = 0)
example_model
}
#> stan_glmer
#>  family:       binomial [logit]
#>  formula:      cbind(incidence, size - incidence) ~ size + period + (1 | herd)
#>  observations: 56
#> ------
#>             Median MAD_SD
#> (Intercept) -1.5    0.6  
#> size         0.0    0.0  
#> period2     -1.0    0.3  
#> period3     -1.1    0.4  
#> period4     -1.6    0.5  
#> 
#> Error terms:
#>  Groups Name        Std.Dev.
#>  herd   (Intercept) 0.79    
#> Num. levels: herd 15 
#> 
#> ------
#> * For help interpreting the printed output see ?print.stanreg
#> * For info on the priors used see ?prior_summary.stanreg
```
