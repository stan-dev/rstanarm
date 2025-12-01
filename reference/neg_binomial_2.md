# Family function for negative binomial GLMs

Specifies the information required to fit a Negative Binomial GLM in a
similar way to
[`negative.binomial`](https://rdrr.io/pkg/MASS/man/negative.binomial.html).
However, here the overdispersion parameter `theta` is not specified by
the user and always estimated (really the *reciprocal* of the dispersion
parameter is estimated). A call to this function can be passed to the
`family` argument of
[`stan_glm`](https://mc-stan.org/rstanarm/reference/stan_glm.md) or
[`stan_glmer`](https://mc-stan.org/rstanarm/reference/stan_glmer.md) to
estimate a Negative Binomial model. Alternatively, the
[`stan_glm.nb`](https://mc-stan.org/rstanarm/reference/stan_glm.md) and
[`stan_glmer.nb`](https://mc-stan.org/rstanarm/reference/stan_glmer.md)
wrapper functions may be used, which call `neg_binomial_2` internally.

## Usage

``` r
neg_binomial_2(link = "log")
```

## Arguments

- link:

  The same as for [`poisson`](https://rdrr.io/r/stats/family.html),
  typically a character vector of length one among `"log"`,
  `"identity"`, and `"sqrt"`.

## Value

An object of class [`family`](https://rdrr.io/r/stats/family.html) very
similar to that of [`poisson`](https://rdrr.io/r/stats/family.html) but
with a different family name.

## Examples

``` r
if (.Platform$OS.type != "windows" || .Platform$r_arch != "i386")
stan_glm(Days ~ Sex/(Age + Eth*Lrn), data = MASS::quine, seed = 123,
         family = neg_binomial_2, QR = TRUE, algorithm = "optimizing") 
#> stan_glm
#>  family:       neg_binomial_2 [log]
#>  formula:      Days ~ Sex/(Age + Eth * Lrn)
#>  observations: 146
#>  predictors:   14
#> ------
#>                 Median MAD_SD
#> (Intercept)      3.1    0.3  
#> SexM            -0.5    0.4  
#> SexF:AgeF1      -0.8    0.3  
#> SexM:AgeF1      -0.7    0.3  
#> SexF:AgeF2      -0.6    0.4  
#> SexM:AgeF2       0.6    0.3  
#> SexF:AgeF3      -0.4    0.4  
#> SexM:AgeF3       1.1    0.4  
#> SexF:EthN       -0.1    0.3  
#> SexM:EthN       -0.7    0.3  
#> SexF:LrnSL       1.0    0.3  
#> SexM:LrnSL       0.2    0.4  
#> SexF:EthN:LrnSL -1.4    0.4  
#> SexM:EthN:LrnSL  0.8    0.5  
#> 
#> Auxiliary parameter(s):
#>                       Median MAD_SD
#> reciprocal_dispersion 1.4    0.2   
#> 
#> ------
#> * For help interpreting the printed output see ?print.stanreg
#> * For info on the priors used see ?prior_summary.stanreg
                
# or, equivalently, call stan_glm.nb() without specifying the family
```
