# Extract the posterior sample

For models fit using MCMC (`algorithm="sampling"`), the posterior sample
—the post-warmup draws from the posterior distribution— can be extracted
from a fitted model object as a matrix, data frame, or array. The
`as.matrix` and `as.data.frame` methods merge all chains together,
whereas the `as.array` method keeps the chains separate. For models fit
using optimization (`"optimizing"`) or variational inference
(`"meanfield"` or `"fullrank"`), there is no posterior sample but rather
a matrix (or data frame) of 1000 draws from either the asymptotic
multivariate Gaussian sampling distribution of the parameters or the
variational approximation to the posterior distribution.

## Usage

``` r
# S3 method for class 'stanreg'
as.matrix(x, ..., pars = NULL, regex_pars = NULL)

# S3 method for class 'stanreg'
as.array(x, ..., pars = NULL, regex_pars = NULL)

# S3 method for class 'stanreg'
as.data.frame(x, ..., pars = NULL, regex_pars = NULL)
```

## Arguments

- x:

  A fitted model object returned by one of the rstanarm modeling
  functions. See
  [`stanreg-objects`](https://mc-stan.org/rstanarm/reference/stanreg-objects.md).

- ...:

  Ignored.

- pars:

  An optional character vector of parameter names.

- regex_pars:

  An optional character vector of [regular
  expressions](https://rdrr.io/r/base/grep.html) to use for parameter
  selection. `regex_pars` can be used in place of `pars` or in addition
  to `pars`. Currently, all functions that accept a `regex_pars`
  argument ignore it for models fit using optimization.

## Value

A matrix, data.frame, or array, the dimensions of which depend on `pars`
and `regex_pars`, as well as the model and estimation algorithm (see the
Description section above).

## See also

[`stanreg-draws-formats`](https://mc-stan.org/rstanarm/reference/stanreg-draws-formats.md),
[`stanreg-methods`](https://mc-stan.org/rstanarm/reference/stanreg-methods.md)

## Examples

``` r
if (.Platform$OS.type != "windows" || .Platform$r_arch != "i386") {
# \donttest{
if (!exists("example_model")) example(example_model)
# Extract posterior sample after MCMC
draws <- as.matrix(example_model)
print(dim(draws))

# For example, we can see that the median of the draws for the intercept 
# is the same as the point estimate rstanarm uses
print(median(draws[, "(Intercept)"]))
print(example_model$coefficients[["(Intercept)"]])

# The as.array method keeps the chains separate
draws_array <- as.array(example_model)
print(dim(draws_array)) # iterations x chains x parameters

# Extract draws from asymptotic Gaussian sampling distribution 
# after optimization
fit <- stan_glm(mpg ~ wt, data = mtcars, algorithm = "optimizing")
draws <- as.data.frame(fit)
print(colnames(draws))
print(nrow(draws)) # 1000 draws are taken

# Extract draws from variational approximation to the posterior distribution
fit2 <- update(fit, algorithm = "meanfield")
draws <- as.data.frame(fit2, pars = "wt")
print(colnames(draws))
print(nrow(draws)) # 1000 draws are taken
# }
}
#> 
#> exmpl_> if (.Platform$OS.type != "windows" || .Platform$r_arch != "i386") {
#> exmpl_+ example_model <- 
#> exmpl_+   stan_glmer(cbind(incidence, size - incidence) ~ size + period + (1|herd),
#> exmpl_+              data = lme4::cbpp, family = binomial, QR = TRUE,
#> exmpl_+              # this next line is only to keep the example small in size!
#> exmpl_+              chains = 2, cores = 1, seed = 12345, iter = 1000, refresh = 0)
#> exmpl_+ example_model
#> exmpl_+ }
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
#>  herd   (Intercept) 0.76    
#> Num. levels: herd 15 
#> 
#> ------
#> * For help interpreting the printed output see ?print.stanreg
#> * For info on the priors used see ?prior_summary.stanreg
#> [1] 1000   21
#> [1] -1.515377
#> [1] -1.515377
#> [1] 500   2  21
#> [1] "(Intercept)" "wt"          "sigma"      
#> [1] 1000
#> Chain 1: ------------------------------------------------------------
#> Chain 1: EXPERIMENTAL ALGORITHM:
#> Chain 1:   This procedure has not been thoroughly tested and may be unstable
#> Chain 1:   or buggy. The interface is subject to change.
#> Chain 1: ------------------------------------------------------------
#> Chain 1: 
#> Chain 1: 
#> Chain 1: 
#> Chain 1: Gradient evaluation took 2.3e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.23 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: Begin eta adaptation.
#> Chain 1: Iteration:   1 / 250 [  0%]  (Adaptation)
#> Chain 1: Iteration:  50 / 250 [ 20%]  (Adaptation)
#> Chain 1: Iteration: 100 / 250 [ 40%]  (Adaptation)
#> Chain 1: Iteration: 150 / 250 [ 60%]  (Adaptation)
#> Chain 1: Iteration: 200 / 250 [ 80%]  (Adaptation)
#> Chain 1: Success! Found best value [eta = 1] earlier than expected.
#> Chain 1: 
#> Chain 1: Begin stochastic gradient ascent.
#> Chain 1:   iter             ELBO   delta_ELBO_mean   delta_ELBO_med   notes 
#> Chain 1:    100         -199.711             1.000            1.000
#> Chain 1:    200         -145.781             0.685            1.000
#> Chain 1:    300         -134.939             0.483            0.370
#> Chain 1:    400         -124.247             0.384            0.370
#> Chain 1:    500         -114.013             0.325            0.090
#> Chain 1:    600          -97.593             0.299            0.168
#> Chain 1:    700          -89.716             0.269            0.090
#> Chain 1:    800          -89.115             0.236            0.090
#> Chain 1:    900          -89.614             0.210            0.088
#> Chain 1:   1000          -89.450             0.190            0.088
#> Chain 1:   1100          -89.124             0.090            0.086
#> Chain 1:   1200          -89.136             0.053            0.080
#> Chain 1:   1300          -89.459             0.045            0.007
#> Chain 1:   1400          -89.339             0.037            0.006
#> Chain 1:   1500          -89.409             0.028            0.004
#> Chain 1:   1600          -89.057             0.012            0.004
#> Chain 1:   1700          -89.342             0.003            0.004
#> Chain 1:   1800          -88.968             0.003            0.004
#> Chain 1:   1900          -89.288             0.003            0.004
#> Chain 1:   2000          -89.399             0.003            0.004
#> Chain 1:   2100          -89.904             0.003            0.004
#> Chain 1:   2200          -89.471             0.003            0.004
#> Chain 1:   2300          -89.648             0.003            0.004
#> Chain 1:   2400          -89.044             0.004            0.004
#> Chain 1:   2500          -89.182             0.004            0.004
#> Chain 1:   2600          -89.249             0.003            0.004
#> Chain 1:   2700          -89.307             0.003            0.004
#> Chain 1:   2800          -89.656             0.003            0.004
#> Chain 1:   2900          -89.055             0.003            0.004
#> Chain 1:   3000          -89.374             0.004            0.004
#> Chain 1:   3100          -89.347             0.003            0.004
#> Chain 1:   3200          -89.044             0.003            0.003
#> Chain 1:   3300          -89.282             0.003            0.003
#> Chain 1:   3400          -89.180             0.002            0.003
#> Chain 1:   3500          -89.445             0.003            0.003
#> Chain 1:   3600          -88.930             0.003            0.003
#> Chain 1:   3700          -89.128             0.003            0.003
#> Chain 1:   3800          -89.273             0.003            0.003
#> Chain 1:   3900          -89.221             0.002            0.003
#> Chain 1:   4000          -88.955             0.002            0.003
#> Chain 1:   4100          -89.087             0.002            0.003
#> Chain 1:   4200          -89.442             0.003            0.003
#> Chain 1:   4300          -89.457             0.002            0.002
#> Chain 1:   4400          -89.148             0.003            0.003
#> Chain 1:   4500          -89.370             0.002            0.002
#> Chain 1:   4600          -88.966             0.002            0.002
#> Chain 1:   4700          -89.159             0.002            0.002
#> Chain 1:   4800          -89.471             0.003            0.003
#> Chain 1:   4900          -89.269             0.003            0.003
#> Chain 1:   5000          -89.090             0.003            0.002
#> Chain 1:   5100          -89.190             0.003            0.002
#> Chain 1:   5200          -89.215             0.002            0.002
#> Chain 1:   5300          -89.220             0.002            0.002
#> Chain 1:   5400          -89.496             0.002            0.002
#> Chain 1:   5500          -89.219             0.002            0.002
#> Chain 1:   5600          -89.212             0.002            0.002
#> Chain 1:   5700          -89.269             0.002            0.002
#> Chain 1:   5800          -89.702             0.002            0.002
#> Chain 1:   5900          -89.541             0.002            0.002
#> Chain 1:   6000          -89.517             0.002            0.001
#> Chain 1:   6100          -89.470             0.001            0.001
#> Chain 1:   6200          -89.060             0.002            0.002
#> Chain 1:   6300          -89.467             0.002            0.003
#> Chain 1:   6400          -89.174             0.002            0.003
#> Chain 1:   6500          -89.106             0.002            0.002
#> Chain 1:   6600          -89.083             0.002            0.002
#> Chain 1:   6700          -89.312             0.002            0.003
#> Chain 1:   6800          -89.392             0.002            0.002
#> Chain 1:   6900          -89.270             0.002            0.001
#> Chain 1:   7000          -88.979             0.002            0.003
#> Chain 1:   7100          -89.230             0.002            0.003
#> Chain 1:   7200          -88.998             0.002            0.003
#> Chain 1:   7300          -89.365             0.002            0.003
#> Chain 1:   7400          -89.208             0.002            0.003
#> Chain 1:   7500          -89.108             0.002            0.003
#> Chain 1:   7600          -89.129             0.002            0.003
#> Chain 1:   7700          -88.939             0.002            0.002
#> Chain 1:   7800          -89.534             0.003            0.003
#> Chain 1:   7900          -88.906             0.003            0.003
#> Chain 1:   8000          -89.087             0.003            0.003
#> Chain 1:   8100          -88.991             0.003            0.002
#> Chain 1:   8200          -88.983             0.003            0.002
#> Chain 1:   8300          -89.323             0.003            0.002
#> Chain 1:   8400          -89.016             0.003            0.002
#> Chain 1:   8500          -89.293             0.003            0.003
#> Chain 1:   8600          -89.358             0.003            0.003
#> Chain 1:   8700          -89.161             0.003            0.003
#> Chain 1:   8800          -88.894             0.003            0.003
#> Chain 1:   8900          -89.349             0.002            0.003
#> Chain 1:   9000          -89.346             0.002            0.003
#> Chain 1:   9100          -89.402             0.002            0.003
#> Chain 1:   9200          -89.223             0.002            0.003
#> Chain 1:   9300          -89.061             0.002            0.002
#> Chain 1:   9400          -88.958             0.002            0.002
#> Chain 1:   9500          -89.476             0.002            0.002
#> Chain 1:   9600          -89.143             0.003            0.002
#> Chain 1:   9700          -89.275             0.002            0.002
#> Chain 1:   9800          -88.954             0.003            0.002
#> Chain 1:   9900          -89.765             0.003            0.002
#> Chain 1:   10000          -89.301             0.003            0.004
#> Chain 1: Informational Message: The maximum number of iterations is reached! The algorithm may not have converged.
#> Chain 1: This variational approximation is not guaranteed to be meaningful.
#> Chain 1: 
#> Chain 1: Drawing a sample of size 1000 from the approximate posterior... 
#> Chain 1: COMPLETED.
#> Warning: Pareto k diagnostic value is 0.9. Resampling is unreliable. Increasing the number of draws or decreasing tol_rel_obj may help.
#> [1] "wt"
#> [1] 1000
```
