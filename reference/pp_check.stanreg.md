# Graphical posterior predictive checks

Interface to the
[PPC](https://mc-stan.org/bayesplot/reference/PPC-overview.html)
(posterior predictive checking) module in the
[bayesplot](https://mc-stan.org/bayesplot/reference/bayesplot-package.html)
package, providing various plots comparing the observed outcome variable
\\y\\ to simulated datasets \\y^{rep}\\ from the posterior predictive
distribution. The `pp_check` method for
[stanreg-objects](https://mc-stan.org/rstanarm/reference/stanreg-objects.md)
prepares the arguments required for the specified bayesplot PPC plotting
function and then calls that function. It is also straightforward to use
the functions from the bayesplot package directly rather than via the
`pp_check` method. Examples of both are given below.

## Usage

``` r
# S3 method for class 'stanreg'
pp_check(object, plotfun = "dens_overlay", nreps = NULL, seed = NULL, ...)
```

## Arguments

- object:

  A fitted model object returned by one of the rstanarm modeling
  functions. See
  [`stanreg-objects`](https://mc-stan.org/rstanarm/reference/stanreg-objects.md).

- plotfun:

  A character string naming the bayesplot
  [PPC](https://mc-stan.org/bayesplot/reference/PPC-overview.html)
  function to use. The default is to call
  [`ppc_dens_overlay`](https://mc-stan.org/bayesplot/reference/PPC-distributions.html).
  `plotfun` can be specified either as the full name of a bayesplot
  plotting function (e.g. `"ppc_hist"`) or can be abbreviated to the
  part of the name following the `"ppc_"` prefix (e.g. `"hist"`). To get
  the names of all available PPC functions see
  [`available_ppc`](https://mc-stan.org/bayesplot/reference/available_ppc.html).

- nreps:

  The number of \\y^{rep}\\ datasets to generate from the [posterior
  predictive
  distribution](https://mc-stan.org/rstanarm/reference/posterior_predict.stanreg.md)
  and show in the plots. The default depends on `plotfun`. For functions
  that plot each `yrep` dataset separately (e.g. `ppc_hist`), `nreps`
  defaults to a small value to make the plots readable. For functions
  that overlay many `yrep` datasets (e.g., `ppc_dens_overlay`) a larger
  number is used by default, and for other functions (e.g. `ppc_stat`)
  the default is to set `nreps` equal to the posterior sample size.

- seed:

  An optional [`seed`](https://rdrr.io/r/base/Random.html) to pass to
  [`posterior_predict`](https://mc-stan.org/rstanarm/reference/posterior_predict.stanreg.md).

- ...:

  Additonal arguments passed to the
  [bayesplot](https://mc-stan.org/bayesplot/reference/bayesplot-package.html)
  function called. For many plotting functions `...` is optional,
  however for functions that require a `group` or `x` argument, these
  arguments should be specified in `...`. If specifying `group` and/or
  `x`, they can be provided as either strings naming variables (in which
  case they are searched for in the model frame) or as vectors
  containing the actual values of the variables. See the **Examples**
  section, below.

## Value

`pp_check` returns a ggplot object that can be further customized using
the ggplot2 package.

## Note

For binomial data, plots of \\y\\ and \\y^{rep}\\ show the proportion of
'successes' rather than the raw count. Also for binomial models see
[`ppc_error_binned`](https://mc-stan.org/bayesplot/reference/PPC-errors.html)
for binned residual plots.

## References

Gelman, A., Carlin, J. B., Stern, H. S., Dunson, D. B., Vehtari, A., and
Rubin, D. B. (2013). *Bayesian Data Analysis.* Chapman & Hall/CRC Press,
London, third edition. (Ch. 6)

Gabry, J. , Simpson, D. , Vehtari, A. , Betancourt, M. and Gelman, A.
(2019), Visualization in Bayesian workflow. *J. R. Stat. Soc. A*, 182:
389-402. doi:10.1111/rssa.12378, [arXiv
preprint](https://arxiv.org/abs/1709.01449), [code on
GitHub](https://github.com/jgabry/bayes-vis-paper))

## See also

- The vignettes in the bayesplot package for many examples. Examples of
  posterior predictive checks can also be found in the rstanarm
  vignettes and demos.

- [`PPC-overview`](https://mc-stan.org/bayesplot/reference/PPC-overview.html)
  (bayesplot) for links to the documentation for all the available
  plotting functions.

- [`posterior_predict`](https://mc-stan.org/rstanarm/reference/posterior_predict.stanreg.md)
  for drawing from the posterior predictive distribution.

- [`color_scheme_set`](https://mc-stan.org/bayesplot/reference/bayesplot-colors.html)
  to change the color scheme of the plots.

## Examples

``` r
if (.Platform$OS.type != "windows" || .Platform$r_arch != "i386") {
fit <- stan_glmer(
  mpg ~ wt + am + (1|cyl),
  data = mtcars,
  iter = 400, # iter and chains small just to keep example quick
  chains = 2,
  refresh = 0
)

# Compare distribution of y to distributions of multiple yrep datasets
pp_check(fit)
pp_check(fit, plotfun = "boxplot", nreps = 10, notch = FALSE)
pp_check(fit, plotfun = "hist", nreps = 3)

# \donttest{
# Same plot (up to RNG noise) using bayesplot package directly
bayesplot::ppc_hist(y = mtcars$mpg, yrep = posterior_predict(fit, draws = 3))

# Check histograms of test statistics by level of grouping variable 'cyl'
pp_check(fit, plotfun = "stat_grouped", stat = "median", group = "cyl")

# Defining a custom test statistic
q25 <- function(y) quantile(y, probs = 0.25)
pp_check(fit, plotfun = "stat_grouped", stat = "q25", group = "cyl")

# Scatterplot of two test statistics
pp_check(fit, plotfun = "stat_2d", stat = c("mean", "sd"))

# Scatterplot of y vs. average yrep
pp_check(fit, plotfun = "scatter_avg") # y vs. average yrep
# Same plot (up to RNG noise) using bayesplot package directly
bayesplot::ppc_scatter_avg(y = mtcars$mpg, yrep = posterior_predict(fit))

# Scatterplots of y vs. several individual yrep datasets
pp_check(fit, plotfun = "scatter", nreps = 3)

# Same plot (up to RNG noise) using bayesplot package directly
bayesplot::ppc_scatter(y = mtcars$mpg, yrep = posterior_predict(fit, draws = 3))

# yrep intervals with y points overlaid
# by default 1:length(y) used on x-axis but can also specify an x variable
pp_check(fit, plotfun = "intervals")
pp_check(fit, plotfun = "intervals", x = "wt") + ggplot2::xlab("wt")

# Same plot (up to RNG noise) using bayesplot package directly
bayesplot::ppc_intervals(y = mtcars$mpg, yrep = posterior_predict(fit),
                         x = mtcars$wt) + ggplot2::xlab("wt")

# predictive errors
pp_check(fit, plotfun = "error_hist", nreps = 6)
pp_check(fit, plotfun = "error_scatter_avg_vs_x", x = "wt") +
  ggplot2::xlab("wt")

# Example of a PPC for ordinal models (stan_polr)
fit2 <- stan_polr(tobgp ~ agegp, data = esoph, method = "probit",
                  prior = R2(0.2, "mean"), init_r = 0.1,
                  refresh = 0)
pp_check(fit2, plotfun = "bars", nreps = 500, prob = 0.5)
pp_check(fit2, plotfun = "bars_grouped", group = esoph$agegp,
         nreps = 500, prob = 0.5)
# }
}
#> Warning: There were 3 divergent transitions after warmup. See
#> https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#> to find out why this is a problem and how to eliminate them.
#> Warning: Examine the pairs() plot to diagnose sampling problems
#> Warning: The largest R-hat is 1.07, indicating chains have not mixed.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#r-hat
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
#> Warning: Markov chains did not converge! Do not analyze results!
#> Error in get(as.character(FUN), mode = "function", envir = envir): object 'q25' of mode 'function' was not found
```
