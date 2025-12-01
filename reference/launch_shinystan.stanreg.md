# Using the ShinyStan GUI with rstanarm models

The ShinyStan interface provides visual and numerical summaries of model
parameters and convergence diagnostics.

## Usage

``` r
# S3 method for class 'stanreg'
launch_shinystan(
  object,
  ppd = TRUE,
  seed = 1234,
  model_name = NULL,
  note = NULL,
  rstudio = getOption("shinystan.rstudio"),
  ...
)
```

## Arguments

- object:

  A fitted model object returned by one of the rstanarm modeling
  functions. See
  [`stanreg-objects`](https://mc-stan.org/rstanarm/reference/stanreg-objects.md).

- ppd:

  Should rstanarm draw from the posterior predictive distribution before
  launching ShinyStan? The default is `TRUE`, although for very large
  objects it can be convenient to set it to `FALSE` as drawing from the
  posterior predictive distribution can be time consuming. If `ppd` is
  `TRUE` then graphical posterior predictive checks are available when
  ShinyStan is launched.

- seed:

  Passed to
  [pp_check](https://mc-stan.org/rstanarm/reference/pp_check.stanreg.md)
  if `ppd` is `TRUE`.

- model_name, note:

  Optional arguments passed to
  [`as.shinystan`](https://mc-stan.org/shinystan/reference/as.shinystan.html).

- rstudio:

  Only relevant for 'RStudio' users. The default (`FALSE`) is to launch
  the app in the user's default web browser rather than the pop-up
  Viewer provided by 'RStudio'. Users can change the default to `TRUE`
  by setting the global option `options(shinystan.rstudio = TRUE)`.

- ...:

  Optional arguments passed to
  [`runApp`](https://rdrr.io/pkg/shiny/man/runApp.html).

## Details

The
[`launch_shinystan`](https://mc-stan.org/shinystan/reference/launch_shinystan.html)
function will accept a
[`stanreg`](https://mc-stan.org/rstanarm/reference/stanreg-objects.md)
object as input. Currently, almost any model fit using one of rstanarm's
model-fitting functions can be used with ShinyStan. The only exception
is that ShinyStan does not currently support rstanarm models fit using
`algorithm='optimizing'`. See the
[shinystan](https://mc-stan.org/shinystan/reference/shinystan-package.html)
package documentation for more information.

## Faster launch times

For some rstanarm models ShinyStan may take a very long time to launch.
If this is the case with one of your models you may be able to speed up
`launch_shinystan` in one of several ways:

- Prevent ShinyStan from preparing graphical posterior predictive
  checks::

  When used with a
  [`stanreg`](https://mc-stan.org/rstanarm/reference/stanreg-objects.md)
  object (rstanarm model object) ShinyStan will draw from the posterior
  predictive distribution and prepare graphical posterior predictive
  checks before launching. That way when you go to the PPcheck page the
  plots are immediately available. This can be time consuming for models
  fit to very large datasets and you can prevent this behavior by
  creating a shinystan object before calling `launch_shinystan`. To do
  this use
  [`as.shinystan`](https://mc-stan.org/shinystan/reference/as.shinystan.html)
  with optional argument `ppd` set to `FALSE` (see the Examples section
  below). When you then launch ShinyStan and go to the PPcheck page the
  plots will no longer be automatically generated and you will be
  presented with the standard interface requiring you to first specify
  the appropriate \\y\\ and \\yrep\\, which can be done for many but not
  all rstanarm models.

- Use a shinystan object::

  Even if you don't want to prevent ShinyStan from preparing graphical
  posterior predictive checks, first creating a shinystan object using
  [`as.shinystan`](https://mc-stan.org/shinystan/reference/as.shinystan.html)
  can reduce *future* launch times. That is, `launch_shinystan(sso)`
  will be faster than `launch_shinystan(fit)`, where `sso` is a
  shinystan object and `fit` is a stanreg object. It still may take some
  time for `as.shinystan` to create `sso` initially, but each time you
  subsequently call `launch_shinystan(sso)` it will reuse `sso` instead
  of internally creating a shinystan object every time. See the Examples
  section below.

## References

Gabry, J. , Simpson, D. , Vehtari, A. , Betancourt, M. and Gelman, A.
(2019), Visualization in Bayesian workflow. *J. R. Stat. Soc. A*, 182:
389-402. doi:10.1111/rssa.12378, [arXiv
preprint](https://arxiv.org/abs/1709.01449), [code on
GitHub](https://github.com/jgabry/bayes-vis-paper))

Muth, C., Oravecz, Z., and Gabry, J. (2018) User-friendly Bayesian
regression modeling: A tutorial with rstanarm and shinystan. *The
Quantitative Methods for Psychology*. 14(2), 99–119.
<https://www.tqmp.org/RegularArticles/vol14-2/p099/p099.pdf>

## Examples

``` r
if (.Platform$OS.type != "windows" || .Platform$r_arch != "i386") {
# \dontrun{
if (!exists("example_model")) example(example_model) 

# Launch the ShinyStan app without saving the resulting shinystan object
if (interactive()) launch_shinystan(example_model)

# Launch the ShinyStan app (saving resulting shinystan object as sso)
if (interactive()) sso <- launch_shinystan(example_model)

# First create shinystan object then call launch_shinystan
sso <- shinystan::as.shinystan(example_model)
if (interactive()) launch_shinystan(sso)

# Prevent ShinyStan from preparing graphical posterior predictive checks that
# can be time consuming. example_model is small enough that it won't matter
# much here but in general this can help speed up launch_shinystan
sso <- shinystan::as.shinystan(example_model, ppd = FALSE)
if (interactive()) launch_shinystan(sso)
# }
}
#> 
#> Hang on... preparing graphical posterior predictive checks for rstanarm model.
#> See help('shinystan', 'rstanarm') for how to disable this feature.
#> Note: in most cases the default test statistic 'mean' is too weak to detect anything of interest.
```
