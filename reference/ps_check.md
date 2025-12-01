# Graphical checks of the estimated survival function

This function plots the estimated marginal survival function based on
draws from the posterior predictive distribution of the fitted joint
model, and then overlays the Kaplan-Meier curve based on the observed
data.

## Usage

``` r
ps_check(
  object,
  check = "survival",
  limits = c("ci", "none"),
  draws = NULL,
  seed = NULL,
  xlab = NULL,
  ylab = NULL,
  ci_geom_args = NULL,
  ...
)
```

## Arguments

- object:

  A fitted model object returned by the
  [`stan_jm`](https://mc-stan.org/rstanarm/reference/stan_jm.md)
  modelling function. See
  [`stanreg-objects`](https://mc-stan.org/rstanarm/reference/stanreg-objects.md).

- check:

  The type of plot to show. Currently only "survival" is allowed, which
  compares the estimated marginal survival function under the joint
  model to the estimated Kaplan-Meier curve based on the observed data.

- limits:

  A quoted character string specifying the type of limits to include in
  the plot. Can be one of: `"ci"` for the Bayesian posterior uncertainty
  interval (often known as a credible interval); or `"none"` for no
  interval limits.

- draws:

  An integer indicating the number of MCMC draws to use to to estimate
  the survival function. The default and maximum number of draws is the
  size of the posterior sample.

- seed:

  An optional [`seed`](https://rdrr.io/r/base/Random.html) to use.

- xlab, ylab:

  An optional axis label passed to
  [`labs`](https://ggplot2.tidyverse.org/reference/labs.html).

- ci_geom_args:

  Optional arguments passed to
  [`geom_ribbon`](https://ggplot2.tidyverse.org/reference/geom_ribbon.html)
  and used to control features of the plotted interval limits. They
  should be supplied as a named list.

- ...:

  Optional arguments passed to
  [`geom_line`](https://ggplot2.tidyverse.org/reference/geom_path.html)
  and used to control features of the plotted trajectory.

## Value

A ggplot object that can be further customized using the ggplot2
package.

## See also

[`posterior_survfit`](https://mc-stan.org/rstanarm/reference/posterior_survfit.md)
for the estimated marginal or subject-specific survival function based
on draws of the model parameters from the posterior distribution,
[`posterior_predict`](https://mc-stan.org/rstanarm/reference/posterior_predict.stanreg.md)
for drawing from the posterior predictive distribution for the
longitudinal submodel, and
[`pp_check`](https://mc-stan.org/rstanarm/reference/pp_check.stanreg.md)
for graphical checks of the longitudinal submodel.

## Examples

``` r
if (.Platform$OS.type != "windows" || .Platform$r_arch != "i386") {
# \donttest{
if (!exists("example_jm")) example(example_jm)
# Compare estimated survival function to Kaplan-Meier curve
ps <- ps_check(example_jm)
ps + 
 ggplot2::scale_color_manual(values = c("red", "black")) + # change colors
 ggplot2::scale_size_manual(values = c(0.5, 3)) + # change line sizes 
 ggplot2::scale_fill_manual(values = c(NA, NA)) # remove fill
# }
}
```
