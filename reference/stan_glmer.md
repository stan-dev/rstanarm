# Bayesian generalized linear models with group-specific terms via Stan

Bayesian inference for GLMs with group-specific coefficients that have
unknown covariance matrices with flexible priors.

## Usage

``` r
stan_glmer(
  formula,
  data = NULL,
  family = gaussian,
  subset,
  weights,
  na.action = getOption("na.action", "na.omit"),
  offset,
  contrasts = NULL,
  ...,
  prior = default_prior_coef(family),
  prior_intercept = default_prior_intercept(family),
  prior_aux = exponential(autoscale = TRUE),
  prior_covariance = decov(),
  prior_PD = FALSE,
  algorithm = c("sampling", "meanfield", "fullrank"),
  adapt_delta = NULL,
  QR = FALSE,
  sparse = FALSE
)

stan_lmer(
  formula,
  data = NULL,
  subset,
  weights,
  na.action = getOption("na.action", "na.omit"),
  offset,
  contrasts = NULL,
  ...,
  prior = default_prior_coef(family),
  prior_intercept = default_prior_intercept(family),
  prior_aux = exponential(autoscale = TRUE),
  prior_covariance = decov(),
  prior_PD = FALSE,
  algorithm = c("sampling", "meanfield", "fullrank"),
  adapt_delta = NULL,
  QR = FALSE
)

stan_glmer.nb(
  formula,
  data = NULL,
  subset,
  weights,
  na.action = getOption("na.action", "na.omit"),
  offset,
  contrasts = NULL,
  link = "log",
  ...,
  prior = default_prior_coef(family),
  prior_intercept = default_prior_intercept(family),
  prior_aux = exponential(autoscale = TRUE),
  prior_covariance = decov(),
  prior_PD = FALSE,
  algorithm = c("sampling", "meanfield", "fullrank"),
  adapt_delta = NULL,
  QR = FALSE
)
```

## Arguments

- formula, data:

  Same as for [`glmer`](https://rdrr.io/pkg/lme4/man/glmer.html). *We
  strongly advise against omitting the `data` argument*. Unless `data`
  is specified (and is a data frame) many post-estimation functions
  (including `update`, `loo`, `kfold`) are not guaranteed to work
  properly.

- family:

  Same as for [`glmer`](https://rdrr.io/pkg/lme4/man/glmer.html) except
  it is also possible to use `family=mgcv::betar` to estimate a Beta
  regression with `stan_glmer`.

- subset, weights, offset:

  Same as [`glm`](https://rdrr.io/r/stats/glm.html).

- na.action, contrasts:

  Same as [`glm`](https://rdrr.io/r/stats/glm.html), but rarely
  specified.

- ...:

  For `stan_glmer`, further arguments passed to `sampling` (e.g. `iter`,
  `chains`, `cores`, etc.) or to `vb` (if `algorithm` is `"meanfield"`
  or `"fullrank"`). For `stan_lmer` and `stan_glmer.nb`, `...` should
  also contain all relevant arguments to pass to `stan_glmer` (except
  `family`).

- prior:

  The prior distribution for the (non-hierarchical) regression
  coefficients.

  The default priors are described in the vignette [*Prior Distributions
  for rstanarm
  Models*](https://mc-stan.org/rstanarm/articles/priors.html). If not
  using the default, `prior` should be a call to one of the various
  functions provided by rstanarm for specifying priors. The subset of
  these functions that can be used for the prior on the coefficients can
  be grouped into several "families":

  |                                 |                                 |
  |---------------------------------|---------------------------------|
  | **Family**                      | **Functions**                   |
  | *Student t family*              | `normal`, `student_t`, `cauchy` |
  | *Hierarchical shrinkage family* | `hs`, `hs_plus`                 |
  | *Laplace family*                | `laplace`, `lasso`              |
  | *Product normal family*         | `product_normal`                |

  See the [priors help
  page](https://mc-stan.org/rstanarm/reference/priors.md) for details on
  the families and how to specify the arguments for all of the functions
  in the table above. To omit a prior —i.e., to use a flat (improper)
  uniform prior— `prior` can be set to `NULL`, although this is rarely a
  good idea.

  **Note:** Unless `QR=TRUE`, if `prior` is from the Student t family or
  Laplace family, and if the `autoscale` argument to the function used
  to specify the prior (e.g.
  [`normal`](https://mc-stan.org/rstanarm/reference/priors.md)) is left
  at its default and recommended value of `TRUE`, then the default or
  user-specified prior scale(s) may be adjusted internally based on the
  scales of the predictors. See the [priors help
  page](https://mc-stan.org/rstanarm/reference/priors.md) and the *Prior
  Distributions* vignette for details on the rescaling and the
  [`prior_summary`](https://mc-stan.org/rstanarm/reference/prior_summary.stanreg.md)
  function for a summary of the priors used for a particular model.

- prior_intercept:

  The prior distribution for the intercept (after centering all
  predictors, see note below).

  The default prior is described in the vignette [*Prior Distributions
  for rstanarm
  Models*](https://mc-stan.org/rstanarm/articles/priors.html). If not
  using the default, `prior_intercept` can be a call to `normal`,
  `student_t` or `cauchy`. See the [priors help
  page](https://mc-stan.org/rstanarm/reference/priors.md) for details on
  these functions. To omit a prior on the intercept —i.e., to use a flat
  (improper) uniform prior— `prior_intercept` can be set to `NULL`.

  **Note:** If using a dense representation of the design matrix —i.e.,
  if the `sparse` argument is left at its default value of `FALSE`— then
  the prior distribution for the intercept is set so it applies to the
  value *when all predictors are centered* (you don't need to manually
  center them). This is explained further in \[Prior Distributions for
  rstanarm Models\](https://mc-stan.org/rstanarm/articles/priors.html)
  If you prefer to specify a prior on the intercept without the
  predictors being auto-centered, then you have to omit the intercept
  from the [`formula`](https://rdrr.io/r/stats/formula.html) and include
  a column of ones as a predictor, in which case some element of `prior`
  specifies the prior on it, rather than `prior_intercept`. Regardless
  of how `prior_intercept` is specified, the reported *estimates* of the
  intercept always correspond to a parameterization without centered
  predictors (i.e., same as in `glm`).

- prior_aux:

  The prior distribution for the "auxiliary" parameter (if applicable).
  The "auxiliary" parameter refers to a different parameter depending on
  the `family`. For Gaussian models `prior_aux` controls `"sigma"`, the
  error standard deviation. For negative binomial models `prior_aux`
  controls `"reciprocal_dispersion"`, which is similar to the `"size"`
  parameter of [`rnbinom`](https://rdrr.io/r/stats/NegBinomial.html):
  smaller values of `"reciprocal_dispersion"` correspond to greater
  dispersion. For gamma models `prior_aux` sets the prior on to the
  `"shape"` parameter (see e.g.,
  [`rgamma`](https://rdrr.io/r/stats/GammaDist.html)), and for
  inverse-Gaussian models it is the so-called `"lambda"` parameter
  (which is essentially the reciprocal of a scale parameter). Binomial
  and Poisson models do not have auxiliary parameters.

  The default prior is described in the vignette [*Prior Distributions
  for rstanarm
  Models*](https://mc-stan.org/rstanarm/articles/priors.html). If not
  using the default, `prior_aux` can be a call to `exponential` to use
  an exponential distribution, or `normal`, `student_t` or `cauchy`,
  which results in a half-normal, half-t, or half-Cauchy prior. See
  [`priors`](https://mc-stan.org/rstanarm/reference/priors.md) for
  details on these functions. To omit a prior —i.e., to use a flat
  (improper) uniform prior— set `prior_aux` to `NULL`.

- prior_covariance:

  Cannot be `NULL`; see
  [`decov`](https://mc-stan.org/rstanarm/reference/priors.md) for more
  information about the default arguments.

- prior_PD:

  A logical scalar (defaulting to `FALSE`) indicating whether to draw
  from the prior predictive distribution instead of conditioning on the
  outcome.

- algorithm:

  A string (possibly abbreviated) indicating the estimation approach to
  use. Can be `"sampling"` for MCMC (the default), `"optimizing"` for
  optimization, `"meanfield"` for variational inference with independent
  normal distributions, or `"fullrank"` for variational inference with a
  multivariate normal distribution. See
  [`rstanarm-package`](https://mc-stan.org/rstanarm/reference/rstanarm-package.md)
  for more details on the estimation algorithms. NOTE: not all fitting
  functions support all four algorithms.

- adapt_delta:

  Only relevant if `algorithm="sampling"`. See the
  [adapt_delta](https://mc-stan.org/rstanarm/reference/adapt_delta.md)
  help page for details.

- QR:

  A logical scalar defaulting to `FALSE`, but if `TRUE` applies a scaled
  [`qr`](https://rdrr.io/r/base/qr.html) decomposition to the design
  matrix. The transformation does not change the likelihood of the data
  but is recommended for computational reasons when there are multiple
  predictors. See the
  [QR-argument](https://mc-stan.org/rstanarm/reference/QR-argument.md)
  documentation page for details on how rstanarm does the transformation
  and important information about how to interpret the prior
  distributions of the model parameters when using `QR=TRUE`.

- sparse:

  A logical scalar (defaulting to `FALSE`) indicating whether to use a
  sparse representation of the design (X) matrix. If `TRUE`, the the
  design matrix is not centered (since that would destroy the sparsity)
  and likewise it is not possible to specify both `QR = TRUE` and
  `sparse = TRUE`. Depending on how many zeros there are in the design
  matrix, setting `sparse = TRUE` may make the code run faster and can
  consume much less RAM.

- link:

  For `stan_glmer.nb` only, the link function to use. See
  [`neg_binomial_2`](https://mc-stan.org/rstanarm/reference/neg_binomial_2.md).

## Value

A [stanreg](https://mc-stan.org/rstanarm/reference/stanreg-objects.md)
object is returned for `stan_glmer, stan_lmer, stan_glmer.nb`.

A list with classes `stanreg`, `glm`, `lm`, and `lmerMod`. The
conventions for the parameter names are the same as in the lme4 package
with the addition that the standard deviation of the errors is called
`sigma` and the variance-covariance matrix of the group-specific
deviations from the common parameters is called `Sigma`, even if this
variance-covariance matrix only has one row and one column (in which
case it is just the group-level variance).

## Details

The `stan_glmer` function is similar in syntax to
[`glmer`](https://rdrr.io/pkg/lme4/man/glmer.html) but rather than
performing (restricted) maximum likelihood estimation of generalized
linear models, Bayesian estimation is performed via MCMC. The Bayesian
model adds priors on the regression coefficients (in the same way as
[`stan_glm`](https://mc-stan.org/rstanarm/reference/stan_glm.md)) and
priors on the terms of a decomposition of the covariance matrices of the
group-specific parameters. See
[`priors`](https://mc-stan.org/rstanarm/reference/priors.md) for more
information about the priors.

The `stan_lmer` function is equivalent to `stan_glmer` with
`family = gaussian(link = "identity")`.

The `stan_glmer.nb` function, which takes the extra argument `link`, is
a wrapper for `stan_glmer` with
`family = `[`neg_binomial_2`](https://mc-stan.org/rstanarm/reference/neg_binomial_2.md)`(link)`.

## References

Gelman, A. and Hill, J. (2007). *Data Analysis Using Regression and
Multilevel/Hierarchical Models.* Cambridge University Press, Cambridge,
UK. (Ch. 11-15)

Muth, C., Oravecz, Z., and Gabry, J. (2018) User-friendly Bayesian
regression modeling: A tutorial with rstanarm and shinystan. *The
Quantitative Methods for Psychology*. 14(2), 99–119.
<https://www.tqmp.org/RegularArticles/vol14-2/p099/p099.pdf>

## See also

[`stanreg-methods`](https://mc-stan.org/rstanarm/reference/stanreg-methods.md)
and [`glmer`](https://rdrr.io/pkg/lme4/man/glmer.html).

The vignette for `stan_glmer` and the *Hierarchical Partial Pooling*
vignette. <https://mc-stan.org/rstanarm/articles/>

## Examples

``` r
if (.Platform$OS.type != "windows" || .Platform$r_arch != "i386") {
# see help(example_model) for details on the model below
if (!exists("example_model")) example(example_model) 
print(example_model, digits = 1)
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
#>  herd   (Intercept) 0.76    
#> Num. levels: herd 15 
#> 
#> ------
#> * For help interpreting the printed output see ?print.stanreg
#> * For info on the priors used see ?prior_summary.stanreg
```
