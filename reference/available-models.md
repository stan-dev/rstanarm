# Modeling functions available in rstanarm

Modeling functions available in rstanarm

## Modeling functions

The model estimating functions are described in greater detail in their
individual help pages and vignettes. Here we provide a very brief
overview:

- [`stan_lm`](https://mc-stan.org/rstanarm/reference/stan_lm.md),
  `stan_aov`, `stan_biglm`:

  Similar to [`lm`](https://rdrr.io/r/stats/lm.html) or
  [`aov`](https://rdrr.io/r/stats/aov.html) but with novel regularizing
  priors on the model parameters that are driven by prior beliefs about
  \\R^2\\, the proportion of variance in the outcome attributable to the
  predictors in a linear model.

- [`stan_glm`](https://mc-stan.org/rstanarm/reference/stan_glm.md),
  `stan_glm.nb`:

  Similar to [`glm`](https://rdrr.io/r/stats/glm.html) but with various
  possible prior distributions for the coefficients and, if applicable,
  a prior distribution for any auxiliary parameter in a Generalized
  Linear Model (GLM) that is characterized by a
  [`family`](https://rdrr.io/r/stats/family.html) object (e.g. the shape
  parameter in Gamma models). It is also possible to estimate a negative
  binomial model in a similar way to the
  [`glm.nb`](https://rdrr.io/pkg/MASS/man/glm.nb.html) function in the
  MASS package.

- [`stan_glmer`](https://mc-stan.org/rstanarm/reference/stan_glmer.md),
  `stan_glmer.nb`, `stan_lmer`:

  Similar to the [`glmer`](https://rdrr.io/pkg/lme4/man/glmer.html),
  [`glmer.nb`](https://rdrr.io/pkg/lme4/man/glmer.nb.html) and
  [`lmer`](https://rdrr.io/pkg/lme4/man/lmer.html) functions in the lme4
  package in that GLMs are augmented to have group-specific terms that
  deviate from the common coefficients according to a mean-zero
  multivariate normal distribution with a highly-structured but unknown
  covariance matrix (for which rstanarm introduces an innovative prior
  distribution). MCMC provides more appropriate estimates of uncertainty
  for models that consist of a mix of common and group-specific
  parameters.

- [`stan_nlmer`](https://mc-stan.org/rstanarm/reference/stan_nlmer.md):

  Similar to [`nlmer`](https://rdrr.io/pkg/lme4/man/nlmer.html) in the
  lme4 package for nonlinear "mixed-effects" models, but the
  group-specific coefficients have flexible priors on their unknown
  covariance matrices.

- [`stan_gamm4`](https://mc-stan.org/rstanarm/reference/stan_gamm4.md):

  Similar to [`gamm4`](https://rdrr.io/pkg/gamm4/man/gamm4.html) in the
  gamm4 package, which augments a GLM (possibly with group-specific
  terms) with nonlinear smooth functions of the predictors to form a
  Generalized Additive Mixed Model (GAMM). Rather than calling
  [`glmer`](https://rdrr.io/pkg/lme4/man/glmer.html) like
  [`gamm4`](https://rdrr.io/pkg/gamm4/man/gamm4.html) does,
  [`stan_gamm4`](https://mc-stan.org/rstanarm/reference/stan_gamm4.md)
  essentially calls
  [`stan_glmer`](https://mc-stan.org/rstanarm/reference/stan_glmer.md),
  which avoids the optimization issues that often crop up with GAMMs and
  provides better estimates for the uncertainty of the parameter
  estimates.

- [`stan_polr`](https://mc-stan.org/rstanarm/reference/stan_polr.md):

  Similar to [`polr`](https://rdrr.io/pkg/MASS/man/polr.html) in the
  MASS package in that it models an ordinal response, but the Bayesian
  model also implies a prior distribution on the unknown cutpoints. Can
  also be used to model binary outcomes, possibly while estimating an
  unknown exponent governing the probability of success.

- [`stan_betareg`](https://mc-stan.org/rstanarm/reference/stan_betareg.md):

  Similar to [`betareg`](https://rdrr.io/pkg/betareg/man/betareg.html)
  in that it models an outcome that is a rate (proportion) but, rather
  than performing maximum likelihood estimation, full Bayesian
  estimation is performed by default, with customizable prior
  distributions for all parameters.

- [`stan_clogit`](https://mc-stan.org/rstanarm/reference/stan_clogit.md):

  Similar to [`clogit`](https://rdrr.io/pkg/survival/man/clogit.html) in
  that it models an binary outcome where the number of successes and
  failures is fixed within each stratum by the research design. There
  are some minor syntactical differences relative to
  [`clogit`](https://rdrr.io/pkg/survival/man/clogit.html) that allow
  `stan_clogit` to accept group-specific terms as in
  [`stan_glmer`](https://mc-stan.org/rstanarm/reference/stan_glmer.md).

- [`stan_mvmer`](https://mc-stan.org/rstanarm/reference/stan_mvmer.md):

  A multivariate form of
  [`stan_glmer`](https://mc-stan.org/rstanarm/reference/stan_glmer.md),
  whereby the user can specify one or more submodels each consisting of
  a GLM with group-specific terms. If more than one submodel is
  specified (i.e. there is more than one outcome variable) then a
  dependence is induced by assuming that the group-specific terms for
  each grouping factor are correlated across submodels.

- [`stan_jm`](https://mc-stan.org/rstanarm/reference/stan_jm.md):

  Estimates shared parameter joint models for longitudinal and
  time-to-event (i.e. survival) data. The joint model can be univariate
  (i.e. one longitudinal outcome) or multivariate (i.e. more than one
  longitudinal outcome). A variety of parameterisations are available
  for linking the longitudinal and event processes (i.e. a variety of
  association structures).

## See also

<https://mc-stan.org/rstanarm/>
