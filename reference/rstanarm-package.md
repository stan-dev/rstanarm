# Applied Regression Modeling via RStan

*Stan Development Team*

The rstanarm package is an appendage to the rstan package that enables
many of the most common applied regression models to be estimated using
Markov Chain Monte Carlo, variational approximations to the posterior
distribution, or optimization. The rstanarm package allows these models
to be specified using the customary R modeling syntax (e.g., like that
of [`glm`](https://rdrr.io/r/stats/glm.html) with a `formula` and a
`data.frame`).

The sections below provide an overview of the modeling functions and
estimation algorithms used by rstanarm.

## Details

The set of models supported by rstanarm is large (and will continue to
grow), but also limited enough so that it is possible to integrate them
tightly with the
[`pp_check`](https://mc-stan.org/rstanarm/reference/pp_check.stanreg.md)
function for graphical posterior predictive checks with
[bayesplot](https://mc-stan.org/bayesplot/reference/bayesplot-package.html)
and the
[`posterior_predict`](https://mc-stan.org/rstanarm/reference/posterior_predict.stanreg.md)
function to easily estimate the effect of specific manipulations of
predictor variables or to predict the outcome in a training set.

The objects returned by the rstanarm modeling functions are called
[`stanreg`](https://mc-stan.org/rstanarm/reference/stanreg-objects.md)
objects. In addition to all of the typical
[`methods`](https://mc-stan.org/rstanarm/reference/stanreg-methods.md)
defined for fitted model objects, stanreg objects can be passed to the
[`loo`](https://mc-stan.org/loo/reference/loo.html) function in the loo
package for model comparison or to the
[`launch_shinystan`](https://mc-stan.org/shinystan/reference/launch_shinystan.html)
function in the shinystan package in order to visualize the posterior
distribution using the ShinyStan graphical user interface. See the
rstanarm vignettes for more details about the entire process.

## Prior distributions

See [priors help page](https://mc-stan.org/rstanarm/reference/priors.md)
and the vignette [*Prior Distributions for rstanarm
Models*](https://mc-stan.org/rstanarm/articles/priors.html) for an
overview of the various choices the user can make for prior
distributions. The package vignettes for the modeling functions also
provide examples of using many of the available priors as well as more
detailed descriptions of some of the novel priors used by rstanarm.

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

## Estimation algorithms

The modeling functions in the rstanarm package take an `algorithm`
argument that can be one of the following:

- **Sampling** (`algorithm="sampling"`):

  Uses Markov Chain Monte Carlo (MCMC) — in particular, Hamiltonian
  Monte Carlo (HMC) with a tuned but diagonal mass matrix — to draw from
  the posterior distribution of the parameters. See `sampling` (rstan)
  for more details. This is the slowest but most reliable of the
  available estimation algorithms and it is **the default and
  recommended algorithm for statistical inference.**

- **Mean-field** (`algorithm="meanfield"`):

  Uses mean-field variational inference to draw from an approximation to
  the posterior distribution. In particular, this algorithm finds the
  set of independent normal distributions in the unconstrained space
  that — when transformed into the constrained space — most closely
  approximate the posterior distribution. Then it draws repeatedly from
  these independent normal distributions and transforms them into the
  constrained space. The entire process is much faster than HMC and
  yields independent draws but **is not recommended for final
  statistical inference**. It can be useful to narrow the set of
  candidate models in large problems, particularly when specifying
  `QR=TRUE` in
  [`stan_glm`](https://mc-stan.org/rstanarm/reference/stan_glm.md),
  [`stan_glmer`](https://mc-stan.org/rstanarm/reference/stan_glmer.md),
  and
  [`stan_gamm4`](https://mc-stan.org/rstanarm/reference/stan_gamm4.md),
  but is **only an approximation to the posterior distribution**.

- **Full-rank** (`algorithm="fullrank"`):

  Uses full-rank variational inference to draw from an approximation to
  the posterior distribution by finding the multivariate normal
  distribution in the unconstrained space that — when transformed into
  the constrained space — most closely approximates the posterior
  distribution. Then it draws repeatedly from this multivariate normal
  distribution and transforms the draws into the constrained space. This
  process is slower than meanfield variational inference but is faster
  than HMC. Although still an approximation to the posterior
  distribution and thus **not recommended for final statistical
  inference**, the approximation is more realistic than that of
  mean-field variational inference because the parameters are not
  assumed to be independent in the unconstrained space. Nevertheless,
  fullrank variational inference is a more difficult optimization
  problem and the algorithm is more prone to non-convergence or
  convergence to a local optimum.

- **Optimizing** (`algorithm="optimizing"`):

  Finds the posterior mode using a C++ implementation of the LBGFS
  algorithm. See `optimizing` for more details. If there is no prior
  information, then this is equivalent to maximum likelihood, in which
  case there is no great reason to use the functions in the rstanarm
  package over the emulated functions in other packages. However, if
  priors are specified, then the estimates are penalized maximum
  likelihood estimates, which may have some redeeming value. Currently,
  optimization is only supported for
  [`stan_glm`](https://mc-stan.org/rstanarm/reference/stan_glm.md).

## References

Bates, D., Maechler, M., Bolker, B., and Walker, S. (2015). Fitting
linear mixed-Effects models using lme4. *Journal of Statistical
Software*. 67(1), 1–48.

Gelman, A., Carlin, J. B., Stern, H. S., Dunson, D. B., Vehtari, A., and
Rubin, D. B. (2013). *Bayesian Data Analysis.* Chapman & Hall/CRC Press,
London, third edition. <https://sites.stat.columbia.edu/gelman/book/>

Gelman, A. and Hill, J. (2007). *Data Analysis Using Regression and
Multilevel/Hierarchical Models.* Cambridge University Press, Cambridge,
UK. <https://sites.stat.columbia.edu/gelman/arm/>

Stan Development Team. *Stan Modeling Language Users Guide and Reference
Manual.* <https://mc-stan.org/users/documentation/>.

Vehtari, A., Gelman, A., and Gabry, J. (2017). Practical Bayesian model
evaluation using leave-one-out cross-validation and WAIC. *Statistics
and Computing*. 27(5), 1413–1432. doi:10.1007/s11222-016-9696-4. arXiv
preprint: <https://arxiv.org/abs/1507.04544>

Yao, Y., Vehtari, A., Simpson, D., and Gelman, A. (2018) Using stacking
to average Bayesian predictive distributions. *Bayesian Analysis*,
advance publication,
[doi:10.1214/17-BA1091](https://doi.org/10.1214/17-BA1091) .

Gabry, J. , Simpson, D. , Vehtari, A. , Betancourt, M. and Gelman, A.
(2019), Visualization in Bayesian workflow. *J. R. Stat. Soc. A*, 182:
389-402. doi:10.1111/rssa.12378, [arXiv
preprint](https://arxiv.org/abs/1709.01449), [code on
GitHub](https://github.com/jgabry/bayes-vis-paper))

Muth, C., Oravecz, Z., and Gabry, J. (2018) User-friendly Bayesian
regression modeling: A tutorial with rstanarm and shinystan. *The
Quantitative Methods for Psychology*. 14(2), 99–119.
<https://www.tqmp.org/RegularArticles/vol14-2/p099/p099.pdf>

## See also

- <https://mc-stan.org/> for more information on the Stan C++ package
  used by rstanarm for model fitting.

- <https://github.com/stan-dev/rstanarm/issues/> to submit a bug report
  or feature request.

- <https://discourse.mc-stan.org> to ask a question about rstanarm on
  the Stan-users forum.

## Author

**Maintainer**: Ben Goodrich <benjamin.goodrich@columbia.edu>

Authors:

- Jonah Gabry <jgabry@gmail.com>

Other contributors:

- Imad Ali \[contributor\]

- Sam Brilleman \[contributor\]

- Jacqueline Buros Novik (R/stan_jm.R) \[contributor\]

- AstraZeneca (R/stan_jm.R) \[contributor\]

- Trustees of Columbia University \[copyright holder\]

- Simon Wood (R/stan_gamm4.R) \[copyright holder\]

- R Core Deveopment Team (R/stan_aov.R) \[copyright holder\]

- Douglas Bates (R/pp_data.R) \[copyright holder\]

- Martin Maechler (R/pp_data.R) \[copyright holder\]

- Ben Bolker (R/pp_data.R) \[copyright holder\]

- Steve Walker (R/pp_data.R) \[copyright holder\]

- Brian Ripley (R/stan_aov.R, R/stan_polr.R) \[copyright holder\]

- William Venables (R/stan_polr.R) \[copyright holder\]

- Paul-Christian Burkner <paul.buerkner@gmail.com> (R/misc.R)
  \[copyright holder\]
