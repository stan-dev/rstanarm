# rstanarm 2.19.3

### Bug fixes

* Allow the vignettes to knit on platforms that do not support version 2 of RMarkdown

# rstanarm 2.19.2

### Bug fixes

* src/Makevars{.win} now uses a more robust way to find StanHeaders

* Fixed bug where `ranef()` and `coef()` methods for `glmer`-style models 
printed the wrong output for certain combinations of varying intercepts
and slopes.

* Fixed a bug where `posterior_predict()` failed for `stan_glmer()` models 
estimated with `family = mgcv::betar`.

* Fixed bug in `bayes_R2()` for bernoulli models. (Thanks to @mcol)

* `loo_R2()` can now be called on the same fitted model object multiple times
with identical (not just up to rng noise) results. (Thanks to @mcol)

### New features and improvements

* New vignette on doing MRP using rstanarm. (Thanks to @lauken13)

* 4x speedup for most GLMs (`stan_glm()`) and GAMs (`stan_gamm4()` without
`random` argument). This comes from using Stan's new compound `_glm` functions
(`normal_id_glm`, `bernoulli_logit_glm`, `poisson_log_glm`,
`neg_binomial_2_log_glm`) under the hood whenever possible. (Thanks 
to @avehtari and @VMatthijs)

* `compare_models()` is deprecated in favor of `loo_compare()` to keep up 
with the loo package ([loo::loo_compare()](http://mc-stan.org/loo/reference/loo_compare))

* The `kfold()` method now has a `cores` argument and parallelizes by fold
rather than by Markov chain (unless otherwise specified), which should be much
more efficient when many cores are available.

* For `stan_glm()` with `algorithm='optimizing'`, Pareto smoothed importance
sampling ([arxiv.org/abs/1507.02646](https://arxiv.org/abs/1507.02646),
[mc-stan.org/loo/reference/psis.html](https://mc-stan.org/loo/reference/psis.html))
is now used to diagnose and improve inference (see
https://avehtari.github.io/RAOS-Examples/BigData/bigdata.html). This also now
means that we can use PSIS-LOO also when `algorithm='optimizing'`. (Thanks 
to @avehtari)

* For `stan_glm()` the `"meanfield"` and `"fullrank"` ADVI algorithms also
include the PSIS diagnostics and adjustments, but so far we have not seen any
example where these would be better than optimzation or MCMC.


# rstanarm 2.18.1

### Bug fixes

* `stan_clogit()` now works even when there are no common predictors
* `prior.info()` works better with models produced by `stan_jm()` and
  `stan_mvmer()`

### New features and improvements

* `stan_glm()` (only) gets a `mean_PPD` argument that when `FALSE`
  avoids drawing from the posterior predictive distribution in the
  Stan code
* `posterior_linpred()` now works even if the model was estimated with
  `algorithm = "optimizing"`

# rstanarm 2.17.4

### Bug fixes

* `stan_jm()` and `stan_mvmer()` now correctly include the intercept in the
longitudinal submodel

### New features and improvements

* Compatible with **loo** package version `>= 2.0`

* `QR = TRUE` no longer ignores the `autoscale` argument and has better behavior when `autoscale = FALSE`

* `posterior_linpred()` now has a draws argument like for `posterior_predict()`

* Dynamic predictions are now supported in `posterior_traj()` for 
`stan_jm` models.

* More options for K-fold CV, including manually specifying the folds or using helper functions to create them for particular model/data combinations.


# rstanarm 2.17.3 

Minor release for build fixes for Solaris and avoiding a test failure

# rstanarm 2.17.2

Lots of good stuff in this release.

### Bug fixes

* `stan_polr()` and `stan_lm()` handle the `K = 1` case better

### Important user-facing improvements

* The prior_aux arguments now defaults to exponential rather than Cauchy. This should be a safer default.

* The Stan programs do not drop any constants and should now be
safe to use with the **bridgesampling** package

* `hs()` and `hs_plus()` priors have new defaults based on a new
paper by Aki Vehtari and Juho Piironen

* `stan_gamm4()` is now more closely based on `mgcv::jagam()`, which may affect
some estimates but the options remain largely the same

* The `product_normal()` prior permits `df = 1`, which is a product of ... one
normal variate

* The build system is more conventional now. It should require less RAM to build
from source but it is slower unless you utilize parallel make and LTO
    
### Big new features

* `stan_jm()` and `stan_mvmer()` contributed by Sam Brilleman

* `bayes_R2()` method to calculate a quantity similar to $R^2$

* `stan_nlmer()`, which is similar to `lme4::nlmer` but watch
out for multimodal posterior distributions

* `stan_clogit()`, which is similar to `survival::clogit` but
accepts lme4-style group-specific terms

* The `mgcv::betar` family is supported for the lme4-like modeling functions,
allowing for beta regressions with lme4-style group terms and / or smooth
nonlinear functions of predictors

# rstanarm 2.15.3

### Bug fixes

* Fix to `stan_glmer()` Bernoulli models with multiple group-specific intercept terms that could result in draws from the wrong posterior distribution

* Fix bug with contrasts in `stan_aov()` (thanks to Henrik Singmann)

* Fix bug with `na.action` in `stan_glmer()` (thanks to Henrik Singmann)

# rstanarm 2.15.1

Minor release with only changes to allow tests to pass on CRAN

# rstanrm 2.14.2

### Bug fixes

* Fix for intercept with identity or square root link functions for the
auxiliary parameter of a beta regression

* Fix for special case where only the intercepts vary by group and a non-default
prior is specified for their standard deviation

* Fix for off-by-one error in some lme4-style models with multiple grouping
terms
  
### New features

* New methods `loo_linpred()`, `loo_pit()`, `loo_predict()`, and `loo_predictive_interval()`

* Support for many more plotfuns in `pp_check()` that are implemented in the **bayesplot** package

* Option to compute latent residuals in `stan_polr()` (Thanks to Nate Sanders)

* The pairs plot now uses the ggplot2 package

# rstanarm 2.14.1

### Bug fixes

  * `VarCorr()` could return duplicates in cases where a `stan_{g}lmer` model used grouping factor level names with spaces
  
  * `The pairs()` function now works with group-specific parameters
  
  * The `stan_gamm4()` function works better now
  
  * Fix a problem with factor levels after estimating a model via `stan_lm()`

### New features

* New model-fitting function(s) `stan_betareg()` (and `stan_betareg.fit()`)
that uses the same likelihoods as those supported by the `betareg()` function in
the **betareg** package (Thanks to Imad Ali)

* New choices for priors on coefficients: `laplace()`, `lasso()`,
`product_normal()`

* The `hs()` and `hs_plus()` priors now have new `global_df` and `global_scale` arguments

* `stan_{g}lmer()` models that only have group-specific intercept shifts are considerably faster now

* Models with Student t priors and low degrees of freedom (that are not 1, 2, or 4) may work better now due to Cornish-Fisher transformations

* Many functions for priors have gained an `autoscale` argument that defaults to
`TRUE` and indicates that rstanarm should make internal changes to the prior
based on the scales of the variables so that they default priors are weakly
informative

* The new `compare_models()` function does more extensive checking that the
models being compared are compatible

### Deprecated arguments

* The `prior_ops` argument to various model fitting functions is deprecated
and replaced by a the `prior_aux` argument for the prior on the auxiliary
parameter of various GLM-like models

# rstanarm 2.13.1

### Bug fixes

  * Fix bug in `reloo()` if data was not specified
  * Fix bug in `pp_validate()` that was only introduced on GitHub
  
### New features

* Uses the new **bayesplot** and **rstantools** R packages

* The new `prior_summary()` function can be used to figure out what priors were actually used

* `stan_gamm4()` is better implemented, can be followed by `plot_nonlinear()`, 
`posterior_predict()` (with newdata), etc.

* Hyperparameters (i.e. covariance matrices in general) for lme4 style models
are now returned by `as.matrix()` and `as.data.frame()`

* `pp_validate()` can now be used if optimization or variational Bayesian
inference was used to estimate the original model


# rstanarm 2.12.1

### Bug fixes

* Fix for bad bug in `posterior_predict()` when factor labels have spaces in lme4-style models

* Fix when weights are used in Poisson models

### New features

* `posterior_linpred()` gains an `XZ` argument to output the design matrix

# rstanarm 2.11.1

### Bug fixes

* Requiring manually specifying offsets when model has an offset and newdata is not NULL
  
### New features

* `stan_biglm()` function that somewhat supports `biglm::biglm`

* `as.array()` method for stanreg objects

# rstanarm 2.10.1

### Bug fixes

  * Works with devtools now

### New features

* `k_threshold` argument to `loo()` to do PSIS-LOO+ 

* `kfold()` for K-fold CV

* Ability to use sparse X matrices (slowly) for many models if memory is an issue

### rstanarm 2.9.0-4

### Bug fixes 

* `posterior_predict()` with newdata now works correctly for ordinal models

* `stan_lm()` now works when intercept is omitted

* `stan_glmer.fit()` no longer permit models with duplicative group-specific terms since they don't make sense and are usually a mistake on the user's part

* `posterior_predict()` with lme4-style models no longer fails if there are
spaces or colons in the levels of the grouping variables

* `posterior_predict()` with ordinal models outputs a character matrix now

### New features

* `pp_validate()` function based on the BayesValidate package by Sam Cook

* `posterior_vs_prior()` function to visualize the effect of conditioning on the data

* Works (again) with R versions back to 3.0.2 (untested though)

# rstanarm 2.9.0-3

### Bug fixes

* Fix problem with models that had group-specific coefficients, which were
mislabled. Although the parameters were estimated correctly, users of previous
versions of rstanarm should run such models again to obtain correct summaries
and posterior predictions. Thanks to someone named Luke for pointing this
problem out on stan-users.

* Vignettes now view correctly on the CRAN webiste thanks to Yihui Xie

* Fix problem with models without intercepts thanks to Paul-Christian Buerkner

* Fix problem with specifying binomial 'size' for posterior_predict using newdata

* Fix problem with lme4-style formulas that use the same grouping factor multiple times

* Fix conclusion in rstanarm vignette thanks to someone named Michael

### New features

* Group-specific design matrices are kept sparse throughout to reduce memory consumption

* The `log_lik()` function now has a `newdata` argument

* New vignette on hierarchical partial pooling
  

# rstanarm 2.9.0-1

Initial CRAN release
