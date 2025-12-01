# Package index

## About rstanarm

These pages provides a summary of the functionality available in
rstanarm.

- [`rstanarm`](https://mc-stan.org/rstanarm/reference/rstanarm-package.md)
  [`rstanarm-package`](https://mc-stan.org/rstanarm/reference/rstanarm-package.md)
  : Applied Regression Modeling via RStan

- [`available-models`](https://mc-stan.org/rstanarm/reference/available-models.md)
  :

  Modeling functions available in rstanarm

- [`available-algorithms`](https://mc-stan.org/rstanarm/reference/available-algorithms.md)
  :

  Estimation algorithms available for rstanarm models

## Fitting models

Functions for model fitting.

- [`stan_betareg()`](https://mc-stan.org/rstanarm/reference/stan_betareg.md)
  [`stan_betareg.fit()`](https://mc-stan.org/rstanarm/reference/stan_betareg.md)
  : Bayesian beta regression models via Stan
- [`stan_biglm()`](https://mc-stan.org/rstanarm/reference/stan_biglm.md)
  [`stan_biglm.fit()`](https://mc-stan.org/rstanarm/reference/stan_biglm.md)
  : Bayesian regularized linear but big models via Stan
- [`stan_clogit()`](https://mc-stan.org/rstanarm/reference/stan_clogit.md)
  : Conditional logistic (clogit) regression models via Stan
- [`stan_gamm4()`](https://mc-stan.org/rstanarm/reference/stan_gamm4.md)
  [`plot_nonlinear()`](https://mc-stan.org/rstanarm/reference/stan_gamm4.md)
  : Bayesian generalized linear additive models with optional
  group-specific terms via Stan
- [`stan_glm()`](https://mc-stan.org/rstanarm/reference/stan_glm.md)
  [`stan_glm.nb()`](https://mc-stan.org/rstanarm/reference/stan_glm.md)
  [`stan_glm.fit()`](https://mc-stan.org/rstanarm/reference/stan_glm.md)
  : Bayesian generalized linear models via Stan
- [`stan_glmer()`](https://mc-stan.org/rstanarm/reference/stan_glmer.md)
  [`stan_lmer()`](https://mc-stan.org/rstanarm/reference/stan_glmer.md)
  [`stan_glmer.nb()`](https://mc-stan.org/rstanarm/reference/stan_glmer.md)
  : Bayesian generalized linear models with group-specific terms via
  Stan
- [`stan_jm()`](https://mc-stan.org/rstanarm/reference/stan_jm.md) :
  Bayesian joint longitudinal and time-to-event models via Stan
- [`stan_aov()`](https://mc-stan.org/rstanarm/reference/stan_lm.md)
  [`stan_lm()`](https://mc-stan.org/rstanarm/reference/stan_lm.md)
  [`stan_lm.wfit()`](https://mc-stan.org/rstanarm/reference/stan_lm.md)
  [`stan_lm.fit()`](https://mc-stan.org/rstanarm/reference/stan_lm.md) :
  Bayesian regularized linear models via Stan
- [`stan_mvmer()`](https://mc-stan.org/rstanarm/reference/stan_mvmer.md)
  : Bayesian multivariate generalized linear models with correlated
  group-specific terms via Stan
- [`stan_nlmer()`](https://mc-stan.org/rstanarm/reference/stan_nlmer.md)
  : Bayesian nonlinear models with group-specific terms via Stan
- [`stan_polr()`](https://mc-stan.org/rstanarm/reference/stan_polr.md)
  [`stan_polr.fit()`](https://mc-stan.org/rstanarm/reference/stan_polr.md)
  : Bayesian ordinal regression models via Stan
- [`normal()`](https://mc-stan.org/rstanarm/reference/priors.md)
  [`student_t()`](https://mc-stan.org/rstanarm/reference/priors.md)
  [`cauchy()`](https://mc-stan.org/rstanarm/reference/priors.md)
  [`hs()`](https://mc-stan.org/rstanarm/reference/priors.md)
  [`hs_plus()`](https://mc-stan.org/rstanarm/reference/priors.md)
  [`laplace()`](https://mc-stan.org/rstanarm/reference/priors.md)
  [`lasso()`](https://mc-stan.org/rstanarm/reference/priors.md)
  [`product_normal()`](https://mc-stan.org/rstanarm/reference/priors.md)
  [`exponential()`](https://mc-stan.org/rstanarm/reference/priors.md)
  [`decov()`](https://mc-stan.org/rstanarm/reference/priors.md)
  [`lkj()`](https://mc-stan.org/rstanarm/reference/priors.md)
  [`dirichlet()`](https://mc-stan.org/rstanarm/reference/priors.md)
  [`R2()`](https://mc-stan.org/rstanarm/reference/priors.md)
  [`default_prior_intercept()`](https://mc-stan.org/rstanarm/reference/priors.md)
  [`default_prior_coef()`](https://mc-stan.org/rstanarm/reference/priors.md)
  : Prior distributions and options

## Methods

Functions to work with fitted model objects.

- [`stanreg-objects`](https://mc-stan.org/rstanarm/reference/stanreg-objects.md)
  : Fitted model objects

- [`as.matrix(`*`<stanreg>`*`)`](https://mc-stan.org/rstanarm/reference/as.matrix.stanreg.md)
  [`as.array(`*`<stanreg>`*`)`](https://mc-stan.org/rstanarm/reference/as.matrix.stanreg.md)
  [`as.data.frame(`*`<stanreg>`*`)`](https://mc-stan.org/rstanarm/reference/as.matrix.stanreg.md)
  : Extract the posterior sample

- [`bayes_R2(`*`<stanreg>`*`)`](https://mc-stan.org/rstanarm/reference/bayes_R2.stanreg.md)
  [`loo_R2(`*`<stanreg>`*`)`](https://mc-stan.org/rstanarm/reference/bayes_R2.stanreg.md)
  : Compute a Bayesian version of R-squared or LOO-adjusted R-squared
  for regression models.

- [`kfold(`*`<stanreg>`*`)`](https://mc-stan.org/rstanarm/reference/kfold.stanreg.md)
  : K-fold cross-validation

- [`launch_shinystan(`*`<stanreg>`*`)`](https://mc-stan.org/rstanarm/reference/launch_shinystan.stanreg.md)
  : Using the ShinyStan GUI with rstanarm models

- [`log_lik(`*`<stanreg>`*`)`](https://mc-stan.org/rstanarm/reference/log_lik.stanreg.md)
  [`log_lik(`*`<stanmvreg>`*`)`](https://mc-stan.org/rstanarm/reference/log_lik.stanreg.md)
  [`log_lik(`*`<stanjm>`*`)`](https://mc-stan.org/rstanarm/reference/log_lik.stanreg.md)
  : Pointwise log-likelihood matrix

- [`loo(`*`<stanreg>`*`)`](https://mc-stan.org/rstanarm/reference/loo.stanreg.md)
  [`waic(`*`<stanreg>`*`)`](https://mc-stan.org/rstanarm/reference/loo.stanreg.md)
  [`loo_compare(`*`<stanreg>`*`)`](https://mc-stan.org/rstanarm/reference/loo.stanreg.md)
  [`loo_compare(`*`<stanreg_list>`*`)`](https://mc-stan.org/rstanarm/reference/loo.stanreg.md)
  [`loo_model_weights(`*`<stanreg_list>`*`)`](https://mc-stan.org/rstanarm/reference/loo.stanreg.md)
  [`compare_models()`](https://mc-stan.org/rstanarm/reference/loo.stanreg.md)
  : Information criteria and cross-validation

- [`loo_predict(`*`<stanreg>`*`)`](https://mc-stan.org/rstanarm/reference/loo_predict.stanreg.md)
  [`loo_linpred(`*`<stanreg>`*`)`](https://mc-stan.org/rstanarm/reference/loo_predict.stanreg.md)
  [`loo_predictive_interval(`*`<stanreg>`*`)`](https://mc-stan.org/rstanarm/reference/loo_predict.stanreg.md)
  : Compute weighted expectations using LOO

- [`pairs(`*`<stanreg>`*`)`](https://mc-stan.org/rstanarm/reference/pairs.stanreg.md)
  : Pairs method for stanreg objects

- [`plot(`*`<stanreg>`*`)`](https://mc-stan.org/rstanarm/reference/plot.stanreg.md)
  : Plot method for stanreg objects

- [`posterior_interval(`*`<stanreg>`*`)`](https://mc-stan.org/rstanarm/reference/posterior_interval.stanreg.md)
  : Posterior uncertainty intervals

- [`posterior_linpred(`*`<stanreg>`*`)`](https://mc-stan.org/rstanarm/reference/posterior_linpred.stanreg.md)
  [`posterior_epred(`*`<stanreg>`*`)`](https://mc-stan.org/rstanarm/reference/posterior_linpred.stanreg.md)
  : Posterior distribution of the (possibly transformed) linear
  predictor

- [`posterior_predict(`*`<stanreg>`*`)`](https://mc-stan.org/rstanarm/reference/posterior_predict.stanreg.md)
  [`posterior_predict(`*`<stanmvreg>`*`)`](https://mc-stan.org/rstanarm/reference/posterior_predict.stanreg.md)
  : Draw from posterior predictive distribution

- [`posterior_vs_prior()`](https://mc-stan.org/rstanarm/reference/posterior_vs_prior.md)
  : Juxtapose prior and posterior

- [`pp_check(`*`<stanreg>`*`)`](https://mc-stan.org/rstanarm/reference/pp_check.stanreg.md)
  : Graphical posterior predictive checks

- [`predict(`*`<stanreg>`*`)`](https://mc-stan.org/rstanarm/reference/predict.stanreg.md)
  : Predict method for stanreg objects

- [`predictive_error(`*`<stanreg>`*`)`](https://mc-stan.org/rstanarm/reference/predictive_error.stanreg.md)
  [`predictive_error(`*`<matrix>`*`)`](https://mc-stan.org/rstanarm/reference/predictive_error.stanreg.md)
  [`predictive_error(`*`<ppd>`*`)`](https://mc-stan.org/rstanarm/reference/predictive_error.stanreg.md)
  : In-sample or out-of-sample predictive errors

- [`predictive_interval(`*`<stanreg>`*`)`](https://mc-stan.org/rstanarm/reference/predictive_interval.stanreg.md)
  [`predictive_interval(`*`<matrix>`*`)`](https://mc-stan.org/rstanarm/reference/predictive_interval.stanreg.md)
  [`predictive_interval(`*`<ppd>`*`)`](https://mc-stan.org/rstanarm/reference/predictive_interval.stanreg.md)
  : Predictive intervals

- [`print(`*`<stanreg>`*`)`](https://mc-stan.org/rstanarm/reference/print.stanreg.md)
  [`print(`*`<stanmvreg>`*`)`](https://mc-stan.org/rstanarm/reference/print.stanreg.md)
  : Print method for stanreg objects

- [`prior_summary(`*`<stanreg>`*`)`](https://mc-stan.org/rstanarm/reference/prior_summary.stanreg.md)
  : Summarize the priors used for an rstanarm model

- [`as_draws(`*`<stanreg>`*`)`](https://mc-stan.org/rstanarm/reference/stanreg-draws-formats.md)
  [`as_draws_matrix(`*`<stanreg>`*`)`](https://mc-stan.org/rstanarm/reference/stanreg-draws-formats.md)
  [`as_draws_array(`*`<stanreg>`*`)`](https://mc-stan.org/rstanarm/reference/stanreg-draws-formats.md)
  [`as_draws_df(`*`<stanreg>`*`)`](https://mc-stan.org/rstanarm/reference/stanreg-draws-formats.md)
  [`as_draws_list(`*`<stanreg>`*`)`](https://mc-stan.org/rstanarm/reference/stanreg-draws-formats.md)
  [`as_draws_rvars(`*`<stanreg>`*`)`](https://mc-stan.org/rstanarm/reference/stanreg-draws-formats.md)
  :

  Create a `draws` object from a `stanreg` object

- [`nobs(`*`<stanmvreg>`*`)`](https://mc-stan.org/rstanarm/reference/stanreg-methods.md)
  [`coef(`*`<stanreg>`*`)`](https://mc-stan.org/rstanarm/reference/stanreg-methods.md)
  [`confint(`*`<stanreg>`*`)`](https://mc-stan.org/rstanarm/reference/stanreg-methods.md)
  [`fitted(`*`<stanreg>`*`)`](https://mc-stan.org/rstanarm/reference/stanreg-methods.md)
  [`nobs(`*`<stanreg>`*`)`](https://mc-stan.org/rstanarm/reference/stanreg-methods.md)
  [`residuals(`*`<stanreg>`*`)`](https://mc-stan.org/rstanarm/reference/stanreg-methods.md)
  [`se(`*`<stanreg>`*`)`](https://mc-stan.org/rstanarm/reference/stanreg-methods.md)
  [`update(`*`<stanreg>`*`)`](https://mc-stan.org/rstanarm/reference/stanreg-methods.md)
  [`vcov(`*`<stanreg>`*`)`](https://mc-stan.org/rstanarm/reference/stanreg-methods.md)
  [`fixef(`*`<stanreg>`*`)`](https://mc-stan.org/rstanarm/reference/stanreg-methods.md)
  [`ngrps(`*`<stanreg>`*`)`](https://mc-stan.org/rstanarm/reference/stanreg-methods.md)
  [`nsamples(`*`<stanreg>`*`)`](https://mc-stan.org/rstanarm/reference/stanreg-methods.md)
  [`ranef(`*`<stanreg>`*`)`](https://mc-stan.org/rstanarm/reference/stanreg-methods.md)
  [`sigma(`*`<stanreg>`*`)`](https://mc-stan.org/rstanarm/reference/stanreg-methods.md)
  [`VarCorr(`*`<stanreg>`*`)`](https://mc-stan.org/rstanarm/reference/stanreg-methods.md)
  : Methods for stanreg objects

- [`summary(`*`<stanreg>`*`)`](https://mc-stan.org/rstanarm/reference/summary.stanreg.md)
  [`print(`*`<summary.stanreg>`*`)`](https://mc-stan.org/rstanarm/reference/summary.stanreg.md)
  [`as.data.frame(`*`<summary.stanreg>`*`)`](https://mc-stan.org/rstanarm/reference/summary.stanreg.md)
  [`summary(`*`<stanmvreg>`*`)`](https://mc-stan.org/rstanarm/reference/summary.stanreg.md)
  [`print(`*`<summary.stanmvreg>`*`)`](https://mc-stan.org/rstanarm/reference/summary.stanreg.md)
  : Summary method for stanreg objects

- [`posterior_survfit()`](https://mc-stan.org/rstanarm/reference/posterior_survfit.md)
  : Estimate subject-specific or standardised survival probabilities

- [`posterior_traj()`](https://mc-stan.org/rstanarm/reference/posterior_traj.md)
  : Estimate the subject-specific or marginal longitudinal trajectory

- [`ps_check()`](https://mc-stan.org/rstanarm/reference/ps_check.md) :
  Graphical checks of the estimated survival function

- [`coef(`*`<stanmvreg>`*`)`](https://mc-stan.org/rstanarm/reference/stanmvreg-methods.md)
  [`fitted(`*`<stanmvreg>`*`)`](https://mc-stan.org/rstanarm/reference/stanmvreg-methods.md)
  [`residuals(`*`<stanmvreg>`*`)`](https://mc-stan.org/rstanarm/reference/stanmvreg-methods.md)
  [`se(`*`<stanmvreg>`*`)`](https://mc-stan.org/rstanarm/reference/stanmvreg-methods.md)
  [`formula(`*`<stanmvreg>`*`)`](https://mc-stan.org/rstanarm/reference/stanmvreg-methods.md)
  [`update(`*`<stanmvreg>`*`)`](https://mc-stan.org/rstanarm/reference/stanmvreg-methods.md)
  [`update(`*`<stanjm>`*`)`](https://mc-stan.org/rstanarm/reference/stanmvreg-methods.md)
  [`fixef(`*`<stanmvreg>`*`)`](https://mc-stan.org/rstanarm/reference/stanmvreg-methods.md)
  [`ngrps(`*`<stanmvreg>`*`)`](https://mc-stan.org/rstanarm/reference/stanmvreg-methods.md)
  [`ranef(`*`<stanmvreg>`*`)`](https://mc-stan.org/rstanarm/reference/stanmvreg-methods.md)
  [`sigma(`*`<stanmvreg>`*`)`](https://mc-stan.org/rstanarm/reference/stanmvreg-methods.md)
  : Methods for stanmvreg objects

- [`plot(`*`<predict.stanjm>`*`)`](https://mc-stan.org/rstanarm/reference/plot.predict.stanjm.md)
  : Plot the estimated subject-specific or marginal longitudinal
  trajectory

- [`plot(`*`<survfit.stanjm>`*`)`](https://mc-stan.org/rstanarm/reference/plot.survfit.stanjm.md)
  [`plot_stack_jm()`](https://mc-stan.org/rstanarm/reference/plot.survfit.stanjm.md)
  : Plot the estimated subject-specific or marginal survival function

## Additional documentation

Misc. other help pages.

- [`rstanarm-datasets`](https://mc-stan.org/rstanarm/reference/rstanarm-datasets.md)
  [`kidiq`](https://mc-stan.org/rstanarm/reference/rstanarm-datasets.md)
  [`roaches`](https://mc-stan.org/rstanarm/reference/rstanarm-datasets.md)
  [`wells`](https://mc-stan.org/rstanarm/reference/rstanarm-datasets.md)
  [`bball1970`](https://mc-stan.org/rstanarm/reference/rstanarm-datasets.md)
  [`bball2006`](https://mc-stan.org/rstanarm/reference/rstanarm-datasets.md)
  [`mortality`](https://mc-stan.org/rstanarm/reference/rstanarm-datasets.md)
  [`tumors`](https://mc-stan.org/rstanarm/reference/rstanarm-datasets.md)
  [`radon`](https://mc-stan.org/rstanarm/reference/rstanarm-datasets.md)
  [`pbcLong`](https://mc-stan.org/rstanarm/reference/rstanarm-datasets.md)
  [`pbcSurv`](https://mc-stan.org/rstanarm/reference/rstanarm-datasets.md)
  : Datasets for rstanarm examples

- [`example_model`](https://mc-stan.org/rstanarm/reference/example_model.md)
  : Example model

- [`example_jm`](https://mc-stan.org/rstanarm/reference/example_jm.md) :
  Example joint longitudinal and time-to-event model

- [`stanreg_list()`](https://mc-stan.org/rstanarm/reference/stanreg_list.md)
  [`stanmvreg_list()`](https://mc-stan.org/rstanarm/reference/stanreg_list.md)
  [`stanjm_list()`](https://mc-stan.org/rstanarm/reference/stanreg_list.md)
  [`print(`*`<stanreg_list>`*`)`](https://mc-stan.org/rstanarm/reference/stanreg_list.md)
  : Create lists of fitted model objects, combine them, or append new
  models to existing lists of models.

- [`adapt_delta`](https://mc-stan.org/rstanarm/reference/adapt_delta.md)
  :

  `adapt_delta`: Target average acceptance probability

- [`QR-argument`](https://mc-stan.org/rstanarm/reference/QR-argument.md)
  :

  The `QR` argument

- [`neg_binomial_2()`](https://mc-stan.org/rstanarm/reference/neg_binomial_2.md)
  : Family function for negative binomial GLMs

- [`prior_options()`](https://mc-stan.org/rstanarm/reference/rstanarm-deprecated.md)
  : Deprecated functions

- [`logit()`](https://mc-stan.org/rstanarm/reference/logit.md)
  [`invlogit()`](https://mc-stan.org/rstanarm/reference/logit.md) :
  Logit and inverse logit
