[<img src="https://raw.githubusercontent.com/stan-dev/logos/master/logo_tm.png" width=100 alt="Stan Logo"/>](http://mc-stan.org)

# rstanarm 

[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/rstanarm?color=blue)](http://cran.r-project.org/package=rstanarm)[![Downloads](http://cranlogs.r-pkg.org/badges/rstanarm?color=blue)](http://cran.rstudio.com/package=rstanarm)

### Bayesian applied regression modeling (arm) via Stan

This is an R package that emulates other R model-fitting functions but uses
[Stan](http://mc-stan.org) (via the **rstan** package) for the back-end
estimation. The primary target audience is people who would be open to Bayesian
inference if using Bayesian software were easier but would use frequentist
software otherwise. 

Fitting models with **rstanarm** is also useful for experienced Bayesian
software users who want to take advantage the pre-compiled Stan programs that
are written by Stan developers and carefully implemented to prioritize numerical
stability and the avoidance of sampling problems.

Click the arrows for more details:
<details><summary>More detail</summary>

The **rstanarm** package is an appendage to the **rstan** package, the R
interface to [Stan](http://mc-stan.org/). **rstanarm** enables many of the most
common applied regression models to be estimated using Markov Chain Monte Carlo,
variational approximations to the posterior distribution, or optimization. The
package allows these models to be specified using the customary R modeling
syntax (e.g., like that of `glm` with a `formula` and `data.frame`).
Additional arguments are provided for specifying prior distributions.

The set of models supported by **rstanarm** is large (and will continue to
grow), but also limited enough so that it is possible to integrate them
tightly with the [`pp_check`](http://mc-stan.org/rstanarm/reference/pp_check.stanreg.html) function for graphical posterior predictive checks using [**bayesplot**](http://mc-stan.org/bayesplot) and the
[`posterior_predict`](http://mc-stan.org/rstanarm/reference/posterior_predict.stanreg.html)
function to easily estimate the effect of specific manipulations of predictor
variables or to predict the outcome in a training set.

The fitted model objects returned by the **rstanarm** modeling functions are
called _stanreg_ objects. In addition to all of the traditional
[methods](http://mc-stan.org/rstanarm/reference/stanreg-methods.html)
defined for fitted model objects, stanreg objects can also be used with the
[**loo**](http://mc-stan.org/rstanarm/reference/loo.stanreg.html) package for
leave-one-out cross-validation, model comparison, and model weighting/averaging
and the [**shinystan**](http://mc-stan.org/rstanarm/reference/shinystan.html) 
package for exploring the posterior distribution and model diagnostics
with a graphical user interface. 

Check out the **rstanarm** [vignettes](http://mc-stan.org/rstanarm/articles/)
for examples and more details about the entire process.
</details>

<details><summary>Modeling functions</summary>

The model estimating functions are described in greater detail in their
individual help pages and vignettes. Here we provide a very brief overview:

* [__`stan_lm`__, __`stan_aov`__,__`stan_biglm`__](http://mc-stan.org/rstanarm/reference/stan_lm.html)

  Similar to  `lm` and `aov` but with novel regularizing priors on the model
  parameters that are driven by prior beliefs about R-squared, the proportion of
  variance in the outcome attributable to the predictors in a linear model.

* [__`stan_glm`__, __`stan_glm.nb`__](http://mc-stan.org/rstanarm/reference/stan_glm.html)

  Similar to `glm` but with various possible prior distributions for the
  coefficients and, if applicable, a prior distribution for any auxiliary
  parameter in a Generalized Linear Model (GLM) that is characterized by a
  `family` object (e.g. the shape parameter in Gamma models). It is also possible
  to estimate a negative binomial model similar to the `glm.nb` function
  in the `MASS` package.

* [__`stan_glmer`__, __`stan_glmer.nb`__, __`stan_lmer`__](http://mc-stan.org/rstanarm/reference/stan_glmer.html)

  Similar to the `glmer`, `glmer.nb`, and `lmer` functions (__lme4__ package) in
  that GLMs are augmented to have group-specific terms that deviate from the
  common coefficients according to a mean-zero multivariate normal distribution
  with a highly-structured but unknown covariance matrix (for which **rstanarm**
  introduces an innovative prior distribution). MCMC provides more appropriate
  estimates of uncertainty for models that consist of a mix of common and
  group-specific parameters.
  
* [__`stan_nlmer`__](http://mc-stan.org/rstanarm/reference/stan_nlmer.html)

  Similar to `nlmer` (__lme4__ package) package for nonlinear "mixed-effects"
  models, but flexible priors can be specified for all parameters in the model, 
  including the unknown covariance matrices for the varying 
  (group-specific) coefficients.

* [__`stan_gamm4`__](http://mc-stan.org/rstanarm/reference/stan_gamm4.html)

  Similar to `gamm4` (__gamm4__ package), which augments a GLM (possibly with
  group-specific terms) with nonlinear smooth functions of the predictors to
  form a Generalized Additive Mixed Model (GAMM). Rather than calling
  `lme4::glmer` like `gamm4` does, `stan_gamm4` essentially calls `stan_glmer`,
  which avoids the optimization issues that often crop up with GAMMs and
  provides better estimates for the uncertainty of the parameter estimates.
 
* [__`stan_polr`__](http://mc-stan.org/rstanarm/reference/stan_polr.html)

  Similar to `polr` (__MASS__ package) in that it models an ordinal response,
  but the Bayesian model also implies a prior distribution on the unknown
  cutpoints. Can also be used to model binary outcomes, possibly while
  estimating an unknown exponent governing the probability of success.
 
* [__`stan_betareg`__](http://mc-stan.org/rstanarm/reference/stan_betareg.html)

  Similar to `betareg` (__betareg__ package) in that it models an outcome that
  is a rate (proportion) but, rather than performing maximum likelihood
  estimation, full Bayesian estimation is performed by default, with
  customizable prior distributions for all parameters.

* [__`stan_clogit`__](http://mc-stan.org/rstanarm/reference/stan_clogit.html)

   Similar to `clogit` (__survival__ package) in that it models an binary outcome
   where the number of successes and failures is fixed within each stratum by
   the research design. There are some minor syntactical differences relative
   to `survival::clogit` that allow `stan_clogit` to accept
   group-specific terms as in `stan_glmer`.

* [__`stan_mvmer`__](http://mc-stan.org/rstanarm/reference/stan_mvmer.html)

   A multivariate form of `stan_glmer`, whereby the user can specify
   one or more submodels each consisting of a GLM with group-specific terms. If
   more than one submodel is specified (i.e. there is more than one outcome
   variable) then a dependence is induced by assuming that the group-specific
   terms for each grouping factor are correlated across submodels.

* [__`stan_jm`__](http://mc-stan.org/rstanarm/reference/stan_jm.html)

   Estimates shared parameter joint models for longitudinal and time-to-event
   (i.e. survival) data. The joint model can be univariate (i.e. one longitudinal
   outcome) or multivariate (i.e. more than one longitudinal outcome). A variety
   of parameterisations are available for linking the longitudinal and event
   processes (i.e. a variety of association structures).

</details>

<details><summary>Estimation algorithms</summary>

The modeling functions in the **rstanarm** package take an `algorithm`
argument that can be one of the following:

* __Sampling__ (`algorithm="sampling"`):
 
 Uses Markov Chain Monte Carlo (MCMC) --- in particular, Stan's implementation
 of Hamiltonian Monte Carlo (HMC) with a tuned but diagonal mass matrix --- 
 to draw from the posterior distribution of the parameters. This is the slowest
 but most reliable of the available estimation algorithms and it is __the
 default and recommended algorithm for statistical inference__.

* __Mean-field__ (`algorithm="meanfield"`):

 Uses mean-field variational inference to draw from an approximation to the
 posterior distribution. In particular, this algorithm finds the set of
 independent normal distributions in the unconstrained space that --- when
 transformed into the constrained space --- most closely approximate the
 posterior distribution. Then it draws repeatedly from these independent
 normal distributions and transforms them into the constrained space. The
 entire process is much faster than HMC and yields independent draws but
 __is not recommended for final statistical inference__. It can be useful to
 narrow the set of candidate models in large problems, particularly when
 specifying `QR=TRUE` in `stan_glm`, `stan_glmer`, and `stan_gamm4`, but is
 __only an approximation to the posterior distribution__.

* __Full-rank__ (`algorithm="fullrank"`):

 Uses full-rank variational inference to draw from an approximation to the
 posterior distribution by finding the multivariate normal distribution in
 the unconstrained space that --- when transformed into the constrained space
 --- most closely approximates the posterior distribution. Then it draws
 repeatedly from this multivariate normal distribution and transforms the
 draws into the constrained space. This process is slower than meanfield
 variational inference but is faster than HMC. Although still an
 approximation to the posterior distribution and thus __not recommended
 for final statistical inference__, the approximation is more realistic than
 that of mean-field variational inference because the parameters are not
 assumed to be independent in the unconstrained space. Nevertheless, fullrank
 variational inference is a more difficult optimization problem and the
 algorithm is more prone to non-convergence or convergence to a local
 optimum.

* __Optimizing__ (`algorithm="optimizing"`):

 Finds the posterior mode using a C++ implementation of the LBGFS algorithm. If
 there is no prior information, then this is equivalent to maximum likelihood,
 in which case there is no great reason to use the functions in the **rstanarm**
 package over the emulated functions in other packages. However, if priors are
 specified, then the estimates are penalized maximum likelihood estimates, which
 may have some redeeming value. Currently, optimization is only supported for
 `stan_glm`.

</details>

---

### Resources

* [mc-stan.org/rstanarm](http://mc-stan.org/rstanarm) (online documentation, vignettes)
* [Ask a question](http://discourse.mc-stan.org) (Stan Forums on Discourse)
* [Open an issue](https://github.com/stan-dev/rstanarm/issues) (GitHub issues for bug reports, feature requests)

### Installation

#### Latest Release

The most recent **rstanarm** release can be installed from CRAN via

```r
install.packages("rstanarm")
```

#### Development Version

To install from GitHub, first make sure that you can install the **rstan**
package and C++ toolchain by following these
[instructions](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started).
Once **rstan** is successfully installed, you can install **rstanarm** from
GitHub using the **devtools** package by executing the following in R:

```r
if (!require(devtools)) {
  install.packages("devtools")
  library(devtools)
}
install_github("stan-dev/rstanarm", build_vignettes = FALSE)
```

You can switch `build_vignettes` to `TRUE` but it takes a lot longer to install and the 
vignettes are already separately available from the 
[Stan website](http://mc-stan.org/rstanarm/articles/index.html) 
and 
[CRAN](https://cran.r-project.org/package=rstanarm/vignettes). 
If installation fails, please let us know by [filing an issue](https://github.com/stan-dev/rstanarm/issues).

### Contributing 

If you are interested in contributing to the development of **rstanarm** please 
see the [developer notes](http://mc-stan.org/rstanarm/dev-notes/index.html) page.
