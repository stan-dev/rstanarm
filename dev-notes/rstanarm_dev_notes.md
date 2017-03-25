---
title: RStanARM Developer Notes
author: Imad Ali <br> March 10, 2017
numbersections: True
---

<!---
Run the following in terminal to render the md doc as a html file
with frost.css:
pandoc --css frost.css rstanarm_dev_notes.md --mathjax --toc --number-sections -o rstanarm_dev_notes.html
-->

# Preliminaries

This note is designed to help a developer contribute to rstanarm, which requires an understanding of how the various components of the package fit together.

Some things to keep in mind going forward:

* rstanarm objects are S3 objects of class `'stanreg'`.
* rstan objects are S4 objects of class `'stanfit'`.

Before working on including a model into rstanarm you should get the general form of the model(s) working using rstan. Include generated quantities such as posterior predictions, `log_lik`, and the mean posterior predictive distribution (`mean_PPD`). This will be useful in debugging the same model that you'll implement in rstanarm.

# Model Fitting

Below we break things down into stuff you need to do in R and stuff you need to do in Stan to get the model working.

## R Stuff

First add the package you plan on emulating in the `DESCRIPTION` file under `suggests:`.

You'll need to write an rstanarm function for each of the models you plan on emulating. These functions should call a .fit workhorse function which will run Stan. Essentially, the workflow will look something like the following:

```
                     ┏━━━━━━ USER ━━━━━━┓
                  model_A.R                  model_B.R
                      ┗━━━━ model.fit ━━━━┛
                                   ┃
                               model.stan
```

Regarding the above diagram, the `model_*.R` and `model.fit` files are in the `R` folder and the `model.stan` file is in the `exec` folder. You should include existing snippets of stan code in your `model.stan` file. These are located in the `inst/chunks` folder.

In the .fit function you should,

* Center the covariates.
* Perform QR decomposition.
* Determine whether the intercept is declared (or whether multiple linear predictors have been declared).
* Extract and process information regarding the priors declared.
* Define the parameters that you want Stan to return as a vector of strings `pars`.
* Put together data as a list that can be called into rstan.

Conditionals should be set up so that the .fit workhorse function can deal with the following algorithm arguments.

* `'optimizing'`
* `'sampling'`
* `'fullrank'` and `'meanfield'` (Variational Bayes methods)

### Tips

* If you added a new argument to your modeling function and are wondering why the argument is not recognized when you run the function (after building the package) it's probably because you haven't NULLified the argument in the `match.call(expand.dots = FALSE)` list.
* If Stan isn't returning parameters that are in the model make sure you've specified them in `pars`.
* If the model's family is using an existing family but the model linear predictor, etc. is different then it might be a good idea to give it an additional class identifier: `class(out) <- c("stanreg", "additional_class")`

## Stan Stuff

In the Stan file you should try to minimize the number of loops, storing n-dimensional objects, redoing calculations that only need to be done once, calculating inverses (e.g. use the precision multinormal instead of multinormal).

Because we center the covariates you have to separate the intercept out of the linear predictor matrix (i.e. `X` should not contain a vector of ones). If `gamma` is the intercept parameter fit using centered predictors and `alpha` is the intercept parameter you want to report then do the following transformation in generated quantities:
```
...
generated quantities {
  real alpha[has_intercept];
  {
    alpha[1] = gamma[1] - dot_product(beta, xbar)
  }
}
```

Don't forget to evaluate `mean_PPD` (the mean of the posterior predictive distribution) in the generated quantities block.

For efficiency posterior predictions and the log likelihood is computed in R.

### Tips

* **Every time you make a change to a Stan file you need to recompile the package**. (Only running `devtools::build()` is not sufficient.) So do the following:  
    1. Move one level up from the `rstanarm/` directory.
    2. From the command line run `R CMD INSTALL --preclean rstanarm`.
    3. In RStudio run `devtools::build()`.
* In the generated quantities block define the variables you want to report globally and use local scopes (i.e. `{}`) to define (and perform calculations on) N-dimensional objects.

## Priors

This varies from model to model.

In the Stan file you should be able to use existing code in `inst/chunks` to apply priors on the intercept (if it exists) and independent priors on the parameters of the predictors.

User `prior_aux` to take care of single scalar parameters in the model (e.g. the spatial autocorrelation coefficient in spatial models)

If you need the user to define a prior distribution that is not currently available then add the function in `R/priors.R`. (Use the existing functions as a guide). Include the appropriate documentation so that prior distribution is defined in the `?priors` help page.

# `R/stanreg.R`

The main things to deal with here are the **coefficients** and the **linear predictors/fitted values**.

The first few conditionals deal with picking up the estimated coefficients. Sometimes `object$family$family` isn't sufficient to pick up on this so you might have to use `is(object, "class_name")` to determine whether the object is of a certain class (in addition to the class "stanreg").

The `linear.predictors` should be an N-dimensional vector of predictions that have not been transformed by the link function. The `fitted.values` are the linear predictors transformed by the link function. (e.g. if `object$family$family == "gaussian"` then the linear predictor and fitted values will be identical since the link function is the identity function.)

Lastly, at the end of `out <- list(...)` you should include any other stuff that you might need for the methods (e.g. spatial models need the spatial weight matrix, stan_betareg needs the info associated with the "z" linear predictor if declared, etc.)

# Methods

Most of the `*.stanreg` methods are in `R/stanreg-methods.R`, but as long as things are done appropriately in the .fit file and in `stanreg.R` all the methods here should work fine.

## `predict`

The main thing here is to make sure predict works appropriately when the user declares new data. As a rough check, the predictions should match the predictions made by the function you're emulating.

Also, if no new data is declared then `predict(fit)` and `fit$fitted.values` should be identical.

## `posterior_predict`

This is a little more involved than the `predict` method. Essentially you need to return and $N \times S$ dimensional matrix where $N$ is the number of observations and $S$ is the number of draws from the posterior distribution. There are two parts to this:

1. Specify `pp_fun`
   * `pp_fun` will call on the posterior prediction function of the form `.pp_*`. So you need to specify the (stochastic) data generating process within `.pp_*`. We use `sapply()` to iterate over the number of draws and compute the fitted values.

2. Specify `pp_args`
   * Include anything you might need for posterior predictions within the `args` list in the `pp_args` function. (Make sure you do any necessary link function transformations here.)

## `posterior_linpred`

## `loo` and `log_lik`

You need to check whether,

1. `loo()` is using the correct log likelihood specified in `log_lik.R`. This is the log likelihood function that corresponds to `object$family` (or some other identifier that you can subset from `object`). If it does then you're done.
2.  If not then you need to specify the appropriate log likelihood to be used in `loo()`.

Getting the loo function to work on a stanreg object can be tricky. It involves creating a log likelihood function for the posterior `llfun` and a set of arguments to be passed through this function `llargs`.

### `llfun`  
The best way to think about this is that you want to  create a $S \times N$ matrix point-wise log likelihood, where $S$ is the number of draws and $N$ is the number of observations (i.e. you're evaluating the log-likelihood of the posterior for each datum and draw from the marginal posterior).

The approach taken with using loo on a stanreg object is to declare a function that iterates over the data, rather than specifying the entire point-wise log likelihood matrix.

### `llargs`  
Within the `llargs` list `data` needs to be a data frame or matrix that can be iterated over $N$ times. `draws` should be a list containing the draws of $S$ dimension. One way to think about it is that data is what you need to iterate over and draws is fixed. ~~This is useful in cases where some variables may be considered as data but you don't actually want to iterate over them, or in cases where you only have one observation and actually need to iterate over the draws (e.g. a multinormal outcome with correlated errors.)~~

## `prior_summary`

The `prior_summary` function is used to report the prior distributions specified on the parameters when the sampler iterates over the target distribution (which is not necessarily identical to what the user declares).

1. Define a `summarize_*_prior` function at the end of the model's .fit file to capture all the prior information. See `stan_glm.fit` for a comprehensive example or `stan_sp.fit` for a simple example.
     * If the user can call `prior_aux` then you need to give this parameter a name in `$prior_aux$aux_name = "prior_aux_name_here"`. (e.g. in spatial models we have `$prior_aux$aux_name = "rho"` and in stan_betareg we have `$prior_aux$aux_name = "phi"`)
2. Call `prior_info <- summarize_*_prior(...)` before you do any model fitting.
3. At end of the `"optimizing"` and `"sampling"` conditionals make sure you `return(structure(stanfit, prior.info = prior_info))`.

If you do this right then everything should work out swimmingly in the `prior_summary.R` file. If it so happens that you've introduced a new prior then you'll need to update the conditional in the relevant `.prior_*_prior` function to pick this information up.

# Documentation and Examples
We use [roxygen](http://r-pkgs.had.co.nz/man.html) for documentation and examples. Some advice follows,

* The title will probably be something like "Bayesian model of awesomeness".
* Following the title you should add a description of the model. Some things to consider are,  
    * What is the model?
    * What are the equation(s) of the model (if it is not already obvious and if they can be stated clearly)?
    * On what parameters can the user specify priors?
* Don't forget to export the function with `@export`.
* Where possible, always use `@template` (and, if relevant, `@tamplateVar`) to pull in the existing templates from the `man-roxygen` folder.
* Document the additional arguments that are not covered by the templates using `@param`
* In `@seealso`,
    * Point the user to the rstanarm vignettes associated with the model.
    * Mention related models (especially if you're implementing multiple models from a single package).
* In `@details` specify,
    * The R package being emulated.
    * What is being done "under-the-hood" at a high-level.
    * What .fit file the model calls.

Note, every time you make a change to the documentation you need to rebuild the documentation (e.g. run `devtools::document()`) to make sure it works. If you want to check that links to other packages work then you'll have to rebuild the package (e.g. run `devtools::build()`).

# Testing

All tests go in `tests/testthat`. Test everything you could possibly think of. If you think it should be a test, then that probably means it should be a test. Arguably, before writing any code you should write (at the very least) some basic tests first.

Make sure you add a test file for the model you're including as `tests/testthat/test_stan_*.R`. Also don't forget to add the relevant tests for the methods associated with your model in the other test files.

For speed, most of the tests should be specified using `algorithm = 'optimizing'`, `QR = TRUE`, and around 200 iterations.

You should have tests for the following,

* All variations of model specification and comparison with the package you're emulating. (e.g. model with a constant, without a constant, etc.)
* Test predict/posterior_predict with and without new data.
* Test that loo and compare work on the stanreg object.

Run a comprehensive test of rstanarm often. Especially if you're altering Stan files. This will help you catch any bugs early on (which means they'll be easier to fix).

This script should be sufficient (at the time of writing) to run all the tests (excluding the vignettes):

```
### script to run all rstanarm tests locally

library(rstanarm)
# library(rstantools)
# library(bayesplot)

remove(list=ls())
### run prerequisite functions
#
example_model <-
  stan_glmer(cbind(incidence, size - incidence) ~ size
  + period + (1|herd),
             data = lme4::cbpp, family = binomial,
             # this next line is only to keep the example small in size!
             chains = 2, cores = 1, seed = 12345, iter = 500)
#
last_dimnames <- function(x) {
  ndim <- length(dim(x))
  dimnames(x)[[ndim]]
}

### run tests
devtools::test()
```

# Vignettes

This should be pretty straightforward if you use the existing vignettes as a template. You should cover the following,

1. Mathematically define the posterior distribution of the model.
2. An example using simulated data.
3. An example using real data.

In both examples above you should,

* Inspect the data.
* Fit a couple of models.
* Run a posterior predictive checks (e.g. `pp_check`).
* Do basic model comparison with `loo`.

Where possible use existing templates in `vignettes/children`. You can include them with:
```
{r, child="children/*.txt"}
```

# An outline of what goes where

A brief description of what (generally) goes into each the various files/folders.

**`/R`**

* Contains all the model functions and corresponding .fit functions (e.g. `stan_glm` and `stan_glm.fit`).  

**`/exec`** and **`/inst/chunks`**

* The `/inst/chunks` folder contains reusable snippits of Stan code. The `/exec` folder contains the Stan models that are used in the relevant R function. For example, `continuous.stan` contains all the models that can be declared by `stan_glm` (as well as some others). You can view the compiled model in R by executing `rstanarm:::stanmodels$continuous`.

**`data`**

* Example data used in the examples/tests.

**`man-roxygen`**

* Templates for documentation.

**`man`**

* Don't edit any documentation here (changes will get overwritten when rebuilding the package).

**`tests/testthat`**

* Tests run using the [testthat](http://r-pkgs.had.co.nz/tests.html) package.

**`R/misc.R`**

* Contains a bunch of helper functions.
