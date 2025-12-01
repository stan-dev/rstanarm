# Information criteria and cross-validation

For models fit using MCMC, compute approximate leave-one-out
cross-validation (LOO, LOOIC) or, less preferably, the Widely Applicable
Information Criterion (WAIC) using the
[loo](https://mc-stan.org/loo/reference/loo-package.html) package. (For
\\K\\-fold cross-validation see
[`kfold.stanreg`](https://mc-stan.org/rstanarm/reference/kfold.stanreg.md).)
Functions for model comparison, and model weighting/averaging are also
provided.

**Note**: these functions are not guaranteed to work properly unless the
`data` argument was specified when the model was fit. Also, as of loo
version `2.0.0` the default number of cores is now only 1, but we
recommend using as many (or close to as many) cores as possible by
setting the `cores` argument or using `options(mc.cores = VALUE)` to set
it for an entire session.

## Usage

``` r
# S3 method for class 'stanreg'
loo(
  x,
  ...,
  cores = getOption("mc.cores", 1),
  save_psis = FALSE,
  k_threshold = NULL,
  r_eff = FALSE
)

# S3 method for class 'stanreg'
waic(x, ...)

# S3 method for class 'stanreg'
loo_compare(x, ..., criterion = c("loo", "kfold", "waic"), detail = FALSE)

# S3 method for class 'stanreg_list'
loo_compare(x, ..., criterion = c("loo", "kfold", "waic"), detail = FALSE)

# S3 method for class 'stanreg_list'
loo_model_weights(x, ..., cores = getOption("mc.cores", 1), k_threshold = NULL)

compare_models(..., loos = list(), detail = FALSE)
```

## Arguments

- x:

  For `loo` and `waic`, a fitted model object returned by one of the
  rstanarm modeling functions. See
  [stanreg-objects](https://mc-stan.org/rstanarm/reference/stanreg-objects.md).

  For the `loo_model_weights` method, `x` should be a "stanreg_list"
  object, which is a list of fitted model objects created by
  [`stanreg_list`](https://mc-stan.org/rstanarm/reference/stanreg_list.md).
  `loo_compare` also allows `x` to be a single stanreg object, with the
  remaining objects passed via `...`, or a single `stanreg_list` object.

- ...:

  For `loo_compare.stanreg`, `...` can contain objects returned by the
  `loo`,
  [`kfold`](https://mc-stan.org/rstanarm/reference/kfold.stanreg.md), or
  `waic` method (see the **Examples** section, below).

  For `loo_model_weights`, `...` should contain arguments (e.g.
  `method`) to pass to the default
  [`loo_model_weights`](https://mc-stan.org/loo/reference/loo_model_weights.html)
  method from the loo package.

- cores, save_psis:

  Passed to [`loo`](https://mc-stan.org/loo/reference/loo.html).

- k_threshold:

  Threshold for flagging estimates of the Pareto shape parameters \\k\\
  estimated by `loo`. See the *How to proceed when `loo` gives warnings*
  section, below, for details.

- r_eff:

  `TRUE` or `FALSE` indicating whether to compute the `r_eff` argument
  to pass to the loo package. If `TRUE`, rstanarm will call
  [`relative_eff`](https://mc-stan.org/loo/reference/relative_eff.html)
  to compute the `r_eff` argument to pass to the loo package. If `FALSE`
  (the default), we avoid computing `r_eff`, which can be very slow.
  `r_eff` measures the amount of autocorrelation in MCMC draws, and is
  used to compute more accurate ESS and MCSE estimates for pointwise and
  total ELPDs. When `r_eff=FALSE`, the reported ESS and MCSE estimates
  may be over-optimistic if the posterior draws are far from
  independent.

- criterion:

  For `loo_compare.stanreg` and `loo_compare.stanreg_list`, should the
  comparison be based on LOO-CV (`criterion="loo"`), K-fold-CV
  (`criterion="kfold"`), or WAIC (`criterion="waic"`). The default is
  LOO-CV. See the **Comparing models** and **Examples** sections below.

- detail:

  For `loo_compare.stanreg` and `loo_compare.stanreg_list`, if `TRUE`
  then extra information about each model (currently just the model
  formulas) will be printed with the output.

- loos:

  a list of objects produced by the `loo` function

## Value

The structure of the objects returned by `loo` and `waic` methods are
documented in detail in the **Value** section in
[`loo`](https://mc-stan.org/loo/reference/loo.html) and
[`waic`](https://mc-stan.org/loo/reference/waic.html) (from the loo
package).

`loo_compare` returns a matrix with class 'compare.loo'. See the
**Comparing models** section below for more details.

## Approximate LOO CV

The `loo` method for stanreg objects provides an interface to the
[loo](https://mc-stan.org/loo/reference/loo-package.html) package for
approximate leave-one-out cross-validation (LOO). The LOO Information
Criterion (LOOIC) has the same purpose as the Akaike Information
Criterion (AIC) that is used by frequentists. Both are intended to
estimate the expected log predictive density (ELPD) for a new dataset.
However, the AIC ignores priors and assumes that the posterior
distribution is multivariate normal, whereas the functions from the loo
package do not make this distributional assumption and integrate over
uncertainty in the parameters. This only assumes that any one
observation can be omitted without having a major effect on the
posterior distribution, which can be judged using the diagnostic plot
provided by the
[`plot.loo`](https://mc-stan.org/loo/reference/pareto-k-diagnostic.html)
method and the warnings provided by the
[`print.loo`](https://mc-stan.org/loo/reference/print.loo.html) method
(see the *How to Use the rstanarm Package* vignette for an example of
this process).

### How to proceed when `loo` gives warnings (k_threshold)

The `k_threshold` argument to the `loo` method for rstanarm models is
provided as a possible remedy when the diagnostics reveal problems
stemming from the posterior's sensitivity to particular observations.
Warnings about Pareto \\k\\ estimates indicate observations for which
the approximation to LOO is problematic (this is described in detail in
Vehtari, Gelman, and Gabry (2017) and the
[loo](https://mc-stan.org/loo/reference/loo-package.html) package
documentation). The `k_threshold` argument can be used to set the \\k\\
value above which an observation is flagged. If `k_threshold` is not
`NULL` and there are \\J\\ observations with \\k\\ estimates above
`k_threshold` then when `loo` is called it will refit the original model
\\J\\ times, each time leaving out one of the \\J\\ problematic
observations. The pointwise contributions of these observations to the
total ELPD are then computed directly and substituted for the previous
estimates from these \\J\\ observations that are stored in the object
created by `loo`. Another option to consider is K-fold cross-validation,
which is documented on a separate page (see
[`kfold`](https://mc-stan.org/rstanarm/reference/kfold.stanreg.md)).

**Note**: in the warning messages issued by `loo` about large Pareto
\\k\\ estimates we recommend setting `k_threshold` to at least \\0.7\\.
There is a theoretical reason, explained in Vehtari, Gelman, and Gabry
(2017), for setting the threshold to the stricter value of \\0.5\\, but
in practice they find that errors in the LOO approximation start to
increase non-negligibly when \\k \> 0.7\\.

## Comparing models

"loo" (or "waic" or "kfold") objects can be passed to the
[`loo_compare`](https://mc-stan.org/loo/reference/loo_compare.html)
function in the loo package to perform model comparison. rstanarm also
provides a `loo_compare.stanreg` method that can be used if the "loo"
(or "waic" or "kfold") object has been added to the fitted model object
(see the **Examples** section below for how to do this). This second
method allows rstanarm to perform some extra checks that can't be done
by the loo package itself (e.g., verifying that all models to be
compared were fit using the same outcome variable).

`loo_compare` will return a matrix with one row per model and columns
containing the ELPD difference and the standard error of the difference.
In the first row of the matrix will be the model with the largest ELPD
(smallest LOOIC) and will contain zeros (there is no difference between
this model and itself). For each of the remaining models the ELPD
difference and SE are reported relative to the model with the best ELPD
(the first row). See the **Details** section at the
[`loo_compare`](https://mc-stan.org/loo/reference/loo_compare.html) page
in the loo package for more information.

## Model weights

The `loo_model_weights` method can be used to compute model weights for
a `"stanreg_list"` object, which is a list of fitted model objects made
with
[`stanreg_list`](https://mc-stan.org/rstanarm/reference/stanreg_list.md).
The end of the **Examples** section has a demonstration. For details see
the
[`loo_model_weights`](https://mc-stan.org/loo/reference/loo_model_weights.html)
documentation in the loo package.

## References

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

## See also

- The [loo package vignettes](https://mc-stan.org/loo/articles/) and
  various [rstanarm vignettes](https://mc-stan.org/rstanarm/articles/)
  for more examples using `loo` and related functions with rstanarm
  models.

- [`pareto-k-diagnostic`](https://mc-stan.org/loo/reference/pareto-k-diagnostic.html)
  in the loo package for more on Pareto \\k\\ diagnostics.

- [`log_lik.stanreg`](https://mc-stan.org/rstanarm/reference/log_lik.stanreg.md)
  to directly access the pointwise log-likelihood matrix.

## Examples

``` r
if (.Platform$OS.type != "windows" || .Platform$r_arch != "i386") {
# \donttest{
fit1 <- stan_glm(mpg ~ wt, data = mtcars, refresh = 0)
fit2 <- stan_glm(mpg ~ wt + cyl, data = mtcars, refresh = 0)

# (for bigger models use as many cores as possible)
loo1 <- loo(fit1, cores = 1)
print(loo1)
loo2 <- loo(fit2, cores = 1)
print(loo2)

# when comparing models the loo objects can be passed to loo_compare
# as individual arguments or as a list of loo objects
loo_compare(loo1, loo2)
loo_compare(list(loo1, loo2))

# if the fitted model objects contain a loo object in the component "loo"
# then the model objects can be passed directly or as a stanreg_list
fit1$loo <- loo1
fit2$loo <- loo2
loo_compare(fit1, fit2)

# if the fitted model objects contain a loo object _and_ a waic or kfold
# object, then the criterion argument determines which of them the comparison
# is based on 
fit1$waic <- waic(fit1)
fit2$waic <- waic(fit2)
loo_compare(fit1, fit2, criterion = "waic")

# the models can also be combined into a stanreg_list object, and more 
# informative model names can be provided to use when printing
model_list <- stanreg_list(fit1, fit2, model_names = c("Fewer predictors", "More predictors"))
loo_compare(model_list)

fit3 <- stan_glm(mpg ~ disp * as.factor(cyl), data = mtcars, refresh = 0)
loo3 <- loo(fit3, cores = 2, k_threshold = 0.7)
loo_compare(loo1, loo2, loo3)

# setting detail=TRUE will also print model formulas if used with
# loo_compare.stanreg or loo_compare.stanreg_list
fit3$loo <- loo3
model_list <- stanreg_list(fit1, fit2, fit3)
loo_compare(model_list, detail=TRUE)

# Computing model weights
#
# if the objects in model_list already have 'loo' components then those
# will be used. otherwise loo will be computed for each model internally
# (in which case the 'cores' argument may also be used and is passed to loo())
loo_model_weights(model_list)  # defaults to method="stacking"
loo_model_weights(model_list,  method = "pseudobma")
loo_model_weights(model_list,  method = "pseudobma", BB = FALSE)

# you can also pass precomputed loo objects directly to loo_model_weights
loo_list <- list(A = loo1, B = loo2, C = loo3) # names optional (affects printing)
loo_model_weights(loo_list)
# }
}
#> 
#> Computed from 4000 by 32 log-likelihood matrix.
#> 
#>          Estimate  SE
#> elpd_loo    -83.4 4.3
#> p_loo         3.2 1.1
#> looic       166.9 8.5
#> ------
#> MCSE of elpd_loo is 0.0.
#> MCSE and ESS estimates assume independent draws (r_eff=1).
#> 
#> All Pareto k estimates are good (k < 0.7).
#> See help('pareto-k-diagnostic') for details.
#> 
#> Computed from 4000 by 32 log-likelihood matrix.
#> 
#>          Estimate  SE
#> elpd_loo    -78.6 4.6
#> p_loo         4.1 1.3
#> looic       157.1 9.2
#> ------
#> MCSE of elpd_loo is 0.1.
#> MCSE and ESS estimates assume independent draws (r_eff=1).
#> 
#> All Pareto k estimates are good (k < 0.7).
#> See help('pareto-k-diagnostic') for details.
#> Warning: 
#> 3 (9.4%) p_waic estimates greater than 0.4. We recommend trying loo instead.
#> Warning: 
#> 3 (9.4%) p_waic estimates greater than 0.4. We recommend trying loo instead.
#> 1 problematic observation(s) found.
#> Model will be refit 1 times.
#> 
#> Fitting model 1 out of 1 (leaving out observation 4)
#> Method: stacking
#> ------
#>   weight
#> A 0.000 
#> B 0.487 
#> C 0.513 
```
