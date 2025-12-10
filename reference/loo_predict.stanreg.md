# Compute weighted expectations using LOO

These functions are wrappers around the
[`E_loo`](https://mc-stan.org/loo/reference/E_loo.html) function (loo
package) that provide compatibility for rstanarm models.

## Usage

``` r
# S3 method for class 'stanreg'
loo_predict(
  object,
  type = c("mean", "var", "quantile"),
  probs = 0.5,
  ...,
  psis_object = NULL
)

# S3 method for class 'stanreg'
loo_linpred(
  object,
  type = c("mean", "var", "quantile"),
  probs = 0.5,
  transform = FALSE,
  ...,
  psis_object = NULL
)

# S3 method for class 'stanreg'
loo_predictive_interval(object, prob = 0.9, ..., psis_object = NULL)
```

## Arguments

- object:

  A fitted model object returned by one of the rstanarm modeling
  functions. See
  [`stanreg-objects`](https://mc-stan.org/rstanarm/reference/stanreg-objects.md).

- type:

  The type of expectation to compute. The options are `"mean"`,
  `"variance"`, `"sd"`, and `"quantile"`.

- probs:

  For computing quantiles, a vector of probabilities.

- ...:

  Currently unused.

- psis_object:

  An object returned by
  [`psis`](https://mc-stan.org/loo/reference/psis.html). If missing then
  `psis` will be run internally, which may be time consuming for models
  fit to very large datasets.

- transform:

  Passed to
  [`posterior_linpred`](https://mc-stan.org/rstanarm/reference/posterior_linpred.stanreg.md).

- prob:

  For `loo_predictive_interval`, a scalar in \\(0,1)\\ indicating the
  desired probability mass to include in the intervals. The default is
  `prob=0.9` (\\90\\% intervals).

## Value

A list with elements `value` and `pareto_k`.

For `loo_predict` and `loo_linpred` the value component is a vector with
one element per observation.

For `loo_predictive_interval` the `value` component is a matrix with one
row per observation and two columns (like
[`predictive_interval`](https://mc-stan.org/rstanarm/reference/predictive_interval.stanreg.md)).
`loo_predictive_interval(..., prob = p)` is equivalent to
`loo_predict(..., type = "quantile", probs = c(a, 1-a))` with
`a = (1 - p)/2`, except it transposes the result and adds informative
column names.

See [`E_loo`](https://mc-stan.org/loo/reference/E_loo.html) and
[`pareto-k-diagnostic`](https://mc-stan.org/loo/reference/pareto-k-diagnostic.html)
for details on the `pareto_k` diagnostic.

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

## Examples

``` r
if (.Platform$OS.type != "windows" || .Platform$r_arch != "i386") {
# \dontrun{
if (!exists("example_model")) example(example_model)

# optionally, log-weights can be pre-computed and reused
psis_result <- loo::psis(log_ratios = -log_lik(example_model))

loo_probs <- loo_linpred(example_model, type = "mean", transform = TRUE, psis_object = psis_result)
str(loo_probs)

loo_pred_var <- loo_predict(example_model, type = "var", psis_object = psis_result)
str(loo_pred_var)

loo_pred_ints <- loo_predictive_interval(example_model, prob = 0.8, psis_object = psis_result)
str(loo_pred_ints)
# }
}
#> Warning: Some Pareto k diagnostic values are too high. See help('pareto-k-diagnostic') for details.
#> Instead of posterior_linpred(..., transform=TRUE) please call posterior_epred(), which provides equivalent functionality.
#> List of 2
#>  $ value   : num [1:56] 0.4294 0.1141 0.0679 0.0993 0.1687 ...
#>  $ pareto_k: num [1:56] 0.853 0.344 0.617 0.201 0.504 ...
#> List of 2
#>  $ value   : num [1:56] 5.062 1.676 0.597 0.499 6.039 ...
#>  $ pareto_k: num [1:56] 0.799 0.344 0.617 0.131 0.829 ...
#> List of 2
#>  $ value   : num [1:56, 1:2] 3 0 0 0 1 0 0 1 0 0 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : chr [1:2] "10%" "90%"
#>  $ pareto_k: num [1:56] 0.791 0.344 0.617 0.131 0.425 ...
```
