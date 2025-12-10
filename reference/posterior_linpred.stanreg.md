# Posterior distribution of the (possibly transformed) linear predictor

Extract the posterior draws of the linear predictor, possibly
transformed by the inverse-link function. This function is occasionally
useful, but it should be used sparingly: inference and model checking
should generally be carried out using the posterior predictive
distribution (i.e., using
[`posterior_predict`](https://mc-stan.org/rstanarm/reference/posterior_predict.stanreg.md)).

## Usage

``` r
# S3 method for class 'stanreg'
posterior_linpred(
  object,
  transform = FALSE,
  newdata = NULL,
  draws = NULL,
  re.form = NULL,
  offset = NULL,
  XZ = FALSE,
  ...
)

# S3 method for class 'stanreg'
posterior_epred(
  object,
  newdata = NULL,
  draws = NULL,
  re.form = NULL,
  offset = NULL,
  XZ = FALSE,
  ...
)
```

## Arguments

- object:

  A fitted model object returned by one of the rstanarm modeling
  functions. See
  [`stanreg-objects`](https://mc-stan.org/rstanarm/reference/stanreg-objects.md).

- transform:

  Should the linear predictor be transformed using the inverse-link
  function? The default is `FALSE`. This argument is still allowed but
  not recommended because the `posterior_epred` function now provides
  the equivalent of `posterior_linpred(..., transform=TRUE)`. See
  **Examples**.

- newdata, draws, re.form, offset:

  Same as for
  [`posterior_predict`](https://mc-stan.org/rstanarm/reference/posterior_predict.stanreg.md).

- XZ:

  If `TRUE` then instead of computing the linear predictor the design
  matrix `X` (or `cbind(X,Z)` for models with group-specific terms)
  constructed from `newdata` is returned. The default is `FALSE`.

- ...:

  Currently ignored.

## Value

The default is to return a `draws` by `nrow(newdata)` matrix of
simulations from the posterior distribution of the (possibly
transformed) linear predictor. The exception is if the argument `XZ` is
set to `TRUE` (see the `XZ` argument description above).

## Details

The `posterior_linpred` function returns the posterior distribution of
the linear predictor, while the `posterior_epred` function returns the
posterior distribution of the conditional expectation. In the special
case of a Gaussian likelihood with an identity link function, these two
concepts are the same. The `posterior_epred` function is a less noisy
way to obtain expectations over the output of
[`posterior_predict`](https://mc-stan.org/rstanarm/reference/posterior_predict.stanreg.md).

## Note

For models estimated with
[`stan_clogit`](https://mc-stan.org/rstanarm/reference/stan_clogit.md),
the number of successes per stratum is ostensibly fixed by the research
design. Thus, when calling `posterior_linpred` with new data and
`transform = TRUE`, the `data.frame` passed to the `newdata` argument
must contain an outcome variable and a stratifying factor, both with the
same name as in the original `data.frame`. Then, the probabilities will
condition on this outcome in the new data.

## See also

[`posterior_predict`](https://mc-stan.org/rstanarm/reference/posterior_predict.stanreg.md)
to draw from the posterior predictive distribution of the outcome, which
is typically preferable.

## Examples

``` r
if (.Platform$OS.type != "windows" || .Platform$r_arch != "i386") {
if (!exists("example_model")) example(example_model)
print(family(example_model))

# linear predictor on log-odds scale
linpred <- posterior_linpred(example_model)
colMeans(linpred)

# probabilities
# same as posterior_linpred(example_model, transform = TRUE)
probs <- posterior_epred(example_model) 
colMeans(probs)

# not conditioning on any group-level parameters
probs2 <- posterior_epred(example_model, re.form = NA)
apply(probs2, 2, median)
}
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#>          1          2          3          4          5          6          7 
#> 0.19321385 0.07963232 0.07058646 0.04248258 0.20038371 0.08293944 0.07416584 
#>          8          9         10         11         12         13         14 
#> 0.20038371 0.08184582 0.07185831 0.04623685 0.19090399 0.07887589 0.07058646 
#>         15         16         17         18         19         20         21 
#> 0.04277052 0.19618560 0.08564462 0.07561027 0.04226464 0.19533469 0.08247894 
#>         22         23         24         25         26         27         28 
#> 0.07292688 0.04623685 0.19424230 0.07887589 0.07058646 0.04248258 0.21139640 
#>         29         30         31         32         33         34         35 
#> 0.19034663 0.07727663 0.07032113 0.04277052 0.20038371 0.08387899 0.07292688 
#>         36         37         38         39         40         41         42 
#> 0.04711961 0.20205369 0.08630690 0.07462140 0.04711961 0.19090399 0.07843710 
#>         43         44         45         46         47         48         49 
#> 0.06977302 0.04248258 0.19933598 0.08501582 0.07355331 0.04754166 0.19751916 
#>         50         51         52         53         54         55         56 
#> 0.07599325 0.06856146 0.04164939 0.19751916 0.08096430 0.07165066 0.04567692 
```
