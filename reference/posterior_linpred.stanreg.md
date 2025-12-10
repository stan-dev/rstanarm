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
#> 0.19279976 0.08131904 0.06966905 0.04371449 0.20133549 0.08448797 0.07581547 
#>          8          9         10         11         12         13         14 
#> 0.20133549 0.08344953 0.07267840 0.04985763 0.19017467 0.08051768 0.06966905 
#>         15         16         17         18         19         20         21 
#> 0.04402754 0.19691580 0.08897117 0.07717166 0.04365934 0.19562395 0.08444267 
#>         22         23         24         25         26         27         28 
#> 0.07349302 0.04985763 0.19500139 0.08051768 0.06966905 0.04371449 0.21312495 
#>         29         30         31         32         33         34         35 
#> 0.19053415 0.07940598 0.06945606 0.04402754 0.20133549 0.08761524 0.07349302 
#>         36         37         38         39         40         41         42 
#> 0.05033955 0.20381612 0.09008201 0.07665394 0.05033955 0.19017467 0.08012851 
#>         43         44         45         46         47         48         49 
#> 0.06800971 0.04371449 0.19996042 0.08856720 0.07389744 0.05058364 0.19895168 
#>         50         51         52         53         54         55         56 
#> 0.07760974 0.06751417 0.04363252 0.19895168 0.08261023 0.07272047 0.04707914 
```
