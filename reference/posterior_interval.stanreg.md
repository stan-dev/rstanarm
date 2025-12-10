# Posterior uncertainty intervals

For models fit using MCMC (`algorithm="sampling"`) or one of the
variational approximations (`"meanfield"` or `"fullrank"`), the
`posterior_interval` function computes Bayesian posterior uncertainty
intervals. These intervals are often referred to as *credible*
intervals, but we use the term *uncertainty* intervals to highlight the
fact that wider intervals correspond to greater uncertainty.

## Usage

``` r
# S3 method for class 'stanreg'
posterior_interval(
  object,
  prob = 0.9,
  type = "central",
  pars = NULL,
  regex_pars = NULL,
  ...
)
```

## Arguments

- object:

  A fitted model object returned by one of the rstanarm modeling
  functions. See
  [`stanreg-objects`](https://mc-stan.org/rstanarm/reference/stanreg-objects.md).

- prob:

  A number \\p \in (0,1)\\ indicating the desired probability mass to
  include in the intervals. The default is to report \\90\\% intervals
  (`prob=0.9`) rather than the traditionally used \\95\\% (see Details).

- type:

  The type of interval to compute. Currently the only option is
  `"central"` (see Details). A central \\100p\\% interval is defined by
  the \\\alpha/2\\ and \\1 - \alpha/2\\ quantiles, where \\\alpha = 1 -
  p\\.

- pars:

  An optional character vector of parameter names.

- regex_pars:

  An optional character vector of [regular
  expressions](https://rdrr.io/r/base/grep.html) to use for parameter
  selection. `regex_pars` can be used in place of `pars` or in addition
  to `pars`. Currently, all functions that accept a `regex_pars`
  argument ignore it for models fit using optimization.

- ...:

  Currently ignored.

## Value

A matrix with two columns and as many rows as model parameters (or the
subset of parameters specified by `pars` and/or `regex_pars`). For a
given value of `prob`, \\p\\, the columns correspond to the lower and
upper \\100p\\% interval limits and have the names \\100\alpha/2\\% and
\\100(1 - \alpha/2)\\%, where \\\alpha = 1-p\\. For example, if
`prob=0.9` is specified (a \\90\\% interval), then the column names will
be `"5%"` and `"95%"`, respectively.

## Details

### Interpretation

Unlike for a frenquentist confidence interval, it is valid to say that,
conditional on the data and model, we believe that with probability
\\p\\ the value of a parameter is in its \\100p\\% posterior interval.
This intuitive interpretation of Bayesian intervals is often erroneously
applied to frequentist confidence intervals. See Morey et al. (2015) for
more details on this issue and the advantages of using Bayesian
posterior uncertainty intervals (also known as credible intervals).

### Default 90% intervals

We default to reporting \\90\\% intervals rather than \\95\\% intervals
for several reasons:

- Computational stability: \\90\\% intervals are more stable than
  \\95\\% intervals (for which each end relies on only \\2.5\\% of the
  posterior draws).

- Relation to Type-S errors (Gelman and Carlin, 2014): \\95\\% of the
  mass in a \\90\\% central interval is above the lower value (and
  \\95\\% is below the upper value). For a parameter \\\theta\\, it is
  therefore easy to see if the posterior probability that \\\theta \>
  0\\ (or \\\theta \< 0\\) is larger or smaller than \\95\\%.

Of course, if \\95\\% intervals are desired they can be computed by
specifying `prob=0.95`.

### Types of intervals

Currently `posterior_interval` only computes central intervals because
other types of intervals are rarely useful for the models that rstanarm
can estimate. Additional possibilities may be provided in future
releases as more models become available.

## References

Gelman, A. and Carlin, J. (2014). Beyond power calculations: assessing
Type S (sign) and Type M (magnitude) errors. *Perspectives on
Psychological Science*. 9(6), 641–51.

Morey, R. D., Hoekstra, R., Rouder, J., Lee, M. D., and Wagenmakers, E.
(2016). The fallacy of placing confidence in confidence intervals.
*Psychonomic Bulletin & Review*. 23(1), 103–123.

## See also

[`confint.stanreg`](https://mc-stan.org/rstanarm/reference/stanreg-methods.md),
which, for models fit using optimization, can be used to compute
traditional confidence intervals.

[`predictive_interval`](https://mc-stan.org/rstanarm/reference/predictive_interval.stanreg.md)
for predictive intervals.

## Examples

``` r
if (.Platform$OS.type != "windows" || .Platform$r_arch != "i386") {
if (!exists("example_model")) example(example_model)
posterior_interval(example_model)
posterior_interval(example_model, regex_pars = "herd")
posterior_interval(example_model, pars = "period2", prob = 0.5)
}
#>               25%        75%
#> period2 -1.185064 -0.7729716
```
