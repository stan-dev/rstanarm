# `adapt_delta`: Target average acceptance probability

Details about the `adapt_delta` argument to rstanarm's modeling
functions.

## Details

For the No-U-Turn Sampler (NUTS), the variant of Hamiltonian Monte Carlo
used used by rstanarm, `adapt_delta` is the target average proposal
acceptance probability during Stan's adaptation period. `adapt_delta` is
ignored by rstanarm if the `algorithm` argument is not set to
`"sampling"`.

The default value of `adapt_delta` is 0.95, except when the prior for
the regression coefficients is
[`R2`](https://mc-stan.org/rstanarm/reference/priors.md),
[`hs`](https://mc-stan.org/rstanarm/reference/priors.md), or
[`hs_plus`](https://mc-stan.org/rstanarm/reference/priors.md), in which
case the default is 0.99.

These defaults are higher (more conservative) than the default of
`adapt_delta=0.8` used in the rstan package, which may result in slower
sampling speeds but will be more robust to posterior distributions with
high curvature.

In general you should not need to change `adapt_delta` unless you see a
warning message about divergent transitions, in which case you can
increase `adapt_delta` from the default to a value *closer* to 1 (e.g.
from 0.95 to 0.99, or from 0.99 to 0.999, etc). The step size used by
the numerical integrator is a function of `adapt_delta` in that
increasing `adapt_delta` will result in a smaller step size and fewer
divergences. Increasing `adapt_delta` will typically result in a slower
sampler, but it will always lead to a more robust sampler.

## References

Stan Development Team. *Stan Modeling Language Users Guide and Reference
Manual.* <https://mc-stan.org/users/documentation/>.

Brief Guide to Stan's Warnings:
<https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup>
