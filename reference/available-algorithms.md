# Estimation algorithms available for rstanarm models

Estimation algorithms available for rstanarm models

## Estimation algorithms

The modeling functions in the rstanarm package take an `algorithm`
argument that can be one of the following:

- **Sampling** (`algorithm="sampling"`):

  Uses Markov Chain Monte Carlo (MCMC) — in particular, Hamiltonian
  Monte Carlo (HMC) with a tuned but diagonal mass matrix — to draw from
  the posterior distribution of the parameters. See `sampling` (rstan)
  for more details. This is the slowest but most reliable of the
  available estimation algorithms and it is **the default and
  recommended algorithm for statistical inference.**

- **Mean-field** (`algorithm="meanfield"`):

  Uses mean-field variational inference to draw from an approximation to
  the posterior distribution. In particular, this algorithm finds the
  set of independent normal distributions in the unconstrained space
  that — when transformed into the constrained space — most closely
  approximate the posterior distribution. Then it draws repeatedly from
  these independent normal distributions and transforms them into the
  constrained space. The entire process is much faster than HMC and
  yields independent draws but **is not recommended for final
  statistical inference**. It can be useful to narrow the set of
  candidate models in large problems, particularly when specifying
  `QR=TRUE` in
  [`stan_glm`](https://mc-stan.org/rstanarm/reference/stan_glm.md),
  [`stan_glmer`](https://mc-stan.org/rstanarm/reference/stan_glmer.md),
  and
  [`stan_gamm4`](https://mc-stan.org/rstanarm/reference/stan_gamm4.md),
  but is **only an approximation to the posterior distribution**.

- **Full-rank** (`algorithm="fullrank"`):

  Uses full-rank variational inference to draw from an approximation to
  the posterior distribution by finding the multivariate normal
  distribution in the unconstrained space that — when transformed into
  the constrained space — most closely approximates the posterior
  distribution. Then it draws repeatedly from this multivariate normal
  distribution and transforms the draws into the constrained space. This
  process is slower than meanfield variational inference but is faster
  than HMC. Although still an approximation to the posterior
  distribution and thus **not recommended for final statistical
  inference**, the approximation is more realistic than that of
  mean-field variational inference because the parameters are not
  assumed to be independent in the unconstrained space. Nevertheless,
  fullrank variational inference is a more difficult optimization
  problem and the algorithm is more prone to non-convergence or
  convergence to a local optimum.

- **Optimizing** (`algorithm="optimizing"`):

  Finds the posterior mode using a C++ implementation of the LBGFS
  algorithm. See `optimizing` for more details. If there is no prior
  information, then this is equivalent to maximum likelihood, in which
  case there is no great reason to use the functions in the rstanarm
  package over the emulated functions in other packages. However, if
  priors are specified, then the estimates are penalized maximum
  likelihood estimates, which may have some redeeming value. Currently,
  optimization is only supported for
  [`stan_glm`](https://mc-stan.org/rstanarm/reference/stan_glm.md).

## See also

<https://mc-stan.org/rstanarm/>
