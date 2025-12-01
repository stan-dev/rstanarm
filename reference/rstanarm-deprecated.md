# Deprecated functions

These functions are deprecated and will be removed in a future release.
The **Arguments** section below provides details on how the
functionality obtained via each of the arguments has been replaced.

## Usage

``` r
prior_options(
  prior_scale_for_dispersion = 5,
  min_prior_scale = 1e-12,
  scaled = TRUE
)
```

## Arguments

- prior_scale_for_dispersion, min_prior_scale, scaled:

  Arguments to deprecated `prior_options` function. The functionality
  provided by the now deprecated `prior_options` function has been
  replaced as follows:

  `prior_scale_for_dispersion`

  :   Instead of using the `prior_scale_for_dispersion` argument to
      `prior_options`, priors for these parameters can now be specified
      directly when calling
      [`stan_glm`](https://mc-stan.org/rstanarm/reference/stan_glm.md)
      (or
      [`stan_glmer`](https://mc-stan.org/rstanarm/reference/stan_glmer.md),
      etc.) using the new `prior_aux` argument.

  `scaled`

  :   Instead of setting `prior_options(scaled=FALSE)`, internal
      rescaling is now toggled using the new `autoscale` arguments to
      [`normal`](https://mc-stan.org/rstanarm/reference/priors.md),
      [`student_t`](https://mc-stan.org/rstanarm/reference/priors.md),
      and [`cauchy`](https://mc-stan.org/rstanarm/reference/priors.md)
      (the other prior distributions do not support 'autoscale').

  `min_prior_scale`

  :   No replacement. `min_prior_scale` (the minimum possible scale
      parameter value that be used for priors) is now fixed to `1e-12`.
