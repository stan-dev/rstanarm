# rstanarm <img src="man/figures/stanlogo.png" align="right" width="120" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/mcdonohue/rstanarm/workflows/R-CMD-check/badge.svg)](https://github.com/mcdonohue/rstanarm/actions)
[![Build Status](https://travis-ci.org/mcdonoue/rstanarm.svg?branch=master)](https://travis-ci.org/mcdonoue/rstanarm)
<!-- badges: end -->

### Latent Time Joint Mixed Effect Models (LTJMM)

This fork of the [rstanarm](https://github.com/stan-dev/rstanarm) package includes the following modifications:

1. **stan_mvmer** models are extended from 5 to 20 longitudinal submodels
2. **stan_ljtmm** extends **stan_mvmer** to accommodate a group-specific (individual-specific) latent time parameters which are shared within group across submodels.

The LTJMM is described in [Li, et al. (2017)](https://doi.org/10.1177/0962280217737566).

### Resources

* [Open an issue](https://github.com/stan-dev/mcdonohue/issues) (GitHub issues for bug reports, feature requests)

### Installation

To install from GitHub, first make sure that you can install the **rstan**
package and C++ toolchain by following these
[instructions](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started).
Once **rstan** is successfully installed, you can install **rstanarm** from
GitHub using the **remotes** package by executing the following in R:

```r
# Change 2 to however many cores you can/want to use to parallelize install
# If you experience crashes or run out RAM during installation, try changing this to 1
Sys.setenv(MAKEFLAGS = "-j2")
Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS" = "true")
remotes::install_github("mcdonohue/rstanarm", INSTALL_opts = "--no-multiarch", force = TRUE)
```

You can switch `build_vignettes` to `TRUE` but it takes a lot longer to install and the 
vignettes are already separately available from the 
[Stan website](https://mc-stan.org/rstanarm/articles/index.html) 
and 
[CRAN](https://cran.r-project.org/package=rstanarm/vignettes). 
If installation fails, please let us know by [filing an issue](https://github.com/mcdonohue/rstanarm/issues).
