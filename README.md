<a href="http://mc-stan.org">
<img src="https://raw.githubusercontent.com/stan-dev/logos/master/logo.png" width=200 alt="Stan Logo"/>
</a>

# rstanarm

[![Build Status](https://travis-ci.org/stan-dev/rstanarm.svg?branch=master)](https://travis-ci.org/stan-dev/rstanarm) 
[![Coverage Status](https://codecov.io/github/stan-dev/rstanarm/coverage.svg?branch=master)](https://codecov.io/github/stan-dev/rstanarm?branch=master) 
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/rstanarm)](http://cran.r-project.org/package=rstanarm)

This is an R package that emulates other R model-fitting functions but uses [Stan](http://mc-stan.org) (via the **rstan** package) 
for the back-end estimation. The primary target audience is people who would be open to Bayesian inference if using Bayesian 
software were easier but would use frequentist software otherwise. That said, this R package often uses posterior medians as 
point estimates.

### Installation

**rstanarm** is not yet on CRAN. To install it, first make sure that you can install the **rstan** package by following these [instructions](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started). Once **rstan** is successfully installed, you can 
install **rstanarm** from GitHub using the **devtools** package by executing the following in R:

```{r}
if (!require(devtools)) {
  install.packages("devtools")
  library(devtools)
}
devtools::install_github("stan-dev/rstanarm", args = "--preclean")
```

Make sure to include the `args = "--preclean"` argument or the package will not install properly. If installation fails, please let us know by [filing an issue](https://github.com/stan-dev/rstanarm/issues).

### Contributing 

If you are interested in contributing to the development of **rstanarm** please see the [Contributing to development](https://github.com/stan-dev/rstanarm/wiki/Contributing-to-development) page of the wiki.
