<a href="http://mc-stan.org">
<img src="https://raw.githubusercontent.com/stan-dev/logos/master/logo.png" width=200 alt="Stan Logo"/>
</a>

# rstanarm

[![Build Status](https://travis-ci.org/stan-dev/rstanarm.svg?branch=master)](https://travis-ci.org/stan-dev/rstanarm)  
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/rstanarm)](http://cran.r-project.org/package=rstanarm)
[![Downloads](http://cranlogs.r-pkg.org/badges/rstanarm?color=brightgreen)](http://cran.rstudio.com/package=rstanarm)


This is an R package that emulates other R model-fitting functions but uses
[Stan](http://mc-stan.org) (via the **rstan** package) for the back-end
estimation. The primary target audience is people who would be open to Bayesian
inference if using Bayesian software were easier but would use frequentist
software otherwise.

## Installation

### Latest Release
The most recent **rstanarm** release can be installed from CRAN via
```{r}
install.packages("rstanarm")
```

### Development Version
To install from GitHub, first make sure that you can install the **rstan**
package and C++ toolchain by following these
[instructions](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started).
Once **rstan** is successfully installed, you can install **rstanarm** from
GitHub using the **devtools** package by executing the following in R:

```{r}
if (!require(devtools)) {
  install.packages("devtools")
  library(devtools)
}
install_github("stan-dev/rstanarm", args = "--preclean", build_vignettes = FALSE)
```

Make sure to include the `args = "--preclean"` and `build_vignettes = FALSE` arguments 
or the development version of package will not install properly. If installation fails, 
please let us know by [filing an issue](https://github.com/stan-dev/rstanarm/issues).

## Getting help

* [File an issue on GitHub](https://github.com/stan-dev/rstanarm/issues)
* [Ask a question on the Stan users Google group](https://groups.google.com/forum/#!forum/stan-users)

## Contributing 

If you are interested in contributing to the development of **rstanarm** please 
see the [Contributing to
development](https://github.com/stan-dev/rstanarm/wiki/Contributing-to-development)
page of the wiki.
