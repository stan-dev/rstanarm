[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/rstanarm?color=blue)](https://cran.r-project.org/package=rstanarm)
[![Downloads](https://cranlogs.r-pkg.org/badges/rstanarm?color=blue)](https://cran.rstudio.com/package=rstanarm)

<br>

<div style="text-align:left">
<span><a href="https://mc-stan.org">
<img src="https://raw.githubusercontent.com/stan-dev/logos/master/logo_tm.png" width=100 alt="Stan Logo"/> </a><h2><strong>rstanarm</strong></h2>
<h4>Bayesian applied regression modeling</h4></span>
</div>

<br>

**rstanarm** is an R package that emulates other R model-fitting functions but
uses Stan (via the [rstan](https://mc-stan.org/rstan/) package)
for the back-end estimation. The primary target audience is people who would be
open to Bayesian inference if using Bayesian software were easier but would use
frequentist software otherwise.

Fitting models with **rstanarm** is also useful for experienced Bayesian
software users who want to take advantage of the pre-compiled Stan programs that
are written by Stan developers and carefully implemented to prioritize numerical
stability and the avoidance of sampling problems.

## Getting started

If you are new to **rstanarm** we recommend starting with the tutorial
[vignettes](https://mc-stan.org/rstanarm/articles/).

## Installation

Install the latest release from **CRAN**

```r
install.packages("rstanarm")
```

Instructions for installing the latest development version from **GitHub** can
be found in the **rstanarm** [Readme](https://github.com/stan-dev/rstanarm#installation).

## Getting help

* [File an issue on GitHub](https://github.com/stan-dev/rstanarm/issues)
* [Ask a question on Discourse](https://discourse.mc-stan.org/)

## Contributing

If you are interested in contributing to the development of **rstanarm** please
see the [Developer notes](https://mc-stan.org/rstanarm/dev-notes/index.html).
