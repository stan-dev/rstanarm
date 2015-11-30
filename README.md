[![Build Status](https://travis-ci.org/stan-dev/rstanarm.svg?branch=master)](https://travis-ci.org/stan-dev/rstanarm)[![Coverage Status](https://img.shields.io/codecov/c/github/stan-dev/rstanarm/master.svg)](https://codecov.io/github/stan-dev/rstanarm?branch=master) [![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/rstanarm)](http://cran.r-project.org/package=rstanarm)


rstanarm
========
This is an R package that emulates other R model-fitting functions but uses (R)Stan for the back-end estimation. The primary target audience is people who would be open to Bayesian inference if using
Bayesian software were easy but would use frequentist software otherwise. That said, this R package
often uses posterior medians as point estimates.

Rules:
  1. The stan\_* (e.g. `stan_glm`) wrapper function should take (almost) all of the same arguments as the function it emulates (e.g. `glm`) but unnecessary arguments can be dropped. The likelihood of the data must be the same in both cases.
  2. The ... is always passed to the estimation function (e.g. `sampling`) so that the user can specify the number of chains, etc. Put the ... in the argument list after the arguments from the emulated function but before any rstanarm-specific arguments like `priors` (unless the ... comes earlier in the emulated function's argument list as in `polr`)
  3. After that, you can add additional arguments that often pertain to the priors used by Stan
  4. The .stan files go in the exec/ subdirectory; try to make them as abstract, numerically stable, fast, etc. as possible at the expense of readability if necessary.
  5. The R/stanmodels.R file just creates a bunch of `stanmodel` objects that are passed to `sampling`, `optimizing`, etc.
  6. The order of the returned parameters should be intercept(s) (if any), regression coefficients, nuisance parameters, generated quantities, anything else. Use the `pars` argument to `stan` if necessary to achieve this ordering.
  7. Do not store any quantities whose number is proportional to the sample size. If necessary, generate such things inside R so that we can aggregate over them incrementally (e.g. importance sampling weighted LOO-CV).
  8. Do generate the (typically scalar, unless outcome is multivariate) sample average posterior predictive distribution of the outcome and call it `mean_PPD`. This is useful for checking.
  9. Utilize user-defined Stan functions wherever possible. Verify they are correct in tests/testthat/test\_stan\_functions.R. These tests should mostly utilize `algorithm="optimizing"` and `NULL` priors so that the results are within numerical tolerance of the function being emulated.
