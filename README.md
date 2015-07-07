rstanarm
========
This is an R package that emulates other R model-fitting functions but uses (R)Stan for the back-end estimation. The primary target audience is people who would be open to Bayesian inference if using
Bayesian software were easy but would use frequentist software otherwise. That said, this R package
often uses posterior means as point estimates.

Rules:
  1. The stan* (e.g. `stan_glm`) wrapper function should take (almost) all of the same arguments as the function it emulates (e.g. `glm`) but ignoring unnecessary arguments is fine
  2. After that, you can add additional arguments that often pertain to the priors used by Stan
  3. The ... is always passed to `stan` so that the user can specify the number of chains, etc.
  4. The .stan files go in the exec/ subdirectory; try to make them as abstract, numerically stable, fast, etc. as possible at the expense of readability if necessary.
  5. The R/stanmodels.R file just creates a bunch of empty internal `stanfit` objects that are implicitly "updated" by the user
  6. The order of the returned parameters should be intercept(s) (if any), regression coefficients, nuisance parameters, other stuff. Use the `pars` argument to `stan` if necessary to achieve this ordering.
  7. Do not store any quantities whose number is proportional to the sample size. If necessary, generate such things inside R so that we can aggregate over them incrementally (e.g. importance sampling weighted LOO-CV).