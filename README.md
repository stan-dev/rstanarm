rstanarm
========
This is an R package that emulates other R model-fitting functions but uses (R)Stan for the back-end estimation.

Rules:
  1. The stan* (e.g. stanglm) wrapper function should take all the same arguments as the function it emulates (e.g. glm)
  2. After that, you can add additional arguments that often pertain to the priors used by Stan
  3. The ... is always passed to stan() so that the user can specify the number of chains, etc.
  4. The .stan files go in the exec/ subdirectory; try to make them as abstract as possible
  5. The R/stanmodels.R file just creates a bunch of empty stanfit objects that are implicitly "updated" by the user
