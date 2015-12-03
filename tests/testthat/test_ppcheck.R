# tests can be run using devtools::test() or manually by loading testthat 
# package and then running the code below possibly with options(mc.cores = 4).

library(rstanarm)
SEED <- 123
ITER <- 10
CHAINS <- 2
CORES <- 1

context("ppcheck")
test_that("ppcheck doesn't throw bad errors", {
  fit <- example_model
  
  expect_silent(p <- ppcheck(fit, check = "dist", overlay = TRUE))
  expect_silent(p <- ppcheck(fit, check = "resid"))
  for (j in 1:2) {
    expect_silent(p <- ppcheck(fit, check = "dist", overlay = FALSE, nreps = j))
    expect_silent(p <- ppcheck(fit, check = "dist", overlay = TRUE, nreps = j))
    expect_silent(p <- ppcheck(fit, check = "resid", nreps = j))
  }
  expect_silent(p <- ppcheck(fit, check = "test"))
  expect_silent(p <- ppcheck(fit, check = "test", test = "sd"))
  expect_silent(p <- ppcheck(fit, check = "test", test = c("mean","sd")))
})

test_that("ppcheck throws appropriate warnings", {
  Ty <- function(x) quantile(x, probs = 0.9)
  if (utils::packageVersion("ggplot2") < "1.0.1.9003") {
    expect_silent(p <- ppcheck(fit, check = "test", test = "Ty"))
  } else {
    expect_warning(p <- ppcheck(fit, check = "test", test = "Ty"), 
                   regexp = "`show_guide` has been deprecated")
  }
  expect_warning(p <- ppcheck(fit, check = "test", nreps = 1), 
                 regexp = "'nreps' is ignored")
})

test_that("ppcheck throws appropriate errors", {
  fit <- stan_glm(mpg ~ wt, data = mtcars, algorithm = "optimizing")
  expect_error(ppcheck(fit), regexp = "algorithm")
  expect_error(ppcheck(rnorm(10)), regexp = "not a stanreg object")
})
