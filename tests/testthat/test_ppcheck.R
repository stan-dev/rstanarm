# tests can be run using devtools::test() or manually by loading testthat 
# package and then running the code below possibly with options(mc.cores = 4).

library(rstanarm)
SEED <- 123
ITER <- 10
CHAINS <- 2
CORES <- 1

fit <- example_model

context("ppcheck")
test_that("ppcheck doesn't throw bad errors", {
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
  Ty <- function(x) quantile(x, probs = 0.9)
  expect_silent(p <- ppcheck(fit, check = "test", test = "Ty"))
})

test_that("ppcheck throws appropriate errors", {
  expect_error(p <- ppcheck(fit, check = "test", test = "10982pqmeaw"), 
               regexp = "not found")
  expect_error(p <- ppcheck(fit, check = "test", test = c("mean", "sd", "var")), 
               regexp = "length 1 or 2")
  
  fito <- stan_glm(mpg ~ wt, data = mtcars, algorithm = "optimizing")
  expect_error(ppcheck(fito), regexp = "algorithm")
  expect_error(ppcheck(rnorm(10)), regexp = "not a stanreg object")
})

test_that("ppcheck throws appropriate warnings", {
  expect_warning(p <- ppcheck(fit, check = "test", nreps = 1), 
                 regexp = "'nreps' is ignored")
})