# tests can be run using devtools::test() or manually by loading testthat 
# package and then running the code below possibly with options(mc.cores = 4).

library(rstanarm)
SEED <- 123
ITER <- 10
CHAINS <- 2
CORES <- 1

fit <- example_model
fit2 <- suppressWarnings(stan_glm(mpg ~ wt, data = mtcars, iter = ITER, 
                                  chains = CHAINS, cores = CORES, seed = SEED))

context("pp_check")
test_that("pp_check doesn't throw bad errors", {
  expect_silent(p <- pp_check(fit, check = "dist", overlay = TRUE, size = 2))
  expect_silent(p <- pp_check(fit, check = "resid"))
  expect_silent(p <- pp_check(fit2, check = "resid", fill = "red", bins = 15))
  expect_silent(p <- pp_check(fit, check = "scatter"))
  expect_silent(p <- pp_check(fit2, check = "scatter", color = "purple"))
  expect_is(p, "ggplot")
  for (j in 1:2) {
    expect_silent(p <- pp_check(fit, check = "dist", overlay = FALSE, nreps = j))
    expect_silent(p <- pp_check(fit, check = "dist", overlay = TRUE, nreps = j))
    expect_silent(p <- pp_check(fit, check = "resid", nreps = j))
    expect_silent(p <- pp_check(fit2, check = "resid", nreps = j))
    expect_silent(p <- pp_check(fit, check = "scat", nreps = j))
    expect_silent(p <- pp_check(fit2, check = "scat", nreps = j))
  }
  expect_silent(p <- pp_check(fit, check = "test"))
  expect_silent(p <- pp_check(fit, check = "test", test = "sd"))
  expect_silent(p <- pp_check(fit, check = "test", test = c("mean","sd")))
  expect_is(p, "ggplot")
})

test_that("pp_check ok for vb", {
  fit3 <- update(fit2, algorithm = "meanfield", iter = 10000)
  expect_silent(p <- pp_check(fit3))
  expect_silent(p <- pp_check(fit3, check = "resid"))
  expect_silent(p <- pp_check(fit3, check = "scat"))
  expect_silent(p <- pp_check(fit3, check = "test"))
})

test_that("pp_check throws appropriate errors", {
  expect_error(p <- pp_check(fit, check = "test", test = "10982pqmeaw"), 
               regexp = "not found")
  expect_error(p <- pp_check(fit, check = "test", test = c("mean", "sd", "var")), 
               regexp = "length 1 or 2")
  
  fito <- stan_glm(mpg ~ wt, data = mtcars, algorithm = "optimizing")
  expect_error(pp_check(fito), regexp = "algorithm")
  expect_error(pp_check(rnorm(10)), regexp = "not a stanreg object")
})

test_that("pp_check throws appropriate warnings", {
  expect_warning(p <- pp_check(fit, check = "test", nreps = 1), 
                 regexp = "'nreps' is ignored")
})
