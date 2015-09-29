# tests can be run using devtools::test() or manually by loading testthat 
# package and then running the code below possibly with options(mc.cores = 4).

library(rstanarm)

context("ORifNULL")
test_that("%ORifNULL% works", {
  `%ORifNULL%` <- rstanarm:::`%ORifNULL%`
  a <- list(NULL, NA, NaN, 1, "a", FALSE, mat.or.vec(5,5))
  b <- 1
  ans <- c(1, a[-1])
  for (j in seq_along(a)) {
    expect_identical(a[[j]] %ORifNULL% b, ans[[j]])
  }
})

context("maybe_broadcast")
test_that("maybe_broadcast works", {
  n <- 5
  x <- list(numeric(0), NULL, 1, c(1,1))
  ans <- list(rep(0,n), rep(0,n), rep(1,n), c(1,1))
  for (j in seq_along(ans)) {
    expect_equal(rstanarm:::maybe_broadcast(x[[j]], n), ans[[j]])  
  }
})

context("set_prior_scale")
test_that("set_prior_scale works", {
  set_prior_scale <- rstanarm:::set_prior_scale
  expect_error(set_prior_scale("a", "b", "c"))
  expect_error(set_prior_scale(1, 1, 1))
  expect_equal(set_prior_scale(NULL, 1, "a"), 1)
  expect_equal(set_prior_scale(NULL, 1, "probit"), dnorm(0) / dlogis(0))
  expect_equal(set_prior_scale(2, 1, "a"), 2)
  expect_equal(set_prior_scale(2, 1, "probit"), 2 * dnorm(0) / dlogis(0))
})

context("validate_parameter_value")
test_that("validate_parameter_value works", {
  validate_parameter_value <- rstanarm:::validate_parameter_value
  expect_error(validate_parameter_value(-1), "should be positive")
  expect_error(validate_parameter_value(0), "should be positive")
  expect_true(validate_parameter_value(NULL))
  expect_true(validate_parameter_value(.01))
  expect_true(validate_parameter_value(.Machine$double.xmax))
})

context("get_x, get_y, get_z")
test_that("get_x, get_y, get_z work properly", {
  fit <- suppressWarnings(stan_glm(mpg ~ wt, data = mtcars, iter = 5, chains = 1))
  x_ans <- cbind("(Intercept)" = 1, wt = mtcars$wt)
  y_ans <- mtcars$mpg
  expect_equal(x_ans, get_x(fit), check.attributes = FALSE)
  expect_equal(y_ans, get_y(fit), check.attributes = FALSE)
  expect_error(get_z(fit), "no applicable method")
  
  fit2 <- suppressWarnings(stan_glmer(mpg ~ wt + (1|cyl), data = mtcars, iter = 5, chains = 1))
  z_ans2 <- model.matrix(mpg ~ -1 + factor(cyl), data = mtcars)
  expect_equal(x_ans, get_x(fit2), check.attributes = FALSE)
  expect_equal(y_ans, get_y(fit2), check.attributes = FALSE)
  expect_equal(z_ans2, get_z(fit2), check.attributes = FALSE)
  

  fit3 <- suppressWarnings(stan_glmer(mpg ~ wt + (1 + wt|cyl), data = mtcars, iter = 5, chains = 1))
  z_ans3 <- mat.or.vec(nr = nrow(mtcars), nc = 6)
  z_ans3[, c(1, 3, 5)] <- model.matrix(mpg ~ 0 + factor(cyl), data = mtcars)
  z_ans3[, c(2, 4, 6)] <- model.matrix(mpg ~ 0 + wt:factor(cyl), data = mtcars)
  expect_equal(x_ans, get_x(fit3), check.attributes = FALSE)
  expect_equal(y_ans, get_y(fit3), check.attributes = FALSE)
  expect_equal(z_ans3, get_z(fit3), check.attributes = FALSE)
})

