# tests can be run using devtools::test() or manually by loading testthat 
# package and then running the code below possibly with options(mc.cores = 4).

library(rstanarm)

context("nlist")
test_that("nlist works", {
  nlist <- rstanarm:::nlist
  a <- 1
  b <- 2
  c <- 3
  val <- list(nlist(a, b, c), 
              nlist(a, b, c = "tornado"), 
              nlist(a = -1, b = -2, c))
  ans <- list(list(a = a, b = b, c = c), 
              list(a = a, b = b, c = "tornado"), 
              list(a = -1, b = -2, c = c))
  expect_identical(val, ans)
})

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
  expect_error(validate_parameter_value("a"), "should be NULL or numeric")
  expect_error(validate_parameter_value(NA), "should be NULL or numeric")
  expect_true(validate_parameter_value(NULL))
  expect_true(validate_parameter_value(.01))
  expect_true(validate_parameter_value(.Machine$double.xmax))
})

context("get_x, get_y, get_z")
test_that("get_x, get_y, get_z work", {
  fit <- suppressWarnings(stan_glm(mpg ~ wt, data = mtcars, iter = 5, chains = 1))
  x_ans <- cbind("(Intercept)" = 1, wt = mtcars$wt)
  y_ans <- mtcars$mpg
  expect_equivalent(x_ans, get_x(fit))
  expect_equivalent(y_ans, get_y(fit))
  expect_error(get_z(fit), "no applicable method")
  
  fit2 <- suppressWarnings(stan_glmer(mpg ~ wt + (1|cyl), data = mtcars, iter = 5, chains = 1))
  z_ans2 <- model.matrix(mpg ~ -1 + factor(cyl), data = mtcars)
  expect_equivalent(x_ans, get_x(fit2))
  expect_equivalent(y_ans, get_y(fit2))
  expect_equivalent(z_ans2, get_z(fit2))
  
  fit3 <- suppressWarnings(stan_glmer(mpg ~ wt + (1 + wt|cyl), data = mtcars, iter = 5, chains = 1))
  z_ans3 <- mat.or.vec(nr = nrow(mtcars), nc = 6)
  z_ans3[, c(1, 3, 5)] <- model.matrix(mpg ~ 0 + factor(cyl), data = mtcars)
  z_ans3[, c(2, 4, 6)] <- model.matrix(mpg ~ 0 + wt:factor(cyl), data = mtcars)
  expect_equivalent(x_ans, get_x(fit3))
  expect_equivalent(y_ans, get_y(fit3))
  expect_equivalent(z_ans3, get_z(fit3))
})

