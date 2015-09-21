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
