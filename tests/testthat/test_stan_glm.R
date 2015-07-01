# tests can be run using devtools::test() or manually by loading testthat 
# package and then running the code below.

library(rstanarm)
set.seed(123)

context("stan_lm")
test_that("stan_lm returns expected result for mtcars example", {
  fit <- stan_lm(mpg ~ wt, data = mtcars, iter = 400, seed = 123)
  val <- cbind(coef(fit), se(fit))
  ans <- summary(lm(mpg ~ wt, data = mtcars))$coefficients[,1:2]
  diff <- abs(val - ans)
  testthat::expect_true(all(diff < 0.1))
})

context("stan_glm")
test_that("stan_glm returns expected result for glm poisson example", {
  # example from help("glm")
  counts <- c(18,17,15,20,10,20,25,13,12)
  outcome <- gl(3,1,9)
  treatment <- gl(3,3)
  fit <- stan_glm(counts ~ outcome + treatment, family = poisson(), 
                  iter = 400, seed = 123)
  val <- cbind(coef(fit), se(fit))
  ans <- summary(glm(counts ~ outcome + treatment, 
                     family = poisson()))$coefficients[,1:2]
  diff <- abs(val - ans)
  testthat::expect_true(all(diff < 0.1))
})

