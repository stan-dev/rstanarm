# tests can be run using devtools::test() or manually by loading testthat 
# package and then running the code below.

library(rstanarm)
set.seed(123)

context("stan_lm")
test_that("stan_lm returns expected result for mtcars example", {
  # example using mtcars dataset
  fit <- stan_lm(mpg ~ wt, data = mtcars, iter = 400, seed = 123)
  val <- cbind(coef(fit), se(fit))
  ans <- summary(lm(mpg ~ wt, data = mtcars))$coefficients[,1:2]
  diff <- abs(val - ans)
  expect_true(all(diff < 0.1))
})
test_that("gaussian(link = 'log') returns expected result for trees example", {
  # example using trees dataset
  fit1 <- stan_glm(Volume ~ log(Girth) + log(Height), data = trees, 
                 family = gaussian(link = "log"), iter = 400, seed = 123)
  fit2 <- stan_lm(log(Volume) ~ log(Girth) + log(Height), data = trees, 
                  iter = 400, seed = 123)
  val1 <- cbind(coef(fit1), se(fit1)); val2 <- cbind(coef(fit2), se(fit2))
  ans <- summary(lm(log(Volume) ~ log(Girth) + log(Height),
                    data = trees))$coefficients[,1:2]
  diff1 <- abs(val1 - ans); diff2 <- abs(val2 - ans)
  expect_true(all(diff1 < 0.1 & diff2 < 0.1))
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
  expect_true(all(diff < 0.1))
})

context("stan_glm")
test_that("stan_glm returns expected result for cars example", {
  # example using cars dataset
  fit <- stan_glm(log(dist) ~ log(speed), data = cars, 
                  family = gaussian(link = "identity"), 
                  iter = 400, seed = 123)
  val <- cbind(coef(fit), se(fit))
  ans <- summary(glm(log(dist) ~ log(speed), data = cars, 
                     family = gaussian(link = "identity")))$coefficients[,1:2]
  diff <- abs(val - ans)
  testthat::expect_true(all(diff < 0.1))
})

