# tests can be run using devtools::test() or manually by loading testthat 
# package and then running the code below possibly with options(mc.cores = 4).

library(rstanarm)
SEED <- 123
threshold <- 0.51

f1 <- function(x) cbind(coef(x), se(x))
f2 <- function(x) summary(x)$coefficients[,1:2]

context("stan_lm")
test_that("stan_lm returns expected result for mtcars example", {
  # example using mtcars dataset
  fit <- stan_lm(mpg ~ ., data = mtcars, prior = R2(location = 0.75), 
                 chains = 2, iter = 400, seed = SEED)
  fit_sigma <- get_posterior_mean(fit$stanfit)["sigma[1]",3]
  lm_sigma <- summary(lm(mpg ~ ., data = mtcars))$sigma
  expect_equal(fit_sigma, lm_sigma, tol = threshold)
})
context("stan_lm")
test_that("stan_lm returns expected result for trees example", {
  # example using trees dataset
  fit1 <- stan_lm(log(Volume) ~ log(Girth) + log(Height), data = trees, 
                  prior = R2(location = 0.9, what = "mean"), 
                  chains = 2, iter = 400, seed = SEED)
  ans <- lm(log(Volume) ~ log(Girth) + log(Height),data = trees)
  expect_equal(coef(fit1), coef(ans), tol = threshold)
})
