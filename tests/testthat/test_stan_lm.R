# tests can be run using devtools::test() or manually by loading testthat 
# package and then running the code below possibly with options(mc.cores = 4).

library(rstanarm)
set.seed(123)

threshold <- 0.51

f1 <- function(x) cbind(coef(x), se(x))
f2 <- function(x) summary(x)$coefficients[,1:2]

context("stan_lm")
test_that("stan_lm returns expected result for mtcars example", {
  # example using mtcars dataset
  fit <- stan_lm(mpg ~ ., data = mtcars, prior = R2(location = 0.75), 
                 iter = 400, seed = 123)
  fit_sigma <- get_posterior_mean(fit$stanfit)["sigma[1]",5]
  lm_sigma <- summary(lm(mpg ~ ., data = mtcars))$sigma
  diff <- abs(lm_sigma - fit_sigma)
  expect_true(all(diff < threshold))
})
context("stan_lm")
test_that("stan_lm returns expected result for trees example", {
  # example using trees dataset
  fit1 <- stan_lm(log(Volume) ~ log(Girth) + log(Height), data = trees, 
                  prior = R2(location = 0.9, what = "mean"), 
                  iter = 400, seed = 123, 
                  control = list(adapt_delta = 0.95, max_treedepth = 12))
  val1 <- f1(fit1)
  ans <- f2(lm(log(Volume) ~ log(Girth) + log(Height),data = trees))
  diff1 <- abs(val1 - ans)
  expect_true(all(diff1 < threshold))
})
