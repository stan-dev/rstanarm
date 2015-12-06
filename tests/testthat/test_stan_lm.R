# tests can be run using devtools::test() or manually by loading testthat 
# package and then running the code below possibly with options(mc.cores = 4).

library(rstanarm)
SEED <- 123
CHAINS <- 2
ITER <- 400
threshold <- 0.21

context("stan_lm")
test_that("stan_lm returns expected result for mtcars example", {
  # example using mtcars dataset
  fit <- stan_lm(mpg ~ ., data = mtcars, prior = R2(location = 0.75), 
                 chains = CHAINS, iter = ITER, seed = SEED)
  fit_sigma <- fit$stan_summary["sigma", "mean"]
  lm_sigma <- summary(lm(mpg ~ ., data = mtcars))$sigma
  expect_equal(fit_sigma, lm_sigma, tol = threshold)
})
test_that("stan_lm returns expected result for trees example", {
  # example using trees dataset
  fit <- stan_lm(log(Volume) ~ log(Girth) + log(Height), data = trees, 
                  prior = R2(location = 0.9, what = "mean"), 
                  chains = CHAINS, iter = ITER, seed = SEED, adapt_delta = 0.999)
  fit_sigma <- fit$stan_summary["sigma", "mean"]
  lm_sigma <- summary(lm(log(Volume) ~ log(Girth) + log(Height),data = trees))$sigma
  expect_equal(fit_sigma, lm_sigma, tol = threshold)
})

context("stan_aov")
test_that("stan_aov returns expected result for npk example", {
  fit <- stan_aov(yield ~ block + N*P*K, data = npk, contrasts = "contr.poly",
           prior = R2(0.5), chains = CHAINS, iter = ITER, seed = SEED)
  fit_sigma <- fit$stan_summary["sigma", "mean"]
  lm_sigma <- summary(lm(yield ~ block + N*P*K, data = npk, contrasts = "contr.poly"))$sigma
  expect_equal(fit_sigma, lm_sigma, tol = threshold)
})


