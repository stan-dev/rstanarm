# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2015, 2016, 2017 Trustees of Columbia University
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

# tests can be run using devtools::test() or manually by loading testthat 
# package and then running the code below possibly with options(mc.cores = 4).

library(rstanarm)
SEED <- 123
CHAINS <- 2
ITER <- 400
threshold <- 0.21
REFRESH <- 0

source(test_path("helpers", "expect_stanreg.R"))
source(test_path("helpers", "SW.R"))

SW(fit <- stan_lm(mpg ~ ., data = mtcars, prior = R2(location = 0.75), 
                  chains = CHAINS, iter = ITER, seed = SEED, refresh = 0))

context("stan_lm")
test_that("stan_lm returns expected result for mtcars example", {
  # example using mtcars dataset
  expect_stanreg(fit)
  
  fit_sigma <- fit$stan_summary["sigma", "mean"]
  lm_sigma <- summary(lm(mpg ~ ., data = mtcars))$sigma
  expect_equal(fit_sigma, lm_sigma, tol = threshold)
})
test_that("stan_lm returns expected result for trees example", {
  # example using trees dataset
  SW(fit <- stan_lm(log(Volume) ~ log(Girth) + log(Height), data = trees, 
                  prior = R2(location = 0.9, what = "mean"), refresh = 0,
                  chains = CHAINS, iter = ITER, seed = SEED, adapt_delta = 0.999))
  expect_stanreg(fit)
  
  fit_sigma <- fit$stan_summary["sigma", "mean"]
  lm_sigma <- summary(lm(log(Volume) ~ log(Girth) + log(Height),data = trees))$sigma
  expect_equal(fit_sigma, lm_sigma, tol = threshold)
})

test_that("stan_lm doesn't break with less common priors", {
  # prior = NULL
  expect_output(fit <- stan_lm(mpg ~ -1 + ., data = mtcars, prior = NULL,
                iter = 10, chains = 1, seed = SEED), regexp = "SAMPLING")
  expect_stanreg(fit)
  
  # prior_intercept = normal()
  fit <- stan_lm(mpg ~ ., data = mtcars, refresh = 0,
                 prior = R2(0.75), prior_intercept = normal(),
                 iter = 10, chains = 1, seed = SEED)
  expect_stanreg(fit)
})

test_that("stan_lm doesn't break with vb algorithms", {
  fit <- stan_lm(mpg ~ ., data = mtcars, 
                 prior = R2(location = 0.75), refresh = 0,
                 algorithm = "meanfield", seed = SEED)
  expect_stanreg(fit)
  
  fit2 <- update(fit, algorithm = "fullrank")
  expect_stanreg(fit2)
})

test_that("stan_lm throws error if only intercept", {
  expect_error(stan_lm(mpg ~ 1, data = mtcars, prior = R2(location = 0.75)), 
               regexp = "not suitable for estimating a mean")
})

test_that("stan_lm throws error if 'location' is a vector", {
  expect_error(stan_lm(mpg ~ ., data = mtcars, prior = R2(location = c(0.25, 0.5))), 
               regexp = "only accepts a single value for 'location'")
})

test_that("stan_lm throws error if N < K", {
  # NOTE: remove this test once N < K is enabled
  expect_error(stan_lm(mpg ~ ., data = mtcars[1:5, ], prior = R2(0.75)), 
               regexp = "more predictors than data points is not yet enabled")
})

test_that("stan_lm throws error if glmer syntax used", {
  expect_error(stan_lm(mpg ~ wt + (1|cyl), data = mtcars, 
                       prior = R2(0.5, "mean")), 
               regexp = "model formula not allowed")
})

context("stan_aov")
test_that("stan_aov returns expected result for npk example", {
  SW(fit <- stan_aov(yield ~ block + N*P*K, data = npk, contrasts = "contr.poly",
                     prior = R2(0.5), chains = CHAINS, iter = ITER, seed = SEED, 
                     refresh = 0))
  expect_stanreg(fit)
  
  fit_sigma <- fit$stan_summary["sigma", "mean"]
  lm_sigma <- summary(lm(yield ~ block + N*P*K, data = npk, contrasts = "contr.poly"))$sigma
  expect_equal(fit_sigma, lm_sigma, tol = threshold)
  expect_output(print(fit), regexp = "stan_aov")
  expect_output(print(fit), regexp = "ANOVA-like table")
})

context("stan_biglm")
test_that("stan_biglm returns stanfit (not stanreg) object ", {
  ols <- lm(mpg ~ wt + qsec + am - 1,
            data = as.data.frame(scale(mtcars, scale = FALSE)))
  b <- coef(ols)
  R <- qr.R(ols$qr)
  SSR <- crossprod(ols$residuals)[1]
  N <- length(ols$fitted.values)
  xbar <- colMeans(mtcars[,c("wt", "qsec", "am")])
  y <- mtcars$mpg
  ybar <- mean(y)
  s_y <- sd(y)
  SW(post <- stan_biglm.fit(b, R, SSR, N, xbar, ybar, s_y, prior = R2(.75),
                           chains = 1, iter = 10, seed = SEED, refresh = 0))
  expect_s4_class(post, "stanfit")
})

test_that("loo/waic for stan_lm works", {
  source(test_path("helpers", "expect_equivalent_loo.R"))
  ll_fun <- rstanarm:::ll_fun
  expect_equivalent_loo(fit)
  expect_identical(ll_fun(fit), rstanarm:::.ll_gaussian_i)
})

context("posterior_predict (stan_lm)")
test_that("posterior_predict compatible with stan_lm", {
  source(test_path("helpers", "check_for_error.R"))
  source(test_path("helpers", "expect_linpred_equal.R"))
  
  check_for_error(fit)
  expect_linpred_equal(fit)
})
