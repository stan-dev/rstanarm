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
library(MASS)

SEED <- 123
ITER <- 100
CHAINS <- 2
CORES <- 1
REFRESH <- 0

threshold <- 0.03

source(test_path("helpers", "expect_stanreg.R"))

context("stan_polr")


f <- tobgp ~ agegp + alcgp
suppressWarnings(capture.output(
  fit1 <- stan_polr(f, data = esoph, method = "logistic", prior_PD = TRUE,
                    prior = R2(location = 0.4, what = "median"),
                    chains = CHAINS, iter = ITER, seed = SEED, refresh = 0),
  fit1vb <- stan_polr(f, data = esoph, method = "loglog",
                      prior = R2(location = 0.4, what = "median"),
                      seed = SEED, algorithm = "fullrank"),
  fit2 <- stan_polr(factor(tobgp == "30+") ~ agegp + alcgp, data = esoph, 
                    prior = R2(location = 0.4), method = "logistic", shape = 2, rate = 2,
                    chains = CHAINS, iter = ITER, seed = SEED, refresh = 0),
  fit2vb <- stan_polr(factor(tobgp == "30+") ~ agegp + alcgp, data = esoph, 
                      method = "probit", seed = SEED, algorithm = "fullrank",
                      prior = NULL, prior_counts = NULL), # test with NULL priors
  fit3 <- stan_polr(factor(tobgp == "30+") ~ agegp + alcgp,
                    data = esoph, prior = R2(location = 0.4),
                    shape = 2, rate = 2, chains = CHAINS, iter = ITER,
                    seed = SEED, refresh = 0)
))

test_that("stan_polr runs for esoph example", {
  expect_stanreg(fit1)
  expect_stanreg(fit2)
  expect_stanreg(fit1vb)
  expect_stanreg(fit2vb)
})

test_that("stan_polr throws error if formula excludes intercept", {
  expect_error(stan_polr(tobgp ~ 0 + agegp + alcgp, data = esoph, 
                         method = "loglog", prior = R2(0.4, "median")), 
               regexp = "an intercept is needed and assumed")
})

test_that("stan_polr throws error if shape,rate specified with >2 outcome levels", {
  expect_error(
    stan_polr(f, data = esoph, method = "loglog", prior = R2(0.4, "median"), shape = 2), 
    "'shape' must be NULL when there are more than 2 outcome categories"
  )
  expect_error(
    stan_polr(f, data = esoph, method = "loglog", prior = R2(0.4, "median"), rate = 2), 
    "'rate' must be NULL when there are more than 2 outcome categories"
  )
})

test_that("gumbel functions ok", {
  # formulas are correct
  # just test a few cases so they're flagged if anything changes by accident
  # maybe should compare to corresponding functions in ordinal package?
  expect_equal(rstanarm:::dgumbel(0), 0.3678794, tol = 0.00001)
  expect_equal(rstanarm:::qgumbel(0), -Inf)
  expect_equal(rstanarm:::qgumbel(0.5), 0.3665129, tol = 0.00001)
  expect_equal(rstanarm:::pgumbel(0.3665129), 0.5, tol = 0.00001)
  expect_equal(rstanarm:::qgumbel(1), Inf)
})

test_that("loo/waic for stan_polr works", {
  source(test_path("helpers", "expect_equivalent_loo.R"))
  ll_fun <- rstanarm:::ll_fun
  expect_equivalent_loo(fit1)
  expect_identical(ll_fun(fit1), rstanarm:::.ll_polr_i)
  
  expect_equivalent_loo(fit2)
  expect_identical(ll_fun(fit2), rstanarm:::.ll_polr_i)
  
  expect_equivalent_loo(fit3)
  expect_identical(ll_fun(fit3), rstanarm:::.ll_polr_i)
})

context("posterior_predict (stan_polr)")
test_that("compatible with stan_polr", {
  source(test_path("helpers", "check_for_error.R"))
  check_for_error(fit1)
  check_for_error(fit2)
  check_for_error(fit3)
})
