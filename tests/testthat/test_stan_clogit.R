# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2017 Trustees of Columbia University
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

# this mostly goes through same code as a logit model so only testing the unique stuff

suppressPackageStartupMessages(library(rstanarm))

SEED <- 123
ITER <- 100
CHAINS <- 2
CORES <- 1
REFRESH <- 0

threshold <- 0.03

context("stan_clogit")

SW(fit <- stan_clogit(case ~ spontaneous + induced, strata = stratum, prior = NULL,
                   data = infert[order(infert$stratum), ], 
                   QR = TRUE, init_r = 0.5,
                   chains = CHAINS, iter = ITER, seed = SEED, refresh = 0))

test_that("stan_clogit is similar to survival::clogit", {
  ref_vals <- c(spontaneous = 1.985876, induced = 1.409012)
  # Account for RNG change in new Stan
  if (utils::packageVersion("StanHeaders") >= "2.36") {
    ref_vals <- c(spontaneous = 2.062676, induced = 1.360712)
  }
  expect_equal(ref_vals, coef(fit), tol = threshold)
})

test_that("stan_clogit runs for infert example", {
  expect_stanreg(fit)
})

test_that("stan_clogit works when y is a factor", {
  d <- infert[order(infert$stratum), ]
  d$case <- factor(d$case, labels = c("A", "B"))
  SW(fit_factor <- stan_clogit(case ~ spontaneous + induced, strata = stratum, prior = NULL,
                        data = infert[order(infert$stratum), ], 
                        QR = TRUE, init_r = 0.5,
                        chains = CHAINS, iter = ITER, seed = SEED, refresh = 0))
  expect_equal(coef(fit_factor), coef(fit))
})

test_that("stan_clogit throws error if data are not sorted", {
  expect_error(update(fit, data = infert), 
               regexp = "Data must be sorted")
})

test_that("loo/waic for stan_clogit works", {
  ll_fun <- rstanarm:::ll_fun
  expect_equivalent_loo(fit)
  expect_identical(ll_fun(fit), rstanarm:::.ll_clogit_i)
})

context("posterior_predict (stan_clogit)")
test_that("compatible with stan_clogit", {
  PPD1 <- posterior_predict(fit)
  PPD2 <- posterior_predict(fit, newdata = infert) # order irrelevant
  expect_identical(rowSums(PPD1), rowSums(PPD2))
  expect_equal(rowSums(PPD1), round(rowSums(
               posterior_linpred(fit, newdata = infert, transform = TRUE))))
})
