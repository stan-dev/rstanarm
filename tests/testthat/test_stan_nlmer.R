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

# tests can be run using devtools::test() or manually by loading testthat 
# package and then running the code below possibly with options(mc.cores = 4).

library(rstanarm)
library(lme4)

SEED <- 12345
ITER <- 100
CHAINS <- 2
CORES <- 2
REFRESH <- 0

threshold <- 0.05

source(test_path("helpers", "expect_stanreg.R"))

context("stan_nlmer")

data("Orange", package = "datasets")
Orange$circumference <- Orange$circumference / 100
Orange$age <- Orange$age / 100
fit <- stan_nlmer(circumference ~ SSlogis(age, Asym, xmid, scal) ~ Asym|Tree, 
                  data = Orange, prior = NULL, cores = CORES, init_r = 1,
                  chains = CHAINS, seed = SEED, refresh = 0, QR = TRUE)
startvec <- c(Asym = 200, xmid = 725, scal = 350) / 100
ml <- nlmer(circumference ~ SSlogis(age, Asym, xmid, scal) ~ Asym|Tree,
            data = Orange, start = startvec)


test_that("stan_nlmer runs for Orange example", {
  expect_stanreg(fit)
})

test_that("stan_nlmer is similar to nlmer on Orange example", {
  expect_equal(fixef(ml), fixef(fit), tol = threshold)
})

test_that("stan_nlmer throws error if formula includes an unknown function", {
  expect_error(stan_nlmer(circumference ~ SSfoo(age, Asym, xmid, scal) ~ Asym|Tree, 
                          data = Orange),
               regexp = "self-starting nonlinear function")
})

test_that("loo/waic for stan_nlmer works", {
  source(test_path("helpers", "expect_equivalent_loo.R"))
  # expect_equivalent_loo(fit)
})

context("posterior_predict (stan_nlmer)")
test_that("compatible with stan_nlmer", {
  source(test_path("helpers", "check_for_error.R"))
  check_for_error(fit)
})
