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


context("pp_check")


library(rstanarm)
SEED <- 123
set.seed(SEED)
ITER <- 10
CHAINS <- 2
REFRESH <- 0

source(test_path("helpers", "SW.R"))
source(test_path("helpers", "expect_gg.R"))

fit <- example_model
SW(fit2 <- stan_glm(mpg ~ wt + am, data = mtcars, iter = ITER, chains = CHAINS,
                    seed = SEED, refresh = 0))


patt <- "rootogram|_bars|vs_x|grouped$|_data$"
ppc_funs_not_grouped <- bayesplot::available_ppc(patt, invert = TRUE)
ppc_funs_grouped <- bayesplot::available_ppc("vs_x|grouped")
ppc_funs_discrete <- bayesplot::available_ppc("rootogram|_bars")


test_that("pp_check.stanreg creates ggplot object", {
  exclude <- c("ppc_bars", 
               "ppc_loo_pit", 
               "ppc_loo_pit_overlay", 
               "ppc_loo_pit_qq",
               "ppc_loo_intervals", 
               "ppc_loo_ribbon", 
               "ppc_rootogram", 
               "ppc_error_binned")
  for (f in ppc_funs_not_grouped) for (j in 1:2) {
    if (!f %in% exclude)
      expect_gg(suppressWarnings(pp_check(fit, plotfun = f, nreps = j)), 
                info = f)
  }
})

test_that("pp_check.stanreg creates ggplot object for grouped functions", {
  for (f in setdiff(ppc_funs_grouped, ppc_funs_discrete)) for (j in 1:2) {
    expect_gg(suppressWarnings(pp_check(fit2, plotfun = f, nreps = j, group = "am", x = "wt")), 
              info = f)
  }
})

test_that("pp_check.stanreg creates ggplot object for count & ordinal outcomes", {
  d <- data.frame(
    counts = c(18,17,15,20,10,20,25,13,12),
    outcome = gl(3,1,9),
    treatment = gl(3,3)
  )
  SW(fit3 <- stan_glm(counts ~ outcome + treatment, data = d, 
                   family = poisson(link="log"),
                   iter = ITER, chains = CHAINS,
                   seed = SEED, refresh = 0))
  expect_gg(pp_check(fit3, plotfun = "rootogram"))
  
  SW(fit4 <- stan_polr(tobgp ~ agegp, data = esoph, method = "probit",
                       prior = R2(0.2, "mean"), init_r = 0.1, 
                       iter = ITER, chains = CHAINS,
                       seed = SEED, refresh = 0))
  expect_gg(pp_check(fit4, plotfun = "bars"))
  expect_gg(pp_check(fit4, plotfun = "bars_grouped", group = "agegp"))
})


test_that("pp_check ok for vb", {
  SW(fit3 <- stan_glm(mpg ~ wt, data = mtcars, algorithm = "meanfield", 
                      seed = SEED, iter = 10000))
  expect_gg(pp_check(fit3))
  expect_gg(pp_check(fit3, plotfun = "error_hist"))
})

# test_that("pp_check binned residual plot works for factors", {
#   ir2 <- iris[-c(1:50), ]
#   ir2$Species <- factor(ir2$Species)
#   SW(fit3 <- stan_glm(Species ~ Petal.Length + Petal.Width + Sepal.Length + Sepal.Width,
#                       data=ir2, family = "binomial", iter = ITER, chains = CHAINS,
#                       seed = SEED, refresh = 0))
#   expect_gg(pp_check(fit3, plotfun = "error_binned"))
# })


# test errors --------------------------------------------------------------
test_that("pp_check throws error if 'stat' arg is bad", {
  expect_error(pp_check(fit, plotfun = "stat", stat = "10982pqmeaw"),
               regexp = "not found")
})
test_that("pp_check throws error if plotfun not found", {
  expect_error(pp_check(fit, plotfun = "9999"), 
               "not a valid PPC function name")
  expect_error(pp_check(fit, plotfun = "mcmc_hist"), 
               "use the 'plot' method")
})
test_that("pp_check throws error if 'group' variable not found", {
  expect_error(pp_check(fit, plotfun = "stat_grouped", group = "herd2"), 
               "not found in model frame")
})
test_that("pp_check throws error for optimizing", {
  SW(fito <- stan_glm(mpg ~ wt, data = mtcars, algorithm = "optimizing", 
                      seed = SEED, refresh = 0))
  expect_error(pp_check(fito), regexp = "algorithm")
})

# test warnings ----------------------------------------------------------
test_that("pp_check throws warning if 'nreps' ignored ", {
  expect_warning(pp_check(fit, plotfun = "stat", nreps = 1),
                 regexp = "'nreps' is ignored")
})
test_that("pp_check throws warning if 'group' or 'x' ignored", {
  expect_warning(pp_check(fit, plotfun = "stat_2d", stat = c("mean", "sd"), group = "herd"),
                 regexp = "ignored: group")
  expect_warning(pp_check(fit, plotfun = "scatter", nreps = 3, group = "herd"),
                 regexp = "ignored: group")
  expect_warning(pp_check(fit, plotfun = "error_hist", x = "herd"),
                 regexp = "ignored: x")
})


# helpers -----------------------------------------------------------------
test_that(".ignore_nreps and .set_nreps work", {
  ignore_nreps <- rstanarm:::.ignore_nreps
  set_nreps <- rstanarm:::.set_nreps
  expect_warning(ignore_nreps(10), "'nreps' is ignored")
  expect_silent(ignore_nreps(NULL))
  
  expect_warning(r <- set_nreps(10, "ppc_stat"), "'nreps' is ignored")
  expect_null(r)
  expect_equal(set_nreps(10, "ppc_hist"), 10)
})

test_that("y coerced to numeric (attributes dropped)", {
  d <- mtcars
  attr(d$mpg, "test") <- "something"
  SW(fit3 <- update(fit2, data = d))
  expect_equal(attr(get_y(fit3), "test"), "something")
  expect_gg(pp_check(fit3, nreps = 3))
})
