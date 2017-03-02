# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2015, 2016 Trustees of Columbia University
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
set.seed(SEED)
ITER <- 10
CHAINS <- 2
REFRESH <- 0

SW <- function(expr) capture.output(suppressWarnings(expr))

fit <- example_model
SW(fit2 <- stan_glm(mpg ~ wt + am, data = mtcars, iter = ITER, chains = CHAINS,
                    seed = SEED, refresh = REFRESH))

expect_gg <- function(x, info = NULL, label = NULL) {
  testthat::expect_is(x, "ggplot", info = info, label = label)
}


context("pp_check")


# test deprecated stuff -----------------------------------------------
test_that("pp_check with deprecated 'check' arg works", {
  expect_warning(p1 <- pp_check(fit, check = "dist"), 
                 "Argument 'check' is deprecated")
  expect_warning(p2 <- pp_check(fit2, check = "dist", overlay = FALSE))
  expect_warning(p3 <- pp_check(fit, check = "resid"))
  expect_warning(p4 <- pp_check(fit2, check = "resid", binwidth = .5))
  expect_warning(p5 <- pp_check(fit, check = "scatter"))
  expect_warning(p6 <- pp_check(fit2, check = "scatter"))
  expect_gg(p1)
  expect_gg(p2)
  expect_gg(p3)
  expect_gg(p4)
  expect_gg(p5)
  expect_gg(p6)
})


# test new pp_check.stanreg  ----------------------------------------------
all_ppc_funs <- grep("^ppc_", getNamespaceExports("bayesplot"), value = TRUE)
ppc_funs_not_grouped <- grep("vs_x|_grouped$", all_ppc_funs, value = TRUE, invert = TRUE)
ppc_funs_grouped <- grep("vs_x|_grouped$", all_ppc_funs, value = TRUE)

test_that("pp_check.stanreg creates ggplot object", {
  for (f in ppc_funs_not_grouped) for (j in 1:2) {
    expect_gg(suppressWarnings(pp_check(fit, plotfun = f, nreps = j)), 
              info = f)
  }
})

test_that("pp_check.stanreg creates ggplot object for grouped functions", {
  for (f in ppc_funs_grouped) for (j in 1:2) {
    expect_gg(suppressWarnings(pp_check(fit2, plotfun = f, nreps = j, group = "am", x = "wt")), 
              info = f)
  }
})


test_that("pp_check ok for vb", {
  SW(fit3 <- stan_glm(mpg ~ wt, data = mtcars, algorithm = "meanfield", 
                      seed = SEED, iter = 10000))
  expect_gg(pp_check(fit3))
  expect_gg(pp_check(fit3, plotfun = "error_hist"))
})

test_that("pp_check binned residual plot works for factors", {
  ir2 <- iris[-c(1:50), ]
  ir2$Species <- factor(ir2$Species)
  SW(fit3 <- stan_glm(Species ~ Petal.Length + Petal.Width + Sepal.Length + Sepal.Width,
                      data=ir2, family = "binomial", iter = ITER, chains = CHAINS,
                      seed = SEED, refresh = REFRESH))
  expect_gg(pp_check(fit3, plotfun = "error_binned"))
})


# test errors --------------------------------------------------------------
test_that("pp_check throws error if 'test' arg is bad", {
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
  SW(fito <- stan_glm(mpg ~ wt, data = mtcars, algorithm = "optimizing", seed = SEED))
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
