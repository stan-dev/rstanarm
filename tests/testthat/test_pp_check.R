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

SW <- suppressWarnings

fit <- example_model
fit2 <- SW(stan_glm(mpg ~ wt, data = mtcars, iter = ITER, chains = CHAINS,
                    seed = SEED, refresh = REFRESH))

expect_gg <- function(x) expect_s3_class(x, "ggplot")


context("pp_check")


# test ggplot object creation -----------------------------------------------
test_that("pp_check creates ggplot objects when it should", {
  expect_gg(pp_check(fit, check = "dist", overlay = TRUE, size = 2))
  expect_gg(pp_check(fit2, check = "dist", overlay = FALSE))

  expect_gg(pp_check(fit, check = "resid"))
  expect_gg(pp_check(fit2, check = "resid", binwidth = .5))

  expect_gg(pp_check(fit, check = "scatter"))
  expect_gg(pp_check(fit2, check = "scatter"))

  for (j in 1:2) {
    expect_gg(pp_check(fit, check = "dist", overlay = FALSE, nreps = j))
    expect_gg(pp_check(fit, check = "dist", overlay = TRUE, nreps = j))
    expect_gg(pp_check(fit, check = "resid", nreps = j))
    expect_gg(pp_check(fit2, check = "resid", nreps = j))
    expect_gg(pp_check(fit, check = "scat", nreps = j))
    expect_gg(pp_check(fit2, check = "scat", nreps = j))
  }

  expect_gg(pp_check(fit, check = "test"))
  expect_gg(pp_check(fit, check = "test", test = "sd"))
  expect_gg(pp_check(fit, check = "test", test = c("mean","sd")))


  # by group
  expect_gg(pp_check(fit, check = "dist", group = "herd"))
  expect_gg(pp_check(fit, check = "scatter", group = "herd"))
  expect_gg(pp_check(fit, check = "test", group = "herd"))
})

test_that("pp_check ok for vb", {
  fit3 <- SW(stan_glm(mpg ~ wt, data = mtcars, iter = ITER,
                      seed = SEED, algorithm = "meanfield", iter = 10000))
  expect_gg(pp_check(fit3))
  expect_gg(pp_check(fit3, check = "resid"))
  expect_gg(pp_check(fit3, check = "scat"))
  expect_gg(pp_check(fit3, check = "test"))
})

test_that("pp_check binned residual plot works for factors", {
  ir2 <- iris[-c(1:50), ]
  ir2$Species <- factor(ir2$Species)
  fit3 <- SW(stan_glm(Species ~ Petal.Length + Petal.Width + Sepal.Length + Sepal.Width,
                      data=ir2, family = "binomial", iter = ITER, chains = CHAINS,
                      seed = SEED, refresh = REFRESH))
  expect_gg(pp_check(fit3, check = "resid"))
})


# test errors --------------------------------------------------------------
test_that("pp_check throws error if 'test' arg is bad", {
  expect_error(pp_check(fit, check = "test", test = "10982pqmeaw"),
               regexp = "not found")
  expect_error(pp_check(fit, check = "test", test = c("mean", "sd", "var")),
               regexp = "length")
})
test_that("pp_check throws error if 'group' variable not found", {
  expect_error(pp_check(fit, group = "herd2"), "not found in model frame")
})
test_that("pp_check throws error for optimizing", {
  fito <- SW(stan_glm(mpg ~ wt, data = mtcars, algorithm = "optimizing", seed = SEED))
  expect_error(pp_check(fito), regexp = "algorithm")
})


# test warnings ----------------------------------------------------------
test_that("pp_check throws warning if 'nreps' ignored ", {
  expect_warning(pp_check(fit, check = "test", nreps = 1),
                 regexp = "'nreps' is ignored")
})
test_that("pp_check throws warning if 'group' ignored", {
  expect_warning(pp_check(fit, check = "test", test = c("mean", "sd"), group = "herd"),
                 regexp = "'group' is ignored")
  expect_warning(pp_check(fit, check = "scatter", nreps = 3, group = "herd"),
                 regexp = "'group' is ignored")
  expect_warning(pp_check(fit, check = "resid", group = "herd"),
                 regexp = "'group' is ignored")
})



# helpers -----------------------------------------------------------------
test_that("ignore_nreps works", {
  ignore_nreps <- rstanarm:::ignore_nreps
  expect_null(ignore_nreps(10, "dist"))
  expect_silent(ignore_nreps(NULL, "test"))
  expect_warning(ignore_nreps(10, "test"), "'nreps' is ignored")
})

test_that("set_group works", {
  set_group <- rstanarm:::set_group
  expect_null(set_group(fit, group = NULL))
  expect_equal(set_group(fit, group = "herd"), model.frame(fit)$herd)
  expect_error(set_group(fit, group = "banana"),
               "variable 'banana' not found in model frame")
})

