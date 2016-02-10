# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2015 Trustees of Columbia University
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

context("pp_check")
test_that("pp_check doesn't throw bad errors", {
  expect_silent(p <- pp_check(fit, check = "dist", overlay = TRUE, size = 2))
  expect_silent(p <- pp_check(fit, check = "resid"))
  expect_silent(p <- pp_check(fit2, check = "resid", fill = "red", bins = 15))
  expect_silent(p <- pp_check(fit, check = "scatter"))
  expect_silent(p <- pp_check(fit2, check = "scatter", color = "purple"))
  expect_is(p, "ggplot")
  for (j in 1:2) {
    expect_silent(p <- pp_check(fit, check = "dist", overlay = FALSE, nreps = j))
    expect_silent(p <- pp_check(fit, check = "dist", overlay = TRUE, nreps = j))
    expect_silent(p <- pp_check(fit, check = "resid", nreps = j))
    expect_silent(p <- pp_check(fit2, check = "resid", nreps = j))
    expect_silent(p <- pp_check(fit, check = "scat", nreps = j))
    expect_silent(p <- pp_check(fit2, check = "scat", nreps = j))
  }
  expect_silent(p <- pp_check(fit, check = "test"))
  expect_silent(p <- pp_check(fit, check = "test", test = "sd"))
  expect_silent(p <- pp_check(fit, check = "test", test = c("mean","sd")))
  expect_is(p, "ggplot")
})

test_that("pp_check ok for vb", {
  fit3 <- SW(update(fit2, algorithm = "meanfield", iter = 10000))
  expect_silent(p <- pp_check(fit3))
  expect_silent(p <- pp_check(fit3, check = "resid"))
  expect_silent(p <- pp_check(fit3, check = "scat"))
  expect_silent(p <- pp_check(fit3, check = "test"))
})

test_that("pp_check throws appropriate errors", {
  expect_error(p <- pp_check(fit, check = "test", test = "10982pqmeaw"), 
               regexp = "not found")
  expect_error(p <- pp_check(fit, check = "test", test = c("mean", "sd", "var")), 
               regexp = "length 1 or 2")
  
  fito <- stan_glm(mpg ~ wt, data = mtcars, algorithm = "optimizing", seed = SEED)
  expect_error(pp_check(fito), regexp = "algorithm")
  expect_error(pp_check(rnorm(10)), regexp = "not a stanreg object")
})

test_that("pp_check throws appropriate warnings", {
  expect_warning(p <- pp_check(fit, check = "test", nreps = 1), 
                 regexp = "'nreps' is ignored")
})

test_that("pp_check binned residual plot ok for factors", {
  ir2 <- iris[-c(1:50), ]
  ir2$Species <- factor(ir2$Species)
  fit3 <- SW(stan_glm(Species ~ Petal.Length + Petal.Width + Sepal.Length + Sepal.Width, 
                   data=ir2, family = "binomial", iter = ITER, chains = CHAINS,
                   seed = SEED, refresh = REFRESH))
  expect_silent(p <- pp_check(fit3, check = "resid"))
})
