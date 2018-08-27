# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2016, 2017 Trustees of Columbia University
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
SEED <- 1234
set.seed(SEED)

context("pp_validate")
test_that("pp_validate throws correct errors", {
  expect_error(pp_validate(example_model$stanfit), "not a stanreg object")
  expect_error(pp_validate(example_model, nreps = 1), "at least 2")
})

test_that("pp_validate runs for very quick example", {
  capture.output(
    fit <- stan_glm(mpg ~ wt, data = mtcars, seed = SEED, refresh = 0, 
                    init_r = 0.1, iter = 500)
  )
  gg <- pp_validate(fit, nreps = 2, seed = SEED)
  expect_s3_class(gg, "ggplot")
})
