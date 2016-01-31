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
SEED <- 12345
ITER <- 10L
CHAINS <- 2L
REFRESH <- ITER

context("family stuff")

test_that("default_prior_params works", {
  defs <- rstanarm:::default_hyperparams
  expect_error(defs("gaussian"), regexp = "'family' should be a family object")
  
  ans <- list(scaled = TRUE, 
              prior_scale_for_dispersion = 5, 
              min_prior_scale = 1e-12)
  
  expect_identical(defs(family = gaussian()), ans)
  expect_identical(defs(Gamma()), ans)
  expect_identical(defs(neg_binomial_2()), ans)
  expect_identical(defs(poisson()), ans)
  expect_identical(defs(binomial()), ans)
  expect_identical(defs(inverse.gaussian()), ans)
  expect_identical(defs(family = t_family()), 
                   c(ans, prior_shape_for_df = 2, prior_rate_for_df = 0.1))
})

test_that("rstanarm_family works", {
  rstanarm_family <- rstanarm:::rstanarm_family
  is.rstanarm_family <- rstanarm:::is.rstanarm_family
  expect_error(rstanarm_family("abcd"))
  expect_error(rstanarm_family("gaussian", link = "abcd"), 
               regexp = "link not recognised")
  expect_error(rstanarm_family("binomial", link = "abcd"), 
               regexp = "link not recognised")
  
  expect_equal(rstanarm_family("binomial", "probit")$family, 
               binomial("probit"))
  expect_equal(rstanarm_family("poisson", "sqrt")$family, 
               poisson("sqrt"))
  
  
  gaus1 <- rstanarm_family("gaussian", "identity", 
                           scaled = FALSE,
                           prior_scale_for_dispersion = 2)
  expect_true(is.rstanarm_family(gaus1))
  expect_equal(gaus1$family, gaussian(link = "identity"))
  expect_equal(gaus1$params$prior_scale_for_dispersion, 2)
  expect_false(gaus1$params$scaled)
  
  gamma1 <- rstanarm_family(family = "Gamma", link = "log",
                            min_prior_scale = 0.1,
                            prior_scale_for_shape = 3, 
                            scaled = FALSE)
  gamma2 <- rstanarm_family("Gamma", "log",
                            min_prior_scale = 0.1,
                            prior_scale_for_shape = 3, 
                            scaled = FALSE)
  expect_equal(gamma1, gamma2)
  expect_true(is.rstanarm_family(gamma1))
  expect_equal(gamma1$family, Gamma(link = "log"))
  expect_equal(gamma1$params$prior_scale_for_dispersion, 3)
  expect_equal(gamma1$params$min_prior_scale, 0.1)
  expect_false(gamma1$params$scaled)
  
  gamma3 <- rstanarm_family("Gamma")
  expect_equal(gamma3$family, Gamma())
  expect_equal(gamma3$params$prior_scale_for_dispersion, 5)
  expect_equal(gamma3$params$min_prior_scale, 1e-12)
  expect_true(gamma3$params$scaled)
  
  t1 <- rstanarm_family(family = "t_family")
  expect_true(is.rstanarm_family(t1))
  expect_equal(t1$family, t_family())
  expect_equal(t1$params$prior_scale_for_dispersion, 5)
  expect_equal(t1$params$prior_shape_for_df, 2)
  expect_equal(t1$params$prior_rate_for_df, 0.1)
  expect_true(t1$params$scaled)
  
  t2 <- rstanarm_family(family = "t_family", link = "log", 
                        prior_rate_for_df = NULL, scaled = FALSE)
  expect_true(is.rstanarm_family(t2))
  expect_equal(t2$family, t_family("log"))
  expect_equal(t2$params$prior_scale_for_dispersion, 5)
  expect_equal(t2$params$prior_shape_for_df, 2)
  expect_equal(t2$params$prior_rate_for_df, NULL)
  expect_false(t2$params$scaled)
})

