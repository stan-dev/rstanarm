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
ITER <- 10
CHAINS <- 2
CORES <- 1

fit <- example_model
fito <- stan_glm(mpg ~ ., data = mtcars, algorithm = "optimizing", seed = SEED)
fitvb <- update(fito, algorithm = "meanfield")

context("plot.stanreg")
test_that("plot method doesn't throw bad errors and creates ggplot objects", {
  expect_default_message <- function(fit, message = "default", ...) {
    expect_message(p <- plot(fit, ...), regexp = message)
  }
  
  plotters1 <- paste0("stan_", c("plot", "trace", "hist", "dens", "ac")) # scat
  plotters2 <- paste0("stan_", c("rhat", "ess", "mcse")) # diag, par
  
  for (j in seq_along(plotters1)) {
    expect_message(p <- plot(fit, plotfun = plotters1[j]), regexp = "default")
    expect_is(p, "ggplot")
    expect_is(p <- plot(fit, plotfun = plotters1[j], pars = "beta"), "ggplot")
    expect_is(p <- plot(fit, plotfun = plotters1[j], pars = "alpha", 
                        regex_pars = "period"), "ggplot")
    expect_is(p <- plot(fit, plotfun = plotters1[j], regex_pars = "period"), 
              "ggplot")
  }
  for (j in seq_along(plotters2)) {
    expect_silent(p <- plot(fit, plotfun = plotters2[j]))
    expect_is(p, "ggplot")
  }
})

test_that("plot.stanreg ok for optimization", {
  expect_silent(plot(fito))
  expect_silent(plot(fito, pars = "alpha"))
  expect_silent(plot(fito, pars = "beta"))
  expect_silent(plot(fito, pars = c("wt", "gear")))
  expect_warning(plot(fito, regex_pars = "wt"), regexp = "'regex_pars' ignored")
  expect_error(plot(fito, "trace"), regexp = "'plotfun'")
})

test_that("plot.stanreg ok for vb", {
  expect_is(plot(fitvb), "ggplot")
  expect_is(plot(fitvb, "hist"), "ggplot")
  expect_is(plot(fitvb, "dens", separate_chains = TRUE), "ggplot")
  samp_only <- c("rhat", "ess", "mcse", "par", "diag", "ac")
  for (j in seq_along(samp_only)) {
    expect_error(plot(fitvb, samp_only[j]), regexp = "MCMC") 
  }
})

context("pairs.stanreg")
test_that("pairs method ok", {
  requireNamespace("rstan")
  requireNamespace("KernSmooth")
  expect_silent(pairs(fit, pars = c("period2", "log-posterior")))
  expect_error(pairs(fito), regexp = "only available for models fit using MCMC")
})
