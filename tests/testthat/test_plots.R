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
ITER <- 10
CHAINS <- 2
CORES <- 1

fit <- example_model
fito <- stan_glm(mpg ~ ., data = mtcars, algorithm = "optimizing", seed = SEED)
fitvb <- update(fito, algorithm = "meanfield")

expect_gg <- function(x) expect_s3_class(x, "ggplot")


# plot.stanreg ------------------------------------------------------------
context("plot.stanreg")
test_that("plot method doesn't throw bad errors and creates ggplot objects", {
  expect_default_message <- function(fit, message = "default", ...) {
    expect_message(p <- plot(fit, ...), regexp = message)
  }
  
  plotters1 <- paste0("stan_", c("plot", "trace", "hist", "dens", "ac")) # scat
  plotters2 <- paste0("stan_", c("rhat", "ess", "mcse")) # diag, par
  
  for (j in seq_along(plotters1)) {
    expect_message(p <- plot(fit, plotfun = plotters1[j]), regexp = "default")
    expect_gg(p)
    expect_gg(plot(fit, plotfun = plotters1[j], pars = "beta"))
    expect_gg(plot(fit, plotfun = plotters1[j], pars = "alpha", regex_pars = "period"))
    expect_gg(plot(fit, plotfun = plotters1[j], regex_pars = "period"))
  }
  for (j in seq_along(plotters2)) {
    expect_gg(plot(fit, plotfun = plotters2[j]))
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
  expect_gg(plot(fitvb))
  expect_gg(plot(fitvb, "hist"))
  expect_gg(plot(fitvb, "dens", separate_chains = TRUE))
  samp_only <- c("rhat", "ess", "mcse", "par", "diag", "ac")
  for (j in seq_along(samp_only)) {
    expect_error(plot(fitvb, samp_only[j]), regexp = "MCMC") 
  }
})


# pairs.stanreg -----------------------------------------------------------
context("pairs.stanreg")
test_that("pairs method ok", {
  requireNamespace("rstan")
  requireNamespace("KernSmooth")
  expect_silent(pairs(fit, pars = c("period2", "log-posterior")))
  expect_error(pairs(fito), regexp = "only available for models fit using MCMC")
})



# posterior_vs_prior ------------------------------------------------------
context("posterior_vs_prior")
test_that("posterior_vs_prior ok", {
  expect_gg(posterior_vs_prior(fit, pars = "beta"))
  expect_gg(posterior_vs_prior(fit, pars = "varying", group_by_parameter = TRUE, 
                               color_by = "vs"))
  expect_gg(posterior_vs_prior(fit, regex_pars = "period", group_by_parameter = FALSE, 
                               color_by = "none", facet_args = list(scales = "free", nrow = 2)))
  
  fit_polr <- stan_polr(tobgp ~ agegp, data = esoph, method = "probit",
                        prior = R2(0.2, "mean"), init_r = 0.1, 
                        seed = SEED, chains = CHAINS, cores = CORES, 
                        iter = 100, refresh = 0)
  expect_gg(posterior_vs_prior(fit_polr))
  expect_gg(posterior_vs_prior(fit_polr, regex_pars = "\\|", group_by_parameter = TRUE, 
                               color_by = "vs"))
})

test_that("posterior_vs_prior throws errors", {
  lmfit <- lm(mpg ~ wt, data = mtcars)
  expect_error(posterior_vs_prior(lmfit), "not a stanreg object")
  expect_error(posterior_vs_prior(fit, prob = 1), "prob < 1")
  expect_error(posterior_vs_prior(fito), 
               "only available for models fit using MCMC")
  expect_error(posterior_vs_prior(fitvb), 
               "only available for models fit using MCMC")
})
