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

library(rstanarm)
SEED <- 123
ITER <- 10
CHAINS <- 2
CORES <- 1

source(test_path("helpers", "SW.R"))
source(test_path("helpers", "expect_gg.R"))

fit <- example_model
SW(fito <- stan_glm(mpg ~ ., data = mtcars, algorithm = "optimizing", seed = SEED, refresh = 0))
SW(fitvb <- update(fito, algorithm = "meanfield"))

# plot.stanreg ------------------------------------------------------------
context("plot.stanreg")
test_that("plot.stanreg errors if chains = 1 but needs multiple", {
  multiple_chain_plots <- c("trace_highlight",
                            "hist_by_chain",
                            "dens_overlay",
                            "violin")
  SW(fit_1chain <- stan_glm(mpg ~ wt, data = mtcars, chains = 1, iter = 100, refresh = 0))
  for (f in multiple_chain_plots) {
    expect_error(plot(fit_1chain, plotfun = f), info = f, 
                 regexp = "requires multiple chains")
  }
})

test_that("other plot.stanreg errors thrown correctly", {
  expect_error(plot(fit, plotfun = "9999"), 
               "not a valid MCMC function name")
  expect_error(plot(fit, plotfun = "ppc_hist"), 
               "use the 'pp_check' method")
  expect_error(plot(fit, plotfun = "stan_diag"), 
               "help('NUTS', 'bayesplot')", fixed = TRUE)
})

test_that("plot.stanreg returns correct object", {
  # ggplot objects
  ggplot_object_plots <- c(
    "intervals", "areas",
    "dens", "dens_overlay", 
    "hist", "hist_by_chain",
    "trace", "trace_highlight",
    "violin", 
    "rhat", "rhat_hist", 
    "neff", "neff_hist", "ess",
    "acf", "acf_bar", "ac"
  )
  for (f in ggplot_object_plots)
    expect_gg(plot(fit, f))
  
  # requires exactly 2 parameters
  expect_gg(plot(fit, "scat", pars = c("period2", "period3")))
  expect_gg(plot(fit, "hex", pars = c("period2", "period3")))
})

test_that("plot method returns correct object for nuts diagnostic plots", {
  # energy plot returns ggplot object
  expect_gg(plot(fit, "nuts_energy"))
  
  # others return gtable objects
  gtable_object_plots <-
    paste0("nuts_",
           c("stepsize", "acceptance", "divergence", "treedepth"))
  for (f in gtable_object_plots)
    expect_s3_class(plot(fit, plotfun = f), "gtable")
})

test_that("plot.stanreg ok for optimization", {
  expect_gg(plot(fito))
  expect_gg(plot(fito, "areas"))
  expect_gg(plot(fito, "dens"))
  expect_gg(plot(fito, "scatter", pars = c("wt", "cyl")))
  expect_gg(plot(fito, pars = c("alpha", "beta")))
  
  expect_warning(plot(fito, regex_pars = "wt"), 
                 regexp = "'regex_pars' ignored")
  expect_error(plot(fito, "trace"), 
               regexp = "only available for models fit using MCMC")
  expect_error(plot(fito, "nuts_acceptance"), 
               regexp = "only available for models fit using MCMC")
  expect_error(plot(fito, "rhat_hist"), 
               regexp = "only available for models fit using MCMC")
})

test_that("plot.stanreg ok for vb", {
  expect_gg(plot(fitvb))
  expect_gg(plot(fitvb, "areas"))
  expect_gg(plot(fitvb, "dens"))
  expect_gg(plot(fitvb, "scatter", pars = c("wt", "cyl")))
  expect_gg(plot(fitvb, pars = c("alpha", "beta")))
  
  expect_error(plot(fitvb, "trace"), 
               regexp = "only available for models fit using MCMC")
  expect_error(plot(fitvb, "nuts_acceptance"), 
               regexp = "only available for models fit using MCMC")
  expect_error(plot(fitvb, "rhat_hist"), 
               regexp = "only available for models fit using MCMC")
  expect_error(plot(fitvb, "mcmc_neff"), 
               regexp = "only available for models fit using MCMC")
})


# pairs.stanreg -----------------------------------------------------------
context("pairs.stanreg")
test_that("pairs method ok", {
  expect_silent(pairs(fit, pars = c("period2", "log-posterior")))
  expect_silent(pairs(fit, pars = "b[(Intercept) herd:15]", regex_pars = "Sigma"))
  expect_silent(pairs(fit, pars = "b[(Intercept) herd:15]", regex_pars = "Sigma", 
                      condition = pairs_condition(nuts = "lp__")))
  expect_error(pairs(fitvb), regexp = "only available for models fit using MCMC")
  expect_error(pairs(fito), regexp = "only available for models fit using MCMC")
})



# posterior_vs_prior ------------------------------------------------------
context("posterior_vs_prior")
test_that("posterior_vs_prior ok", {
  SW(p1 <- posterior_vs_prior(fit, pars = "beta"))
  expect_gg(p1)
  
  SW(p2 <- posterior_vs_prior(fit, pars = "varying", group_by_parameter = TRUE, 
                              color_by = "vs"))
  expect_gg(p2)
  SW(p3 <- posterior_vs_prior(fit, regex_pars = "period", 
                              group_by_parameter = FALSE, 
                              color_by = "none", 
                              facet_args = list(scales = "free", nrow = 2)))
  expect_gg(p3)
  
  SW(fit_polr <- stan_polr(tobgp ~ agegp, data = esoph, method = "probit",
                          prior = R2(0.2, "mean"), init_r = 0.1, 
                          seed = SEED, chains = CHAINS, cores = CORES, 
                          iter = 100, refresh = 0))
  SW(p4 <- posterior_vs_prior(fit_polr))
  SW(p5 <- posterior_vs_prior(fit_polr, regex_pars = "\\|", 
                              group_by_parameter = TRUE, 
                              color_by = "vs"))
  expect_gg(p4)
  expect_gg(p5)
})

test_that("posterior_vs_prior throws errors", {
  lmfit <- lm(mpg ~ wt, data = mtcars)
  expect_error(posterior_vs_prior(lmfit), "no applicable method")
  expect_error(posterior_vs_prior(fit, prob = 1), "prob < 1")
  expect_error(posterior_vs_prior(fito), 
               "only available for models fit using MCMC")
  expect_error(posterior_vs_prior(fitvb), 
               "only available for models fit using MCMC")
})
