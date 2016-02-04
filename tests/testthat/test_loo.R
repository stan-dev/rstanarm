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

library(rstanarm)
library(loo)
options(loo.cores = 2)
SEED <- 1234
set.seed(SEED)
CHAINS <- 2
ITER <- 40 # small iter for speed but large enough for psis
REFRESH <- 0

SW <- suppressWarnings

# These tests just check that the loo.stanreg method (which calls loo.function
# method) results are identical to the loo.matrix results. Since for these tests 
# the log-likelihood matrix is computed using the log-likelihood function, the 
# only thing these tests really do is make sure that loo.stanreg and all the 
# log-likelihood functions don't return any errors and whatnot (it does not check
# that the results returned by loo are actually correct). 

context("loo and waic")

expect_identical_loo <- function(fit) {
  expect_identical(SW(loo(fit)), SW(loo(log_lik(fit))))
  expect_equal(waic(fit), waic(log_lik(fit)))
}
mcmc_only_error <- function(fit) {
  msg <- "only available for models fit using MCMC"
  expect_error(loo(fit), regexp = msg)
  expect_error(waic(fit), regexp = msg)
}

test_that("loo & waic throw error for non mcmc models", {
  fito <- stan_glm(mpg ~ wt, data = mtcars, algorithm = "optimizing", seed = SEED)
  fitvb1 <- update(fito, algorithm = "meanfield")
  fitvb2 <- update(fito, algorithm = "fullrank")
  mcmc_only_error(fito)
  mcmc_only_error(fitvb1)
  mcmc_only_error(fitvb2)
})

test_that("loo/waic for stan_glm works", {
  # gaussian
  fit_gaus <- SW(stan_glm(mpg ~ wt, data = mtcars, chains = CHAINS, iter = ITER, 
                          seed = SEED, refresh = REFRESH))
  expect_identical_loo(fit_gaus)
  
  # binomial
  dat <- data.frame(ldose = rep(0:5, 2), sex = factor(rep(c("M", "F"), c(6, 6))))
  numdead <- c(1, 4, 9, 13, 18, 20, 0, 2, 6, 10, 12, 16)
  SF <- cbind(numdead, numalive = 20-numdead)
  fit_binom <- SW(stan_glm(SF ~ sex*ldose, data = dat, family = binomial, 
                           chains = CHAINS, iter = ITER, seed = SEED, 
                           refresh = REFRESH))
  dead <- rbinom(length(numdead), 1, prob = 0.5)
  fit_binom2 <- SW(update(fit_binom, formula = factor(dead) ~ .))
  expect_identical_loo(fit_binom)
  expect_identical_loo(fit_binom2)
  
  # poisson 
  d.AD <- data.frame(treatment = gl(3,3), outcome =  gl(3,1,9), 
                     counts = c(18,17,15,20,10,20,25,13,12))
  fit_pois <- SW(stan_glm(counts ~ outcome + treatment, data = d.AD, family = poisson,
                          chains = CHAINS, iter = ITER, seed = SEED, refresh = REFRESH))
  expect_identical_loo(fit_pois)
  
  # negative binomial
  fit_negbin <- SW(update(fit_pois, family = neg_binomial_2))
  expect_identical_loo(fit_negbin)
  
  # gamma
  clotting <- data.frame(log_u = log(c(5,10,15,20,30,40,60,80,100)),
                         lot1 = c(118,58,42,35,27,25,21,19,18),
                         lot2 = c(69,35,26,21,18,16,13,12,12))
  fit_gamma <- SW(stan_glm(lot1 ~ log_u, data = clotting, family = Gamma, 
                           chains = CHAINS, iter = ITER, seed = SEED, refresh = REFRESH))
  expect_identical_loo(fit_gamma)
  
  # inverse gaussian
  fit_igaus <- SW(update(fit_gamma, family = inverse.gaussian))
  expect_identical_loo(fit_igaus)
})

test_that("loo/waic for stan_polr works", {
  # logistic
  fit_logistic <- SW(stan_polr(tobgp ~ agegp, data = esoph, prior = R2(0.2, "mean"),  
                               init_r = 0.1, chains = CHAINS, iter = ITER, 
                               seed = SEED, refresh = REFRESH))
  expect_identical_loo(fit_logistic)
})

test_that("loo/waic for stan_lm works", {
  fit_lm <- SW(stan_lm(mpg ~ ., data = mtcars, prior = R2(0.75), 
                       chains = CHAINS, iter = ITER, seed = SEED, refresh = REFRESH))
  expect_identical_loo(fit_lm)
})

test_that("loo/waic for stan_glmer works", {
  # gaussian
  fit_glmer1 <- SW(stan_glmer(mpg ~ wt + (1|cyl) + (1+wt|gear), data = mtcars, 
                              chains = CHAINS, iter = ITER, seed = SEED, 
                              refresh = REFRESH))
  expect_identical_loo(fit_glmer1)
  
  # binomial
  expect_identical_loo(example_model)
})


context("loo and waic helpers")
test_that("ll_fun works", {
  ll_fun <- rstanarm:::ll_fun
  expect_identical(ll_fun(gaussian(link = "log")), rstanarm:::.ll_gaussian_i)
  expect_identical(ll_fun(binomial()), rstanarm:::.ll_binomial_i)
  expect_identical(ll_fun(poisson()), rstanarm:::.ll_poisson_i)
  expect_identical(ll_fun("logistic"), rstanarm:::.ll_polr_i)
  expect_error(ll_fun(example_model), "must be a family or a character string")
})
test_that(".weighted works", {
  f <- rstanarm:::.weighted
  expect_equal(f(2, NULL), 2)
  expect_equal(f(2, 3), 6)
  expect_equal(f(8, 0.25), 2)
  expect_error(f(2), "missing, with no default")
})
