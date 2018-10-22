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
SEED <- 12345
set.seed(SEED)
ITER <- 25
CHAINS <- 2

context("stan_bym2")
# source(paste0("tests/testthat/",(file.path("helpers", "expect_stanreg.R"))))
# source(paste0("tests/testthat/",(file.path("helpers", "SW.R"))))
source(file.path("helpers", "expect_stanreg.R"))
source(file.path("helpers", "SW.R"))

data("lattice", package = "rstanarm")

# Convert a spatial polygon to an N-by-N weight matrix
sp2weightmatrix <- function(spatialpolygon) {
  spdep::nb2mat(spdep::poly2nb(spatialpolygon, queen = TRUE), style = "B", zero.policy = TRUE)
}
W <- sp2weightmatrix(grid_sim15)
spatial_data <- grid_sim15@data

SW(fit_gauss <- stan_bym2(y_gauss ~ 1 + x, data = spatial_data, family = gaussian(link = "identity"),
                          prior_intercept = normal(0,1), prior = normal(0,1), prior_aux = normal(0,1),
                          prior_structured = normal(0,1), prior_mixing = beta(0.5,0.5),
                          W = W, iter = ITER, chains = CHAINS))
SW(fit_binom <- stan_bym2(y_binom ~ 1 + x, trials = spatial_data$trials, data = spatial_data,
                          prior_intercept = normal(0,1), prior = normal(0,1),
                          prior_structured = normal(0,1), prior_mixing = beta(0.5,0.5),
                          family = binomial(link = "logit"),
                          W = W, iter = ITER, chains = CHAINS))
SW(fit_pois <- stan_bym2(y_pois ~ 1 + x, data = spatial_data,
                         prior_intercept = normal(0,1), prior = normal(0,1),
                         prior_structured = normal(0,1), prior_mixing = beta(0.5,0.5),
                         family = poisson(link = "log"),
                         W = W, iter = ITER, chains = CHAINS))
SW(fit_nb2 <- stan_bym2(y_pois ~ 1 + x, data = spatial_data,
                        prior_intercept = normal(0,1), prior = normal(0,1), prior_aux = normal(0,1),
                        prior_structured = normal(0,1), prior_mixing = beta(0.5,0.5),
                        family = neg_binomial_2(link = "log"),
                        W = W, iter = ITER, chains = CHAINS))
SW(fit_gamma <- stan_bym2(y_gamma ~ 1 + x, data = spatial_data, family = Gamma(link = "log"),
                          prior_intercept = normal(0,1), prior = normal(0,1), prior_aux = normal(0,1),
                          prior_structured = normal(0,1), prior_mixing = beta(0.5,0.5),
                          W = W, iter = ITER, chains = CHAINS))

# # test family/link combinations
# test_that("family = 'gaussian' works", {
#   SW(fit_gauss_ident <- stan_bym2(y_gauss ~ 1 + x, data = spatial_data, family = gaussian(link = "identity"),
#                                    W = W, iter = 10, chains = CHAINS))
#   SW(fit_gauss_log <- stan_bym2(y_gauss ~ 1 + x, data = spatial_data, family = gaussian(link = "log"),
#                                  W = W, iter = 10, chains = CHAINS))
#   SW(fit_gauss_inv <- stan_bym2(y_gauss ~ 1 + x, data = spatial_data, family = gaussian(link = "inverse"),
#                                  W = W, iter = 10, chains = CHAINS))
#   expect_stanreg(fit_gauss_ident)
#   expect_stanreg(fit_gauss_log)
#   expect_stanreg(fit_gauss_inv)
# })
# test_that("family = 'binomial' works", {
#   SW(fit_binom_logit <- stan_bym2(y_binom ~ 1 + x, trials = spatial_data$trials, data = spatial_data,
#                                    family = binomial(link = "logit"),
#                                    W = W, iter = 10, chains = CHAINS))
#   SW(fit_binom_probit <- stan_bym2(y_binom ~ 1 + x, trials = spatial_data$trials, data = spatial_data,
#                                     family = binomial(link = "probit"),
#                                     W = W, iter = 10, chains = CHAINS))
#   SW(fit_binom_cauchit <- stan_bym2(y_binom ~ 1 + x, trials = spatial_data$trials, data = spatial_data,
#                                      family = binomial(link = "cauchit"),
#                                      W = W, iter = 10, chains = CHAINS))
#   SW(fit_binom_log <- stan_bym2(y_binom ~ 1 + x, trials = spatial_data$trials, data = spatial_data,
#                                  family = binomial(link = "log"),
#                                  W = W, iter = 10, chains = CHAINS))
#   SW(fit_binom_cloglog <- stan_bym2(y_binom ~ 1 + x, trials = spatial_data$trials, data = spatial_data,
#                                      family = binomial(link = "cloglog"),
#                                      W = W, iter = 10, chains = CHAINS))
#   expect_stanreg(fit_binom_logit)
#   expect_stanreg(fit_binom_probit)
#   expect_stanreg(fit_binom_cauchit)
#   expect_stanreg(fit_binom_log)
#   expect_stanreg(fit_binom_cloglog)
#   
# })
# test_that("family = 'poisson' works", {
#   SW(fit_pois_log <- stan_bym2(y_pois ~ 1 + x, data = spatial_data, family = poisson(link = "log"),
#                                 W = W, iter = 10, chains = CHAINS))
#   SW(fit_pois_ident <- stan_bym2(y_pois ~ 1 + x, data = spatial_data, family = poisson(link = "identity"),
#                                   W = W, iter = 10, chains = CHAINS))
#   SW(fit_pois_sqrt <- stan_bym2(y_pois ~ 1 + x, data = spatial_data, family = poisson(link = "sqrt"),
#                                  W = W, iter = 10, chains = CHAINS))
#   expect_stanreg(fit_pois_log)
#   expect_stanreg(fit_pois_ident)
#   expect_stanreg(fit_pois_sqrt)
# })
# 
# # the Gamma and neg_binomial_2 likelihoods are buggy with identity/inverse links
# test_that("family = 'neg_binomial_2' works", {
#   SW(fit_nb2_log <- stan_bym2(y_pois ~ 1 + x, data = spatial_data, family = neg_binomial_2(link = "log"),
#                                W = W, iter = 10, chains = CHAINS))
#   SW(fit_nb2_ident <- stan_bym2(y_pois ~ 1 + x, data = spatial_data, family = neg_binomial_2(link = "identity"),
#                                  W = W, iter = 10, chains = CHAINS))
#   SW(fit_nb2_sqrt <- stan_bym2(y_pois ~ 1 + x, data = spatial_data, family = neg_binomial_2(link = "sqrt"),
#                                 W = W, iter = 10, chains = CHAINS))
#   expect_stanreg(fit_nb2_log)
#   expect_stanreg(fit_nb2_ident)
#   expect_stanreg(fit_nb2_sqrt)
# })
# 
# test_that("family = 'Gamma' works", {
#   SW(fit_gamma_log <- stan_bym2(y_gamma ~ 1 + x, data = spatial_data, family = Gamma(link = "identity"),
#                                  W = W, iter = 10, chains = CHAINS))
#   SW(fit_gamma_ident <- stan_bym2(y_gamma ~ 1 + x, data = spatial_data, family = Gamma(link = "log"),
#                                    W = W, iter = 10, chains = CHAINS))
#   SW(fit_gamma_sqrt <- stan_bym2(y_gamma + 30 ~ 1 + x, data = spatial_data, family = Gamma(link = "inverse"),
#                                   W = W, iter = 10, chains = CHAINS))
#   expect_stanreg(fit_gamma_log)
#   expect_stanreg(fit_gamma_ident)
#   expect_stanreg(fit_gamma_sqrt)
# })

# test QR 
test_that("QR errors when number of predictors is <= 1", {
  expect_error(
    stan_bym2(y_gauss ~ x, data = spatial_data, family = gaussian(), W = W, seed = SEED, QR = TRUE),
    "'QR' can only be specified when there are multiple predictors"
  )
})

test_that("QR works when number of x predictors is > 1", {
  SW(stan_bym <- stan_bym2(y_gauss ~ 1 + x + I(x^2), data = spatial_data, family = gaussian(),
                             W = W, iter = ITER, chains = CHAINS, QR = TRUE))
  expect_stanreg(stan_bym)
})

test_that("stan_bym2 errors with algorithm = 'optimizing'", {
  expect_error(stan_bym2(y_gamma ~ 1 + x, data = spatial_data, family = gaussian(),
                          W = W, iter = ITER, chains = CHAINS, algorithm = "optimizing"),
               "'arg' should be one of “sampling”, “meanfield”, “fullrank”")
})

test_that("loo/waic for stan_bym2 works", {
  loo(fit_gauss)
  loo(fit_binom)
  loo(fit_pois)
  loo(fit_nb2)
  loo(fit_gamma)
})

test_that("posterior_predict works for stan_bym2", {
  preds_gauss <- posterior_predict(fit_gauss)
  preds_binom <- posterior_predict(fit_binom)
  preds_pois <- posterior_predict(fit_pois)
  preds_nb2 <- posterior_predict(fit_nb2)
  preds_gamma <- posterior_predict(fit_gamma)
})

test_that("predict works for stan_bym2", {
  new_dat <- data.frame(x = rnorm(nrow(W), 2, 1))
  predict(fit_gauss)
  predict(fit_binom)
  predict(fit_pois)
  predict(fit_nb2)
  predict(fit_gamma)
  predict(fit_gauss, newdata = new_dat)
  predict(fit_binom, newdata = new_dat)
  predict(fit_pois, newdata = new_dat)
  predict(fit_nb2, newdata = new_dat)
  predict(fit_gamma, newdata = new_dat)
})

test_that("predict errors if nrow(newdata) < number of spatial regions", {
  expect_error(predict(fit_gauss, newdata = data.frame(x = rnorm(10, 2, 1))),
               "'newdata' is less than the number of spatial regions.")
})
