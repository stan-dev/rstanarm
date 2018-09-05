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
ITER <- 50
CHAINS <- 2

context("stan_besag")
# for line-byline testing only
# source(paste0("tests/testthat/",(file.path("helpers", "expect_stanreg.R"))))
# source(paste0("tests/testthat/",(file.path("helpers", "SW.R"))))
# for full package testing
source(file.path("helpers", "expect_stanreg.R"))
source(file.path("helpers", "SW.R"))

data("lattice", package = "rstanarm")

# Convert a spatial polygon to an N-by-N weight matrix
sp2weightmatrix <- function(spatialpolygon) {
  spdep::nb2mat(spdep::poly2nb(spatialpolygon, queen = TRUE), style = "B", zero.policy = TRUE)
}
W <- sp2weightmatrix(grid_sim15)
spatial_data <- grid_sim15@data

SW(fit_gauss <- stan_besag(y_gauss ~ 1 + x, data = spatial_data, family = gaussian(link = "identity"),
                           prior_intercept = normal(0,1), prior = normal(0,1), prior_aux = normal(0,1),
                           W = W, iter = 100, chains = 4))
SW(fit_binom <- stan_besag(y_binom ~ 1 + x, trials = spatial_data$trials, data = spatial_data,
                           prior_intercept = normal(0,1), prior = normal(0,1),
                           family = binomial(link = "logit"),
                           W = W, iter = 100, chains = 4))
SW(fit_pois <- stan_besag(y_pois ~ 1 + x, data = spatial_data,
                          prior_intercept = normal(0,1), prior = normal(0,1),
                          family = poisson(link = "log"),
                          W = W, iter = 500, chains = 4))
SW(fit_nb2 <- stan_besag(y_pois ~ 1 + x, data = spatial_data,
                         prior_intercept = normal(0,1), prior = normal(0,1), prior_aux = normal(0,1),
                         family = neg_binomial_2(link = "log"),
                         W = W, iter = 500, chains = 4))
SW(fit_gamma <- stan_besag(y_gamma ~ 1 + x, data = spatial_data, family = Gamma(link = "log"),
                           prior_intercept = normal(0,1), prior = normal(0,1), prior_aux = normal(0,1),
                           W = W, iter = 500, chains = 4))

# compare answers with INLA (NB2 reciprocal_dispersion param fails!)
# test_that("stan_besag estimates match INLA", {
#   inla_gauss_est <- c(-0.0475, 0.4329, 1/0.8277, 1/0.9830)
#   inla_binom_est <- c(-1.9444, 0.3340, 1/0.4229)
#   inla_pois_est <- c(0.7801, 0.4366, 1/0.3448)
#   inla_nb_est <- c(0.8053, 0.4362, 1/0.3574, 48.6003)
#   inla_gamma_est <- c(0.8119, -1.3407, 1/0.559, 1/1.106)
#   besag_gauss <- unname(fit_gauss$stan_summary[1:4,"mean"])
#   besag_binom <- unname(fit_binom$stan_summary[1:3,"mean"])
#   besag_pois <- unname(fit_pois$stan_summary[1:3,"mean"])
#   besag_nb2 <- unname(fit_nb2$stan_summary[1:4,"mean"])
#   besag_gamma <- unname(fit_gamma$stan_summary[1:4,"mean"])
#   expect_equal(besag_gauss, inla_gauss_est, tol = 0.2)
#   expect_equal(besag_binom, inla_binom_est, tol = 0.2)
#   expect_equal(besag_pois, inla_pois_est, tol = 0.2)
#   expect_equal(besag_nb2, inla_nb_est, tol = 0.2)
#   expect_equal(besag_gamma, inla_gamma_est, tol = 0.2)
# })

# test family/link combinations
test_that("family = 'gaussian' works", {
  SW(fit_gauss_ident <- stan_besag(y_gauss ~ 1 + x, data = spatial_data, family = gaussian(link = "identity"),
                                   W = W, iter = ITER, chains = CHAINS))
  SW(fit_gauss_log <- stan_besag(y_gauss ~ 1 + x, data = spatial_data, family = gaussian(link = "log"),
                                 W = W, iter = ITER, chains = CHAINS))
  SW(fit_gauss_inv <- stan_besag(y_gauss ~ 1 + x, data = spatial_data, family = gaussian(link = "inverse"),
                                 W = W, iter = ITER, chains = CHAINS))
  expect_stanreg(fit_gauss_ident)
  expect_stanreg(fit_gauss_log)
  expect_stanreg(fit_gauss_inv)
})
test_that("family = 'binomial' works", {
  SW(fit_binom_logit <- stan_besag(y_binom ~ 1 + x, trials = spatial_data$trials, data = spatial_data,
                                   family = binomial(link = "logit"),
                                   W = W, iter = ITER, chains = CHAINS))
  SW(fit_binom_probit <- stan_besag(y_binom ~ 1 + x, trials = spatial_data$trials, data = spatial_data,
                                    family = binomial(link = "probit"),
                                   W = W, iter = ITER, chains = CHAINS))
  SW(fit_binom_cauchit <- stan_besag(y_binom ~ 1 + x, trials = spatial_data$trials, data = spatial_data,
                                     family = binomial(link = "cauchit"),
                                   W = W, iter = ITER, chains = CHAINS))
  SW(fit_binom_log <- stan_besag(y_binom ~ 1 + x, trials = spatial_data$trials, data = spatial_data,
                                 family = binomial(link = "log"),
                                   W = W, iter = ITER, chains = CHAINS))
  SW(fit_binom_cloglog <- stan_besag(y_binom ~ 1 + x, trials = spatial_data$trials, data = spatial_data,
                                     family = binomial(link = "cloglog"),
                                   W = W, iter = ITER, chains = CHAINS))
  expect_stanreg(fit_binom_logit)
  expect_stanreg(fit_binom_probit)
  expect_stanreg(fit_binom_cauchit)
  expect_stanreg(fit_binom_log)
  expect_stanreg(fit_binom_cloglog)
  
})
test_that("family = 'poisson' works", {
  SW(fit_pois_log <- stan_besag(y_pois ~ 1 + x, data = spatial_data, family = poisson(link = "log"),
                                 W = W, iter = ITER, chains = CHAINS))
  SW(fit_pois_ident <- stan_besag(y_pois ~ 1 + x, data = spatial_data, family = poisson(link = "identity"),
                                W = W, iter = ITER, chains = CHAINS))
  SW(fit_pois_sqrt <- stan_besag(y_pois ~ 1 + x, data = spatial_data, family = poisson(link = "sqrt"),
                                W = W, iter = ITER, chains = CHAINS))
  expect_stanreg(fit_pois_log)
  expect_stanreg(fit_pois_ident)
  expect_stanreg(fit_pois_sqrt)
})

# the Gamma and neg_binomial_2 likelihoods are buggy with identity/inverse links
test_that("family = 'neg_binomial_2' works", {
  SW(fit_nb2_log <- stan_besag(y_pois ~ 1 + x, data = spatial_data, family = neg_binomial_2(link = "log"),
                                W = W, iter = ITER, chains = CHAINS))
  SW(fit_nb2_ident <- stan_besag(y_pois ~ 1 + x, data = spatial_data, family = neg_binomial_2(link = "identity"),
                                  W = W, iter = ITER, chains = CHAINS))
  SW(fit_nb2_sqrt <- stan_besag(y_pois ~ 1 + x, data = spatial_data, family = neg_binomial_2(link = "sqrt"),
                                 W = W, iter = ITER, chains = CHAINS))
  expect_stanreg(fit_nb2_log)
  expect_stanreg(fit_nb2_ident)
  expect_stanreg(fit_nb2_sqrt)
})

test_that("family = 'Gamma' works", {
  SW(fit_gamma_log <- stan_besag(y_gamma ~ 1 + x, data = spatial_data, family = Gamma(link = "identity"),
                               W = W, iter = ITER, chains = CHAINS))
  SW(fit_gamma_ident <- stan_besag(y_gamma ~ 1 + x, data = spatial_data, family = Gamma(link = "log"),
                                 W = W, iter = ITER, chains = CHAINS))
  SW(fit_gamma_sqrt <- stan_besag(y_gamma + 30 ~ 1 + x, data = spatial_data, family = Gamma(link = "inverse"),
                                W = W, iter = ITER, chains = CHAINS))
  expect_stanreg(fit_gamma_log)
  expect_stanreg(fit_gamma_ident)
  expect_stanreg(fit_gamma_sqrt)
})

# test QR 
test_that("QR errors when number of predictors is <= 1", {
  expect_error(
    stan_besag(y_gauss ~ x, data = spatial_data, W = W, family = gaussian(), seed = SEED, QR = TRUE),
    "'QR' can only be specified when there are multiple predictors"
  )
})

test_that("QR works when number of x predictors is > 1", {
  SW(fit_besag <- stan_besag(y_gauss ~ 1 + x + I(x^2), data = spatial_data, family = gaussian(),
                                   W = W, iter = ITER, chains = CHAINS, QR = TRUE))
  expect_stanreg(fit_besag)
})

test_that("stan_besag errors with algorithm = 'optimizing'", {
  expect_error(stan_besag(y_gamma ~ 1 + x, data = spatial_data, family = gaussian(),
                          W = W, iter = ITER, chains = CHAINS, algorithm = "optimizing"),
               "'arg' should be one of “sampling”, “meanfield”, “fullrank”")
})

# test_that("loo/waic for stan_besag works", {
#   loo(fit_gauss)
#   loo(fit_binom)
#   loo(fit_pois)
#   loo(fit_nb2)
#   loo(fit_gamma)
# })

test_that("posterior_predict works for stan_besag", {
  preds_gauss <- posterior_predict(fit_gauss)
  preds_binom <- posterior_predict(fit_binom)
  preds_pois <- posterior_predict(fit_pois)
  preds_nb2 <- posterior_predict(fit_nb2)
  preds_gamma <- posterior_predict(fit_gamma)
})

test_that("predict works for stan_besag", {
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

