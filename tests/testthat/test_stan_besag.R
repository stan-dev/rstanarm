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

if (.Platform$OS.type != "windows" && require(betareg)) {
  library(rstanarm)
  SEED <- 12345
  set.seed(SEED)
  ITER <- 100
  CHAINS <- 2
  
  context("stan_besag")
  
  source(file.path("helpers", "expect_stanreg.R"))
  source(file.path("helpers", "SW.R"))
  
  data("lattice10", package = "rstanarm")
  
  # Convert a spatial polygon to an N-by-N weight matrix
  sp2weightmatrix <- function(spatialpolygon) {
    spdep::nb2mat(spdep::poly2nb(spatialpolygon, queen = TRUE), style = "B", zero.policy = TRUE)
  }
  W <- sp2weightmatrix(grid_sim)
  x <- rnorm(nrow(W), 2, 1)
  spatial_data <- data.frame(x, phi = grid_sim@data$gmrf)
  spatial_data$y_gauss <- rnorm(nrow(W), 0 + 0.4 * x + spatial_data$phi, 1)
  spatial_data$y_gauss <- rnorm(nrow(W), 0 + 0.4 * x + spatial_data$phi, 1)
  spatial_data$y_pois <- rpois(nrow(W), exp(-1 + 1.4 * x + spatial_data$phi))
  spatial_data$trials <- rep(10, nrow(W))
  spatial_data$y_binom <- rbinom(nrow(W), spatial_data$trials, binomial(link = "logit")$linkinv(1 + 0.4 * x + spatial_data$phi))
  
  fit_besag <- stan_besag(y_gauss ~ 1 + x, data = spatial_data, W = W, iter = 1e3, chains = 4,
                          prior_intercept = normal(0,1), prior = normal(0,1), prior_rho = normal(0,1),
                          prior_sigma = normal(0,1))
  
  # test family/link combinations
  test_that("family = 'gaussian' works", {
    SW(fit_gauss_ident <- stan_besag(y_gauss ~ 1 + x, data = spatial_data, family = gaussian(),
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
                                  W = W, iter = 10, chains = CHAINS))
    expect_stanreg(fit_pois_log)
    expect_stanreg(fit_pois_ident)
    expect_stanreg(fit_pois_sqrt)
  })
  
  # test QR 
  test_that("QR errors when number of predictors is <= 1", {
    expect_error(
      stan_besag(y_gauss ~ x, data = spatial_data, family = gaussian(), seed = SEED, QR = TRUE),
      "'QR' can only be specified when there are multiple predictors"
    )
  })
  
  test_that("QR works when number of x and/or z predictors is >= 1", {
    SW(fit_besag <- stan_besag(y_gauss ~ 1 + x + I(x^2), data = spatial_data, family = gaussian(),
                                     W = W, iter = ITER, chains = CHAINS, QR = TRUE))
    expect_stanreg(fit_besag)
  })
  
  test_that("stan_besag works with QR = TRUE and algorithm = 'optimizing'", {

  })
  
  test_that("loo/waic for stan_besag works", {
    SW(fit_gauss <- stan_besag(y_gauss ~ 1 + x + I(x^2), data = spatial_data, family = gaussian(),
                               W = W, iter = ITER, chains = CHAINS, QR = TRUE))
    SW(fit_binom <- stan_besag(y_binom ~ 1 + x, trials = spatial_data$trials, data = spatial_data,
                                     family = binomial(link = "logit"),
                                     W = W, iter = ITER, chains = CHAINS))
    SW(fit_pois <- stan_besag(y_pois ~ 1 + x, data = spatial_data,
                               family = poisson(link = "log"),
                               W = W, iter = ITER, chains = CHAINS))
    loo(fit_gauss)
    loo(fit_binom)
    loo(fit_pois)
  })
  
  test_that("posterior_predict works for stan_besag", {
    SW(fit_gauss <- stan_besag(y_gauss ~ 1 + x + I(x^2), data = spatial_data, family = gaussian(),
                               W = W, iter = ITER, chains = CHAINS, QR = TRUE))
    SW(fit_binom <- stan_besag(y_binom ~ 1 + x, trials = spatial_data$trials, data = spatial_data,
                               family = binomial(link = "logit"),
                               W = W, iter = ITER, chains = CHAINS))
    SW(fit_pois <- stan_besag(y_pois ~ 1 + x, data = spatial_data,
                              family = poisson(link = "log"),
                              W = W, iter = ITER, chains = CHAINS))
    preds_gauss <- posterior_predict(fit_gauss)
    preds_binom <- posterior_predict(fit_binom)
    preds_pois <- posterior_predict(fit_pois)
  })
  
  test_that("predict works for stan_besag", {
    SW(fit_besag <- stan_besag(y_gauss ~ 1 + x + I(x^2), data = spatial_data, family = gaussian(),
                               W = W, iter = ITER, chains = CHAINS, QR = TRUE))
    pred_besag <- predict(fit_besag)
    pred_new_besag <- predict(fit_besag, newdata = data.frame(x = rnorm(100, 2, 1)))
  })
  
  test_that("predict errors if nrow(newdata) < number of spatial regions", {
    SW(fit_besag <- stan_besag(y_gauss ~ 1 + x + I(x^2), data = spatial_data, family = gaussian(),
                               W = W, iter = ITER, chains = CHAINS, QR = TRUE))
    pred_besag <- predict(fit_besag)
    expect_error(pred_new_besag <- predict(fit_besag, newdata = data.frame(x = rnorm(10, 2, 1))),
                 "'newdata' is less than the number of spatial regions.")
  })
  
}