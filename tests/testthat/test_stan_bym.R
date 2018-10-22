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

context("stan_bym")
# source(paste0("tests/testthat/",(file.path("helpers", "expect_stanreg.R"))))
# source(paste0("tests/testthat/",(file.path("helpers", "SW.R"))))
source(file.path("helpers", "expect_stanreg.R"))
source(file.path("helpers", "SW.R"))

data("lattice", package = "rstanarm")

# Convert a spatial polygon to an N-by-N weight matrix
sp2weightmatrix <- function(spatialpolygon) {
  spdep::nb2mat(spdep::poly2nb(spatialpolygon, queen = TRUE), style = "B", zero.policy = TRUE)
}
adj <- sp2weightmatrix(grid_sim15)
spatial_data <- grid_sim15@data

SW(fit_gauss <- stan_bym(y_gauss ~ 1 + x, data = spatial_data, family = gaussian(link = "identity"),
                         prior_intercept = normal(0,1), prior = normal(0,1), prior_aux = normal(0,1),
                         prior_unstructured = normal(0,1), prior_structured = normal(0,1),
                         W = adj, iter = ITER, chains = CHAINS))
SW(fit_binom <- stan_bym(y_binom ~ 1 + x, trials = spatial_data$trials, data = spatial_data,
                         prior_intercept = normal(0,1), prior = normal(0,1),
                         prior_unstructured = normal(0,1), prior_structured = normal(0,1),
                         family = binomial(link = "logit"),
                         W = adj, iter = ITER, chains = CHAINS))
SW(fit_pois <- stan_bym(y_pois ~ 1 + x, data = spatial_data,
                        prior_intercept = normal(0,1), prior = normal(0,1),
                        prior_unstructured = normal(0,1), prior_structured = normal(0,1),
                        family = poisson(link = "log"),
                        W = adj, iter = ITER, chains = CHAINS))
SW(fit_nb2 <- stan_bym(y_pois ~ 1 + x, data = spatial_data,
                       prior_intercept = normal(0,1), prior = normal(0,1), prior_aux = normal(0,1),
                       prior_unstructured = normal(0,1), prior_structured = normal(0,1),
                       family = neg_binomial_2(link = "log"),
                       W = adj, iter = ITER, chains = CHAINS))
SW(fit_gamma <- stan_bym(y_gamma ~ 1 + x, data = spatial_data, family = Gamma(link = "log"),
                         prior_intercept = normal(0,1), prior = normal(0,1), prior_aux = normal(0,1),
                         prior_unstructured = normal(0,1), prior_structured = normal(0,1),
                         W = adj, iter = ITER, chains = CHAINS))

# test QR 
test_that("QR errors when number of predictors is <= 1", {
  expect_error(
    stan_bym(y_gauss ~ x, data = spatial_data, family = gaussian(), W = adj, seed = SEED, QR = TRUE),
    "'QR' can only be specified when there are multiple predictors"
  )
})

test_that("QR works when number of x predictors is > 1", {
  SW(fit_bym <- stan_bym(y_gauss ~ 1 + x + I(x^2), data = spatial_data, family = gaussian(),
                             W = adj, iter = ITER, chains = CHAINS, QR = TRUE))
  expect_stanreg(fit_bym)
})

test_that("stan_bym errors with algorithm = 'optimizing'", {
  expect_error(stan_bym(y_gamma ~ 1 + x, data = spatial_data, family = gaussian(),
                          W = adj, iter = ITER, chains = CHAINS, algorithm = "optimizing"),
               "'arg' should be one of “sampling”, “meanfield”, “fullrank”")
})

test_that("loo/waic for stan_bym works", {
  loo(fit_gauss)
  loo(fit_binom)
  loo(fit_pois)
  loo(fit_nb2)
  loo(fit_gamma)
})

test_that("posterior_predict works for stan_bym", {
  preds_gauss <- posterior_predict(fit_gauss)
  preds_binom <- posterior_predict(fit_binom)
  preds_pois <- posterior_predict(fit_pois)
  preds_nb2 <- posterior_predict(fit_nb2)
  preds_gamma <- posterior_predict(fit_gamma)
})

test_that("predict works for stan_bym", {
  new_dat <- data.frame(x = rnorm(nrow(adj), 2, 1))
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
