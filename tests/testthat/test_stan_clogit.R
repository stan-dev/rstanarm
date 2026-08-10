# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2017 Trustees of Columbia University
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

# this mostly goes through same code as a logit model so only testing the unique stuff

suppressPackageStartupMessages(library(rstanarm))

SEED <- 123
ITER <- 100
CHAINS <- 2
CORES <- 1
REFRESH <- 0

threshold <- 0.03

context("stan_clogit")

SW(fit <- stan_clogit(case ~ spontaneous + induced, strata = stratum, prior = NULL,
                   data = infert[order(infert$stratum), ], 
                   QR = TRUE, init_r = 0.5,
                   chains = CHAINS, iter = ITER, seed = SEED, refresh = 0))

test_that("stan_clogit is similar to survival::clogit", {
  ref_vals <- c(spontaneous = 1.985876, induced = 1.409012)
  # Account for RNG change in new Stan
  if (utils::packageVersion("StanHeaders") >= "2.36") {
    ref_vals <- c(spontaneous = 2.062676, induced = 1.360712)
  }
  expect_equal(ref_vals, coef(fit), tol = threshold)
})

test_that("stan_clogit runs for infert example", {
  expect_stanreg(fit)
})

test_that("stan_clogit works when y is a factor", {
  d <- infert[order(infert$stratum), ]
  d$case <- factor(d$case, labels = c("A", "B"))
  SW(fit_factor <- stan_clogit(case ~ spontaneous + induced, strata = stratum, prior = NULL,
                        data = infert[order(infert$stratum), ], 
                        QR = TRUE, init_r = 0.5,
                        chains = CHAINS, iter = ITER, seed = SEED, refresh = 0))
  expect_equal(coef(fit_factor), coef(fit))
})

test_that("stan_clogit throws error if data are not sorted", {
  expect_error(update(fit, data = infert), 
               regexp = "Data must be sorted")
})

test_that("loo/waic for stan_clogit works", {
  ll_fun <- rstanarm:::ll_fun
  expect_equivalent_loo(fit)
  expect_identical(ll_fun(fit), rstanarm:::.ll_clogit_i)
})

SW(fit_mer <- stan_clogit(case ~ spontaneous + induced + (1 | education),
                          strata = stratum,
                          data = infert[order(infert$stratum), ],
                          QR = TRUE, init_r = 0.5,
                          chains = CHAINS, iter = ITER, seed = SEED, refresh = 0))

# linear predictor built by matching coefficients to columns by name
clogit_eta_ref <- function(object) {
  mat <- as.matrix(object)
  x <- get_x(object)
  eta <- tcrossprod(mat[, colnames(x), drop = FALSE], x)
  b <- grep("^b\\[", colnames(mat))
  if (length(b))
    eta <- eta + tcrossprod(mat[, b, drop = FALSE], as.matrix(get_z(object)))
  eta
}

# one conditional log-likelihood term per stratum; every stratum in infert has
# exactly one case, so the denominator reduces to a log-sum-exp
clogit_ll_ref <- function(object) {
  eta <- clogit_eta_ref(object)
  y <- as.vector(get_y(object))
  g <- droplevels(as.factor(model.frame(object)[, "(weights)"]))
  vapply(levels(g), FUN.VALUE = numeric(nrow(eta)), FUN = function(s) {
    j <- which(g == s)
    e <- eta[, j, drop = FALSE]
    mx <- apply(e, 1, max)
    eta[, j[y[j] == 1]] - (mx + log(rowSums(exp(e - mx))))
  })
}

test_that("stan_clogit with group-specific terms has no intercept column", {
  expect_false("(Intercept)" %in% colnames(get_x(fit_mer)))
  expect_identical(colnames(get_x(fit_mer)), colnames(fit_mer$x))
  expect_identical(colnames(model.matrix(fit_mer)), colnames(fit_mer$x))
})

test_that("log_lik matches coefficients to columns for stan_clogit with group terms", {
  ll <- log_lik(fit_mer)
  g <- droplevels(as.factor(model.frame(fit_mer)[, "(weights)"]))
  expect_equal(ncol(ll), nlevels(g))
  expect_equivalent(ll, clogit_ll_ref(fit_mer))
})

test_that("posterior_linpred matches coefficients to columns for stan_clogit with group terms", {
  expect_equivalent(posterior_linpred(fit_mer), clogit_eta_ref(fit_mer))
  expect_equivalent(posterior_linpred(fit_mer, newdata = infert[order(infert$stratum), ]),
                    clogit_eta_ref(fit_mer))
})

test_that("loo for stan_clogit with group terms gives a sensible p_loo", {
  SW(loo_mer <- loo(fit_mer))
  expect_lt(loo_mer$estimates["p_loo", "Estimate"], 10)
})

context("posterior_predict (stan_clogit)")
test_that("compatible with stan_clogit", {
  PPD1 <- posterior_predict(fit)
  PPD2 <- posterior_predict(fit, newdata = infert) # order irrelevant
  expect_identical(rowSums(PPD1), rowSums(PPD2))
  expect_equal(rowSums(PPD1), round(rowSums(
               posterior_linpred(fit, newdata = infert, transform = TRUE))))
})
