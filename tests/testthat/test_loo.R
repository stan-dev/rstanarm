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

library(rstanarm)
LOO.CORES <- ifelse(.Platform$OS.type == "windows", 1, 2)
SEED <- 1234L
set.seed(SEED)
CHAINS <- 2
ITER <- 40 # small iter for speed but large enough for psis
REFRESH <- 0

source(test_path("helpers", "SW.R"))

# loo and waic ------------------------------------------------------------
context("loo and waic")

# These tests just check that the loo.stanreg method (which calls loo.function
# method) results are identical to the loo.matrix results. Since for these tests
# the log-likelihood matrix is computed using the log-likelihood function, the
# only thing these tests really do is make sure that loo.stanreg and all the
# log-likelihood functions don't return any errors and whatnot (it does not
# check that the results returned by loo are actually correct).

expect_equivalent_loo <- function(fit) {
  l <- suppressWarnings(loo(fit, cores = LOO.CORES))
  w <- suppressWarnings(waic(fit))
  expect_s3_class(l, "loo")
  expect_s3_class(w, "loo")
  expect_s3_class(w, "waic")

  att_names <- c("names", "dims", "class", "model_name", "discrete", "yhash", "formula")
  expect_named(attributes(l), att_names)
  expect_named(attributes(w), att_names)

  discrete <- attr(l, "discrete")
  expect_true(!is.na(discrete) && is.logical(discrete))

  llik <- log_lik(fit)
  r <- loo::relative_eff(exp(llik), chain_id = rstanarm:::chain_id_for_loo(fit))
  l2 <- suppressWarnings(loo(llik, r_eff = r, cores = LOO.CORES))
  expect_equal(l$estimates, l2$estimates)
  expect_equivalent(w, suppressWarnings(waic(log_lik(fit))))
}

test_that("loo & waic do something for non mcmc models", {
  SW(fito <- stan_glm(mpg ~ wt, data = mtcars, algorithm = "optimizing",
                      seed = 1234L, prior_intercept = NULL,
                      prior = NULL, prior_aux = NULL))
  capture.output(SW(fitvb1 <- update(fito, algorithm = "meanfield", iter = ITER)))
  capture.output(SW(fitvb2 <- update(fito, algorithm = "fullrank", iter = ITER)))
  expect_true("importance_sampling_loo" %in% class(loo(fito)))
  expect_true("importance_sampling_loo" %in% class(loo(fitvb1)))
  expect_true("importance_sampling_loo" %in% class(loo(fitvb2)))
})

test_that("loo errors if model has weights", {
  SW(
    fit <- stan_glm(mpg ~ wt, data = mtcars,
                    weights = rep_len(c(1,2), nrow(mtcars)),
                    seed = SEED, refresh = 0, iter = 50)
  )
  expect_error(loo(fit), "not supported")
  expect_error(loo(fit), "'kfold'")
})

# loo with refitting ------------------------------------------------------
context("loo then refitting")

test_that("loo issues errors/warnings", {
  expect_warning(loo(example_model, cores = LOO.CORES, k_threshold = 2),
                 "Setting 'k_threshold' > 1 is not recommended")
  expect_error(loo(example_model, k_threshold = -1),
               "'k_threshold' < 0 not allowed.")
  expect_error(loo(example_model, k_threshold = 1:2),
               "'k_threshold' must be a single numeric value")

  expect_warning(rstanarm:::recommend_kfold(5), "Found 5")
  expect_warning(rstanarm:::recommend_kfold(5), "10-fold")
  expect_warning(rstanarm:::recommend_reloo(7), "Found 7")
})

test_that("loo with k_threshold works", {
#  fit <- SW(stan_glm(mpg ~ wt, prior = normal(0, 500), data = mtcars[25:32,],
#                     seed = 12345, iter = 300, chains = 4, cores = 1,
#                     refresh = 0))
#  expect_warning(loo_x <- loo(fit, k_threshold = 0.5),
#                 "We recommend calling 'loo' again")
#  expect_message(rstanarm:::reloo(fit, loo_x, obs = 1:10, refit = FALSE),
#                 "Model will be refit 10 times")
#  expect_output(SW(rstanarm:::reloo(fit, loo_x, obs = 1, refit = TRUE)),
#                "Elapsed Time")

#  # test that no errors from binomial model because it's trickier to get the
#  # data right internally in reloo (matrix outcome)
#  loo_x <- loo(example_model)
#  expect_message(SW(rstanarm:::reloo(example_model, loo_x, obs = 1)),
#                 "Model will be refit 1 times")
})

test_that("loo with k_threshold works for edge case(s)", {
  # without 'data' argument
  y <- mtcars$mpg
  x <- rexp(length(y))
  SW(fit <- stan_glm(y ~ 1, refresh = 0, iter = 200))
  expect_message(
    SW(res <- loo(fit, k_threshold = 0.1, cores = LOO.CORES)), # low k_threshold to make sure reloo is triggered
    "problematic observation\\(s\\) found"
  )
  expect_s3_class(res, "loo")
})



# kfold -------------------------------------------------------------------
context("kfold")

test_that("kfold does not throw an error for non mcmc models", {
  SW(
    fito <- stan_glm(mpg ~ wt, data = mtcars, algorithm = "optimizing",
                     seed = 1234L, refresh = 0)
  )
  expect_true("kfold" %in% class(kfold(fito)))
})

test_that("kfold throws error if K <= 1 or K > N", {
  expect_error(kfold(example_model, K = 1), "K > 1", fixed = TRUE)
  expect_error(kfold(example_model, K = 1e5), "K <= nobs(x)", fixed = TRUE)
})

test_that("kfold throws error if folds arg is bad", {
  expect_error(kfold(example_model, K = 2, folds = 1:100), "length(folds) == N is not TRUE", fixed = TRUE)
  expect_error(kfold(example_model, K = 2, folds = 1:2), "length(folds) == N is not TRUE", fixed = TRUE)
  expect_error(kfold(example_model, K = 2, folds = seq(1,100, length.out = 56)), "all(folds == as.integer(folds)) is not TRUE", fixed = TRUE)
})

test_that("kfold throws error if model has weights", {
  SW(
    fit <- stan_glm(mpg ~ wt, data = mtcars,
                    iter = ITER, chains = CHAINS, refresh = 0,
                    weights = runif(nrow(mtcars), 0.5, 1.5))
  )
  expect_error(kfold(fit), "not currently available for models fit using weights")
})

test_that("kfold works on some examples", {
  mtcars2 <- mtcars
  mtcars2$wt[1] <- NA # make sure kfold works if NAs are dropped from original data
  SW(
    fit_gaus <- stan_glm(mpg ~ wt, data = mtcars2, refresh = 0,
                         chains = 1, iter = 200)
  )
  SW(kf <- kfold(fit_gaus, 4))
  SW(kf2 <- kfold(example_model, 2))

  expect_named(kf, c("estimates", "pointwise", "elpd_kfold", "se_elpd_kfold"))
  expect_named(kf2, c("estimates", "pointwise", "elpd_kfold", "se_elpd_kfold"))
  expect_named(attributes(kf),  c("names", "class", "K", "dims", "model_name", "discrete", "yhash", "formula"))
  expect_named(attributes(kf2), c("names", "class", "K", "dims", "model_name", "discrete", "yhash", "formula"))
  expect_s3_class(kf, c("kfold", "loo"))
  expect_s3_class(kf2, c("kfold", "loo"))

  SW(kf <- kfold(fit_gaus, K = 2, save_fits = TRUE))

  expect_true("fits" %in% names(kf))
  expect_s3_class(kf$fits[[1, "fit"]], "stanreg")
  expect_type(kf$fits[[2, "omitted"]], "integer")
  expect_length(kf$fits[[2, "omitted"]], 16)
})

# loo_compare ----------------------------------------------------------
test_that("loo_compare throws correct errors", {
  SW(capture.output({
    mtcars$mpg <- as.integer(mtcars$mpg)
    fit1 <- stan_glm(mpg ~ wt, data = mtcars, iter = 5, chains = 2, refresh = 0)
    fit2 <- update(fit1, data = mtcars[-1, ])
    fit3 <- update(fit1, formula. = log(mpg) ~ .)
    fit4 <- update(fit1, family = poisson("log"))

    l1 <- loo(fit1, cores = LOO.CORES)
    l2 <- loo(fit2, cores = LOO.CORES)
    l3 <- loo(fit3, cores = LOO.CORES)
    l4 <- loo(fit4, cores = LOO.CORES)

    w1 <- waic(fit1)
    k1 <- kfold(fit1, K = 3)
  }))


  # this uses loo::loo_compare
  expect_error(loo_compare(l1, l2),
               "Not all models have the same number of data points")
  expect_error(loo_compare(list(l4, l2, l3)),
               "Not all models have the same number of data points")
  
  # using loo_compare.stanreg (can do extra checks)
  fit1$loo <- l1
  fit2$loo <- l2
  fit3$loo <- l3
  fit4$loo <- l4
  
  expect_error(loo_compare(fit1, fit2), "Not all models have the same number of data points")
  expect_warning(loo_compare(fit1, fit3), "Not all models have the same y variable")
  expect_error(loo_compare(fit1, fit4),
               "Discrete and continuous observation models can't be compared")


  expect_error(loo_compare(l1, fit1),
               "All inputs should have class 'loo'")
  expect_error(loo_compare(l1),
               "requires at least two models")
})

test_that("loo_compare works", {
  suppressWarnings(capture.output({
    mtcars$mpg <- as.integer(mtcars$mpg)
    fit1 <- stan_glm(mpg ~ wt, data = mtcars, iter = 40, chains = 2, refresh = 0)
    fit2 <- update(fit1, formula. = . ~ . + cyl)
    fit3 <- update(fit2, formula. = . ~ . + gear)
    fit4 <- update(fit1, family = "poisson")
    fit5 <- update(fit1, family = "neg_binomial_2")

    fit1$loo <- loo(fit1, cores = LOO.CORES)
    fit2$loo <- loo(fit2, cores = LOO.CORES)
    fit3$loo <- loo(fit3, cores = LOO.CORES)
    fit4$loo <- loo(fit4, cores = LOO.CORES)
    fit5$loo <- loo(fit5, cores = LOO.CORES)

    k1 <- kfold(fit1, K = 2)
    k2 <- kfold(fit2, K = 2)
    k3 <- kfold(fit3, K = 3)
    k4 <- kfold(fit4, K = 2)
    k5 <- kfold(fit5, K = 2)
  }))

  expect_false(attr(fit1$loo, "discrete"))
  expect_false(attr(fit2$loo, "discrete"))
  expect_false(attr(fit3$loo, "discrete"))
  expect_true(attr(fit4$loo, "discrete"))
  expect_true(attr(fit5$loo, "discrete"))

  comp1 <- loo_compare(fit1, fit2)
  comp2 <- loo_compare(fit1, fit2, fit3)
  expect_s3_class(comp1, "compare.loo")
  expect_s3_class(comp2, "compare.loo")
  expect_equal(comp1[, "elpd_diff"], loo_compare(list(fit1$loo, fit2$loo))[, "elpd_diff"])
  expect_equal(comp2[, "elpd_diff"], loo_compare(list(fit1$loo, fit2$loo, fit3$loo))[, "elpd_diff"])
  
  comp1_detail <- loo_compare(fit1, fit2, detail=TRUE)
  expect_output(print(comp1_detail), "Model formulas")
  
  # equivalent to stanreg_list method
  expect_equivalent(comp2, loo_compare(stanreg_list(fit1, fit2, fit3)))

  # for kfold
  comp3 <- loo_compare(k1, k2, k3)
  expect_true(attr(k4, "discrete"))
  expect_true(attr(k5, "discrete"))
  expect_s3_class(loo_compare(k4, k5), "compare.loo")
})


# helpers -----------------------------------------------------------------
context("loo and waic helpers")

test_that("kfold_and_reloo_data works", {
  f <- rstanarm:::kfold_and_reloo_data
  d <- f(example_model)
  expect_identical(d, lme4::cbpp[, colnames(d)])

  # if 'data' arg not originally specified when fitting the model
  y <- rnorm(40)
  SW(fit <- stan_glm(y ~ 1, iter = ITER, chains = CHAINS, refresh = 0))
  expect_equivalent(f(fit), model.frame(fit))

  # if 'subset' arg specified when fitting the model
  SW(fit2 <- stan_glm(mpg ~ wt, data = mtcars, subset = gear != 5, iter = ITER,
                      chains = CHAINS, refresh = 0))
  expect_equivalent(f(fit2), subset(mtcars[mtcars$gear != 5, c("mpg", "wt")]))
})

test_that(".weighted works", {
  f <- rstanarm:::.weighted
  expect_equal(f(2, NULL), 2)
  expect_equal(f(2, 3), 6)
  expect_equal(f(8, 0.25), 2)
  expect_error(f(2), "missing, with no default")
})
