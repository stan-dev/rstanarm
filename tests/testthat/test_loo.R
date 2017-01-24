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

library(rstanarm)
options(loo.cores = 2)
SEED <- 1234
set.seed(SEED)
CHAINS <- 2
ITER <- 40 # small iter for speed but large enough for psis
REFRESH <- 0

SW <- function(expr) capture.output(suppressWarnings(expr))

# loo and waic ------------------------------------------------------------
context("loo and waic")

# These tests just check that the loo.stanreg method (which calls loo.function
# method) results are identical to the loo.matrix results. Since for these tests
# the log-likelihood matrix is computed using the log-likelihood function, the
# only thing these tests really do is make sure that loo.stanreg and all the
# log-likelihood functions don't return any errors and whatnot (it does not
# check that the results returned by loo are actually correct).

expect_equivalent_loo <- function(fit) {
  l <- suppressWarnings(loo(fit))
  w <- suppressWarnings(waic(fit))
  expect_s3_class(l, "loo")
  expect_s3_class(w, "loo")
  expect_s3_class(w, "waic")
  
  att_names <- c("names", "log_lik_dim", "class", "name", "discrete", "yhash")
  expect_named(attributes(l), att_names)
  expect_named(attributes(w), att_names)
  
  discrete <- attr(l, "discrete")
  expect_true(!is.na(discrete) && is.logical(discrete))
  
  expect_equivalent(l, suppressWarnings(loo(log_lik(fit))))
  expect_equivalent(w, suppressWarnings(waic(log_lik(fit))))
}
mcmc_only_error <- function(fit) {
  msg <- "only available for models fit using MCMC"
  expect_error(loo(fit), regexp = msg)
  expect_error(waic(fit), regexp = msg)
}

test_that("loo & waic throw error for non mcmc models", {
  SW(fito <- stan_glm(mpg ~ wt, data = mtcars, algorithm = "optimizing",
                     seed = SEED))
  SW(fitvb1 <- update(fito, algorithm = "meanfield"))
  SW(fitvb2 <- update(fito, algorithm = "fullrank"))
  mcmc_only_error(fito)
  mcmc_only_error(fitvb1)
  mcmc_only_error(fitvb2)
})

test_that("loo errors if model has weights", {
  SW(
    fit <- stan_glm(mpg ~ wt, data = mtcars, 
                    weights = rep_len(c(1,2), nrow(mtcars)),
                    seed = SEED, refresh = REFRESH, iter = 50)
  )
  expect_error(loo(fit), "not supported")
  expect_error(loo(fit), "'kfold'")
})

test_that("loo/waic for stan_glm works", {
  # gaussian
  SW(
    fit_gaus <- stan_glm(mpg ~ wt, data = mtcars, 
                         chains = CHAINS, iter = ITER,
                         seed = SEED, refresh = REFRESH)
  )
  expect_equivalent_loo(fit_gaus)
  expect_identical(ll_fun(fit_gaus), rstanarm:::.ll_gaussian_i)

  # binomial
  dat <- data.frame(ldose = rep(0:5, 2),
                    sex = factor(rep(c("M", "F"), c(6, 6))))
  numdead <- c(1, 4, 9, 13, 18, 20, 0, 2, 6, 10, 12, 16)
  SF <- cbind(numdead, numalive = 20-numdead)
  SW(
    fit_binom <- stan_glm(SF ~ sex*ldose, data = dat, family = binomial,
                          chains = CHAINS, iter = ITER, seed = SEED,
                          refresh = REFRESH)
  )
  dead <- rbinom(length(numdead), 1, prob = 0.5)
  SW(fit_binom2 <- update(fit_binom, formula = factor(dead) ~ .))
  expect_equivalent_loo(fit_binom)
  expect_equivalent_loo(fit_binom2)
  expect_identical(ll_fun(fit_binom), rstanarm:::.ll_binomial_i)
  expect_identical(ll_fun(fit_binom2), rstanarm:::.ll_binomial_i)

  # poisson
  d.AD <- data.frame(treatment = gl(3,3), outcome =  gl(3,1,9),
                     counts = c(18,17,15,20,10,20,25,13,12))
  SW(fit_pois <- stan_glm(counts ~ outcome + treatment, data = d.AD,
                          family = poisson, chains = CHAINS, iter = ITER,
                          seed = SEED, refresh = REFRESH))
  expect_equivalent_loo(fit_pois)
  expect_identical(ll_fun(fit_pois), rstanarm:::.ll_poisson_i)

  # negative binomial
  SW(fit_negbin <- update(fit_pois, family = neg_binomial_2))
  expect_equivalent_loo(fit_negbin)
  expect_identical(ll_fun(fit_negbin), rstanarm:::.ll_neg_binomial_2_i)

  # gamma
  clotting <- data.frame(log_u = log(c(5,10,15,20,30,40,60,80,100)),
                         lot1 = c(118,58,42,35,27,25,21,19,18),
                         lot2 = c(69,35,26,21,18,16,13,12,12))
  SW(fit_gamma <- stan_glm(lot1 ~ log_u, data = clotting, family = Gamma,
                           chains = CHAINS, iter = ITER, seed = SEED,
                           refresh = REFRESH))
  expect_equivalent_loo(fit_gamma)
  expect_identical(ll_fun(fit_gamma), rstanarm:::.ll_Gamma_i)

  # inverse gaussian
  SW(fit_igaus <- update(fit_gamma, family = inverse.gaussian))
  expect_equivalent_loo(fit_igaus)
  expect_identical(ll_fun(fit_igaus), rstanarm:::.ll_inverse.gaussian_i)
})

test_that("loo/waic for stan_polr works", {
  SW(fit_ord_logistic <- stan_polr(tobgp ~ agegp, data = esoph,
                               prior = R2(0.2, "mean"), init_r = 0.1,
                               chains = CHAINS, iter = ITER, seed = SEED,
                               refresh = REFRESH))
  expect_equivalent_loo(fit_ord_logistic)
  expect_identical(ll_fun(fit_ord_logistic), rstanarm:::.ll_polr_i)

  SW(fit_probit <- stan_polr(factor(tobgp == "30+") ~ agegp + alcgp,
                             data = esoph, prior = R2(location = 0.4),
                             method = "probit", chains = CHAINS, iter = ITER,
                             seed = SEED, refresh = REFRESH))
  expect_equivalent_loo(fit_probit)
  expect_identical(ll_fun(fit_probit), rstanarm:::.ll_binomial_i)

  SW(fit_scobit <- stan_polr(factor(tobgp == "30+") ~ agegp + alcgp,
                             data = esoph, prior = R2(location = 0.4),
                             shape = 2, rate = 2, chains = CHAINS, iter = ITER,
                             seed = SEED, refresh = REFRESH))
  expect_equivalent_loo(fit_scobit)
  expect_identical(ll_fun(fit_scobit), rstanarm:::.ll_polr_i)
})

test_that("loo/waic for stan_lm works", {
  SW(fit_lm <- stan_lm(mpg ~ ., data = mtcars, prior = R2(0.75),
                       chains = CHAINS, iter = ITER, seed = SEED,
                       refresh = REFRESH))
  expect_equivalent_loo(fit_lm)
  expect_identical(ll_fun(fit_lm), rstanarm:::.ll_gaussian_i)
})

test_that("loo/waic for stan_glmer works", {
  # gaussian
  SW(fit_glmer1 <- stan_glmer(mpg ~ wt + (1|cyl) + (1+wt|gear), data = mtcars,
                              chains = CHAINS, iter = ITER, seed = SEED,
                              refresh = REFRESH))
  expect_equivalent_loo(fit_glmer1)
  expect_identical(ll_fun(fit_glmer1), rstanarm:::.ll_gaussian_i)

  # binomial
  expect_equivalent_loo(example_model)
  expect_identical(ll_fun(example_model), rstanarm:::.ll_binomial_i)
})

test_that("loo/waic for stan_betareg works", {
  data("GasolineYield", package = "betareg")
  SW(fit_logit <- stan_betareg(yield ~ batch + temp | temp, data = GasolineYield,
                               link = "logit",
                               chains = CHAINS, iter = ITER,
                               seed = SEED, refresh = REFRESH))
  expect_equivalent_loo(fit_logit)
  expect_identical(ll_fun(fit_logit), rstanarm:::.ll_beta_i)
})


# loo with refitting ------------------------------------------------------
context("loo then refitting")

test_that("loo issues errors/warnings", {
  expect_warning(loo(example_model, k_threshold = 2),
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
  
  # test that no errors from binomial model because it's trickier to get the
  # data right internally in reloo (matrix outcome)
#  loo_x <- loo(example_model)
#  expect_message(SW(rstanarm:::reloo(example_model, loo_x, obs = 1)), 
#                 "Model will be refit 1 times")
})

test_that("loo with k_threshold works for edge case(s)", {
  # without 'data' argument
  y <- mtcars$mpg
  SW(fit <- stan_glm(y ~ 1, refresh = 0, iter = 200))
  expect_message(
    SW(res <- loo(fit, k_threshold = 0.1)), # low k_threshold to make sure reloo is triggered
    "problematic observation\\(s\\) found"
  )
  expect_s3_class(res, "loo")
})



# kfold -------------------------------------------------------------------
context("kfold")

test_that("kfold throws error for non mcmc models", {
  SW(
    fito <- stan_glm(mpg ~ wt, data = mtcars, algorithm = "optimizing",seed = SEED)
  )
  expect_error(kfold(fito), "MCMC")
})

test_that("kfold throws error if K <= 1 or K > N", {
  expect_error(kfold(example_model, K = 1), "K > 1", fixed = TRUE)
  expect_error(kfold(example_model, K = 1e5), "K <= nobs(x)", fixed = TRUE)
})

test_that("kfold throws error if model has weights", {
  SW(
    fit <- stan_glm(mpg ~ wt, data = mtcars, 
                    iter = ITER, chains = CHAINS, refresh = REFRESH,
                    weights = runif(nrow(mtcars), 0.5, 1.5))
  )
  expect_error(kfold(fit), "not currently available for models fit using weights")
})

test_that("kfold works on some examples", {
  mtcars2 <- mtcars
  mtcars2$wt[1] <- NA # make sure kfold works if NAs are dropped from original data
  SW(
    fit_gaus <- stan_glm(mpg ~ wt, data = mtcars2, seed = 12345, refresh = 0, 
                         chains = 1, iter = 100)
  )
  SW(kf <- kfold(fit_gaus, 4))
  SW(kf2 <- kfold(example_model, 2))
  
  expect_named(attributes(kf), c("names", "class", "K", "name", "discrete", "yhash"))
  expect_s3_class(kf, c("kfold", "loo"))
  expect_identical(invisible(print(kf)), kf)
  expect_output(print(kf), "4-fold cross-validation")

  expect_s3_class(kf2, c("kfold", "loo"))
  expect_identical(invisible(print(kf2)), kf2)
  expect_output(print(kf2), "2-fold cross-validation")
})



# compare_models ----------------------------------------------------------
test_that("compare_models throws correct errors", {
  SW(capture.output({
    mtcars$mpg <- as.integer(mtcars$mpg)
    fit1 <- stan_glm(mpg ~ wt, data = mtcars, iter = 40, chains = 2, refresh = -1)
    fit2 <- update(fit1, data = mtcars[-1, ])
    fit3 <- update(fit1, formula. = log(mpg) ~ .)
    fit4 <- update(fit1, family = poisson("log"))
    
    l1 <- loo(fit1)
    l2 <- loo(fit2)
    l3 <- loo(fit3)
    l4 <- loo(fit4)
  
    w1 <- waic(fit1)
    k1 <- kfold(fit1, K = 3)
  }))
  
  
  expect_error(compare_models(l1, l2), 
               "Not all models have the same y variable")
  expect_error(compare_models(l1, l3), 
               "Not all models have the same y variable")
  expect_error(compare_models(loos = list(l4, l2, l3)), 
               "Not all models have the same y variable")
  expect_error(compare_models(l1, l4), 
               "Discrete and continuous observation models can't be compared")
  
  
  expect_error(compare_models(l1, fit1), 
               "All objects must have class 'loo'")
  expect_error(compare_models(l1, k1),
               "Can't mix objects computed using 'loo', 'waic', and 'kfold'.")
  expect_error(compare_models(k1, w1, k1, w1), 
               "Can't mix objects computed using 'loo', 'waic', and 'kfold'.")
  
  expect_error(compare_models(l1, loos = list(l2, l3)), 
               "'...' and 'loos' can't both be specified")
  expect_error(compare_models(l1), 
               "At least two objects are required for model comparison")
})

test_that("compare_models works", {
  SW(capture.output({
    mtcars$mpg <- as.integer(mtcars$mpg)
    fit1 <- stan_glm(mpg ~ wt, data = mtcars, iter = 40, chains = 2, refresh = -1)
    fit2 <- update(fit1, formula. = . ~ . + cyl)
    fit3 <- update(fit2, formula. = . ~ . + gear)
    fit4 <- update(fit1, family = "poisson")
    fit5 <- update(fit1, family = "neg_binomial_2")
    
    l1 <- loo(fit1)
    l2 <- loo(fit2)
    l3 <- loo(fit3)
    l4 <- loo(fit4)
    l5 <- loo(fit5)
    
    k1 <- kfold(fit1, K = 2)
    k2 <- kfold(fit2, K = 2)
    k3 <- kfold(fit3, K = 3)
    k4 <- kfold(fit4, K = 2)
    k5 <- kfold(fit5, K = 2)
  }))
  
  expect_false(attr(l1, "discrete"))
  expect_false(attr(l2, "discrete"))
  expect_false(attr(l3, "discrete"))
  
  comp1 <- compare_models(l1, l2)
  comp2 <- compare_models(l1, l2, l3)
  expect_named(comp1, c("elpd_diff", "se"))
  expect_true(is.matrix(comp2))
  expect_equal(ncol(comp2), 6)
  expect_s3_class(comp1, "compare.loo")
  expect_s3_class(comp2, "compare.loo")
  expect_identical(comp1, compare_models(loos = list(l1, l2)))
  expect_identical(comp2, compare_models(loos = list(l1, l2, l3)))
  
  comp3 <- compare_models(k1, k2, k3)
  expect_equal(ncol(comp3), 2)
  expect_s3_class(comp3, "compare.loo")
  
  expect_true(attr(l4, "discrete"))
  expect_true(attr(l5, "discrete"))
  expect_silent(comp4 <- compare_models(l4, l5))
  expect_silent(compare_models(loos = list(l4, l5)))
  expect_s3_class(comp4, "compare.loo")
  
  expect_true(attr(k4, "discrete"))
  expect_true(attr(k5, "discrete"))
  expect_s3_class(compare_models(k4, k5), "compare.loo")
})


# helpers -----------------------------------------------------------------
context("loo and waic helpers")

test_that("kfold_and_reloo_data works", {
  f <- rstanarm:::kfold_and_reloo_data
  d <- f(example_model)
  mf <- model.frame(example_model)
  expect_equal(colnames(d), c("incidence", "size", "period", "herd"))
  expect_equal(colnames(mf), c("cbind(incidence, size - incidence)", "size", "period", "herd"))
  mf2 <- data.frame(incidence = mf[, 1][, 1], mf[, -1])
  expect_equal(d, mf2)
  
  # if 'data' arg not originally specified when fitting the model
  y <- rnorm(40)
  SW(fit <- stan_glm(y ~ 1, iter = ITER, chains = CHAINS, refresh = REFRESH))
  expect_equivalent(f(fit), model.frame(fit))
})

test_that(".weighted works", {
  f <- rstanarm:::.weighted
  expect_equal(f(2, NULL), 2)
  expect_equal(f(2, 3), 6)
  expect_equal(f(8, 0.25), 2)
  expect_error(f(2), "missing, with no default")
})
