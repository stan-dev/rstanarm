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
library(lme4)
SEED <- 123
set.seed(SEED)
ITER <- 100
CHAINS <- 2
REFRESH <- 0

SW <- suppressWarnings

test_that("posterior_predict returns object with correct classes", {
  expect_s3_class(posterior_predict(example_model), 
                  c("ppd", "matrix"))
})

# Error messages ----------------------------------------------------------
context("posterior_predict (error messages)")
test_that("posterior_predict errors if not a stanreg object", {
  expect_error(posterior_predict(example_model$stanfit), "no applicable method")
  expect_error(posterior_predict(summary(example_model)), "no applicable method")
})
test_that("posterior_predict does not error if model fit using optimization", {
  fit1 <- stan_glm(mpg ~ wt, data = mtcars, algorithm = "optimizing", 
                   seed = SEED, refresh = 0)
  expect_silent(posterior_predict(fit1))
  expect_silent(posterior_linpred(fit1))
})
test_that("posterior_predict errors if NAs in newdata", {
  nd <- model.frame(example_model)
  nd$period[1] <- NA
  expect_error(posterior_predict(example_model, newdata = nd), 
               regexp = "NAs are not allowed in 'newdata'")
  expect_error(posterior_linpred(example_model, newdata = nd), 
               regexp = "NAs are not allowed in 'newdata'")
})
test_that("posterior_predict errors if draws > posterior sample size", {
  expect_error(posterior_predict(example_model, draws = 1e6), 
               regexp = "'draws' should be <= posterior sample size")
})

# VB ----------------------------------------------------------------------
context("posterior_predict ok for vb")
test_that("errors for optimizing and silent for vb", {
  fit1 <- stan_glm(mpg ~ wt + cyl + am, data = mtcars, algorithm = "meanfield", 
                   seed = SEED, refresh = 0)
  fit2 <- update(fit1, algorithm = "fullrank", refresh = 0)
  expect_silent(posterior_predict(fit1))
  expect_silent(posterior_predict(fit2))
  expect_silent(posterior_linpred(fit1))
  expect_silent(posterior_linpred(fit2))
})


# MCMC --------------------------------------------------------------------

test_that("edge cases for posterior_predict work correctly", {
  dims <- c(nrow(as.matrix(example_model)), nrow(lme4::cbpp))
  expect_identical(posterior_predict(example_model, re.form = NA, seed = SEED),
                   posterior_predict(example_model, re.form = ~0, seed = SEED))
  expect_identical(posterior_linpred(example_model, re.form = NA),
                   posterior_linpred(example_model, re.form = ~0))
  expect_identical(posterior_predict(example_model, seed = SEED),
                   posterior_predict(example_model, newdata = lme4::cbpp, seed = SEED))
  expect_identical(posterior_linpred(example_model),
                   posterior_linpred(example_model, newdata = lme4::cbpp))
  expect_error(posterior_predict(example_model, re.form = ~1))
  expect_error(posterior_predict(example_model, re.form = ~(1|foo)))
  expect_error(posterior_linpred(example_model, re.form = ~1))
  expect_error(posterior_linpred(example_model, re.form = ~(1|foo)))
})

test_that("lme4 tests work similarly", {
  # loosely following predict tests from lme4
  
  sfit <- example_model
  nd <- lme4::cbpp
  
  p1 <- posterior_predict(sfit, seed = SEED)
  p1b <- posterior_predict(sfit, newdata = nd, seed = SEED)
  expect_equal(p1, p1b)
  
  p2 <- posterior_predict(sfit, re.form = NA, seed = SEED)
  expect_equal(ncol(p2), nrow(nd))
  
  
  nd2 <- with(nd, expand.grid(period = unique(period), 
                              herd = unique(herd), 
                              size = 20))
  nd2$incidence <- 0
  
  p3 <- posterior_predict(sfit, nd2, seed = SEED)
  p4 <- expect_silent(posterior_predict(sfit, nd2, re.form = NA, seed = SEED))
  p5 <- posterior_predict(sfit, nd2, re.form = ~(1|herd), seed = SEED)
  expect_equal(p3, p5)
  
  # new levels
  nd3 <- rbind(nd2, data.frame(period = as.character(1:4), 
                               herd = rep("new",4), 
                               size = 20, incidence = 0))

  p6 <- posterior_predict(sfit, nd3, allow.new.levels = TRUE, seed = SEED)
  expect_equal(colMeans(p3), colMeans(p6[, 1:ncol(p3)]), tol = 0.05)
  expect_equal(apply(p3, 2, sd), apply(p6[, 1:ncol(p3)], 2, sd), tol = 0.05)
  
  # multiple groups
  lfit <- lmer(diameter ~ (1|plate) + (1|sample), Penicillin)
  sfit <- SW(stan_lmer(diameter ~ (1|plate) + (1|sample), data = Penicillin, 
                    iter = 400, chains = CHAINS, seed = SEED, refresh = 0))
 
  nd <- with(Penicillin, expand.grid(plate=unique(plate), sample=unique(sample)))
  p1 <- posterior_predict(sfit, re.form = NA, seed = SEED)
  p2 <- posterior_predict(sfit, nd, seed = SEED)
  p3 <- posterior_predict(sfit, nd, re.form = NA, seed = SEED)
  p4 <- posterior_predict(sfit, nd, re.form=~(1|plate)+(~1|sample), seed = SEED)
  p4b <- posterior_predict(sfit, nd, re.form=~(1|sample)+(~1|plate), seed = SEED)
  expect_equal(p2,p4,p4b)
  p5 <- posterior_predict(sfit, nd, re.form=~(1|plate), seed = SEED)
})


# spaces in factor levels -------------------------------------------------
context("posterior_linpred/predict with spaces in factor levels")

test_that("posterior_linpred not sensitive to spaces in factor levels", {
  df <- data.frame(
    y = rnorm(10), 
    fac_nospace = gl(2, 5, labels = c("levelone", "leveltwo")), 
    char_nospace = rep(c("levelone", "leveltwo"), each = 5),
    fac_space = gl(2, 5, labels = c("level one", "level two")), 
    char_space = rep(c("level one", "level two"), each = 5),
    fac_mix = gl(2, 5, labels = c("level one", "leveltwo")), 
    char_mix = rep(c("level one", "leveltwo"), each = 5),
    int = rep(1:2, each = 5)
  )
  SW(capture.output(
    fit1 <- stan_lmer(y ~ (1 | fac_nospace), data = df, seed = 123, 
                      chains = 2, iter = 25, refresh = 0),
    fit2 <- update(fit1, formula. = . ~ (1 | char_nospace)),
    fit3 <- update(fit1, formula. = . ~ (1 | fac_space)),
    fit4 <- update(fit1, formula. = . ~ (1 | char_space)),
    fit5 <- update(fit1, formula. = . ~ (1 | fac_mix)),
    fit6 <- update(fit1, formula. = . ~ (1 | char_mix)),
    fit7 <- update(fit1, formula. = . ~ (1 | int))
  ))
  
  # not adding a new level
  nd1 <- df[c(1, 10), ]
  ans1 <- posterior_linpred(fit1, newdata = nd1)
  expect_equal(ans1, posterior_linpred(fit2, newdata = nd1))
  expect_equal(ans1, posterior_linpred(fit3, newdata = nd1))
  expect_equal(ans1, posterior_linpred(fit4, newdata = nd1))
  expect_equal(ans1, posterior_linpred(fit5, newdata = nd1))
  expect_equal(ans1, posterior_linpred(fit6, newdata = nd1))
  expect_equal(ans1, posterior_linpred(fit7, newdata = nd1))
  
  # adding new levels
  nd2 <- data.frame(
    fac_nospace = gl(4, 1, labels = c("levelone", "leveltwo", "levelthree", "levelfour")),
    char_nospace = c("levelone", "leveltwo", "levelthree", "levelfour"),
    fac_space = gl(4, 1, labels = c("level one", "level two", "level three", "level four")),
    char_space = c("level one", "level two", "level three", "level four"), 
    fac_mix = gl(4, 1, labels = c("level one", "leveltwo", "level three", "levelfour")),
    char_mix = c("level one", "leveltwo", "level three", "levelfour"),
    int = 1:4
  )
  ans2 <- posterior_linpred(fit1, newdata = nd2)
  # should be same as ans1 except for cols 3:4 with new levels
  expect_equal(ans2[, 1:2], ans1, check.attributes = FALSE)
  expect_equal(ans2, posterior_linpred(fit2, newdata = nd2))
  expect_equal(ans2, posterior_linpred(fit3, newdata = nd2))
  expect_equal(ans2, posterior_linpred(fit4, newdata = nd2))
  expect_equal(ans2, posterior_linpred(fit5, newdata = nd2))
  expect_equal(ans2, posterior_linpred(fit6, newdata = nd2))
  expect_equal(ans2, posterior_linpred(fit7, newdata = nd2))
})

test_that("posterior_linpred with spaces in factor levels ok with complicated formula", {
  d <- mtcars
  d$cyl_fac <- factor(d$cyl, labels = c("cyl 4", "cyl 6", "cyl 8"))
  d$gear_fac <- factor(d$gear, labels = c("gear 3", "gear 4", "gear 5"))
  
  SW(capture.output(
    fit1 <- stan_lmer(mpg ~ (1 + wt|cyl/gear), data = d,
                      iter = 50, chains = 1, seed = 123, refresh = 0),
    fit2 <- update(fit1, formula. = . ~ (1 + wt|cyl_fac/gear_fac))
  ))
  expect_equal(posterior_linpred(fit1), posterior_linpred(fit2))
  
  # no new levels, all orig levels present in newdata
  nd1 <- data.frame(wt = 2, cyl = d$cyl, gear = d$gear)
  nd2 <- data.frame(wt = 2, cyl_fac = d$cyl_fac, gear_fac = d$gear_fac)
  expect_equal(posterior_linpred(fit1, newdata = nd1), 
               posterior_linpred(fit2, newdata = nd2))
  
  # no new levels, subset of orig levels present in newdata
  nd3 <- data.frame(wt = 2, cyl = 4, gear = 3)
  nd4 <- data.frame(wt = 2, cyl_fac = "cyl 4", gear_fac = factor(3, labels = "gear 3"))
  expect_equal(posterior_linpred(fit1, newdata = nd3), 
               posterior_linpred(fit2, newdata = nd4))
  
  # with new levels
  nd5 <- data.frame(wt = 2, cyl = 98, gear = 99)
  nd6 <- data.frame(wt = 2, cyl_fac = "new cyl", gear_fac = "new gear")
  expect_equal(posterior_linpred(fit1, newdata = nd5), 
               posterior_linpred(fit2, newdata = nd6))
})


# helper functions --------------------------------------------------------
context("posterior_predict helper functions")
test_that("pp_binomial_trials works", {
  ppbt <- rstanarm:::pp_binomial_trials
  
  # binomial
  expect_equal(ppbt(example_model), cbpp$size)
  expect_equal(ppbt(example_model, newdata = cbpp[1:5, ]), cbpp[1:5, "size"])
  
  # bernoulli
  fit <- SW(stan_glm(I(mpg > 25) ~ wt, data = mtcars, family = binomial, 
                     iter = ITER, refresh = 0, chains = CHAINS, 
                     seed = SEED))
  expect_equal(ppbt(fit), rep(1, nrow(mtcars)))
  # expect_equal(ppbt(fit, newdata = mtcars[1:5, ]), rep(1, 5))
})

