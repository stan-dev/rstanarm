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

# These tests just make sure that posterior_predict doesn't throw errors and
# that result has correct dimensions
check_for_error <- function(fit, data = NULL, offset = NULL) {
  nsims <- nrow(as.data.frame(fit))
  mf <- if (!is.null(data)) 
    data else model.frame(fit)
  if (identical(deparse(substitute(fit)), "example_model"))
    mf <- lme4::cbpp
  
  expect_silent(yrep1 <- posterior_predict(fit))
  expect_silent(lin1 <- posterior_linpred(fit))
  expect_silent(posterior_linpred(fit, transform = TRUE))
  expect_equal(dim(yrep1), c(nsims, nobs(fit)))
  expect_equal(dim(lin1), c(nsims, nobs(fit)))

  expect_silent(yrep2 <- posterior_predict(fit, draws = 1))
  expect_equal(dim(yrep2), c(1, nobs(fit)))
  
  offs <- if (!is.null(offset)) offset[1] else offset
  expect_silent(yrep3 <- posterior_predict(fit, newdata = mf[1,], offset = offs))
  expect_silent(lin3 <- posterior_linpred(fit, newdata = mf[1,], offset = offs))
  expect_equal(dim(yrep3), c(nsims, 1))
  expect_equal(dim(lin3), c(nsims, 1))
  
  expect_silent(yrep4 <- posterior_predict(fit, draws = 2, newdata = mf[1,], offset = offs))
  expect_equal(dim(yrep4), c(2, 1))
  
  offs <- if (!is.null(offset)) offset[1:5] else offset
  expect_silent(yrep5 <- posterior_predict(fit, newdata = mf[1:5,], offset = offs))
  expect_silent(lin5 <- posterior_linpred(fit, newdata = mf[1:5,], offset = offs))
  expect_equal(dim(yrep5), c(nsims, 5))
  expect_equal(dim(lin5), c(nsims, 5))
  
  expect_silent(yrep6 <- posterior_predict(fit, draws = 3, newdata = mf[1:5,], offset = offs))
  expect_equal(dim(yrep6), c(3, 5))
  
  expect_error(posterior_predict(fit, draws = nsims + 1), 
               regexep = "posterior sample size is only")
}

expect_linpred_equal <- function(object, tol = 0.1) {
  linpred <- posterior_linpred(object)
  expect_equal(apply(linpred, 2, median), object$linear.predictors, 
               tolerance = tol, 
               check.attributes = FALSE)
}

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
test_that("posterior_predict errors if model fit using optimization", {
  fit1 <- stan_glm(mpg ~ wt + cyl + am, data = mtcars, algorithm = "optimizing", 
                   seed = SEED)
  expect_error(posterior_predict(fit1), regexp = "optimizing")
  expect_error(posterior_linpred(fit1), regexp = "optimizing")
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
                   seed = SEED)
  fit2 <- update(fit1, algorithm = "fullrank")
  expect_silent(posterior_predict(fit1))
  expect_silent(posterior_predict(fit2))
  expect_silent(posterior_linpred(fit1))
  expect_silent(posterior_linpred(fit2))
})


# MCMC --------------------------------------------------------------------
context("posterior_predict (stan_lm)")
test_that("posterior_predict compatible with stan_lm", {
  fit <- SW(stan_lm(mpg ~ wt + cyl + am, data = mtcars, prior = R2(log(0.5), what = "log"),
                 iter = ITER, chains = CHAINS,  seed = SEED, refresh = REFRESH))
  check_for_error(fit)
  expect_linpred_equal(fit)
})

context("posterior_predict (stan_glm)")
test_that("compatible with gaussian glm", {
  fit <- SW(stan_glm(mpg ~ wt, data = mtcars, 
                     iter = ITER, chains = CHAINS, seed = SEED, refresh = REFRESH))
  check_for_error(fit)
  expect_linpred_equal(fit)
})
test_that("compatible with glm with offset", {
  mtcars2 <- mtcars
  mtcars2$offs <- runif(nrow(mtcars))
  fit <- SW(stan_glm(mpg ~ wt, data = mtcars2, offset = offs,
                     iter = ITER, chains = CHAINS, seed = SEED, refresh = REFRESH))
  fit2 <- SW(stan_glm(mpg ~ wt + offset(offs), data = mtcars2,
                      iter = ITER, chains = CHAINS, seed = SEED, refresh = REFRESH))
  
  expect_warning(posterior_predict(fit, newdata = mtcars[1:5, ]), 
                 "offset")
  check_for_error(fit, data = mtcars2, offset = mtcars2$offs)
  check_for_error(fit2, data = mtcars2, offset = mtcars2$offs)
  expect_linpred_equal(fit)
  expect_linpred_equal(fit2)
})
test_that("compatible with poisson & negbin glm", {
  counts <- c(18,17,15,20,10,20,25,13,12)
  outcome <- gl(3,1,9)
  treatment <- gl(3,3)
  fit <- SW(stan_glm(counts ~ outcome + treatment, family = poisson(), 
                     iter = ITER, chains = CHAINS, seed = SEED, 
                     refresh = REFRESH))
  fitnb <- SW(update(fit, family = neg_binomial_2))
  check_for_error(fit)
  check_for_error(fitnb)
  expect_linpred_equal(fit)
  expect_linpred_equal(fitnb)
})
test_that("posterior_predict compatible with gamma & inverse.gaussian glm", {
  clotting <- data.frame(log_u = log(c(5,10,15,20,30,40,60,80,100)),
                         lot1 = c(118,58,42,35,27,25,21,19,18),
                         lot2 = c(69,35,26,21,18,16,13,12,12))
  fit <- SW(stan_glm(lot1 ~ log_u, data = clotting, family = Gamma, 
                  chains = CHAINS, iter = ITER,  seed = SEED, refresh = REFRESH))
  fit_igaus <- SW(update(fit, family = inverse.gaussian))
  
  check_for_error(fit)
  check_for_error(fit_igaus)
  expect_linpred_equal(fit)
  expect_linpred_equal(fit_igaus)
})

context("posterior_predict (stan_polr)")
test_that("compatible with stan_polr", {
  fit <- SW(stan_polr(tobgp ~ agegp + alcgp, data = esoph, prior = R2(location = 0.4),
                   iter = ITER, chains = CHAINS, seed = SEED, refresh = REFRESH))
  
  esoph$tobgp_fac <- factor(esoph$tobgp == "30+")
  fit_binary <- SW(stan_polr(tobgp_fac ~ agegp + alcgp, 
                             data = esoph, prior = R2(location = 0.4), 
                             chains = CHAINS, iter = ITER, 
                             seed = SEED, refresh = REFRESH))
  fit_binary_scobit <- SW(update(fit_binary, shape = 2, rate = 2))
  
  check_for_error(fit)
  check_for_error(fit_binary)
  check_for_error(fit_binary_scobit)
})

context("posterior_predict (stan_gamm4)")
test_that("stan_gamm4 returns expected result for sleepstudy example", {
  fit <- SW(stan_gamm4(Reaction / 10 ~ s(Days), data = sleepstudy,
                       random = ~(1|Subject), chains = CHAINS, iter = ITER, 
                       seed = SEED, refresh = REFRESH))
  expect_silent(yrep1 <- posterior_predict(fit))
  # expect_equal(dim(yrep1), c(nrow(as.data.frame(fit)), nobs(fit)))
  expect_silent(yrep2 <- posterior_predict(fit, draws = 1))
  # expect_equal(dim(yrep2), c(1, nobs(fit)))
  expect_silent(posterior_predict(fit, newdata = sleepstudy))
})


context("posterior_predict (stan_(g)lmer)")
test_that("compatible with stan_lmer", {
  fit <- SW(stan_lmer(mpg ~ wt + (1|cyl) + (1 + wt|gear), data = mtcars, 
                      prior = normal(0,1), iter = ITER, chains = CHAINS,
                      seed = SEED, refresh = REFRESH))
  check_for_error(fit)
  expect_linpred_equal(fit)
})
test_that("compatible with stan_glmer (binomial)", {
  check_for_error(example_model)
  expect_linpred_equal(example_model)
  predprob <- posterior_linpred(example_model, transform = TRUE)
  expect_true(all(predprob > 0) && all(predprob < 1))
})
test_that("compatible with stan_(g)lmer with transformation in formula", {
  d <- mtcars
  d$cyl <- as.factor(d$cyl)
  args <- list(formula = mpg ~ log1p(wt) + (1|cyl) + (1|gear), data = d, 
               iter = ITER, chains = CHAINS,  seed = SEED, refresh = REFRESH)
  fit1 <- SW(do.call("stan_lmer", args))
  fit2 <- SW(do.call("stan_glmer", args))
  nd <- d[6:10, ]
  nd$wt <- runif(5)
  expect_silent(posterior_predict(fit1))
  expect_silent(posterior_predict(fit2))
  expect_silent(posterior_predict(fit1, newdata = nd))
  expect_silent(posterior_predict(fit2, newdata = nd))
  
  expect_silent(posterior_linpred(fit1))
  expect_silent(posterior_linpred(fit2))
  expect_silent(posterior_linpred(fit1, newdata = nd))
  expect_silent(posterior_linpred(fit2, newdata = nd))
})

test_that("compatible with stan_lmer with offset", {
  offs <- rnorm(nrow(mtcars))
  fit <- SW(stan_lmer(mpg ~ wt + (1|cyl) + (1 + wt|gear), data = mtcars, 
                      prior = normal(0,1), iter = ITER, chains = CHAINS,
                      seed = SEED, refresh = REFRESH, offset = offs))
  
  expect_warning(posterior_predict(fit, newdata = mtcars[1:2, ], offset = offs),
                 "STATS")
  check_for_error(fit, offset = offs)
})

context("posterior_predict (stan_betareg)")
test_that("compatible with stan_betareg with z", {
  data("GasolineYield", package = "betareg")
  fit <- SW(stan_betareg(yield ~ pressure + temp | temp, data = GasolineYield,
                         iter = ITER*5, chains = 2*CHAINS, seed = SEED, 
                         refresh = REFRESH))
  check_for_error(fit)
  expect_linpred_equal(fit)
})
test_that("compatible with stan_betareg without z", {
  data("GasolineYield", package = "betareg")
  fit <- SW(stan_betareg(yield ~ temp, data = GasolineYield, 
                     iter = ITER, chains = CHAINS, seed = SEED, refresh = REFRESH))
  check_for_error(fit)
  expect_linpred_equal(fit)
})
test_that("compatible with betareg with offset", {
  GasolineYield2 <- GasolineYield
  GasolineYield2$offs <- runif(nrow(GasolineYield2))
  fit <- SW(stan_betareg(yield ~ temp, data = GasolineYield2, offset = offs,
                     iter = ITER*5, chains = CHAINS, seed = SEED, refresh = REFRESH))
  fit2 <- SW(stan_betareg(yield ~ temp + offset(offs), data = GasolineYield2,
                      iter = ITER*5, chains = CHAINS, seed = SEED, refresh = REFRESH))
  
  expect_warning(posterior_predict(fit, newdata = GasolineYield), 
                 "offset")
  check_for_error(fit, data = GasolineYield2, offset = GasolineYield2$offs)
  check_for_error(fit2, data = GasolineYield2, offset = GasolineYield2$offs)
  expect_linpred_equal(fit)
  expect_linpred_equal(fit2)
})

# compare to lme4 ---------------------------------------------------------
context("posterior_predict (compare to lme4)")
test_that("posterior_predict close to predict.merMod for gaussian", {
  mod1 <- as.formula(mpg ~ wt + (1|cyl) + (1|gear))
  mod2 <- as.formula(mpg ~ log1p(wt) + I(disp/100) + (1|cyl))
  mod3 <- as.formula(mpg ~ wt + (1|cyl) + (1 + wt|gear))
  mod4 <- as.formula(log(mpg) ~ wt + (1 + wt|cyl) + (1 + wt + am|gear))
  
  lfit1 <- lmer(mod1, data = mtcars)
  sfit1 <- stan_glmer(mod1, data = mtcars, iter = 400,
                      chains = CHAINS, seed = SEED, refresh = REFRESH)
  lfit2 <- update(lfit1, formula = mod2)
  sfit2 <- update(sfit1, formula = mod2)
  lfit3 <- update(lfit1, formula = mod3)
  sfit3 <- update(sfit1, formula = mod3)
  lfit4 <- update(lfit1, formula = mod4)
  sfit4 <- update(sfit1, formula = mod4)
  
  nd <- nd2 <- mtcars[1:5, ]
  nd2$cyl[2] <- 5 # add new levels
  nd3 <- nd2
  nd3$gear[2] <- 7
  nd3$gear[5] <- 1
  
  tol <- 0.3
  for (j in 1:4) {
    expect_equal(
      colMeans(posterior_predict(get(paste0("sfit", j)), newdata = nd, seed = SEED)),
      unname(predict(get(paste0("lfit", j)), newdata = nd)),
      tol = tol)
    expect_equal(
      colMeans(posterior_predict(get(paste0("sfit", j)), newdata = nd2, seed = SEED,
                                 allow.new.levels = TRUE)),
      unname(predict(get(paste0("lfit", j)), newdata = nd2, allow.new.levels = TRUE)),
      tol = tol)
    expect_equal(
      colMeans(posterior_predict(get(paste0("sfit", j)), newdata = nd3, seed = SEED,
                                 allow.new.levels = TRUE)),
      unname(predict(get(paste0("lfit", j)), newdata = nd3, allow.new.levels = TRUE)),
      tol = tol)
  }
})

test_that("posterior_predict close to predict.merMod for binomial", {
  d <- nd <- lme4::cbpp
  sfit <- example_model
  lfit <- glmer(formula(example_model), data = d, family = "binomial")
  levels(nd$herd) <- c(levels(nd$herd), "99")
  nd$herd[1:2] <- "99"
  lpred <- simulate(lfit, newdata = nd, re.form = NULL, allow.new.levels = TRUE,
                    nsim = 500, seed = SEED)
  for (j in 1:ncol(lpred)) {
    lpred[, j] <- lpred[, j][, 1] / rowSums(lpred[, j])
  }
  lpred <- t(as.matrix(lpred))
  spred <- posterior_predict(sfit, draws = 500, newdata = nd, 
                             seed = SEED)
  spred <- sweep(spred, 2, rowSums(get_y(sfit)), "/")
  expect_equal(colMeans(spred), unname(colMeans(lpred)),
               tol = .125)
})

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
                    iter = 400, chains = CHAINS, seed = SEED, refresh = REFRESH))
 
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
    fit1 <- stan_lmer(y ~ (1 | fac_nospace), data = df, seed = 123, chains = 2, iter = 25),
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
  expect_equal(ans2[, 1:2], ans1) # should be same as ans1 except for cols 3:4 with new levels
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
                      iter = 50, chains = 1, seed = 123),
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
                     iter = ITER, refresh = REFRESH, chains = CHAINS, 
                     seed = SEED))
  expect_equal(ppbt(fit), rep(1, nrow(mtcars)))
  expect_equal(ppbt(fit, newdata = mtcars[1:5, ]), rep(1, 5))
})

