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

# tests can be run using devtools::test() or manually by loading testthat 
# package and then running the code below possibly with options(mc.cores = 4).

library(rstanarm)
library(lme4)
SEED <- 123
set.seed(SEED)
ITER <- 10
CHAINS <- 2
REFRESH <- 0

SW <- suppressWarnings

# These tests just make sure that posterior_predict doesn't throw errors and
# that result has correct dimensions
check_for_error <- function(fit) {
  nsims <- nrow(as.data.frame(fit))
  mf <- model.frame(fit)
  if (identical(deparse(substitute(fit)), "example_model"))
    mf <- lme4::cbpp
  
  
  expect_silent(yrep1 <- posterior_predict(fit))
  expect_equal(dim(yrep1), c(nsims, nobs(fit)))

  expect_silent(yrep2 <- posterior_predict(fit, draws = 1))
  expect_equal(dim(yrep2), c(1, nobs(fit)))
  
  expect_silent(yrep3 <- posterior_predict(fit, newdata = mf[1,]))
  expect_equal(dim(yrep3), c(nsims, 1))
  
  expect_silent(yrep4 <- posterior_predict(fit, draws = 2, newdata = mf[1,]))
  expect_equal(dim(yrep4), c(2, 1))
  
  expect_silent(yrep5 <- posterior_predict(fit, newdata = mf[1:5,]))
  expect_equal(dim(yrep5), c(nsims, 5))
  
  expect_silent(yrep6 <- posterior_predict(fit, draws = 3, newdata = mf[1:5,]))
  expect_equal(dim(yrep6), c(3, 5))
  
  expect_error(posterior_predict(fit, draws = nsims + 1), 
               regexep = "posterior sample size is only")
}

context("posterior_predict (stan_lm)")
test_that("posterior_predict compatible with stan_lm", {
  fit <- SW(stan_lm(mpg ~ wt + cyl + am, data = mtcars, prior = R2(log(0.5), what = "log"),
                 iter = ITER, chains = CHAINS,  seed = SEED, refresh = REFRESH))
  check_for_error(fit)
})

context("posterior_predict (stan_glm)")
test_that("compatible with gaussian glm", {
  fit <- SW(stan_glm(mpg ~ wt, data = mtcars, 
                     iter = ITER, chains = CHAINS, seed = SEED, refresh = REFRESH))
  check_for_error(fit)
  fit_off <- SW(update(fit, offset = runif(nrow(mtcars))))
  check_for_error(fit)
})
test_that("compatible with poisson & negbin glm", {
  counts <- c(18,17,15,20,10,20,25,13,12)
  outcome <- gl(3,1,9)
  treatment <- gl(3,3)
  fit <- SW(stan_glm(counts ~ outcome + treatment, family = poisson(), 
                     iter = ITER, chains = CHAINS, seed = SEED, refresh = REFRESH))
  fitnb <- SW(update(fit, family = neg_binomial_2))
  check_for_error(fit)
  check_for_error(fitnb)
})
test_that("posterior_predict compatible with gamma & inverse.gaussian glm", {
  clotting <- data.frame(log_u = log(c(5,10,15,20,30,40,60,80,100)),
                         lot1 = c(118,58,42,35,27,25,21,19,18),
                         lot2 = c(69,35,26,21,18,16,13,12,12))
  fit <- SW(stan_glm(lot1 ~ log_u, data = clotting, family = Gamma, 
                  chains = CHAINS, iter = ITER,  seed = SEED, refresh = REFRESH))
  check_for_error(fit)
  
  # inverse gaussian
  fit_igaus <- SW(update(fit, family = inverse.gaussian))
  check_for_error(fit_igaus)
})

context("posterior_predict (stan_polr)")
test_that("compatible with stan_polr", {
  fit <- SW(stan_polr(tobgp ~ agegp + alcgp, data = esoph, prior = R2(location = 0.4),
                   iter = ITER, chains = CHAINS, seed = SEED, refresh = REFRESH))
  check_for_error(fit)
})

context("posterior_predict (stan_(g)lmer)")
test_that("compatible with stan_lmer", {
  fit <- SW(stan_lmer(mpg ~ wt + (1|cyl) + (1 + wt|gear), data = mtcars, 
                      prior = normal(0,1), iter = ITER, chains = CHAINS,
                      seed = SEED, refresh = REFRESH))
  check_for_error(fit)
})
test_that("compatible with stan_glmer (binomial)", {
  check_for_error(example_model)
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
})


context("posterior_predict (optimizing and vb)")
test_that("errors for optimizing and silent for vb", {
  fit1 <- stan_glm(mpg ~ wt + cyl + am, data = mtcars, algorithm = "optimizing", 
                   seed = SEED)
  fit2 <- update(fit1, algorithm = "meanfield")
  fit3 <- update(fit1, algorithm = "fullrank")
  expect_error(posterior_predict(fit1), regexp = "optimizing")
  expect_silent(posterior_predict(fit2))
  expect_silent(posterior_predict(fit3))
})


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
      predict(get(paste0("sfit", j))),
      unname(predict(get(paste0("lfit", j)))),
      tol = tol)
    expect_equal(
      predict(get(paste0("sfit", j)), newdata = nd),
      predict(get(paste0("lfit", j)), newdata = nd),
      tol = tol)
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
               tol = .1)
})

test_that("edge cases for posterior_predict work correctly", {
  dims <- c(nrow(as.matrix(example_model)), nrow(lme4::cbpp))
  expect_identical(posterior_predict(example_model, re.form = NA, seed = SEED),
                   posterior_predict(example_model, re.form = ~0, seed = SEED))
  expect_identical(posterior_predict(example_model, seed = SEED),
                   posterior_predict(example_model, newdata = lme4::cbpp, seed = SEED))
  expect_error(posterior_predict(example_model, re.form = ~1))
  expect_error(posterior_predict(example_model, re.form = ~(1|foo)))
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

