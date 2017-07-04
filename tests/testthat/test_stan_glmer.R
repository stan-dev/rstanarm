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
stopifnot(require(lme4))
# stopifnot(require(gamm4))
stopifnot(require(HSAUR3))
ITER <- 400
CHAINS <- 2
SEED <- 123
REFRESH <- ITER
set.seed(SEED)
if (interactive()) options(mc.cores = parallel::detectCores())
FIXEF_tol <- 0.05
RANEF_tol <- 0.20 

source(file.path("helpers", "expect_stanreg.R"))
source(file.path("helpers", "SW.R"))

SW(fit <- stan_lmer(Reaction / 10 ~ Days + (Days | Subject), 
                    data = sleepstudy, refresh = REFRESH,
                    init_r = 0.05, chains = CHAINS, iter = ITER, seed = SEED))

context("stan_glmer")
test_that("draws from stan_glmer (gaussian) same as from stan_lmer", {
  SW(fit1 <- stan_glmer(mpg ~ wt + (1|cyl), data = mtcars, iter = 10, chains = 1, seed = SEED))
  SW(fit2 <- stan_lmer(mpg ~ wt + (1|cyl), data = mtcars, iter = 10, chains = 1, seed = SEED))
  expect_identical(as.matrix(fit1), as.matrix(fit2))
})
test_that("stan_glmer returns expected result for binomial cbpp example", {
  links <- c("logit", "probit", "cauchit", "log", "cloglog")
  # for (i in seq_along(links)) {
  i <- 1L # it seems only logit gives results similar to glmer with same link 
    fmla <- cbind(incidence, size - incidence) ~ period + (1 | herd)
    SW(fit <- stan_glmer(fmla, data = cbpp, family = binomial(links[i]),
                      chains = CHAINS, iter = ITER, seed = SEED, refresh = REFRESH))
    expect_stanreg(fit)
    
    ans <- glmer(fmla, data = cbpp, family = binomial(links[i]))
    expect_equal(fixef(fit), fixef(ans), tol = FIXEF_tol)
    expect_equal(ranef(fit), ranef(ans), tol = RANEF_tol)
    expect_equal(ngrps(fit), ngrps(ans))
  # }
})

context("stan_glmer.nb")
test_that("stan_glmer.nb ok", {
  dd <- expand.grid(f1 = factor(1:3),
                    f2 = LETTERS[1:2], g=1:9, rep=1:15,
                    KEEP.OUT.ATTRS=FALSE)
  mu <- 5*(-4 + with(dd, as.integer(f1) + 4*as.numeric(f2)))
  dd$y <- rnbinom(nrow(dd), mu = mu, size = 0.5)
  fmla <- as.formula(y ~ f1*f2 + (1|g))
  SW(fit <- stan_glmer.nb(formula = fmla, data = dd, init_r = 1,
                          iter = ITER, seed = SEED, algorithm = "meanfield"))
  expect_stanreg(fit)
  
  ans <- glmer.nb(formula = fmla, data = dd)
  # ans is messed up
  # expect_equal(fixef(fit), fixef(ans), tol = FIXEF_tol)
  # expect_equal(ranef(fit), ranef(ans), tol = RANEF_tol)
  expect_equal(ngrps(fit), ngrps(ans))
})

context("stan_lmer")
test_that("stan_lmer returns expected result for slepstudy example", {
  fmla <- formula(fit)
  expect_stanreg(fit)
  
  ans <- lmer(fmla, data = sleepstudy)
  expect_equal(fixef(fit), fixef(ans), tol = FIXEF_tol)
  # expect_equal(ranef(fit), ranef(ans), tol = RANEF_tol)
  expect_equal(ngrps(fit), ngrps(ans))
})
test_that("stan_lmer returns expected result for Penicillin example", {
  fmla <- as.formula(diameter ~ (1|plate) + (1|sample))
  SW(fit <- stan_lmer(fmla, data = Penicillin, chains = CHAINS, iter = ITER, 
                      seed = SEED, refresh = REFRESH, sparse = TRUE))
  expect_stanreg(fit)
  
  ans <- lmer(fmla, data = Penicillin)
  expect_equal(fixef(fit), fixef(ans), tol = FIXEF_tol)
  expect_equal(ranef(fit), ranef(ans), tol = RANEF_tol)
  expect_identical(ngrps(fit), ngrps(ans))
})
test_that("stan_lmer ok if global intercept forced to 0", {
  SW(fit <- stan_lmer(mpg ~ 0 + (1|cyl), data = mtcars, iter = 10, seed = SEED))
  expect_stanreg(fit)
})
test_that("stan_lmer returns an error when multiple group-specific terms are specified", {
  expect_error(
    stan_lmer(Reaction / 10 ~ Days + (Days | Subject) + (1|Subject), data = sleepstudy), 
    "formulas with duplicate group-specific terms"
  )
})
test_that("stan_lmer returns an error when 'family' specified", {
  expect_error(
    stan_lmer(Reaction / 10 ~ Days + (Days | Subject), family = "gaussian", data = sleepstudy),
    "'family' should not be specified"
  )
})


context("stan_gamm4")
source(file.path("helpers", "expect_gg.R"))
test_that("stan_gamm4 returns stanreg object", {
  sleepstudy$y <- sleepstudy$Reaction / 10
  SW(fit <- stan_gamm4(y ~ s(Days), data = sleepstudy, sparse = TRUE,
                       random = ~(1|Subject), chains = CHAINS, iter = ITER, 
                       seed = SEED, refresh = REFRESH))
  expect_stanreg(fit)
  # ans <- gamm4(Reaction / 10 ~ s(Days), data = sleepstudy, 
  #              random = ~(1|Subject))$mer
  # expect_equal(fixef(fit)[-1], fixef(ans)[-1], tol = FIXEF_tol, check.attributes = FALSE)
  # expect_equal(ranef(fit), ranef(ans), tol = RANEF_tol)
  # expect_identical(ngrps(fit), ngrps(ans))
  
  p1 <- plot_nonlinear(fit)
  p2 <- plot_nonlinear(fit, smooths = "s(Days)")
  expect_gg(p1)
  expect_gg(p2)
})

test_that("loo/waic for stan_glmer works", {
  ll_fun <- rstanarm:::ll_fun
  source(file.path("helpers", "expect_equivalent_loo.R"))
  
  # gaussian
  expect_equivalent_loo(fit)
  expect_identical(ll_fun(fit), rstanarm:::.ll_gaussian_i)
  
  # binomial
  expect_equivalent_loo(example_model)
  expect_identical(ll_fun(example_model), rstanarm:::.ll_binomial_i)
})
SW <- suppressWarnings
context("posterior_predict (stan_gamm4)")
test_that("stan_gamm4 returns expected result for sleepstudy example", {
  sleepstudy$y <- sleepstudy$Reaction / 10
  fit <- SW(stan_gamm4(y ~ s(Days), data = sleepstudy,
                       random = ~(1|Subject), chains = CHAINS, iter = ITER, 
                       seed = SEED, refresh = REFRESH))
  expect_silent(yrep1 <- posterior_predict(fit))
  # expect_equal(dim(yrep1), c(nrow(as.data.frame(fit)), nobs(fit)))
  expect_silent(yrep2 <- posterior_predict(fit, draws = 1))
  # expect_equal(dim(yrep2), c(1, nobs(fit)))
  expect_silent(posterior_predict(fit, newdata = sleepstudy))
})


context("posterior_predict (stan_(g)lmer)")
source(file.path("helpers", "check_for_error.R"))
source(file.path("helpers", "expect_linpred_equal.R"))
test_that("compatible with stan_lmer", {
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
      tol = tol, check.attributes = FALSE)
    expect_equal(
      colMeans(posterior_predict(get(paste0("sfit", j)), newdata = nd2, seed = SEED,
                                 allow.new.levels = TRUE)),
      unname(predict(get(paste0("lfit", j)), newdata = nd2, allow.new.levels = TRUE)),
      tol = tol, check.attributes = FALSE)
    expect_equal(
      colMeans(posterior_predict(get(paste0("sfit", j)), newdata = nd3, seed = SEED,
                                 allow.new.levels = TRUE)),
      unname(predict(get(paste0("lfit", j)), newdata = nd3, allow.new.levels = TRUE)),
      tol = tol, check.attributes = FALSE)
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
               tol = .125, check.attributes = FALSE)
})
