# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2015, 2016 Trustees of Columbia University
# Copyright (C) 2017 Sam Brilleman
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
ITER <- 1000
CHAINS <- 1
SEED <- 12345
REFRESH <- 0L
set.seed(SEED)
if (interactive()) 
  options(mc.cores = parallel::detectCores())

TOLSCALES <- list(
  lmer_fixef = 0.25,  # how many SEs can stan_jm fixefs be from lmer fixefs
  lmer_ranef = 0.05, # how many SDs can stan_jm ranefs be from lmer ranefs
  glmer_fixef = 0.3, # how many SEs can stan_jm fixefs be from glmer fixefs
  glmer_ranef = 0.1 # how many SDs can stan_jm ranefs be from glmer ranefs
)

source(test_path("helpers", "expect_matrix.R"))
source(test_path("helpers", "expect_stanreg.R"))
source(test_path("helpers", "expect_stanmvreg.R"))
source(test_path("helpers", "expect_survfit.R"))
source(test_path("helpers", "expect_ppd.R"))
source(test_path("helpers", "expect_identical_sorted_stanmats.R"))
source(test_path("helpers", "SW.R"))
source(test_path("helpers", "get_tols.R"))
source(test_path("helpers", "recover_pars.R"))

context("stan_mvmer")

#----  Data (for non-Gaussian families)

pbcLong$ybern <- as.integer(pbcLong$logBili >= mean(pbcLong$logBili))
pbcLong$ybino <- as.integer(rpois(nrow(pbcLong), 5))
pbcLong$ypois <- as.integer(pbcLong$albumin)
pbcLong$ynbin <- as.integer(rnbinom(nrow(pbcLong), 3, .3))
pbcLong$ygamm <- as.numeric(pbcLong$platelet / 10)
pbcLong$xbern <- as.numeric(pbcLong$platelet / 100)
pbcLong$xpois <- as.numeric(pbcLong$platelet / 100)
pbcLong$xgamm <- as.numeric(pbcLong$logBili)

#----  Models

# univariate GLM
fm1 <- logBili ~ year + (year | id)
o<-SW(m1 <- stan_mvmer(fm1, pbcLong, iter = 100, chains = 1, seed = SEED))

# multivariate GLM
fm2 <- list(logBili ~ year + (year | id), albumin ~ year + (year | id))
o<-SW(m2 <- stan_mvmer(fm2, pbcLong, iter = 100, chains = 1, seed = SEED))

#----  Tests for stan_mvmer arguments

test_that("formula argument works", {
  SW(m991 <- update(m1, formula. = list(fm1)))
  expect_identical(as.matrix(m1), as.matrix(m991)) # fm as list
})

test_that("data argument works", {
  SW(m991 <- update(m1, data = list(pbcLong)))
  SW(m992 <- update(m2, data = list(pbcLong, pbcLong)))
  expect_identical(as.matrix(m1), as.matrix(m991)) # data as list
  expect_identical(as.matrix(m2), as.matrix(m992))
})

test_that("family argument works", {
  
  expect_output(ret <- update(m1, family = "gaussian"))
  expect_output(ret <- update(m1, family = gaussian))
  expect_output(ret <- update(m1, family = gaussian(link = identity)))
  
  expect_output(ret <- update(m1, formula. = ybern ~ ., family = binomial))
  expect_output(ret <- update(m1, formula. = ypois ~ ., family = poisson))
  expect_output(ret <- update(m1, formula. = ypois ~ ., family = neg_binomial_2))
  expect_output(ret <- update(m1, formula. = ygamm ~ ., family = Gamma, init = 0))
  expect_output(ret <- update(m1, formula. = ygamm ~ ., family = inverse.gaussian, init = 0))
  
  expect_error(ret <- update(m1, formula. = ybino ~ ., family = binomial))
  
  # multivariate model with combinations of family
  expect_output(ret <- update(m2, formula. = list(~ ., ybern ~ .), 
                              family = list(gaussian, binomial)))
})

test_that("prior_PD argument works", {
  expect_output(update(m1, prior_PD = TRUE))
})

test_that("adapt_delta argument works", {
  expect_output(update(m1, adapt_delta = NULL))
  expect_output(update(m1, adapt_delta = 0.8))
})

test_that("error message occurs for arguments not implemented", {
  expect_error(update(m1, weights = 1:10), "not yet implemented")
  expect_error(update(m1, QR = TRUE), "not yet implemented")
  expect_error(update(m1, sparse = TRUE), "not yet implemented")
})

#----  Check models with multiple grouping factors

test_that("multiple grouping factors are ok", {
  
  tmpdat <- pbcLong
  tmpdat$practice <- cut(pbcLong$id, c(0,10,20,30,40))
  
  tmpfm1 <- logBili ~ year + (year | id) + (1 | practice)
  SW(ok_mod1 <- update(m1, formula. = tmpfm1, data = tmpdat, iter = 1, refresh = 0, init = 0))
  expect_stanmvreg(ok_mod1)
  
  tmpfm2 <- list(
    logBili ~ year + (year | id) + (1 | practice),
    albumin ~ year + (year | id))
  SW(ok_mod2 <- update(m2, formula. = tmpfm2, data = tmpdat, iter = 1, refresh = 0, init = 0))
  expect_stanmvreg(ok_mod2)
  
  tmpfm3 <- list(
    logBili ~ year + (year | id) + (1 | practice),
    albumin ~ year + (year | id) + (1 | practice))
  SW(ok_mod3 <- update(m2, formula. = tmpfm3, data = tmpdat, iter = 1, refresh = 0, init = 0))
  expect_stanmvreg(ok_mod3)
  
  # check reordering grouping factors is ok
  # NB it seems these comparisons must be made using init = 0 and one iteration,
  # probably because the order of the parameters passed to Stan affects the 
  # sequence of MCMC samples even when the same seed is used. An alternative
  # would be to test equality of the stanmat colMeans with specified tolerance?
  tmpfm4 <- list(
    logBili ~ year + (1 | practice) + (year | id),
    albumin ~ year + (year | id))
  SW(ok_mod4 <- update(ok_mod2, formula. = tmpfm4))
  expect_identical_sorted_stanmats(ok_mod2, ok_mod4)

  tmpfm5 <- list(
    logBili ~ year + (1 | practice) + (year | id),
    albumin ~ year + (year | id) + (1 | practice))
  SW(ok_mod5 <- update(ok_mod3, formula. = tmpfm5))
  expect_identical_sorted_stanmats(ok_mod3, ok_mod5)
  
  tmpfm6 <- list(
    logBili ~ year + (1 | practice) + (year | id),
    albumin ~ year + (1 | practice) + (year | id))
  SW(ok_mod6 <- update(ok_mod3, formula. = tmpfm6))
  expect_identical_sorted_stanmats(ok_mod3, ok_mod6)
})

#----  Compare estimates: univariate stan_mvmer vs stan_glmer

if (interactive()) {
  compare_glmer <- function(fmLong, fam = gaussian, ...) {
    SW(y1 <- stan_glmer(fmLong, pbcLong, fam, iter = 1000, chains = CHAINS, seed = SEED))
    SW(y2 <- stan_mvmer(fmLong, pbcLong, fam, iter = 1000, chains = CHAINS, seed = SEED, ...))
    tols <- get_tols(y1, tolscales = TOLSCALES)
    pars <- recover_pars(y1)
    pars2 <- recover_pars(y2)
    for (i in names(tols$fixef))
      expect_equal(pars$fixef[[i]], pars2$fixef[[i]], tol = tols$fixef[[i]])     
    for (i in names(tols$ranef))
      expect_equal(pars$ranef[[i]], pars2$ranef[[i]], tol = tols$ranef[[i]])
    expect_equal(colMeans(log_lik(y1)), 
                 colMeans(log_lik(y2)), tol = 0.15)
    nd <- pbcLong[stats::complete.cases(pbcLong), , drop = FALSE]
    expect_equal(colMeans(log_lik(y1, newdata = nd)), 
                 colMeans(log_lik(y2, newdata = nd)), tol = 0.15)
  }
  test_that("coefs same for stan_jm and stan_lmer/coxph", {
    # fails in many cases
    # compare_glmer(logBili ~ year + (1 | id), gaussian)
    })
  # fails in some cases
  # test_that("coefs same for stan_jm and stan_glmer, bernoulli", {
  #   compare_glmer(ybern ~ year + xbern + (1 | id), binomial)})
  test_that("coefs same for stan_jm and stan_glmer, poisson", {
    compare_glmer(ypois ~ year + xpois + (1 | id), poisson, init = 0)})
  test_that("coefs same for stan_jm and stan_glmer, negative binomial", {
    compare_glmer(ynbin ~ year + xpois + (1 | id), neg_binomial_2)})
  test_that("coefs same for stan_jm and stan_glmer, Gamma", {
    compare_glmer(ygamm ~ year + xgamm + (1 | id), Gamma(log))})
#  test_that("coefs same for stan_jm and stan_glmer, inverse gaussian", {
#    compare_glmer(ygamm ~ year + xgamm + (1 | id), inverse.gaussian)})  
}

#----  Check methods and post-estimation functions

tmpdat <- pbcLong
tmpdat$practice <- cut(pbcLong$id, c(0,10,20,30,40))

o<-SW(f1 <- update(m1, formula. = list(logBili ~ year + (year | id)), data = tmpdat))
o<-SW(f2 <- update(f1, formula. = list(logBili ~ year + (year | id) + (1 | practice))))
o<-SW(f3 <- update(m2, formula. = list(logBili ~ year + (year | id) + (1 | practice),
                                       albumin ~ year + (year | id)), data = tmpdat))
o<-SW(f4 <- update(f3, formula. = list(logBili ~ year + (year | id) + (1 | practice),
                                       albumin ~ year + (year | id) + (1 | practice))))
o<-SW(f5 <- update(f3, formula. = list(logBili ~ year + (year | id) + (1 | practice),
                                       ybern ~ year + (year | id) + (1 | practice)),
                   family = list(gaussian, binomial)))

for (j in 1:5) {
  mod <- get(paste0("f", j))
  cat("Checking model:", paste0("f", j), "\n")

  expect_error(posterior_traj(mod), "stanjm")
  expect_error(posterior_survfit(mod), "stanjm")
     
  test_that("posterior_predict works with estimation data", {
    pp <- posterior_predict(mod, m = 1)
    expect_ppd(pp)
    if (mod$n_markers > 1L) {
      pp <- posterior_predict(mod, m = 2)
      expect_ppd(pp)
    }
  })  
  test_that("log_lik works with estimation data", {
    ll <- log_lik(mod)
    expect_matrix(ll)
    expect_identical(ll, log_lik(mod, m = 1))
    if (mod$n_markers > 1L)
      expect_matrix(log_lik(mod, m = 2))
  })   
  
  nd <- tmpdat[tmpdat$id == 2,]
  test_that("posterior_predict works with new data (one individual)", {
    pp <- posterior_predict(mod, m = 1, newdata = nd)
    expect_ppd(pp)
    if (mod$n_markers > 1L) {
      pp <- posterior_predict(mod, m = 2, newdata = nd)
      expect_ppd(pp)
    }
  })     
  test_that("log_lik works with new data (one individual)", {
    ll <- log_lik(mod, newdata = nd)
    expect_matrix(ll)
    expect_identical(ll, log_lik(mod, m = 1, newdata = nd))
    if (mod$n_markers > 1L)
      expect_matrix(log_lik(mod, m = 2, newdata = nd))
    # log_lik is only designed for one submodel at a time so passing
    # newdata as a list should generate an error in validate_newdata
    expect_error(log_lik(mod, newdata = list(nd)), "data frame") 
  }) 
  
  nd <- tmpdat[tmpdat$id %in% c(1,2),]
  test_that("posterior_predict works with new data (multiple individuals)", {
    pp <- posterior_predict(mod, m = 1, newdata = nd)
    expect_ppd(pp)
    if (mod$n_markers > 1L) {
      pp <- posterior_predict(mod, m = 2, newdata = nd)
      expect_ppd(pp)
    }
  })
  test_that("log_lik works with estimation data", {
    expect_matrix(log_lik(mod, newdata = nd))
    if (mod$n_markers > 1L)
      expect_matrix(log_lik(mod, m = 2, newdata = nd))
  }) 
  
  test_that("loo and waic work", {
    l <- suppressWarnings(loo(mod))
    w <- suppressWarnings(waic(mod))
    expect_s3_class(l, "loo")
    expect_s3_class(w, "loo")
    expect_s3_class(w, "waic")
    att_names <- c('names', 'dims', 'class', 'model_name', 'discrete', 'yhash', 'formula')
    expect_named(attributes(l), att_names)
    expect_named(attributes(w), att_names)
  })
  
  test_that("extraction methods work", {
    M <- mod$n_markers
    fe <- fixef(mod)
    re <- ranef(mod)
    ce <- coef(mod)
    mf <- model.frame(mod)
    tt <- terms(mod)
    fm <- formula(mod)
    fam <- family(mod)
    sig <- sigma(mod)
    expect_is(fe, "list"); expect_identical(length(fe), M)
    expect_is(re, "list"); expect_identical(length(re), M)
    expect_is(ce, "list"); expect_identical(length(re), M)
    expect_is(mf, "list"); expect_identical(length(mf), M); lapply(mf, function(x) expect_is(x, "data.frame"))
    expect_is(tt, "list"); expect_identical(length(tt), M); lapply(tt, function(x) expect_is(x, "terms"))
    expect_is(fm, "list"); expect_identical(length(fm), M); lapply(fm, function(x) expect_is(x, "formula"))
    expect_is(fam,"list"); expect_identical(length(fam),M); lapply(fam, function(x) expect_is(x, "family"))
    expect_is(sig, "numeric");
  })
  
  test_that("these extraction methods are currently disallowed", {
    expect_error(se(mod), "Not currently implemented")
    expect_error(fitted(mod), "Not currently implemented")
    expect_error(residuals(mod), "Not currently implemented")
  })
}

