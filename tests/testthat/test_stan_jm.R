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
stopifnot(require(lme4))
stopifnot(require(survival))
ITER <- 1000
CHAINS <- if (interactive()) 1 else 1
SEED <- 12345
REFRESH <- ITER
set.seed(SEED)
if (interactive()) options(mc.cores = parallel::detectCores())
TOLSCALE_lmer_fixef <- 0.2  # how many SEs can stan_jm fixefs be from lmer fixefs
TOLSCALE_lmer_ranef <- 0.05 # how many SDs can stan_jm ranefs be from lmer ranefs
TOLSCALE_glmer_fixef <- 0.3 # how many SEs can stan_jm fixefs be from glmer fixefs
TOLSCALE_glmer_ranef <- 0.1 # how many SDs can stan_jm ranefs be from glmer ranefs
TOLSCALE_event <- 0.2 # how many SEs can stan_jm fixefs be from coxph fixefs
FIXEF_tol <- 0.02
RANEF_tol <- 0.05
EVENT_tol <- 0.05

expect_matrix  <- function(x) expect_identical(class(x), "matrix")
expect_stanjm  <- function(x) expect_s3_class(x, "stanjm")
expect_survfit <- function(x) expect_s3_class(x, "survfit.stanjm")
expect_ppd     <- function(x) expect_s3_class(x, "ppd")
expect_stanreg <- function(x) expect_s3_class(x, "stanreg")
SW <- function(expr) capture.output(suppressWarnings(expr))

context("stan_jm")

#--------  Models

examplejm1 <- 
  stan_jm(
    formulaLong = logBili ~ year + (1 | id), 
    dataLong = pbcLong,
    formulaEvent = Surv(futimeYears, death) ~ sex + trt, 
    dataEvent = pbcSurv,
    time_var = "year", iter = 1,
    chains = CHAINS, seed = SEED)

examplejm2 <- 
  stan_jm(
    formulaLong = logBili ~ year + (year | id), 
    dataLong = pbcLong,
    formulaEvent = Surv(futimeYears, death) ~ sex + trt, 
    dataEvent = pbcSurv,
    time_var = "year", iter = 1,
    chains = CHAINS, seed = SEED)

examplejm3 <- 
  stan_jm(
    formulaLong = list(
      logBili ~ year + (year | id),
      albumin ~ year + (1 | id)),
    dataLong = pbcLong,
    formulaEvent = Surv(futimeYears, death) ~ sex + trt, 
    dataEvent = pbcSurv,
    time_var = "year", iter = 1,
    chains = CHAINS, seed = SEED)

examplejm4 <- 
  stan_jm(
    formulaLong = list(
      logBili ~ year + (year | id),
      albumin ~ year + (year | id)),
    dataLong = pbcLong,
    formulaEvent = Surv(futimeYears, death) ~ sex + trt, 
    dataEvent = pbcSurv,
    time_var = "year", iter = 1,
    chains = CHAINS, seed = SEED)

#--------  Tests

test_that("formula argument works", {
  
})


test_that("data argument works", {
  
})


test_that("id_var argument works", {
  
  # Models with a single grouping factor
  expect_output(update(examplejm1, id_var = "id"))
  expect_output(expect_warning(update(examplejm1, id_var = "year"), 
                               "are not the same; 'id_var' will be ignored"))
  
  # Models with more than one grouping factor
  tmpdat <- pbcLong
  tmpdat$practice <- cut(pbcLong$id, c(0,10,20,30,40))
  tmpfm <- logBili ~ year + (1 | id) + (1 | practice)
  ok_mod <- update(examplejm1, formulaLong. = tmpfm, dataLong = tmpdat, id_var = "id", init = 0)
  expect_stanjm(ok_mod)
  expect_error(update(ok_mod, id_var = NULL), "'id_var' must be specified")
  expect_error(update(ok_mod, id_var = "year"), "'id_var' must be included as a grouping factor")
  expect_error(update(ok_mod, id_var = "practice"), "'id_var' must correspond to the lowest level of clustering")
})


test_that("family argument works", {
  
})


test_that("assoc argument works", {
  
  # Univariate joint models
  
  expect_output(ret <- update(examplejm2, assoc = NULL))
  expect_output(ret <- update(examplejm2, assoc = "null"))
  expect_output(ret <- update(examplejm2, assoc = "etavalue"))
  expect_output(ret <- update(examplejm2, assoc = "etaslope"))
  expect_output(ret <- update(examplejm2, assoc = "muvalue"))
  expect_output(ret <- update(examplejm2, assoc = "muslope"))
  expect_output(ret <- update(examplejm2, assoc = c("etavalue", "etaslope"))) 
  expect_output(ret <- update(examplejm2, assoc = c("etavalue", "muslope"))) 
  expect_output(ret <- update(examplejm2, assoc = c("muvalue", "etaslope"))) 
  expect_output(ret <- update(examplejm2, assoc = c("muvalue", "muslope")))
 
  expect_error(ret <- update(examplejm2, assoc = c("etavalue", "muvalue")), "cannot be specified together")
  expect_error(ret <- update(examplejm2, assoc = c("etaslope", "muslope")), "cannot be specified together")
  
  expect_output(ret <- update(examplejm2, assoc = "shared_b"))
  expect_output(ret <- update(examplejm2, assoc = "shared_b(1)"))
  expect_output(ret <- update(examplejm2, assoc = "shared_b(2)"))
  expect_output(ret <- update(examplejm2, assoc = "shared_b(1:2)"))
  expect_output(ret <- update(examplejm2, assoc = "shared_b(1,2)"))
  
  expect_output(ret <- update(examplejm2, assoc = "shared_coef"))
  expect_output(ret <- update(examplejm2, assoc = "shared_coef(1)"))
  expect_output(ret <- update(examplejm2, assoc = "shared_coef(2)"))
  expect_output(ret <- update(examplejm2, assoc = "shared_coef(1:2)"))
  expect_output(ret <- update(examplejm2, assoc = "shared_coef(1,2)"))
  
  expect_error(ret <- update(examplejm2, assoc = "shared_b(10)"), "greater than the number of")
  expect_error(ret <- update(examplejm2, assoc = "shared_coef(10)"), "greater than the number of")
  expect_error(ret <- update(examplejm2, assoc = c("shared_b(1)", "shared_coef(1)")), "should not be specified in both")
  expect_error(ret <- update(examplejm2, assoc = c("shared_b", "shared_coef")), "should not be specified in both")
  
  expect_output(ret <- update(examplejm2, assoc = list(NULL)))
  expect_output(ret <- update(examplejm2, assoc = list("null")))
  expect_output(ret <- update(examplejm2, assoc = list("etavalue")))
  expect_output(ret <- update(examplejm2, assoc = list("etaslope")))
  expect_output(ret <- update(examplejm2, assoc = list("muvalue")))
  expect_output(ret <- update(examplejm2, assoc = list("muslope")))
  expect_output(ret <- update(examplejm2, assoc = list(c("etavalue", "etaslope")))) 
  expect_output(ret <- update(examplejm2, assoc = list(c("etavalue", "muslope")))) 
  expect_output(ret <- update(examplejm2, assoc = list(c("muvalue", "etaslope")))) 
  expect_output(ret <- update(examplejm2, assoc = list(c("muvalue", "muslope"))))  
  
  expect_error(ret <- update(examplejm2, assoc = NA), "'assoc' should be") 
  expect_error(ret <- update(examplejm2, assoc = 123), "'assoc' should be") 
  expect_error(ret <- update(examplejm2, assoc = c(1,2,3)), "'assoc' should be") 
  
  expect_error(ret <- update(examplejm2, assoc = c("wrong")), "unsupported association type") 
  expect_error(ret <- update(examplejm2, assoc = list("wrong")), "unsupported association type") 
  
  expect_error(ret <- update(examplejm2, assoc = list(NULL, NULL)), "incorrect length") 
  expect_error(ret <- update(examplejm2, assoc = list("etavalue", "etavalue")), "incorrect length") 
  expect_error(ret <- update(examplejm2, assoc = list(c("etavalue", "etaslope"), "etavalue")), "incorrect length") 
  
  # Multivariate joint models
  
  expect_output(ret <- update(examplejm3, assoc = "etavalue"))
  expect_output(ret <- update(examplejm3, assoc = "etaslope"))
  expect_output(ret <- update(examplejm3, assoc = "muvalue"))
  expect_output(ret <- update(examplejm3, assoc = "muslope"))
  expect_output(ret <- update(examplejm3, assoc = c("etavalue", "etaslope"))) 
  expect_output(ret <- update(examplejm3, assoc = c("etavalue", "muslope"))) 
  expect_output(ret <- update(examplejm3, assoc = c("muvalue", "etaslope"))) 
  expect_output(ret <- update(examplejm3, assoc = c("muvalue", "muslope")))
  
  expect_output(ret <- update(examplejm3, assoc = list("etavalue")))
  expect_output(ret <- update(examplejm3, assoc = list("etavalue", "etavalue")))
  expect_output(ret <- update(examplejm3, assoc = list(c("etavalue", "etaslope"), "etavalue")))
  expect_output(ret <- update(examplejm3, assoc = list("etavalue", c("etavalue", "etaslope"))))
  expect_output(ret <- update(examplejm3, assoc = list(c("etavalue", "etaslope"), c("muvalue", "muslope"))))
  
  expect_error(ret <- update(examplejm3, assoc = list("wrong", "etavalue")), "unsupported association type")
  expect_error(ret <- update(examplejm3, assoc = list("null", "etavalue", "etaslope")), "incorrect length")
  expect_error(ret <- update(examplejm3, assoc = data.frame("etavalue", "etaslope")), "'assoc' should be") 
  
})


test_that("basehaz argument works", {
  
  expect_output(update(examplejm1, basehaz = "weibull"))
  expect_output(update(examplejm1, basehaz = "bs"))
  expect_output(update(examplejm1, basehaz = "piecewise"))
  
  expect_output(update(examplejm1, basehaz = "bs", basehaz_ops = list(df = 5)))
  expect_output(update(examplejm1, basehaz = "bs", basehaz_ops = list(knots = c(1,3,5))))
  expect_output(update(examplejm1, basehaz = "piecewise", basehaz_ops = list(df = 5)))
  expect_output(update(examplejm1, basehaz = "piecewise", basehaz_ops = list(knots = c(1,3,5))))
  
  expect_output(expect_warning(update(examplejm1, basehaz = "weibull", basehaz_ops = list(df = 1)), "'df' will be ignored"))
  expect_output(expect_warning(update(examplejm1, basehaz = "weibull", basehaz_ops = list(knots = 1)), "'knots' will be ignored"))
  
  expect_output(update(examplejm1, basehaz = "piecewise", basehaz_ops = list(knots = c(1,3,5))))
  
  expect_error(update(examplejm1, basehaz = "bs", basehaz_ops = list(df = 1)), "must be at least 3")
  expect_error(update(examplejm1, basehaz = "bs", basehaz_ops = list(knots = -1)), "'knots' must be non-negative")
  expect_error(update(examplejm1, basehaz = "piecewise", basehaz_ops = list(knots = -1)), "'knots' must be non-negative")
  expect_error(update(examplejm1, basehaz = "piecewise", basehaz_ops = list(knots = c(1,2,50))), "cannot be greater than the largest event time")
  
})


test_that("quadnodes argument works", {
  
  expect_output(update(examplejm1, quadnodes = 7))
  expect_output(update(examplejm1, quadnodes = 11))
  expect_output(update(examplejm1, quadnodes = 15))
  
  expect_error(update(examplejm1, quadnodes = 1), "'quadnodes' must be either 7, 11 or 15")
  expect_error(update(examplejm1, quadnodes = c(1,2)), "should be a numeric vector of length 1")
  expect_error(update(examplejm1, quadnodes = "wrong"), "should be a numeric vector of length 1")
  
})


test_that("weights argument works", {
  
  idvec0 <- pbcSurv[["id"]]
  idvec1 <- head(idvec0)            # missing IDs
  idvec2 <- rep(idvec0, each = 2)   # repeated IDs
  idvec3 <- c(idvec0, 9998, 9999)    # extra IDs not in model
  
  wts0 <- data.frame(id = idvec0, weights = rep_len(c(1,2), length(idvec0)))
  wts1 <- data.frame(id = idvec1, weights = rep_len(c(1,2), length(idvec1)))
  wts2 <- data.frame(id = idvec2, weights = rep_len(c(1,2), length(idvec2)))
  wts3 <- data.frame(id = idvec0, weights = rep_len(c(1,2), length(idvec0)),
                     junkcol = idvec0)
  wts4 <- data.frame(id = idvec0, weights = rep_len(c("word"), length(idvec0)))
  wts5 <- data.frame(id = idvec0, weights = rep_len(c(NA), length(idvec0)))
  wts6 <- data.frame(id = idvec0, weights = rep_len(c(-1, 1), length(idvec0)))
  wts7 <- data.frame(id = idvec3, weights = rep_len(c(1,2), length(idvec3)))
  
  expect_output(update(examplejm1, weights = wts0, iter = 5))
  expect_output(update(examplejm1, weights = wts7, iter = 5)) # ok to supply extra IDs in weights
  
  expect_error(update(examplejm1, weights = as.matrix(wts0)), "should be a data frame")
  expect_error(update(examplejm1, weights = wts1), "do not have weights supplied")
  expect_error(update(examplejm1, weights = wts2), "should only have one row")
  expect_error(update(examplejm1, weights = wts3), "should be a data frame with two columns")
  expect_error(update(examplejm1, weights = wts4), "weights supplied must be numeric")
  expect_error(update(examplejm1, weights = wts5), "weights supplied must be numeric")
  expect_error(update(examplejm1, weights = wts6), "Negative weights are not allowed")
  
})

test_that("init argument works", {
  expect_output(update(examplejm1, init = "model_based"))
  expect_output(update(examplejm1, init = "0"))
  expect_output(update(examplejm1, init = 0))
  expect_output(update(examplejm1, init = "random"))
})

test_that("prior_PD argument works", {
  expect_output(update(examplejm1, prior_PD = TRUE))
})

test_that("adapt_delta argument works", {
  expect_output(update(examplejm1, adapt_delta = NULL))
  expect_output(update(examplejm1, adapt_delta = 0.8))
})

test_that("max_treedepth argument works", {
  expect_output(update(examplejm1, max_treedepth = NULL))
  expect_output(update(examplejm1, max_treedepth = 5))
  expect_output(update(examplejm1, max_treedepth = 5L))
})

test_that("error message occurs for arguments not implemented", {
  expect_error(update(examplejm1, offset = 1:10), "not yet implemented")
  expect_error(update(examplejm1, QR = TRUE), "not yet implemented")
  expect_error(update(examplejm1, sparse = TRUE), "not yet implemented")
  expect_error(update(examplejm1, dataAssoc = data.frame(a=1:10)), "not yet implemented")
})

#--------  Check parameter estimates with NULL assoc against stan_glmer and coxph

if (interactive()) {

  compare_lmer <- function(fm) {
    y1 <- stan_lmer(fm, pbcLong, iter = 2000, chains = CHAINS, seed = SEED)
    s1 <- coxph(Surv(futimeYears, death) ~ sex + trt, data = pbcSurv)
    j1 <- stan_jm(fm, pbcLong, Surv(futimeYears, death) ~ sex + trt, pbcSurv,
                  time_var = "year", assoc = NULL, 
                  iter = 2000, chains = CHAINS, seed = SEED) 
    ranef_sds <- attr(VarCorr(y1)[[1]], "stddev")
    tols <- TOLSCALE_lmer_ranef * ranef_sds
    for (i in 1:length(tols))
      expect_equal(ranef(y1)[[1]][[i]], ranef(j1)$Long1[[1]][[i]], tol = tols[[i]])     
    fixef_ses <- y1$ses
    tols <- TOLSCALE_lmer_fixef * fixef_ses
    if ("(Intercept)" %in% names(tols)) # use weaker tolerance for intercept
      tols[["(Intercept)"]] <- 2 * tols[["(Intercept)"]] 
    for (i in 1:length(fixef(y1)))
      expect_equal(fixef(y1)[[i]], fixef(j1)$Long1[[i]], tol = tols[[i]])
    event_ses <- summary(s1)$coefficients[, "se(coef)"]
    tols <- TOLSCALE_event * event_ses
    if ("(Intercept)" %in% names(tols)) # use weaker tolerance for intercept
      tols[["(Intercept)"]] <- 2 * tols[["(Intercept)"]] 
    for (i in 1:length(coef(s1)))
      expect_equal(coef(s1)[[i]], fixef(j1)$Event[[i+1]], tol = tols[[i]])        
  }  
    
  compare_glmer <- function(fm, fam) {
    y1 <- stan_glmer(fm, pbcLong, fam, iter = 2000, chains = CHAINS, seed = SEED)
    s1 <- coxph(Surv(futimeYears, death) ~ sex + trt, data = pbcSurv)
    j1 <- stan_jm(fm, pbcLong, Surv(futimeYears, death) ~ sex + trt, pbcSurv,
                  time_var = "year", assoc = NULL, family = fam,
                  iter = 2000, chains = CHAINS, seed = SEED)
    ranef_sds <- attr(VarCorr(y1)[[1]], "stddev")
    tols <- TOLSCALE_glmer_ranef * ranef_sds
    for (i in 1:length(tols))
      expect_equal(ranef(y1)[[1]][[i]], ranef(j1)$Long1[[1]][[i]], tol = tols[[i]])     
    fixef_ses <- y1$ses
    tols <- TOLSCALE_glmer_fixef * fixef_ses
    if ("(Intercept)" %in% names(tols)) # use weaker tolerance for intercept
      tols[["(Intercept)"]] <- 2 * tols[["(Intercept)"]]     
    for (i in 1:length(fixef(y1)))
      expect_equal(fixef(y1)[[i]], fixef(j1)$Long1[[i]], tol = tols[[i]])
    event_ses <- summary(s1)$coefficients[, "se(coef)"]
    tols <- TOLSCALE_event * event_ses
    if ("(Intercept)" %in% names(tols)) # use weaker tolerance for intercept
      tols[["(Intercept)"]] <- 2 * tols[["(Intercept)"]] 
    for (i in 1:length(coef(s1)))
      expect_equal(coef(s1)[[i]], fixef(j1)$Event[[i+1]], tol = tols[[i]])     
  } 
  
  pbcLong$ybino_successes <- as.integer(cut(pbcLong$logBili, 6))
  pbcLong$ybino_failures  <- 6 - as.integer(cut(pbcLong$logBili, 6)) 
  pbcLong$ybino_trials    <- rep_len(6, length(pbcLong$ybino_successes))
  pbcLong$ybino_prop      <- pbcLong$ybino_successes / pbcLong$ybino_trials
  pbcLong$ybern <- as.integer(pbcLong$logBili >= mean(pbcLong$logBili))
  pbcLong$ypois <- as.integer(pbcLong$albumin)
  pbcLong$ygamm <- as.integer(1 + (2 * pbcLong$platelet / 100))
  pbcLong$xbern <- as.numeric(pbcLong$platelet / 100)
  pbcLong$xpois <- as.numeric(pbcLong$platelet / 100)
  pbcLong$xgamm <- as.numeric(pbcLong$logBili)
 
  test_that("coefs same for stan_jm and stan_lmer/coxph", {
    compare_lmer(logBili ~ year + (1 | id))})
  test_that("coefs same for stan_jm and stan_glmer, binomial as cbind(success, failure)", {
    compare_glmer(cbind(ybino_successes, ybino_failures) ~ year + xbern + (1 | id), binomial)})
  test_that("coefs same for stan_jm and stan_glmer, binomial as cbind(success, trials-success)", {
    compare_glmer(cbind(ybino_successes, ybino_trials - ybino_successes) ~ year + xbern + (1 | id), binomial)})
  test_that("coefs same for stan_jm and stan_glmer, bernoulli", {
    compare_glmer(ybern ~ year + xbern + (1 | id), binomial)})
  test_that("coefs same for stan_jm and stan_glmer, poisson", {
    compare_glmer(ypois ~ year + xpois + (1 | id), poisson)})
  test_that("coefs same for stan_jm and stan_glmer, negative binomial", {
    compare_glmer(ypois ~ year + xpois + (1 | id), neg_binomial_2)})
  test_that("coefs same for stan_jm and stan_glmer, Gamma", {
    compare_glmer(ygamm ~ year + xgamm + (1 | id), Gamma)})
  test_that("coefs same for stan_jm and stan_glmer, inverse gaussian", {
    compare_glmer(ygamm ~ year + xgamm + (1 | id), inverse.gaussian)})

  test_that("coefs same for stan_jm and stan_glmer, binomial as prop as outcome and trials as weights", {
    jm_weights <- data.frame(id = pbcLong$id, weights = pbcLong$ybino_trials)
    y1 <- stan_glmer(ybino_prop ~ year + (1 | id), pbcLong, weights = ybino_trials, 
                     iter = 2000, chains = CHAINS, seed = SEED)
    s1 <- coxph(Surv(futimeYears, death) ~ sex + trt, data = pbcSurv)
    j1 <- stan_jm(ybino_prop ~ year + (1 | id), data = pbcLong, weights = jm_weights,
                  formulaEvent = Surv(futimeYears, death) ~ sex + trt,
                  dataLong = pbcLong, dataEvent = pbcSurv,
                  time_var = "year", assoc = NULL,
                  iter = 2000, chains = CHAINS, seed = SEED)        
    expect_equal(ranef(y1), ranef(j1)$Long1, tol = RANEF_tol)
    tols <- TOLSCALE_lmer * y1$ses
    for (i in 1:length(fixef(y1)))
      expect_equal(fixef(y1)[[i]], fixef(j1)$Long1[[i]], tol = tols[[i]])
    tols <- TOLSCALE_event * summary(s1)$coefficients[, "se(coef)"]
    for (i in 1:length(coef(s1)))
      expect_equal(coef(s1)[[i]], fixef(j1)$Event[[i+1]], tol = tols[[i]])       
  })  
  
}

#--------  Check stanjm prediction functions work with various formula specifications

if (interactive()) {
  
f1 <- stan_jm(formulaLong = logBili ~ year + (1 | id), 
              dataLong = pbcLong,
              formulaEvent = Surv(futimeYears, death) ~ sex + trt, 
              dataEvent = pbcSurv,
              time_var = "year",
              # this next line is only to keep the example small in size!
              chains = 1, cores = 1, seed = 12345, iter = 10)

f2 <- stan_jm(formulaLong = exp(logBili) ~ year + (1 | id), 
              dataLong = pbcLong,
              formulaEvent = Surv(futimeYears, death) ~ sex + trt, 
              dataEvent = pbcSurv,
              time_var = "year",
              # this next line is only to keep the example small in size!
              chains = 1, cores = 1, seed = 12345, iter = 10)

f3 <- stan_jm(formulaLong = logBili ~ poly(year, degree = 2) + (1 | id), 
              dataLong = pbcLong,
              formulaEvent = Surv(futimeYears, death) ~ sex + trt, 
              dataEvent = pbcSurv,
              time_var = "year",
              # this next line is only to keep the example small in size!
              chains = 1, cores = 1, seed = 12345, iter = 10)

f4 <- stan_jm(formulaLong = exp(logBili) ~ poly(year, degree = 2) + (1 | id), 
              dataLong = pbcLong,
              formulaEvent = Surv(futimeYears, death) ~ sex + trt, 
              dataEvent = pbcSurv,
              time_var = "year",
              # this next line is only to keep the example small in size!
              chains = 1, cores = 1, seed = 12345, iter = 10)

pbcLong$trials <- rpois(nrow(pbcLong), 6)
pbcLong$succ <- rbinom(nrow(pbcLong), pbcLong$trials, .7)
pbcLong$fail <- pbcLong$trials - pbcLong$succ
f5 <- stan_jm(formulaLong = cbind(succ, fail) ~ poly(year, degree = 2) + (1 | id), 
              dataLong = pbcLong,
              formulaEvent = Surv(futimeYears, death) ~ sex + trt, 
              dataEvent = pbcSurv,
              time_var = "year", family = binomial,
              # this next line is only to keep the example small in size!
              chains = 1, cores = 1, seed = 12345, iter = 10, init = 0)

for (j in 1:4) {
  mod <- get(paste0("f", j))
  
  test_that("log_lik works with estimation data", {
    ll <- log_lik(mod)
    expect_matrix(ll)
  })
  test_that("posterior_survfit works with estimation data", {
    ps <- posterior_survfit(mod)
    expect_survfit(ps)
  })
  test_that("posterior_predict works with estimation data", {
    pp <- posterior_predict(mod, m = 1)
    expect_ppd(pp)
  }) 
  
  ndL <- pbcLong[pbcLong$id == 2,]
  ndE <- pbcSurv[pbcSurv$id == 2,]
  test_that("log_lik works with new data (one individual)", {
    ll <- log_lik(mod, newdataLong = ndL, newdataEvent = ndE)
    expect_matrix(ll)
  })
  test_that("posterior_survfit works with new data (one individual)", {
    ps <- posterior_survfit(mod, newdataLong = ndL, newdataEvent = ndE)
    expect_survfit(ps)
  })  
  test_that("posterior_predict works with new data (one individual)", {
    pp <- posterior_predict(mod, m = 1, newdataLong = ndL, newdataEvent = ndE)
    expect_ppd(pp)
  })  
  
  ndL <- pbcLong[pbcLong$id %in% c(1,2),]
  ndE <- pbcSurv[pbcSurv$id %in% c(1,2),]
  test_that("log_lik works with new data (multiple individuals)", {
    ll <- log_lik(mod, newdataLong = ndL, newdataEvent = ndE)
    expect_matrix(ll)
  })
  test_that("posterior_survfit works with new data (multiple individuals)", {
    ps <- posterior_survfit(mod, newdataLong = ndL, newdataEvent = ndE)
    expect_survfit(ps)
  })
  test_that("posterior_predict works with new data (multiple individuals)", {
    pp <- posterior_predict(mod, m = 1, newdataLong = ndL, newdataEvent = ndE)
    expect_ppd(pp)
  })
  
}

}



