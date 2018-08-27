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
library(survival)
ITER <- 1000
CHAINS <- 1
SEED <- 12345
REFRESH <- 0L
set.seed(SEED)
if (interactive()) 
  options(mc.cores = parallel::detectCores())

TOLSCALES <- list(
  lmer_fixef = 0.25, # how many SEs can stan_jm fixefs be from lmer fixefs
  lmer_ranef = 0.05, # how many SDs can stan_jm ranefs be from lmer ranefs
  glmer_fixef = 0.5, # how many SEs can stan_jm fixefs be from glmer fixefs
  glmer_ranef = 0.1, # how many SDs can stan_jm ranefs be from glmer ranefs
  event = 0.3        # how many SEs can stan_jm fixefs be from coxph fixefs
)

source(test_path("helpers", "expect_matrix.R"))
source(test_path("helpers", "expect_stanreg.R"))
source(test_path("helpers", "expect_stanmvreg.R"))
source(test_path("helpers", "expect_survfit.R"))
source(test_path("helpers", "expect_ppd.R"))
source(test_path("helpers", "expect_equivalent_loo.R"))
source(test_path("helpers", "SW.R"))
# SW <- function(expr) eval(expr)
source(test_path("helpers", "get_tols.R"))
source(test_path("helpers", "recover_pars.R"))

context("stan_jm")

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

# univariate joint model
fmLong1 <- logBili ~ year + (year | id)
fmSurv1 <- Surv(futimeYears, death) ~ sex + trt
jm1 <- stan_jm(
  fmLong1, pbcLong, fmSurv1, pbcSurv, time_var = "year", 
  iter = 1, refresh = 0, chains = 1, seed = SEED)

# multivariate joint model
fmLong2 <- list(
  logBili ~ year + (year | id),
  albumin ~ year + (year | id))
fmSurv2 <- Surv(futimeYears, death) ~ sex + trt
jm2 <- stan_jm(
  fmLong2, pbcLong, fmSurv2, pbcSurv, time_var = "year", 
  iter = 1, refresh = 0, chains = 1, seed = SEED)

#----  Tests for stan_jm arguments

test_that("formula argument works", {
  expect_identical(as.matrix(jm1), as.matrix(update(jm1, formulaLong. = list(fmLong1)))) # fm as list
})

test_that("data argument works", {
  expect_identical(as.matrix(jm1), 
                   as.matrix(update(jm1, dataLong = list(pbcLong)))) # data as list
  expect_identical(as.matrix(jm2), 
                   as.matrix(update(jm2, dataLong = list(pbcLong, pbcLong))))
})

test_that("id_var argument works", {
  
  # Models with a single grouping factor
  expect_output(update(jm1, id_var = "id"))
  expect_output(expect_warning(update(jm1, id_var = "year"), 
                               "are not the same; 'id_var' will be ignored"))
  
  # Models with more than one grouping factor
  tmpdat <- pbcLong
  tmpdat$practice <- cut(pbcLong$id, c(0,10,20,30,40))
  tmpfm <- logBili ~ year + (year | id) + (1 | practice)
  ok_mod <- update(jm1, formulaLong. = tmpfm, dataLong = tmpdat, id_var = "id", init = 0)
  expect_stanmvreg(ok_mod)
  expect_error(update(ok_mod, id_var = NULL), "'id_var' must be specified")
  expect_error(update(ok_mod, id_var = "year"), "'id_var' must be included as a grouping factor")
})

test_that("family argument works", {
  
  expect_output(ret <- update(jm1, family = "gaussian"))
  expect_output(ret <- update(jm1, family = gaussian))
  expect_output(ret <- update(jm1, family = gaussian(link = identity)))
  
  expect_output(ret <- update(jm1, formulaLong. = ypois ~ ., family = poisson, init = 0))
  expect_output(ret <- update(jm1, formulaLong. = ynbin ~ ., family = neg_binomial_2))
  #expect_output(ret <- update(jm1, formulaLong. = ygamm ~ ., family = Gamma))
  #expect_output(ret <- update(jm1, formulaLong. = ygamm ~ ., family = inverse.gaussian))

  expect_error(ret <- update(jm1, formulaLong. = ybino ~ ., family = binomial))
})

test_that("assoc argument works", {
  
  # NB: muslope, shared_b, and shared_coef have been temporarily 
  # disallowed, and will be reinstated in a future release

  expect_error(ret <- update(jm1, assoc = "muslope"), "temporarily disallowed")
  expect_error(ret <- update(jm1, assoc = "shared_b"), "temporarily disallowed")
  expect_error(ret <- update(jm1, assoc = "shared_coef"), "temporarily disallowed")
      
  # Univariate joint models
  
  expect_output(ret <- update(jm1, assoc = NULL))
  expect_output(ret <- update(jm1, assoc = "null"))
  expect_output(ret <- update(jm1, assoc = "etavalue"))
  expect_output(ret <- update(jm1, assoc = "muvalue"))
  expect_output(ret <- update(jm1, assoc = "etaslope"))
  #expect_output(ret <- update(jm1, assoc = "muslope"))
  expect_output(ret <- update(jm1, assoc = "etaauc"))
  expect_output(ret <- update(jm1, assoc = "muauc"))
  expect_output(ret <- update(jm1, assoc = c("etavalue", "etaslope"))) 
  #expect_output(ret <- update(jm1, assoc = c("etavalue", "muslope"))) 
  expect_output(ret <- update(jm1, assoc = c("etavalue", "etaauc"))) 
  expect_output(ret <- update(jm1, assoc = c("etavalue", "muauc"))) 
  expect_output(ret <- update(jm1, assoc = c("muvalue", "etaslope"))) 
  #expect_output(ret <- update(jm1, assoc = c("muvalue", "muslope")))
  expect_output(ret <- update(jm1, assoc = c("muvalue", "etaauc"))) 
  expect_output(ret <- update(jm1, assoc = c("muvalue", "muauc"))) 
  
  expect_error(ret <- update(jm1, assoc = c("etavalue", "muvalue")), "cannot be specified together")
  #expect_error(ret <- update(jm1, assoc = c("etaslope", "muslope")), "cannot be specified together")
  expect_error(ret <- update(jm1, assoc = c("etaauc", "muauc")), "cannot be specified together")
  
  #expect_output(ret <- update(jm1, assoc = "shared_b"))
  #expect_output(ret <- update(jm1, assoc = "shared_b(1)"))
  #expect_output(ret <- update(jm1, assoc = "shared_b(2)"))
  #expect_output(ret <- update(jm1, assoc = "shared_b(1:2)"))
  #expect_output(ret <- update(jm1, assoc = "shared_b(1,2)"))
  
  #expect_output(ret <- update(jm1, assoc = "shared_coef"))
  #expect_output(ret <- update(jm1, assoc = "shared_coef(1)"))
  #expect_output(ret <- update(jm1, assoc = "shared_coef(2)"))
  #expect_output(ret <- update(jm1, assoc = "shared_coef(1:2)"))
  #expect_output(ret <- update(jm1, assoc = "shared_coef(1,2)"))
  
  #expect_error(ret <- update(jm1, assoc = "shared_b(10)"), "greater than the number of")
  #expect_error(ret <- update(jm1, assoc = "shared_coef(10)"), "greater than the number of")
  #expect_error(ret <- update(jm1, assoc = c("shared_b(1)", "shared_coef(1)")), "should not be specified in both")
  #expect_error(ret <- update(jm1, assoc = c("shared_b", "shared_coef")), "should not be specified in both")
  
  expect_output(ret <- update(jm1, assoc = list(NULL)))
  expect_output(ret <- update(jm1, assoc = list("null")))
  expect_output(ret <- update(jm1, assoc = list("etavalue")))
  expect_output(ret <- update(jm1, assoc = list("muvalue")))
  expect_output(ret <- update(jm1, assoc = list("etaslope")))
  #expect_output(ret <- update(jm1, assoc = list("muslope")))
  expect_output(ret <- update(jm1, assoc = list("etaauc")))
  expect_output(ret <- update(jm1, assoc = list("muauc")))
  expect_output(ret <- update(jm1, assoc = list(c("etavalue", "etaslope")))) 
  #expect_output(ret <- update(jm1, assoc = list(c("etavalue", "muslope")))) 
  expect_output(ret <- update(jm1, assoc = list(c("muvalue", "etaslope")))) 
  #expect_output(ret <- update(jm1, assoc = list(c("muvalue", "muslope"))))  
  
  expect_error(ret <- update(jm1, assoc = NA), "'assoc' should be") 
  expect_error(ret <- update(jm1, assoc = 123), "'assoc' should be") 
  expect_error(ret <- update(jm1, assoc = c(1,2,3)), "'assoc' should be") 
  
  expect_error(ret <- update(jm1, assoc = c("wrong")), "unsupported association type") 
  expect_error(ret <- update(jm1, assoc = list("wrong")), "unsupported association type") 
  
  expect_error(ret <- update(jm1, assoc = list(NULL, NULL)), "incorrect length") 
  expect_error(ret <- update(jm1, assoc = list("etavalue", "etavalue")), "incorrect length") 
  expect_error(ret <- update(jm1, assoc = list(c("etavalue", "etaslope"), "etavalue")), "incorrect length") 
  
  # Multivariate joint models
  
  expect_output(ret <- update(jm2, assoc = "etavalue"))
  expect_output(ret <- update(jm2, assoc = "muvalue"))
  expect_output(ret <- update(jm2, assoc = "etaslope"))
  #expect_output(ret <- update(jm2, assoc = "muslope"))
  expect_output(ret <- update(jm2, assoc = "etaauc"))
  expect_output(ret <- update(jm2, assoc = "muauc"))
  expect_output(ret <- update(jm2, assoc = c("etavalue", "etaslope"))) 
  expect_output(ret <- update(jm2, assoc = c("etavalue", "etaauc"))) 
  expect_output(ret <- update(jm2, assoc = c("etaslope", "etaauc"))) 

  expect_output(ret <- update(jm2, assoc = list("etavalue")))
  expect_output(ret <- update(jm2, assoc = list("etavalue", "etavalue")))
  expect_output(ret <- update(jm2, assoc = list(c("etavalue", "etaslope"), "etavalue")))
  expect_output(ret <- update(jm2, assoc = list("etavalue", c("etavalue", "etaslope"))))
  expect_output(ret <- update(jm2, assoc = list(c("etavalue", "etaslope"), c("muvalue", "muauc"))))
  
  expect_error(ret <- update(jm2, assoc = list("wrong", "etavalue")), "unsupported association type")
  expect_error(ret <- update(jm2, assoc = list("null", "etavalue", "etaslope")), "incorrect length")
  expect_error(ret <- update(jm2, assoc = data.frame("etavalue", "etaslope")), "'assoc' should be") 
  
})

test_that("basehaz argument works", {
  
  expect_output(update(jm1, basehaz = "weibull"))
  expect_output(update(jm1, basehaz = "bs"))
  expect_output(update(jm1, basehaz = "piecewise"))
  
  expect_output(update(jm1, basehaz = "bs", basehaz_ops = list(df = 5)))
  expect_output(update(jm1, basehaz = "bs", basehaz_ops = list(knots = c(1,3,5))))
  expect_output(update(jm1, basehaz = "piecewise", basehaz_ops = list(df = 5)))
  expect_output(update(jm1, basehaz = "piecewise", basehaz_ops = list(knots = c(1,3,5))))
  
  expect_output(expect_warning(update(jm1, basehaz = "weibull", basehaz_ops = list(df = 1)), "'df' will be ignored"))
  expect_output(expect_warning(update(jm1, basehaz = "weibull", basehaz_ops = list(knots = 1)), "'knots' will be ignored"))
  
  expect_output(update(jm1, basehaz = "piecewise", basehaz_ops = list(knots = c(1,3,5))))
  
  expect_error(update(jm1, basehaz = "bs", basehaz_ops = list(df = 1)), "must be at least 3")
  expect_error(update(jm1, basehaz = "bs", basehaz_ops = list(knots = -1)), "'knots' must be non-negative")
  expect_error(update(jm1, basehaz = "piecewise", basehaz_ops = list(knots = -1)), "'knots' must be non-negative")
  expect_error(update(jm1, basehaz = "piecewise", basehaz_ops = list(knots = c(1,2,50))), "cannot be greater than the largest event time")
  
})

test_that("qnodes argument works", {
  
  expect_output(update(jm1, qnodes = 7))
  expect_output(update(jm1, qnodes = 11))
  expect_output(update(jm1, qnodes = 15))
  
  expect_error(update(jm1, qnodes = 1), "'qnodes' must be either 7, 11 or 15")
  expect_error(update(jm1, qnodes = c(1,2)), "should be a numeric vector of length 1")
  expect_error(update(jm1, qnodes = "wrong"), "should be a numeric vector of length 1")
  
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
  
  expect_error(update(jm1, weights = wts0, iter = 5), "not yet implemented")
  
  #expect_output(update(jm1, weights = wts0, iter = 5))
  #expect_output(update(jm1, weights = wts7, iter = 5)) # ok to supply extra IDs in weights
  
  #expect_error(update(jm1, weights = as.matrix(wts0)), "should be a data frame")
  #expect_error(update(jm1, weights = wts1), "do not have weights supplied")
  #expect_error(update(jm1, weights = wts2), "should only have one row")
  #expect_error(update(jm1, weights = wts3), "should be a data frame with two columns")
  #expect_error(update(jm1, weights = wts4), "weights supplied must be numeric")
  #expect_error(update(jm1, weights = wts5), "weights supplied must be numeric")
  #expect_error(update(jm1, weights = wts6), "Negative weights are not allowed")
  
})

test_that("init argument works", {
  expect_output(update(jm1, init = "prefit"))
  expect_output(update(jm1, init = "0"))
  expect_output(update(jm1, init = 0))
  expect_output(update(jm1, init = "random"))
})

test_that("prior_PD argument works", {
  expect_output(update(jm1, prior_PD = TRUE))
})

test_that("adapt_delta argument works", {
  expect_output(update(jm1, adapt_delta = NULL))
  expect_output(update(jm1, adapt_delta = 0.8))
  expect_output(update(jm1, control = list(adapt_delta = NULL)))
  expect_output(update(jm1, control = list(adapt_delta = 0.8)))
})

test_that("max_treedepth argument works", {
  expect_output(update(jm1, max_treedepth = NULL))
  expect_output(update(jm1, max_treedepth = 5))
  expect_output(update(jm1, control = list(max_treedepth = NULL)))
  expect_output(update(jm1, control = list(max_treedepth = 5)))
})

test_that("error message occurs for arguments not implemented", {
  expect_error(update(jm1, QR = TRUE), "not yet implemented")
  expect_error(update(jm1, sparse = TRUE), "not yet implemented")
})

#----  Compare parameter estimates: stan_jm(assoc = NULL) vs stan_glmer/coxph

compare_glmer <- function(fmLong, fam = gaussian, ...) {
  require(survival)
  fmSurv <- Surv(futimeYears, death) ~ sex + trt
  y1 <- stan_glmer(fmLong, pbcLong, fam, iter = 1000, chains = CHAINS, seed = SEED)
  s1 <- coxph(fmSurv, data = pbcSurv)
  j1 <- stan_jm(fmLong, pbcLong, fmSurv, pbcSurv, time_var = "year", family = fam, 
                assoc = NULL, iter = 1000, chains = CHAINS, seed = SEED, ...) 
  tols <- get_tols(y1, s1, tolscales = TOLSCALES)
  pars <- recover_pars(y1, s1)
  parsjm <- recover_pars(j1)
  for (i in names(tols$fixef))
    expect_equal(pars$fixef[[i]], parsjm$fixef[[i]], tol = tols$fixef[[i]], info = fam)     
  for (i in names(tols$ranef))
    expect_equal(pars$ranef[[i]], parsjm$ranef[[i]], tol = tols$ranef[[i]], info = fam)
  for (i in names(tols$event))
    expect_equal(pars$event[[i]], parsjm$event[[i]], tol = tols$event[[i]], info = fam)
}

# test_that("coefs same for stan_jm and stan_lmer/coxph", {
#   compare_glmer(logBili ~ year + (1 | id), gaussian)})
# test_that("coefs same for stan_jm and stan_glmer, bernoulli", {
#   compare_glmer(ybern ~ year + xbern + (1 | id), binomial)})
# test_that("coefs same for stan_jm and stan_glmer, poisson", {
#   compare_glmer(ypois ~ year + xpois + (1 | id), poisson, init = 0)})
# test_that("coefs same for stan_jm and stan_glmer, negative binomial", {
#   compare_glmer(ynbin ~ year + xpois + (1 | id), neg_binomial_2)})
# test_that("coefs same for stan_jm and stan_glmer, Gamma", {
#   compare_glmer(ygamm ~ year + xgamm + (1 | id), Gamma(log))})
#test_that("coefs same for stan_jm and stan_glmer, inverse gaussian", {
#  compare_glmer(ygamm ~ year + xgamm + (1 | id), inverse.gaussian)})  

#--------  Check (post-)estimation functions work with various model specifications

# No functions in formula
o<-SW(f1 <- stan_jm(formulaLong = logBili ~ year + (year | id), 
                    dataLong = pbcLong,
                    formulaEvent = Surv(futimeYears, death) ~ sex + trt, 
                    dataEvent = pbcSurv,
                    time_var = "year",
                    # this next line is only to keep the example small in size!
                    chains = 1, cores = 1, seed = 12345, iter = 10))

# Functions on LHS of formula
o<-SW(f2 <- update(f1, formulaLong. = exp(logBili) ~ year + (year | id)))

# Functions on RHS of formula
# o<-SW(f3 <- update(f1, formulaLong. = logBili ~ poly(year, degree = 2) + (poly(year, degree = 2) | id)))

# Functions on LHS and RHS of formula
# o<-SW(f4 <- update(f1, formulaLong. = exp(logBili) ~ poly(year, degree = 2) + (poly(year, degree = 2) | id)))

# Intercept only event submodel
o<-SW(f5 <- update(f1, formulaEvent. = Surv(futimeYears, death) ~ 1))

# Different baseline hazards
o<-SW(f7 <- update(f1, basehaz = "weibull"))
o<-SW(f8 <- update(f1, basehaz = "bs"))
#o<-SW(f9 <- update(f1, basehaz = "piecewise")) # posterior_survfit not yet implemented for piecewise

# Different association structures
o<-SW(f10 <- update(f1, assoc = NULL))
o<-SW(f11 <- update(f1, assoc = "etavalue"))
o<-SW(f12 <- update(f1, assoc = "etaslope"))
o<-SW(f13 <- update(f1, assoc = "etaauc"))
o<-SW(f14 <- update(f1, assoc = "muvalue"))
#o<-SW(f15 <- update(f1, assoc = "muslope"))
o<-SW(f16 <- update(f1, assoc = "muauc"))
o<-SW(f17 <- update(f1, assoc = c("etavalue", "etaslope")))
o<-SW(f18 <- update(f1, assoc = c("etavalue", "etaauc")))

# Different association structures with intercept only submodel
o<-SW(f19 <- update(f5, assoc = NULL))
o<-SW(f20 <- update(f5, assoc = "etavalue"))
o<-SW(f21 <- update(f5, assoc = "etaslope"))
o<-SW(f22 <- update(f5, assoc = "etaauc", iter = 100))
o<-SW(f23 <- update(f5, assoc = "muvalue"))
#o<-SW(f24 <- update(f5, assoc = "muslope"))
o<-SW(f25 <- update(f5, assoc = "muauc"))
o<-SW(f26 <- update(f5, assoc = c("etavalue", "etaslope")))
o<-SW(f27 <- update(f5, assoc = c("etavalue", "etaauc")))  

# Shared random effect association structures
#o<-SW(f28 <- update(f1, assoc = c("shared_b")))
#o<-SW(f29 <- update(f1, assoc = c("shared_coef")))

# Multivariate models
o<-SW(f31 <- stan_jm(formulaLong = list(logBili ~ year + (year | id),
                                        albumin ~ sex + trt + year + (1 | id)),
                     dataLong = pbcLong,
                     formulaEvent = Surv(futimeYears, death) ~ sex + trt, 
                     dataEvent = pbcSurv,
                     time_var = "year",
                     # this next line is only to keep the example small in size!
                     chains = 1, cores = 1, seed = 12345, iter = 50))
o<-SW(f32 <- update(f31, assoc = list("etaslope", c("etavalue", "etaauc"))))  

# New data for predictions
ndL1 <- pbcLong[pbcLong$id == 2,]
ndE1 <- pbcSurv[pbcSurv$id == 2,]
ndL2 <- pbcLong[pbcLong$id %in% c(1,2),]
ndE2 <- pbcSurv[pbcSurv$id %in% c(1,2),]

# Test the models
for (j in c(1:30)) {
  mod <- try(get(paste0("f", j)), silent = TRUE)
  if (class(mod)[1L] == "try-error") {
    cat("Model not found:", paste0("f", j), "\n")
  } else {
    cat("Checking model:", paste0("f", j), "\n")
    
    test_that("log_lik works with estimation data", {
      ll <- log_lik(mod)
      expect_matrix(ll)
      expect_error(log_lik(mod, m = 1), "should not be specified")
    })
    
    test_that("log_lik works with new data (one individual)", {
      ll <- log_lik(mod, newdataLong = ndL1, newdataEvent = ndE1)
      expect_matrix(ll)
    })
    
    test_that("log_lik works with new data (multiple individuals)", {
      ll <- log_lik(mod, newdataLong = ndL2, newdataEvent = ndE2)
      expect_matrix(ll)
    })
    
    test_that("loo and waic work", {
      expect_equivalent_loo(mod)
    })
    
    test_that("posterior_predict works with estimation data", {
      pp <- posterior_predict(mod, seed = SEED)
      expect_ppd(pp)
      expect_identical(pp, posterior_predict(mod, m = 1, seed = SEED))
      if (mod$n_markers > 1L) {
        pp <- posterior_predict(mod, m = 2)
        expect_ppd(pp)
      }  
    }) 
    
    test_that("posterior_predict works with new data (one individual)", {
      pp <- posterior_predict(mod, newdata = ndL1, seed = SEED)
      expect_ppd(pp)
      expect_identical(pp, posterior_predict(mod, m = 1, newdata = ndL1, seed = SEED))
      if (mod$n_markers > 1L) {
        pp <- posterior_predict(mod, m = 2, newdata = ndL1)
        expect_ppd(pp)
      }
      expect_error(posterior_predict(mod, newdataLong = ndL1), "should not be specified")
      expect_error(posterior_predict(mod, newdataEvent = ndE1), "should not be specified")
    })  
    
    test_that("posterior_predict works with new data (multiple individuals)", {
      pp <- posterior_predict(mod, newdata = ndL2, seed = SEED)
      expect_ppd(pp)
      expect_identical(pp, posterior_predict(mod, m = 1, newdata = ndL2, seed = SEED))
      if (mod$n_markers > 1L) {
        pp <- posterior_predict(mod, m = 2, newdata = ndL2)
        expect_ppd(pp)
      }            
    })
    
    test_that("posterior_traj works with estimation data", {
      pp <- posterior_traj(mod)
      expect_s3_class(pp, "predict.stanjm")
      if (mod$n_markers > 1L) {
        pp <- posterior_traj(mod, m = 2)
        expect_s3_class(pp, "predict.stanjm")
      }  
    }) 
    
    test_that("posterior_traj works with new data (one individual)", {
      pp <- posterior_traj(mod, newdataLong = ndL1, dynamic = FALSE)
      expect_s3_class(pp, "predict.stanjm")
      if (mod$n_markers > 1L) {
        pp <- posterior_traj(mod, m = 2, newdataLong = ndL1, dynamic = FALSE)
        expect_s3_class(pp, "predict.stanjm")
      }
    })  
    
    test_that("posterior_traj works with new data (multiple individuals)", {
      pp <- posterior_traj(mod, newdataLong = ndL2, dynamic = FALSE)
      expect_s3_class(pp, "predict.stanjm")
      if (mod$n_markers > 1L) {
        pp <- posterior_traj(mod, m = 2, newdataLong = ndL2, dynamic = FALSE)
        expect_s3_class(pp, "predict.stanjm")
      }            
    })        
    
    test_that("posterior_survfit works with estimation data", {
      SW(ps <- posterior_survfit(mod))
      expect_survfit(ps)
    })
    
    test_that("posterior_survfit works with new data (one individual)", {
      SW(ps <- posterior_survfit(mod, newdataLong = ndL1, newdataEvent = ndE1))
      expect_survfit(ps)
    })  
    
    test_that("posterior_survfit works with new data (multiple individuals)", {
      SW(ps <- posterior_survfit(mod, newdataLong = ndL2, newdataEvent = ndE2))
      expect_survfit(ps)
    })
  }
}
