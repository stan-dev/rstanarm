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
CHAINS <- 1
SEED <- 12345
REFRESH <- ITER
set.seed(SEED)
if (interactive()) 
  options(mc.cores = parallel::detectCores())

TOLSCALES <- list(
  lmer_fixef = 0.25,  # how many SEs can stan_jm fixefs be from lmer fixefs
  lmer_ranef = 0.05, # how many SDs can stan_jm ranefs be from lmer ranefs
  glmer_fixef = 0.3, # how many SEs can stan_jm fixefs be from glmer fixefs
  glmer_ranef = 0.1, # how many SDs can stan_jm ranefs be from glmer ranefs
  event = 0.2        # how many SEs can stan_jm fixefs be from coxph fixefs
)

expect_matrix  <- function(x) expect_identical(class(x), "matrix")
expect_survfit <- function(x) expect_s3_class(x, "survfit.stanmvreg")
expect_ppd     <- function(x) expect_s3_class(x, "ppd")
source(file.path("helpers", "expect_stanreg.R"))
source(file.path("helpers", "expect_stanmvreg.R"))
source(file.path("helpers", "SW.R"))
source(file.path("helpers", "get_tols.R"))
source(file.path("helpers", "recover_pars.R"))

context("stan_jm")

#----  Data (for non-Gaussian families)

pbcLong$ybern <- as.integer(pbcLong$logBili >= mean(pbcLong$logBili))
pbcLong$ybino <- as.integer(rpois(nrow(pbcLong), 5))
pbcLong$ypois <- as.integer(pbcLong$albumin)
pbcLong$ygamm <- as.integer(1 + (2 * pbcLong$platelet / 100))
pbcLong$xbern <- as.numeric(pbcLong$platelet / 100)
pbcLong$xpois <- as.numeric(pbcLong$platelet / 100)
pbcLong$xgamm <- as.numeric(pbcLong$logBili)

#----  Models

# univariate joint model
fmLong <- logBili ~ year + (year | id)
fmSurv <- Surv(futimeYears, death) ~ sex + trt
jm1 <- stan_jm(
  fmLong, pbcLong, fmSurv, pbcSurv, time_var = "year", 
  iter = 1, chains = 1, seed = SEED)

# multivariate joint model
jm2 <- update(
  jm1, formulaLong. = list(
    logBili ~ year + (year | id), 
    albumin ~ year + (year | id)))

#----  Tests for stan_jm arguments

test_that("formula argument works", {
  expect_identical(jm1, update(jm1, formulaLong. = list(fmLong))) # fm as list
})

test_that("data argument works", {
  expect_identical(jm1, update(jm1, dataLong = list(pbcLong))) # data as list
  expect_identical(jm2, update(jm2, dataLong = list(pbcLong, pbcLong)))
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
  expect_error(update(ok_mod, id_var = "practice"), "'id_var' must correspond to the lowest level of clustering")
})

test_that("family argument works", {
  
  expect_output(ret <- update(jm1, family = "gaussian"))
  expect_output(ret <- update(jm1, family = gaussian))
  expect_output(ret <- update(jm1, family = gaussian(link = identity)))
  
  expect_output(ret <- update(jm1, formulaLong. = ybern ~ ., family = binomial))
  expect_output(ret <- update(jm1, formulaLong. = ypois ~ ., family = poisson))
  expect_output(ret <- update(jm1, formulaLong. = ypois ~ ., family = neg_binomial_2))
  expect_output(ret <- update(jm1, formulaLong. = ygamm ~ ., family = Gamma))
  expect_output(ret <- update(jm1, formulaLong. = ygamm ~ ., family = inverse.gaussian))
  
  expect_error(ret <- update(jm1, formulaLong. = ybino ~ ., family = binomial))
})

test_that("assoc argument works", {
  
  # Univariate joint models
  
  expect_output(ret <- update(jm1, assoc = NULL))
  expect_output(ret <- update(jm1, assoc = "null"))
  expect_output(ret <- update(jm1, assoc = "etavalue"))
  expect_output(ret <- update(jm1, assoc = "etaslope"))
  expect_output(ret <- update(jm1, assoc = "muvalue"))
  expect_output(ret <- update(jm1, assoc = "muslope"))
  expect_output(ret <- update(jm1, assoc = c("etavalue", "etaslope"))) 
  expect_output(ret <- update(jm1, assoc = c("etavalue", "muslope"))) 
  expect_output(ret <- update(jm1, assoc = c("muvalue", "etaslope"))) 
  expect_output(ret <- update(jm1, assoc = c("muvalue", "muslope")))
 
  expect_error(ret <- update(jm1, assoc = c("etavalue", "muvalue")), "cannot be specified together")
  expect_error(ret <- update(jm1, assoc = c("etaslope", "muslope")), "cannot be specified together")
  
  expect_output(ret <- update(jm1, assoc = "shared_b"))
  expect_output(ret <- update(jm1, assoc = "shared_b(1)"))
  expect_output(ret <- update(jm1, assoc = "shared_b(2)"))
  expect_output(ret <- update(jm1, assoc = "shared_b(1:2)"))
  expect_output(ret <- update(jm1, assoc = "shared_b(1,2)"))
  
  expect_output(ret <- update(jm1, assoc = "shared_coef"))
  expect_output(ret <- update(jm1, assoc = "shared_coef(1)"))
  expect_output(ret <- update(jm1, assoc = "shared_coef(2)"))
  expect_output(ret <- update(jm1, assoc = "shared_coef(1:2)"))
  expect_output(ret <- update(jm1, assoc = "shared_coef(1,2)"))
  
  expect_error(ret <- update(jm1, assoc = "shared_b(10)"), "greater than the number of")
  expect_error(ret <- update(jm1, assoc = "shared_coef(10)"), "greater than the number of")
  expect_error(ret <- update(jm1, assoc = c("shared_b(1)", "shared_coef(1)")), "should not be specified in both")
  expect_error(ret <- update(jm1, assoc = c("shared_b", "shared_coef")), "should not be specified in both")
  
  expect_output(ret <- update(jm1, assoc = list(NULL)))
  expect_output(ret <- update(jm1, assoc = list("null")))
  expect_output(ret <- update(jm1, assoc = list("etavalue")))
  expect_output(ret <- update(jm1, assoc = list("etaslope")))
  expect_output(ret <- update(jm1, assoc = list("muvalue")))
  expect_output(ret <- update(jm1, assoc = list("muslope")))
  expect_output(ret <- update(jm1, assoc = list(c("etavalue", "etaslope")))) 
  expect_output(ret <- update(jm1, assoc = list(c("etavalue", "muslope")))) 
  expect_output(ret <- update(jm1, assoc = list(c("muvalue", "etaslope")))) 
  expect_output(ret <- update(jm1, assoc = list(c("muvalue", "muslope"))))  
  
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
  expect_output(ret <- update(jm2, assoc = "etaslope"))
  expect_output(ret <- update(jm2, assoc = "muvalue"))
  expect_output(ret <- update(jm2, assoc = "muslope"))
  expect_output(ret <- update(jm2, assoc = c("etavalue", "etaslope"))) 
  expect_output(ret <- update(jm2, assoc = c("etavalue", "muslope"))) 
  expect_output(ret <- update(jm2, assoc = c("muvalue", "etaslope"))) 
  expect_output(ret <- update(jm2, assoc = c("muvalue", "muslope")))
  
  expect_output(ret <- update(jm2, assoc = list("etavalue")))
  expect_output(ret <- update(jm2, assoc = list("etavalue", "etavalue")))
  expect_output(ret <- update(jm2, assoc = list(c("etavalue", "etaslope"), "etavalue")))
  expect_output(ret <- update(jm2, assoc = list("etavalue", c("etavalue", "etaslope"))))
  expect_output(ret <- update(jm2, assoc = list(c("etavalue", "etaslope"), c("muvalue", "muslope"))))
  
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

test_that("quadnodes argument works", {
  
  expect_output(update(jm1, quadnodes = 7))
  expect_output(update(jm1, quadnodes = 11))
  expect_output(update(jm1, quadnodes = 15))
  
  expect_error(update(jm1, quadnodes = 1), "'quadnodes' must be either 7, 11 or 15")
  expect_error(update(jm1, quadnodes = c(1,2)), "should be a numeric vector of length 1")
  expect_error(update(jm1, quadnodes = "wrong"), "should be a numeric vector of length 1")
  
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
  
  expect_output(update(jm1, weights = wts0, iter = 5))
  expect_output(update(jm1, weights = wts7, iter = 5)) # ok to supply extra IDs in weights
  
  expect_error(update(jm1, weights = as.matrix(wts0)), "should be a data frame")
  expect_error(update(jm1, weights = wts1), "do not have weights supplied")
  expect_error(update(jm1, weights = wts2), "should only have one row")
  expect_error(update(jm1, weights = wts3), "should be a data frame with two columns")
  expect_error(update(jm1, weights = wts4), "weights supplied must be numeric")
  expect_error(update(jm1, weights = wts5), "weights supplied must be numeric")
  expect_error(update(jm1, weights = wts6), "Negative weights are not allowed")
  
})

test_that("init argument works", {
  expect_output(update(jm1, init = "model_based"))
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
})

test_that("max_treedepth argument works", {
  expect_output(update(jm1, max_treedepth = NULL))
  expect_output(update(jm1, max_treedepth = 5))
  expect_output(update(jm1, max_treedepth = 5L))
})

test_that("error message occurs for arguments not implemented", {
  expect_error(update(jm1, offset = 1:10), "not yet implemented")
  expect_error(update(jm1, QR = TRUE), "not yet implemented")
  expect_error(update(jm1, sparse = TRUE), "not yet implemented")
  expect_error(update(jm1, dataAssoc = data.frame(a=1:10)), "not yet implemented")
})

#----  Compare parameter estimates: stan_jm(assoc = NULL) vs stan_glmer/coxph

if (interactive()) {
  compare_glmer <- function(fmLong, fam = gaussian, ...) {
    fmSurv <- Surv(futimeYears, death) ~ sex + trt
    y1 <- stan_glmer(fmLong, pbcLong, fam, iter = 1000, chains = CHAINS, seed = SEED)
    s1 <- coxph(fmSurv, data = pbcSurv)
    j1 <- stan_jm(fmLong, pbcLong, fmSurv, pbcSurv, time_var = "year", family = fam, 
                  assoc = NULL, iter = 1000, chains = CHAINS, seed = SEED, ...) 
    tols <- get_tols(y1, s1, tolscales = TOLSCALES)
    pars <- recover_pars(y1, s1)
    parsjm <- recover_pars(j1)
    for (i in names(tols$fixef))
      expect_equal(pars$fixef[[i]], parsjm$fixef[[i]], tol = tols$fixef[[i]])     
    for (i in names(tols$ranef))
      expect_equal(pars$ranef[[i]], parsjm$ranef[[i]], tol = tols$ranef[[i]])
    for (i in names(tols$event))
      expect_equal(pars$event[[i]], parsjm$event[[i]], tol = tols$event[[i]])
  }
  test_that("coefs same for stan_jm and stan_lmer/coxph", {
    compare_glmer(logBili ~ year + (1 | id), gaussian)})
  test_that("coefs same for stan_jm and stan_glmer, bernoulli", {
    compare_glmer(ybern ~ year + xbern + (1 | id), binomial)})
  test_that("coefs same for stan_jm and stan_glmer, poisson", {
    compare_glmer(ypois ~ year + xpois + (1 | id), poisson, init = 0)})
  test_that("coefs same for stan_jm and stan_glmer, negative binomial", {
    compare_glmer(ypois ~ year + xpois + (1 | id), neg_binomial_2)})
  test_that("coefs same for stan_jm and stan_glmer, Gamma", {
    compare_glmer(ygamm ~ year + xgamm + (1 | id), Gamma)})
  #test_that("coefs same for stan_jm and stan_glmer, inverse gaussian", {
  #  compare_glmer(ygamm ~ year + xgamm + (1 | id), inverse.gaussian)})  
}

#--------  Check (post-)estimation functions work with various model specifications

if (interactive()) {

  # No functions in formula
  f1 <- stan_jm(formulaLong = logBili ~ year + (1 | id), 
                dataLong = pbcLong,
                formulaEvent = Surv(futimeYears, death) ~ sex + trt, 
                dataEvent = pbcSurv,
                time_var = "year",
                # this next line is only to keep the example small in size!
                chains = 1, cores = 1, seed = 12345, iter = 10)
  
  # Functions on LHS of formula
  f2 <- update(f1, formulaLong. = exp(logBili) ~ year + (1 | id))
  
  # Functions on RHS of formula
  f3 <- update(f1, formulaLong. = logBili ~ poly(year, degree = 2) + (1 | id))
  
  # Functions on LHS and RHS of formula
  f4 <- update(f1, formulaLong. = exp(logBili) ~ poly(year, degree = 2) + (1 | id))
  
  # Intercept only event submodel
  f5 <- update(f1, formulaEvent. = Surv(futimeYears, death) ~ 1)
    
  # Binomial outcome on LHS of formula
  pbcLong$trials <- rpois(nrow(pbcLong), 6)
  pbcLong$succ <- rbinom(nrow(pbcLong), pbcLong$trials, .7)
  pbcLong$fail <- pbcLong$trials - pbcLong$succ
  f6 <- update(f1, formulaLong. = cbind(succ, fail) ~ poly(year, degree = 2) + (1 | id), 
               family = binomial, init = "prefit_vb", iter = 1000)

  # Different baseline hazards
  f7 <- update(f1, basehaz = "weibull")
  f8 <- update(f1, basehaz = "bs")
#  f9 <- update(f1, basehaz = "piecewise") # posterior_survfit not yet implemented for piecewise
  
  # Different association structures
  f10 <- update(f1, assoc = NULL)
  f11 <- update(f1, assoc = "etavalue")
  f12 <- update(f1, assoc = "etaslope")
  f13 <- update(f1, assoc = "etaauc")
  f14 <- update(f1, assoc = "muvalue")
  f15 <- update(f1, assoc = "muslope")
  f16 <- update(f1, assoc = "muauc")
  f17 <- update(f1, assoc = c("etavalue", "etaslope")) 
  f18 <- update(f1, assoc = c("muvalue", "muslope"))

  # Different association structures with intercept only submodel
  f19 <- update(f5, assoc = NULL)
  f20 <- update(f5, assoc = "etavalue")
  f21 <- update(f5, assoc = "etaslope")
  f22 <- update(f5, assoc = "etaauc")
  f23 <- update(f5, assoc = "muvalue")
  f24 <- update(f5, assoc = "muslope")
  f25 <- update(f5, assoc = "muauc")
  f26 <- update(f5, assoc = c("etavalue", "etaslope")) 
  f27 <- update(f5, assoc = c("muvalue", "muslope"))   
  
  # Shared random effect association structures
  f28 <- update(f1, assoc = c("shared_b"))
  f29 <- update(f1, assoc = c("shared_coef"))
  
  # Test the models
  for (j in 1:29) {
    tryCatch({
      mod <- get(paste0("f", j))
      cat("Checking model:", paste0("f", j), "\n")
      
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
    }, error = function(e)
      cat(" Failed for model", paste0("f", j), " due to error:\n", paste(e)))
  }
  
}



