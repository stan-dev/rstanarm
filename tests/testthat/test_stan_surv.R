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

library(testthat)
library(rstanarm)
library(survival)
library(simsurv)
ITER    <- 1000
CHAINS  <- 1
SEED    <- 12345
REFRESH <- 0L
set.seed(SEED)
if (interactive())
  options(mc.cores = parallel::detectCores())

TOLSCALES <- list(
  hr_fixef  = 0.5, # how many SEs can stan_surv HRs be from coxph/stpm2 HRs
  tde_fixef = 0.5  # how many SEs can stan_surv tde HRs be from coxph/stpm2 tde HRs
)

source(test_path("helpers", "expect_matrix.R"))
source(test_path("helpers", "expect_stanreg.R"))
source(test_path("helpers", "expect_stanmvreg.R"))
source(test_path("helpers", "expect_survfit_surv.R"))
source(test_path("helpers", "expect_ppd.R"))
source(test_path("helpers", "expect_equivalent_loo.R"))
source(test_path("helpers", "SW.R"))
source(test_path("helpers", "get_tols_surv.R"))
source(test_path("helpers", "recover_pars_surv.R"))

eo <- function(...) { expect_output (...) }
ee <- function(...) { expect_error  (...) }
ew <- function(...) { expect_warning(...) }
es <- function(...) { expect_stanreg(...) }
up <- function(...) { update(...) }

#-----------------------------  Models -----------------------------------

#--- Time fixed covariates, time fixed coefficients

cov1 <- data.frame(id = 1:50,
                   x1 = stats::rbinom(50, 1, 0.5),
                   x2 = stats::rnorm (50, -1, 0.5))
dat1 <- simsurv(lambdas = 0.1,
                gammas  = 1.5,
                betas   = c(x1 = -0.5, x2 = -0.3),
                x       = cov1,
                maxt    = 5)
dat1 <- merge(dat1, cov1)
fm1  <- Surv(eventtime, status) ~ x1 + x2
o<-SW(testmod <- stan_surv(fm1, dat1, chains = 1, refresh = 0L, iter = 50, basehaz = "ms"))

# mod1a <- stan_surv(fm1, dat1, chains = 1, refresh = 0L, iter = 1000, basehaz = "ms")
# mod1b <- stan_surv(fm1, dat1, chains = 1, refresh = 0L, iter = 1000, basehaz = "bs")
# mod1c <- stan_surv(fm1, dat1, chains = 1, refresh = 0L, iter = 1000, basehaz = "exp")
# mod1d <- stan_surv(fm1, dat1, chains = 1, refresh = 0L, iter = 1000, basehaz = "weibull")
# mod1e <- stan_surv(fm1, dat1, chains = 1, refresh = 0L, iter = 1000, basehaz = "gompertz")


#--------------------------  Arguments -----------------------------------

test_that("prior_PD argument works", {
  es(up(testmod, prior_PD = TRUE))
})

test_that("adapt_delta argument works", {
  es(up(testmod, adapt_delta = NULL))
  es(up(testmod, adapt_delta = 0.8))
  es(up(testmod, control = list(adapt_delta = NULL)))
  es(up(testmod, control = list(adapt_delta = 0.8)))
})

test_that("init argument works", {
  es(up(testmod, init = "prefit"))
  es(up(testmod, init = "0"))
  es(up(testmod, init = 0))
  es(up(testmod, init = "random"))
})

test_that("qnodes argument works", {
  es(up(testmod, qnodes = 7,  basehaz = "bs"))
  es(up(testmod, qnodes = 11, basehaz = "bs"))
  es(up(testmod, qnodes = 15, basehaz = "bs"))
  
  ew(up(testmod, qnodes = 1),       "is being ignored")
  ew(up(testmod, qnodes = c(1,2)),  "is being ignored")
  ew(up(testmod, qnodes = "wrong"), "is being ignored")
  
  ee(up(testmod, qnodes = 1,       basehaz = "bs"), "7, 11 or 15")
  ee(up(testmod, qnodes = c(1,2),  basehaz = "bs"), "numeric vector of length 1")
  ee(up(testmod, qnodes = "wrong", basehaz = "bs"), "numeric vector of length 1")
})

test_that("basehaz argument works", {

  es(up(testmod, basehaz = "exp"))
  es(up(testmod, basehaz = "weibull"))
  es(up(testmod, basehaz = "gompertz"))
  es(up(testmod, basehaz = "ms"))
  es(up(testmod, basehaz = "bs"))

  dfl <- list(df = 5)
  knl <- list(knots = c(1,3,5))
  es(up(testmod, basehaz = "ms", basehaz_ops = dfl))
  es(up(testmod, basehaz = "ms", basehaz_ops = knl))
  es(up(testmod, basehaz = "bs", basehaz_ops = dfl))
  es(up(testmod, basehaz = "bs", basehaz_ops = knl))

  ee(up(testmod, basehaz_ops = list(junk = 3)), "can only include")

  ee(up(testmod, basehaz_ops = list(df = 1)),            "cannot be negative")
  ee(up(testmod, basehaz_ops = list(knots = -1)),        "earliest entry time")
  ee(up(testmod, basehaz_ops = list(knots = c(1,2,50))), "latest event time")

})

test_that("prior arguments work", {
  es(up(testmod, prior = normal()))
  es(up(testmod, prior = student_t()))
  es(up(testmod, prior = cauchy()))
  es(up(testmod, prior = hs()))
  es(up(testmod, prior = hs_plus()))
  es(up(testmod, prior = lasso()))
  es(up(testmod, prior = laplace()))
 
  es(up(testmod, prior_intercept = normal()))
  es(up(testmod, prior_intercept = student_t()))
  es(up(testmod, prior_intercept = cauchy()))
  
  es(up(testmod, prior_aux = normal()))
  es(up(testmod, prior_aux = student_t()))
  es(up(testmod, prior_aux = cauchy()))
  
  es(up(testmod, prior_smooth = exponential()))
  es(up(testmod, prior_smooth = normal()))
  es(up(testmod, prior_smooth = student_t()))
  es(up(testmod, prior_smooth = cauchy()))
  
  ee(up(testmod, prior_intercept = lasso()), "prior distribution")
  ee(up(testmod, prior_aux       = lasso()), "prior distribution")
  ee(up(testmod, prior_smooth    = lasso()), "prior distribution")
})


#----  Compare parameter estimates: stan_surv vs coxph

  compare_surv <- function(data, basehaz = "weibull", ...) {
    require(survival)
    fm    <- Surv(eventtime, status) ~ X1 + X2
    surv1 <- coxph(fm, data)
    stan1 <- stan_surv(formula = fm,
                       data    = data,
                       basehaz = basehaz,
                       iter    = 1000,
                       refresh = 0L,
                       chains  = CHAINS,
                       seed    = SEED, ...)
    tols <- get_tols(surv1, tolscales = TOLSCALES)
    pars_surv <- recover_pars(surv1)
    pars_stan <- recover_pars(stan1)
    for (i in names(tols$fixef))
      expect_equal(pars_surv$fixef[[i]],
                   pars_stan$fixef[[i]],
                   tol = tols$fixef[[i]],
                   info = basehaz)
  }
  
  #---- weibull data
  
  set.seed(543634)
  covs <- data.frame(id = 1:300,
                     X1 = rbinom(300, 1, 0.3),
                     X2 = rnorm (300, 2, 2.0))
  dat <- simsurv(dist    = "weibull",
                 lambdas = 0.1,
                 gammas  = 1,
                 betas   = c(X1 = 0.3, X2 = -0.5),
                 x       = covs)
  dat <- merge(dat, covs)
  
  compare_surv(data = dat, basehaz = "exp")
  
  #---- weibull data

  set.seed(543634)
  covs <- data.frame(id = 1:300,
                     X1 = rbinom(300, 1, 0.3),
                     X2 = rnorm (300, 2, 2.0))
  dat <- simsurv(dist    = "weibull",
                 lambdas = 0.1,
                 gammas  = 1.3,
                 betas   = c(X1 = 0.3, X2 = -0.5),
                 x       = covs)
  dat <- merge(dat, covs)

  compare_surv(data = dat, basehaz = "weibull")
  compare_surv(data = dat, basehaz = "ms")
  compare_surv(data = dat, basehaz = "bs")
  
  #---- gompertz data

  set.seed(45357)
  covs <- data.frame(id = 1:300,
                     X1 = rbinom(300, 1, 0.3),
                     X2 = rnorm (300, 2, 2.0))
  dat <- simsurv(dist    = "gompertz",
                 lambdas = 0.1,
                 gammas  = 0.05,
                 betas   = c(X1 = -0.6, X2 = -0.4),
                 x       = covs)
  dat <- merge(dat, covs)

  compare_surv(data = dat, basehaz = "gompertz")


#----  Compare parameter estimates: stan_surv vs icenReg
  
  #---- interval censored weibull data
  
  library(icenReg)
  set.seed(321)
  sim_data <- simIC_weib(n     = 5000, 
                         b1    = 0.3, 
                         b2    = -0.3, 
                         model = 'ph', 
                         shape = 2, 
                         scale = 2, 
                         inspections   = 6, 
                         inspectLength = 1)
  sim_data$l[sim_data$l == 0] <- -Inf # left limit = 0 is actually left censoring
  fm <- Surv(l, u, type = 'interval2') ~ x1 + x2
  ic_icen  <- ic_par(fm, data = sim_data)
  ic_stan  <- stan_surv(fm, data = sim_data, basehaz = "weibull")
  truepars <- c('x1' = 0.3, 'x2' = -0.3, 'weibull-shape' = 2)
  stanpars <- fixef(ic_stan)
  ll_icen  <- ic_icen$llk
  ll_stan  <- mean(rowSums(log_lik(ic_stan)))
  expect_equal(stanpars[['x1']], 
               truepars[['x1']], 
               tol  = 0.01, 
               info = "compare estimates (x1) with icenReg")
  expect_equal(stanpars[['x2']], 
               truepars[['x2']], 
               tol  = 0.01, 
               info = "compare estimates (x2) with icenReg")
  expect_equal(stanpars[['weibull-shape']], 
               truepars[['weibull-shape']], 
               tol  = 0.1, 
               info = "compare estimates (weibull-shape) with icenReg")  
  expect_equal(ll_icen, 
               ll_stan, 
               tol = 5, 
               info = "compare log lik with icenReg")
  
  
#--------  Check post-estimation functions work

  pbcSurv$t0 <- 0
  pbcSurv$t0[pbcSurv$futimeYears > 2] <- 1 # delayed entry

  # different baseline hazards
  o<-SW(f1  <- stan_surv(Surv(futimeYears, death) ~ sex + trt,
                         data    = pbcSurv,
                         basehaz = "ms",
                         chains  = 1,
                         cores   = 1,
                         iter    = 40,
                         refresh = 0,
                         seed    = 12345))
  o<-SW(f2  <- update(f1, basehaz = "bs"))
  o<-SW(f3  <- update(f1, basehaz = "exp"))
  o<-SW(f4  <- update(f1, basehaz = "weibull"))
  o<-SW(f5  <- update(f1, basehaz = "gompertz"))
  
  # time-dependent effects
  o<-SW(f6  <- update(f1, Surv(futimeYears, death) ~ sex + tde(trt)))
  o<-SW(f7  <- update(f2, Surv(futimeYears, death) ~ sex + tde(trt)))
  o<-SW(f8  <- update(f3, Surv(futimeYears, death) ~ sex + tde(trt)))
  o<-SW(f9  <- update(f4, Surv(futimeYears, death) ~ sex + tde(trt)))
  o<-SW(f10 <- update(f5, Surv(futimeYears, death) ~ sex + tde(trt)))
  
  # start-stop notation (incl. delayed entry)
  o<-SW(f11 <- update(f1, Surv(t0, futimeYears, death) ~ sex + trt))
  o<-SW(f12 <- update(f1, Surv(t0, futimeYears, death) ~ sex + tde(trt)))
  
  # interval censoring
  o<-SW(f13 <- update(f1, Surv(t0, futimeYears, type = "interval2") ~ sex + trt))
  #o<-SW(f14 <- update(f1, Surv(t0, futimeYears, type = "interval2") ~ sex + tde(trt)))
  
  # new data for predictions
  nd1 <- pbcSurv[pbcSurv$id == 2,]
  nd2 <- pbcSurv[pbcSurv$id %in% c(1,2),]

  # test the models
  for (j in c(1:13)) {

    mod <- try(get(paste0("f", j)), silent = TRUE)

    if (class(mod)[1L] == "try-error") {

      cat("Model not found:", paste0("f", j), "\n")

    } else {

      cat("Checking model:", paste0("f", j), "\n")

      test_that("log_lik works with estimation data", {
        ll <- log_lik(mod)
        expect_matrix(ll)
      })

      test_that("log_lik works with new data (one individual)", {
        ll <- log_lik(mod, newdata = nd1)
        expect_matrix(ll)
      })

      test_that("log_lik works with new data (multiple individuals)", {
        ll <- log_lik(mod, newdata = nd2)
        expect_matrix(ll)
      })

      test_that("loo and waic work", {
        expect_equivalent_loo(mod)
      })

      if (mod$ndelayed == 0) # only test if no delayed entry
        test_that("posterior_survfit works with estimation data", {
          SW(ps <- posterior_survfit(mod))
          expect_survfit(ps)
        })

      test_that("posterior_survfit works with new data (one individual)", {
        SW(ps <- posterior_survfit(mod, newdata = nd1))
        expect_survfit(ps)
      })

      test_that("posterior_survfit works with new data (multiple individuals)", {
        SW(ps <- posterior_survfit(mod, newdata = nd2))
        expect_survfit(ps)
      })

    }
  }
