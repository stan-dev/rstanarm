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
  tve_fixef = 0.5  # how many SEs can stan_surv tve HRs be from coxph/stpm2 tve HRs
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

run_sims <- FALSE # if TRUE then long running simulations are run

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
  es(up(testmod, basehaz = "exp-aft"))
  es(up(testmod, basehaz = "weibull-aft"))
  
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
  
  es(up(testmod, prior_aux = dirichlet()))
  es(up(testmod, prior_aux = normal(),      basehaz = "weibull"))
  es(up(testmod, prior_aux = student_t(),   basehaz = "weibull"))
  es(up(testmod, prior_aux = cauchy(),      basehaz = "weibull"))
  es(up(testmod, prior_aux = exponential(), basehaz = "weibull"))
  
  es(up(testmod, prior_smooth = exponential()))
  es(up(testmod, prior_smooth = normal()))
  es(up(testmod, prior_smooth = student_t()))
  es(up(testmod, prior_smooth = cauchy()))
  
  ee(up(testmod, prior_intercept = lasso()), "prior distribution")
  ee(up(testmod, prior_aux       = lasso()), "prior distribution")
  ee(up(testmod, prior_smooth    = lasso()), "prior distribution")
})

test_that("tve function works", {
  
  # single tve call
  es(up(testmod, formula. = 
          Surv(eventtime, status) ~ tve(x1) + x2))
  
  # multiple tve calls
  es(up(testmod, formula. = 
          Surv(eventtime, status) ~ tve(x1) + tve(x2)))
  
  # b-spline and piecewise tve in same model
  es(up(testmod, formula. = 
          Surv(eventtime, status) ~ tve(x1, type = "bs") + tve(x2, type = "pw")))

  # b-spline tve optional arguments
  es(up(testmod, formula. = 
          Surv(eventtime, status) ~ tve(x1, type = "bs", knots = c(1,2)) + x2))
  es(up(testmod, formula. = 
          Surv(eventtime, status) ~ tve(x1, type = "bs", df = 4) + x2))
  es(up(testmod, formula. = 
          Surv(eventtime, status) ~ tve(x1, type = "bs", degree = 2) + x2))
  ee(up(testmod, formula. = 
          Surv(eventtime, status) ~ tve(x1, type = "bs", junk = 2) + x2), 
          "Invalid argument to 'tve' function.")
  
  # piecewise tve optional arguments
  es(up(testmod, formula. = 
          Surv(eventtime, status) ~ tve(x1, type = "pw", knots = c(1,2)) + x2))
  es(up(testmod, formula. = 
          Surv(eventtime, status) ~ tve(x1, type = "pw", df = 4) + x2))
  ee(up(testmod, formula. = 
          Surv(eventtime, status) ~ tve(x1, type = "pw", degree = 2) + x2), 
          "Invalid argument to 'tve' function.")
})


#----  Compare parameter estimates: stan_surv vs coxph

compare_surv <- function(data, basehaz = "weibull", ...) {
  require(survival)
  fm    <- Surv(eventtime, status) ~ X1 + X2
  surv1 <- coxph(fm, data)
  stan1 <- stan_surv(formula = fm,
                     data    = data,
                     basehaz = basehaz,
                     iter    = ITER,
                     refresh = REFRESH,
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

#---- exponential data

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


#----  Compare parameter estimates: stan_surv vs survreg

compare_surv <- function(data, basehaz = "weibull-aft", ...) {
  require(survival)
  fm    <- Surv(eventtime, status) ~ X1 + X2
  dist  <- ifelse(basehaz == "weibull-aft", "weibull", "exponential")
  surv1 <- survreg(fm, data, dist = dist)
  stan1 <- stan_surv(formula = fm,
                     data    = data,
                     basehaz = basehaz,
                     iter    = ITER,
                     refresh = REFRESH,
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

#---- exponential data

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

compare_surv(data = dat, basehaz = "exp-aft")

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

compare_surv(data = dat, basehaz = "weibull-aft")


# COMMENTED OUT TO AVOID ADDING PACKAGES TO SUGGESTS
#  
# #----  Compare parameter estimates: stan_surv vs icenReg (interval censored)
#   
#   #---- simulated interval censored weibull data
#   
#   library(icenReg); set.seed(321)
#   
#   # simulate interval censored data
#   sim_data <- simIC_weib(n     = 5000, 
#                          b1    = 0.3, 
#                          b2    = -0.3, 
#                          model = 'ph', 
#                          shape = 2, 
#                          scale = 2, 
#                          inspections   = 6, 
#                          inspectLength = 1)
#   
#   # lower limit = 0 is actually left censoring (stan_surv doesn't accept 0's)
#   sim_data$l[sim_data$l == 0] <- -Inf 
#   
#   # fit stan model to interval censored data
#   fm <- Surv(l, u, type = 'interval2') ~ x1 + x2
#   ic_stan  <- stan_surv(fm, 
#                         data    = sim_data, 
#                         basehaz = "weibull",
#                         iter    = ITER,
#                         refresh = REFRESH,
#                         chains  = CHAINS,
#                         seed    = SEED)
#   
#   # compare stan estimates to known values from data generating model
#   truepars <- c('x1' = 0.3, 'x2' = -0.3, 'weibull-shape' = 2)
#   stanpars <- fixef(ic_stan)
#   expect_equal(stanpars[['x1']], 
#                truepars[['x1']], 
#                tol  = 0.01, 
#                info = "compare estimates (x1) with icenReg")
#   expect_equal(stanpars[['x2']], 
#                truepars[['x2']], 
#                tol  = 0.01, 
#                info = "compare estimates (x2) with icenReg")
#   expect_equal(stanpars[['weibull-shape']], 
#                truepars[['weibull-shape']], 
#                tol  = 0.1, 
#                info = "compare estimates (weibull-shape) with icenReg")  
# 
#   # fit model using icenReg package & compare log_lik with stan model
#   ic_icen  <- ic_par(fm, data = sim_data)
#   ll_icen  <- ic_icen$llk
#   ll_stan  <- mean(rowSums(log_lik(ic_stan)))
#   expect_equal(ll_icen, 
#                ll_stan, 
#                tol = 5, 
#                info = "compare log lik with icenReg")
# 
#     
# #----  Compare parameter estimates: stan_surv vs phreg (tvc & delayed entry)
#   
#   #---- mortality data: contains a time-varying covariate
#   
#   library(eha); library(dplyr); set.seed(987)
#   
#   # add a time-fixed covariate to the mortality data
#   data(mort); mort <- mort %>% group_by(id) %>% mutate(sesfixed = ses[[1]])
#   
#   # fit models using the time-fixed covariate & compare HR estimates
#   fm <- Surv(enter, exit, event) ~ sesfixed
#   f_weib <- phreg(fm, data = mort)
#   f_stan <- stan_surv(fm, 
#                       data    = mort, 
#                       basehaz = "weibull",
#                       iter    = ITER,
#                       refresh = REFRESH,
#                       chains  = CHAINS,
#                       seed    = SEED)
#   expect_equal(coef(f_weib)['sesfixedupper'],
#                coef(f_stan)['sesfixedupper'],
#                tol = 0.01)
#   
#   # fit models using the time-varying covariate & compare HR estimates
#   fm <- Surv(enter, exit, event) ~ ses
#   v_weib <- phreg(fm, data = mort)
#   v_stan <- stan_surv(fm, 
#                       data    = mort, 
#                       basehaz = "weibull",
#                       iter    = ITER,
#                       refresh = REFRESH,
#                       chains  = CHAINS,
#                       seed    = SEED)
#   expect_equal(coef(v_weib)['sesupper'],
#                coef(v_stan)['sesupper'],
#                tol = 0.01)
#   
#   # stupidity check; to make sure the hazard ratios actually differed 
#   # between the models with the time-fixed and time-varying covariate
#   expect_error(expect_equal(coef(f_weib)['sesfixedupper'][[1]],
#                             coef(v_weib)['sesupper'][[1]],
#                             tol = 0.1), "not equal") 

#---- Check tve models against coxph

#---- piecewise constant

set.seed(1919002)
covs <- data.frame(id = 1:1000,
                   X1 = rbinom(1000, 1, 0.3),
                   X2 = rnorm (1000, 2, 2.0))
dat <- simsurv(dist    = "exponential",
               lambdas = 0.1,
               betas   = c(X1 = 0.3, X2 = -0.3),
               x       = covs,
               tve     = c(X1 = -0.6),
               tvefun  = function(t) as.numeric(t > 10),
               maxt    = 30)
dat <- merge(dat, covs)

fmsurv <- Surv(eventtime, status) ~ X1 + tt(X1) + X2
o<-SW(surv1 <- coxph(fmsurv, dat, tt = function(x, t, ...) { x * as.numeric(t > 10) }))

fmstan <- Surv(eventtime, status) ~ tve(X1, type = "pw", knots = c(10)) + X2
o<-SW(stan1 <- stan_surv(fmstan, dat, chains = 1, refresh = 0L, iter = 1000, basehaz = "exp"))

tols <- get_tols(surv1, tolscales = TOLSCALES)
pars_surv <- recover_pars(surv1)
pars_stan <- recover_pars(stan1)
for (i in names(tols$fixef))
  expect_equal(pars_surv$fixef[[i]],
               pars_stan$fixef[[i]],
               tol = tols$fixef[[i]],
               info = "compare_estimates_tve_pw")


#--------  Check post-estimation functions work

pbcSurv$t0 <- 0
pbcSurv$t0[pbcSurv$futimeYears > 2] <- 1 # delayed entry

pbcSurv$t1 <- pbcSurv$futimeYears - 1 # lower limit for interval censoring
pbcSurv$t1[pbcSurv$t1 <= 0] <- -Inf   # left censoring

pbcSurv$clinic <- cut(pbcSurv$id, breaks = c(0,10,20,30,40), labels = FALSE)

# different baseline hazards
o<-SW(f1  <- stan_surv(Surv(futimeYears, death) ~ sex + trt,
                       data    = pbcSurv,
                       basehaz = "ms",
                       chains  = 1,
                       iter    = 100,
                       refresh = REFRESH,
                       seed    = SEED))
o<-SW(f2  <- update(f1, basehaz = "bs"))
o<-SW(f3  <- update(f1, basehaz = "exp"))
o<-SW(f4  <- update(f1, basehaz = "weibull"))
o<-SW(f5  <- update(f1, basehaz = "gompertz"))
o<-SW(f6  <- update(f1, basehaz = "exp-aft"))
o<-SW(f7  <- update(f1, basehaz = "weibull-aft"))

# time-varying effects
o<-SW(f8  <- update(f1, Surv(futimeYears, death) ~ sex + tve(trt)))
o<-SW(f9  <- update(f2, Surv(futimeYears, death) ~ sex + tve(trt)))
o<-SW(f10 <- update(f3, Surv(futimeYears, death) ~ sex + tve(trt)))
o<-SW(f11 <- update(f4, Surv(futimeYears, death) ~ sex + tve(trt)))
o<-SW(f12 <- update(f5, Surv(futimeYears, death) ~ sex + tve(trt)))
o<-SW(f13 <- update(f6, Surv(futimeYears, death) ~ sex + tve(trt)))
o<-SW(f14 <- update(f7, Surv(futimeYears, death) ~ sex + tve(trt)))

o<-SW(f15 <- update(f1, Surv(futimeYears, death) ~ sex + tve(trt, type = "pw")))
o<-SW(f16 <- update(f2, Surv(futimeYears, death) ~ sex + tve(trt, type = "pw")))
o<-SW(f17 <- update(f3, Surv(futimeYears, death) ~ sex + tve(trt, type = "pw")))
o<-SW(f18 <- update(f4, Surv(futimeYears, death) ~ sex + tve(trt, type = "pw")))
o<-SW(f19 <- update(f5, Surv(futimeYears, death) ~ sex + tve(trt, type = "pw")))
o<-SW(f20 <- update(f6, Surv(futimeYears, death) ~ sex + tve(trt, type = "pw")))
o<-SW(f21 <- update(f7, Surv(futimeYears, death) ~ sex + tve(trt, type = "pw")))

# start-stop notation (incl. delayed entry)
o<-SW(f22 <- update(f1, Surv(t0, futimeYears, death) ~ sex + trt))
o<-SW(f23 <- update(f1, Surv(t0, futimeYears, death) ~ sex + tve(trt)))
o<-SW(f24 <- update(f6, Surv(t0, futimeYears, death) ~ sex + tve(trt)))
o<-SW(f25 <- update(f6, Surv(t0, futimeYears, death) ~ sex + tve(trt)))

# left and interval censoring
o<-SW(f26 <- update(f1, Surv(t1, futimeYears, type = "interval2") ~ sex + trt))
o<-SW(f27 <- update(f1, Surv(t1, futimeYears, type = "interval2") ~ sex + tve(trt)))
o<-SW(f28 <- update(f6, Surv(t1, futimeYears, type = "interval2") ~ sex + trt))
o<-SW(f29 <- update(f6, Surv(t1, futimeYears, type = "interval2") ~ sex + tve(trt)))

# frailty models
o<-SW(f30 <- update(f1, Surv(futimeYears, death) ~ trt + (trt | clinic)))
o<-SW(f31 <- update(f1, Surv(futimeYears, death) ~ tve(trt) + (1 | clinic)))
o<-SW(f32 <- update(f1, Surv(t0, futimeYears, death) ~ trt + (trt | clinic)))
o<-SW(f33 <- update(f1, Surv(t1, futimeYears, type = "interval2") ~ trt + (trt | clinic)))

# new data for predictions
nd1 <- pbcSurv[pbcSurv$id == 2,]
nd2 <- pbcSurv[pbcSurv$id %in% c(1,2),]

# test the models
for (j in c(30:33)) {
  
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


#--------  Check hazard models with group-specific terms

#--- test estimates for each model type using one simulated dataset

# define a function to simulate a survival dataset
make_data <- function(n = 10,                  # number of patients per site
                      K = 30,                  # number of sites
                      dist    = "exponential", # basehaz for simulation
                      delay   = FALSE,         # induce delayed entry 
                      icens   = FALSE) {       # induce interval censoring 
  
  if (delay && icens)
    stop("'delay' and 'icens' cannot both be TRUE.")
  
  # dimensions
  N <- n * K # total num individuals
  
  # true sd for the random intercepts
  true_sd <- 1
  
  # sample random intercept for each site
  bb <- rnorm(K, 0, true_sd)
  
  # covariate data
  cov <- data.frame(id   = 1:N,
                    site = rep(1:K, each = n),
                    trt  = rbinom(N, 1, 0.5),
                    b    = bb[rep(1:K, each = n)])
  
  # simulate event times
  dat <- simsurv(dist    = dist,
                 lambdas = 0.1,
                 gammas  = switch(dist,
                                  "weibull"  = 1.3,
                                  "gompertz" = 0.05,
                                  NULL),
                 x       = cov,
                 betas   = c(trt = 0.3, b = 1),
                 maxt    = 15)
  
  # create delayed entry
  if (delay) {
    dat[["start"]] <- runif(N, 0, dat[["eventtime"]] / 2)
    dat[["stop"]]  <- dat[["eventtime"]]
  }
  
  # create interval censoring
  if (icens) {
    
    dd <- dat[["status"]] # event indicator
    dat[["lower"]] <- rep(NA, nrow(dat))
    dat[["upper"]] <- rep(NA, nrow(dat))
    
    # construct lower/upper interval cens times for right censored individuals
    dat[dd == 0, "lower"] <- dat[dd == 0, "eventtime"]
    dat[dd == 0, "upper"] <- Inf
    
    # construct lower/upper interval cens times for individuals with events
    dat[dd == 1, "lower"] <- runif(sum(dd == 1), 
                                   dat[dd == 1, "eventtime"] / 2, 
                                   dat[dd == 1, "eventtime"])
    dat[dd == 1, "upper"] <- runif(sum(dd == 1), 
                                   dat[dd == 1, "eventtime"],
                                   dat[dd == 1, "eventtime"] * 1.5)
    
  }
  
  merge(cov, dat)
  
}

# true parameter values used to simulate & corresponding tolerances for tests
true <- c(intercept = log(0.1), trt = 0.3, b_sd = 1)
tols <- c(0.2, 0.1, 0.2)

# function to return the parameter estimates to test
get_ests <- function(mod) {
  c(intercept = fixef(mod)[["(Intercept)"]],
    trt       = fixef(mod)[["trt"]],
    b_sd      = attr(VarCorr(mod)[[1]], "stddev")[[1]])
}

# fit right censored models
set.seed(5434)
dat <- make_data(n = 20, K = 50)
ff  <- Surv(eventtime, status) ~ trt + (1 | site)
m1  <- stan_surv(ff, data = dat, chains = 1, basehaz = "exp")
m2  <- stan_surv(ff, data = dat, chains = 1, basehaz = "weibull")
m3  <- stan_surv(ff, data = dat, chains = 1, basehaz = "gompertz")
m4  <- stan_surv(ff, data = dat, chains = 1, basehaz = "ms")
for (i in 1:3)
  expect_equal(get_ests(m1)[[i]], true[[i]], tol = tols[[i]])
for (i in 1:3)
  expect_equal(get_ests(m2)[[i]], true[[i]], tol = tols[[i]])
for (i in 1:3)
  expect_equal(get_ests(m3)[[i]], true[[i]], tol = tols[[i]])
for (i in 2:3)
  expect_equal(get_ests(m4)[[i]], true[[i]], tol = tols[[i]])

# fit delayed entry models
set.seed(8765)
dat_delay <- make_data(n = 20, K = 50, delay = TRUE)
ffd <- Surv(start, stop, status) ~ trt + (1 | site)
m5  <- stan_surv(ffd, data = dat_delay, chains = 1, refresh = 0, basehaz = "exp")
m6  <- stan_surv(ffd, data = dat_delay, chains = 1, refresh = 0, basehaz = "weibull")
m7  <- stan_surv(ffd, data = dat_delay, chains = 1, refresh = 0, basehaz = "gompertz")
m8  <- stan_surv(ffd, data = dat_delay, chains = 1, refresh = 0, basehaz = "ms")  
for (i in 1:3)
  expect_equal(get_ests(m5)[[i]], true[[i]], tol = tols[[i]])
for (i in 1:3)
  expect_equal(get_ests(m6)[[i]], true[[i]], tol = tols[[i]])
for (i in 1:3)
  expect_equal(get_ests(m7)[[i]], true[[i]], tol = tols[[i]])
for (i in 2:3)
  expect_equal(get_ests(m8)[[i]], true[[i]], tol = tols[[i]])

# fit interval censored models
set.seed(3254)
dat_icens <- make_data(n = 20, K = 50, icens = TRUE)
ffi <- Surv(lower, upper, type = "interval2") ~ trt + (1 | site)
m9  <- stan_surv(ffi, data = dat_icens, chains = 1, refresh = 0, iter = 1000, basehaz = "exp")
m10 <- stan_surv(ffi, data = dat_icens, chains = 1, refresh = 0, iter = 1000, basehaz = "weibull")
m11 <- stan_surv(ffi, data = dat_icens, chains = 1, refresh = 0, iter = 1000, basehaz = "gompertz")
m12 <- stan_surv(ffi, data = dat_icens, chains = 1, refresh = 0, iter = 1000, basehaz = "ms")  
for (i in 1:3)
  expect_equal(get_ests(m9)[[i]], true[[i]], tol = tols[[i]])
for (i in 1:3)
  expect_equal(get_ests(m10)[[i]], true[[i]], tol = tols[[i]])
for (i in 1:3)
  expect_equal(get_ests(m11)[[i]], true[[i]], tol = tols[[i]])
for (i in 2:3)
  expect_equal(get_ests(m12)[[i]], true[[i]], tol = tols[[i]])

#--- previous tests use really weak tolerances to check the 
#    parameter estimates; therefore the next part conducts a full 
#    simulation study to test each model specification and uses a
#    stronger tolerance, checking that relative bias is less than 5%

if (run_sims) {
  
  # number of simulations (for each model specification)
  n_sims <- 2 
  
  # define a function to fit the model to one simulated dataset
  sim_run <- function(n = 10,                  # number of patients per site
                      K = 30,                  # number of sites
                      basehaz = "exp",         # basehaz for analysis
                      dist    = "exponential", # basehaz for simulation
                      delay   = FALSE,         # induce delayed entry 
                      icens   = FALSE,         # induce interval censoring
                      return_relb = FALSE) {         
    
    # simulate data
    dat <- make_data(n = n, K = K, dist = dist, delay = delay, icens = icens)
    
    # define appropriate model formula
    if (delay) {
      ff <- Surv(start, stop, status) ~ trt + (1 | site)
    } else if (icens) {
      ff <- Surv(lower, upper, type = "interval2") ~ trt + (1 | site)
    } else {
      ff <- Surv(eventtime, status) ~ trt + (1 | site)
    }
    
    # fit model
    mod <- stan_surv(formula = ff,
                     data    = dat, 
                     basehaz = basehaz, 
                     chains  = 1,
                     refresh = 0,
                     iter    = 2000)
    
    # true parameters (hard coded here)
    true <- c(intercept = log(0.1),
              trt       = 0.3,
              b_sd      = 1)
    
    # extract parameter estimates
    ests <- c(intercept = fixef(mod)["(Intercept)"],
              trt       = fixef(mod)["trt"],
              b_sd      = attr(VarCorr(mod)[[1]], "stddev")[[1]])
    
    # intercept is irrelevant for spline model
    if (basehaz %in% c("ms", "bs")) {
      true <- true[2:3]
      ests <- ests[2:3]
    }
    
    if (return_relb)
      return(as.vector((ests - true) / true))
    
    list(true = true,
         ests = ests,
         bias = ests - true,
         relb = (ests - true) / true)
  }
  
  # functions to summarise the simulations and check relative bias
  summarise_sims <- function(x) {
    rbind(true = colMeans(do.call(rbind, x["true",])),
          ests = colMeans(do.call(rbind, x["ests",])),
          bias = colMeans(do.call(rbind, x["bias",])),
          relb = colMeans(do.call(rbind, x["relb",])))
  }
  validate_relbias <- function(x, tol = 0.05) {
    relb <- as.vector(summarise_sims(x)["relb",])
    expect_equal(relb, rep(0, length(relb)), tol = tol)
  }
  
  # right censored models
  set.seed(5050)
  sims_exp <- replicate(n_sims, sim_run(basehaz = "exp"))
  validate_relbias(sims_exp)
  
  set.seed(6060)
  sims_weibull <- replicate(n_sims, sim_run(basehaz = "weibull"))
  validate_relbias(sims_weibull)
  
  set.seed(7070)
  sims_gompertz <- replicate(n_sims, sim_run(basehaz = "gompertz"))
  validate_relbias(sims_gompertz)
  
  set.seed(8080)
  sims_ms <- replicate(n_sims, sim_run(basehaz = "ms"))
  validate_relbias(sims_ms)
  
  # delayed entry models
  set.seed(5050)
  sims_exp_d <- replicate(n_sims, sim_run(basehaz = "exp", delay = TRUE))
  validate_relbias(sims_exp_d)
  
  set.seed(6060)
  sims_weibull_d <- replicate(n_sims, sim_run(basehaz = "weibull", delay = TRUE))
  validate_relbias(sims_weibull_d)
  
  set.seed(7070)
  sims_gompertz_d <- replicate(n_sims, sim_run(basehaz = "gompertz", delay = TRUE))
  validate_relbias(sims_gompertz_d)
  
  set.seed(8080)
  sims_ms_d <- replicate(n_sims, sim_run(basehaz = "ms", delay = TRUE))
  validate_relbias(sims_ms_d)
  
  # interval censored models
  set.seed(5050)
  sims_exp_i <- replicate(n_sims, sim_run(basehaz = "exp", icens = TRUE))
  validate_relbias(sims_exp_i)
  
  set.seed(6060)
  sims_weibull_i <- replicate(n_sims, sim_run(basehaz = "weibull", icens = TRUE))
  validate_relbias(sims_weibull_i)
  
  set.seed(7070)
  sims_gompertz_i <- replicate(n_sims, sim_run(basehaz = "gompertz", icens = TRUE))
  validate_relbias(sims_gompertz_i)
  
  set.seed(8080)
  sims_ms_i <- replicate(n_sims, sim_run(basehaz = "ms", icens = TRUE))
  validate_relbias(sims_ms_i)
  
}
