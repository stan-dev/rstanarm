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
REFRESH <- 0L
SEED    <- 12345; set.seed(SEED)
if (interactive())
  options(mc.cores  = parallel::detectCores(),
          loo.cores = parallel::detectCores())

TOLSCALES <- list(
  hr_fixef  = 0.5 # how many SEs can stan_surv HRs be from coxph/stpm2 HRs
)

context("stan_surv")

eo <- function(...) { expect_output (...) }
ee <- function(...) { expect_error  (...) }
ew <- function(...) { expect_warning(...) }
es <- function(...) { expect_stanreg(...) }
up <- function(...) { update(...) }

run_sims <- FALSE # if TRUE then long running simulations are run


#-----------------  Check model fitting arguments work  -----------------------

cov1 <- data.frame(id = 1:50,
                   x1 = stats::rbinom(50, 1, 0.5),
                   x2 = stats::rnorm (50, -1, 0.5))

dat1 <- simsurv(lambdas = 0.1,
                gammas  = 1.5,
                betas   = c(x1 = -0.5, x2 = -0.3),
                x       = cov1,
                maxt    = 5)

dat1$s <- Surv(dat1$eventtime, dat1$status) # abbreviated Surv object

o<-SW(testmod <- stan_surv(formula = s ~ x1 + x2,
                           data    = merge(dat1, cov1),
                           basehaz = "ms",
                           iter    = 20,
                           chains  = CHAINS,
                           refresh = REFRESH,
                           seed    = SEED))

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
  es(up(testmod, formula. = s ~ tve(x1) + x2))
  es(up(testmod, formula. = s ~ tve(x1) + tve(x2)))
  es(up(testmod, formula. = s ~ tve(x1, knots = 1) + tve(x2, knots = 2)))
})  

test_that("tve function works: b-spline optional arguments", {
  es(up(testmod, formula. = s ~ tve(x1, knots = c(1,2)) + x2))
  es(up(testmod, formula. = s ~ tve(x1, df = 4)         + x2))
  es(up(testmod, formula. = s ~ tve(x1, degree = 0)     + x2))
  ee(up(testmod, formula. = s ~ tve(x1, junk = 2)       + x2), "unused")
})


#----------------  Check post-estimation functions work  ----------------------

# use PBC data
pbcSurv$t0 <- 0
pbcSurv$t0[pbcSurv$futimeYears > 2] <- 1 # fake delayed entry
pbcSurv$t1 <- pbcSurv$futimeYears - 1    # fake lower limit for interval censoring
pbcSurv$t1[pbcSurv$t1 <= 0] <- -Inf      # fake left censoring
pbcSurv$site <- cut(pbcSurv$id,          # fake group for frailty models
                    breaks = c(0,10,20,30,40), 
                    labels = FALSE)

# different baseline hazards
o<-SW(f1  <- stan_surv(Surv(futimeYears, death) ~ sex + trt,
                       data    = pbcSurv,
                       basehaz = "ms",
                       chains  = 1,
                       iter    = 20,
                       refresh = REFRESH,
                       seed    = SEED))
o<-SW(f2  <- up(f1, basehaz = "bs"))
o<-SW(f3  <- up(f1, basehaz = "exp"))
o<-SW(f4  <- up(f1, basehaz = "weibull"))
o<-SW(f5  <- up(f1, basehaz = "gompertz"))
o<-SW(f6  <- up(f1, basehaz = "exp-aft"))
o<-SW(f7  <- up(f1, basehaz = "weibull-aft"))

# time-varying effects
o<-SW(f8  <- up(f1, Surv(futimeYears, death) ~ sex + tve(trt)))
o<-SW(f9  <- up(f2, Surv(futimeYears, death) ~ sex + tve(trt)))
o<-SW(f10 <- up(f3, Surv(futimeYears, death) ~ sex + tve(trt)))
o<-SW(f11 <- up(f4, Surv(futimeYears, death) ~ sex + tve(trt)))
o<-SW(f12 <- up(f5, Surv(futimeYears, death) ~ sex + tve(trt)))
o<-SW(f13 <- up(f6, Surv(futimeYears, death) ~ sex + tve(trt)))
o<-SW(f14 <- up(f7, Surv(futimeYears, death) ~ sex + tve(trt)))

o<-SW(f15 <- up(f1, Surv(futimeYears, death) ~ sex + tve(trt, degree = 0)))
o<-SW(f16 <- up(f2, Surv(futimeYears, death) ~ sex + tve(trt, degree = 0)))
o<-SW(f17 <- up(f3, Surv(futimeYears, death) ~ sex + tve(trt, degree = 0)))
o<-SW(f18 <- up(f4, Surv(futimeYears, death) ~ sex + tve(trt, degree = 0)))
o<-SW(f19 <- up(f5, Surv(futimeYears, death) ~ sex + tve(trt, degree = 0)))
o<-SW(f20 <- up(f6, Surv(futimeYears, death) ~ sex + tve(trt, degree = 0)))
o<-SW(f21 <- up(f7, Surv(futimeYears, death) ~ sex + tve(trt, degree = 0)))

# start-stop notation (incl. delayed entry)
o<-SW(f22 <- up(f1, Surv(t0, futimeYears, death) ~ sex + trt))
o<-SW(f23 <- up(f1, Surv(t0, futimeYears, death) ~ sex + tve(trt)))
o<-SW(f24 <- up(f6, Surv(t0, futimeYears, death) ~ sex + tve(trt)))
o<-SW(f25 <- up(f6, Surv(t0, futimeYears, death) ~ sex + tve(trt)))

# left and interval censoring
o<-SW(f26 <- up(f1, Surv(t1, futimeYears, type = "interval2") ~ sex + trt))
o<-SW(f27 <- up(f1, Surv(t1, futimeYears, type = "interval2") ~ sex + tve(trt)))
o<-SW(f28 <- up(f6, Surv(t1, futimeYears, type = "interval2") ~ sex + trt))
o<-SW(f29 <- up(f6, Surv(t1, futimeYears, type = "interval2") ~ sex + tve(trt)))

# frailty models
o<-SW(f30 <- up(f1, Surv(futimeYears, death) ~ trt + (trt | site)))
o<-SW(f31 <- up(f1, Surv(futimeYears, death) ~ tve(trt) + (1 | site)))
o<-SW(f32 <- up(f1, Surv(t0, futimeYears, death) ~ trt + (trt | site)))
o<-SW(f33 <- up(f1, Surv(t1, futimeYears, type = "interval2") ~ trt + (trt | site)))

# new data for predictions
nd1 <- pbcSurv[pbcSurv$id == 2,]
nd2 <- pbcSurv[pbcSurv$id %in% c(1,2),]

# test the models
for (j in c(1:33)) {
  
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
    
    if (mod$ndelayed == 0) # only test if no delayed entry
      test_that("posterior_survfit works with estimation data", {
        SW(ps <- posterior_survfit(mod))
        expect_survfit_surv(ps)
      })
    
    test_that("posterior_survfit works with new data (one individual)", {
      SW(ps <- posterior_survfit(mod, newdata = nd1))
      expect_survfit_surv(ps)
    })
    
    test_that("posterior_survfit works with new data (multiple individuals)", {
      SW(ps <- posterior_survfit(mod, newdata = nd2))
      expect_survfit_surv(ps)
    })
    
  }
}

# test loo for a few models only (too slow to test them all)
for (j in c(1,2,8,26,30)) {
  
  mod <- try(get(paste0("f", j)), silent = TRUE)
  
  if (class(mod)[1L] == "try-error") {
    
    cat("Model not found:", paste0("f", j), "\n")
    
  } else {
    
    cat("Checking loo for model:", paste0("f", j), "\n")
    
    test_that("loo and waic work", {
      loo_try <- try(expect_equivalent_loo(mod), silent = TRUE)
      if (class(loo_try)[1L] == "try-error") {
        # sometimes loo fails with a small number of draws so refit with more
        expect_equivalent_loo(up(mod, iter = 80))
      }
    })
    
  }
}


#----  Check accuracy of log_lik and posterior_survfit: for frailty models  ---

fake_data <-
  data.frame(id        = c(1,2,3,4),
             trt       = c(1,0,1,0),
             age       = c(5,8,2,4),
             site      = c(1,1,2,2),
             eventtime = c(2,4,6,8),
             status    = c(1,1,1,1))

o<-SW(stan1 <- stan_surv(formula = Surv(eventtime, status) ~ trt + age + (1 | site), 
                         data    = fake_data, 
                         basehaz = "weibull",
                         chains  = 1, 
                         refresh = 0L, 
                         iter    = 100,
                         warmup  = 95))

stanmat <- as.matrix(stan1)

stanpars <- list(int   = stanmat[, "(Intercept)"],
                 trt   = stanmat[, "trt"],
                 age   = stanmat[, "age"],
                 shape = stanmat[, "weibull-shape"],
                 site  = list(stanmat[, "b[(Intercept) site:1]"],
                              stanmat[, "b[(Intercept) site:2]"]))

N  <- nrow(fake_data)
S  <- nrow(stanmat)

# define function to calculate log likelihood manually
llfun <- function(i, j, data, pars) { 
  exp_eta_ij <- exp(pars$int[j] + 
                      pars$trt[j] * data$trt[i] + 
                      pars$age[j] * data$age[i] +
                      pars$site[[data$site[i]]][[j]])
  h_ij <- pars$shape[j] * data$eventtime[i] ^ (pars$shape[j] - 1) * exp_eta_ij
  H_ij <- data$eventtime[i] ^ (pars$shape[j]) * exp_eta_ij
  return(data$status[i] * log(h_ij) - H_ij)
}

# define function to calculate survival probability manually
survfun <- function(i, j, t, data, pars) {
  exp_eta_ij <- exp(pars$int[j] + 
                      pars$trt[j] * data$trt[i] + 
                      pars$age[j] * data$age[i] +
                      pars$site[[data$site[i]]][[j]])
  H_ij <- t ^ (pars$shape[j]) * exp_eta_ij
  return(exp(- H_ij))
}

# check log likelihood
L1 <- log_lik(stan1)
L2 <- log_lik(stan1, newdata = fake_data)
L3 <- matrix(NA, S, N) # manually evaluated log likelihood
for (i in 1:N) {
  for (j in 1:S) {
    L3[j,i] <- llfun(i, j, data = fake_data, pars = stanpars)
  }
}
for (i in 1:N) {
  for (j in 1:S) {
    expect_equal(as.vector(L3[j,i]), as.vector(L1[j,i]))
    expect_equal(as.vector(L3[j,i]), as.vector(L2[j,i]))
  }
}

# check survival probability
P1 <- posterior_survfit(stan1, times = 5, extrapolate = FALSE)
P2 <- posterior_survfit(stan1, newdata = fake_data, times = 5, extrapolate = FALSE)
P3 <- matrix(NA, S, N) # manually evaluated survival probability
for (i in 1:N) {
  for (j in 1:S) {
    P3[j,i] <- survfun(i, j, t = 5, data = fake_data, pars = stanpars)
  }
}
for (i in 1:N) {
  expect_equal(median(P3[,i]), P1[i, "median"])
  expect_equal(median(P3[,i]), P2[i, "median"])
}


#----------------  Check parameter estimates: stan vs coxph  -----------------

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
  tols <- get_tols_surv(surv1, tolscales = TOLSCALES)
  pars_surv <- recover_pars_surv(surv1)
  pars_stan <- recover_pars_surv(stan1)
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


#-----------  Check parameter estimates: stan (AFT) vs survreg  ---------------

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
                     seed    = SEED, 
                     ...)
  tols <- get_tols_surv(surv1, tolscales = TOLSCALES)
  pars_surv <- recover_pars_surv(surv1)
  pars_stan <- recover_pars_surv(stan1)
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


#--------  Check parameter estimates: stan (tve) vs coxph (tt)  ---------------

# NB: this only checks piecewise constant hazard ratio

set.seed(SEED)

N <- 1000 # number of individuals to simulate

covs <- data.frame(id = 1:N,
                   X1 = rbinom(N, 1, 0.3),
                   X2 = rnorm (N, 2, 2.0))

dat <- simsurv(dist    = "exponential",
               x       = covs,
               lambdas = c(0.1),
               betas   = c(X1 = 0.3, X2 = -0.3),
               tve     = c(X1 = -0.6),
               tvefun  = function(t) as.numeric(t > 10),
               maxt    = 30)

o<-SW(surv1 <- coxph(
  formula = Surv(eventtime, status) ~ X1 + tt(X1) + X2, 
  data    = merge(dat, covs), 
  tt      = function(x, t, ...) { x * as.numeric(t > 10) }))

o<-SW(stan1 <- stan_surv(
  formula = Surv(eventtime, status) ~ tve(X1, degree = 0, knots = c(10)) + X2, 
  data    = merge(dat, covs),
  basehaz = "exp",
  chains  = CHAINS, 
  refresh = REFRESH, 
  iter    = ITER))

tols <- get_tols_surv(surv1, tolscales = TOLSCALES)

pars_surv <- recover_pars_surv(surv1)
pars_stan <- recover_pars_surv(stan1)

for (i in names(tols$fixef))
  expect_equal(pars_surv$fixef[[i]],
               pars_stan$fixef[[i]],
               tol = tols$fixef[[i]],
               info = "compare_estimates_tve_pw")


# COMMENTED OUT TO AVOID ADDING PACKAGES TO SUGGESTS
#  
# #-------  Compare parameter estimates: stan (icens) vs icenReg  -------------
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
# #----  Compare parameter estimates: stan (tvc & delayed entry) vs phreg  ----
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


#--------  Check parameter estimates: stan (frailty) vs simulated  -----------

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

# simulate datasets
set.seed(SEED)
dat       <- make_data(n = 20, K = 50)
dat_delay <- make_data(n = 20, K = 50, delay = TRUE)
dat_icens <- make_data(n = 20, K = 50, icens = TRUE)

# formulas
ff  <- Surv(eventtime, status)                ~ trt + (1 | site) # right cens
ffd <- Surv(start, stop, status)              ~ trt + (1 | site) # delayed entry
ffi <- Surv(lower, upper, type = "interval2") ~ trt + (1 | site) # interval cens

# fit the starting model
o<-SW(m1  <- stan_surv(formula = ff, 
                       data    = dat, 
                       basehaz = "exp", 
                       iter    = ITER,
                       refresh = REFRESH,
                       chains  = CHAINS,
                       seed    = SEED))

# fit the additional models
o<-SW(m2  <- up(m1, formula. = ff, data = dat, basehaz = "weibull"))
o<-SW(m3  <- up(m1, formula. = ff,  data = dat,       basehaz = "gompertz"))
o<-SW(m4  <- up(m1, formula. = ff,  data = dat,       basehaz = "ms"))
o<-SW(m5  <- up(m1, formula. = ffd, data = dat_delay, basehaz = "exp"))
o<-SW(m6  <- up(m1, formula. = ffd, data = dat_delay, basehaz = "weibull"))
o<-SW(m7  <- up(m1, formula. = ffd, data = dat_delay, basehaz = "gompertz"))
o<-SW(m8  <- up(m1, formula. = ffd, data = dat_delay, basehaz = "ms"))
o<-SW(m9  <- up(m1, formula. = ffi, data = dat_icens, basehaz = "exp"))
o<-SW(m10 <- up(m1, formula. = ffi, data = dat_icens, basehaz = "weibull"))
o<-SW(m11 <- up(m1, formula. = ffi, data = dat_icens, basehaz = "gompertz"))
o<-SW(m12 <- up(m1, formula. = ffi, data = dat_icens, basehaz = "ms"))

# check the estimates against the true parameters
for (j in c(1:12)) {
  modfrail <- get(paste0("m", j))
  for (i in 1:3)
    expect_equal(get_ests(modfrail)[[i]], true[[i]], tol = tols[[i]])
}


#--- previous tests use really weak tolerances to check the 
#    parameter estimates; therefore the next part conducts a full 
#    simulation study to test each model specification and uses a
#    stronger tolerance, checking that relative bias is less than 5%

if (run_sims) {
  
  # number of simulations (for each model specification)
  n_sims <- 200
  
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
    
    
    # check Rhat
    rhats <- summary(mod)[, "Rhat"]
    rhats <- rhats[!names(rhats) %in% c("lp__", "log-posterior")]
    
    converged <- (all(rhats <= 1.1, na.rm = TRUE))

    if (!converged)
      ests <- rep(NA, length(ests)) # set estimates to NA if model didn't converge
    
    if (return_relb)
      return(as.vector((ests - true) / true))
    
    list(true = true,
         ests = ests,
         bias = ests - true,
         relb = (ests - true) / true)
  }
  
  # functions to summarise the simulations and check relative bias
  summarise_sims <- function(x) {
    message("Number of simulations that converged: ", 
            sum(!is.na(do.call(rbind, x["ests",])[,1])))
    rbind(true = colMeans(do.call(rbind, x["true",]), na.rm = TRUE),
          ests = colMeans(do.call(rbind, x["ests",]), na.rm = TRUE),
          bias = colMeans(do.call(rbind, x["bias",]), na.rm = TRUE),
          relb = colMeans(do.call(rbind, x["relb",]), na.rm = TRUE))
  }
  
  validate_relbias <- function(x, tol = 0.05) {
    message("Number of simulations that converged: ", 
            sum(!is.na(do.call(rbind, x["ests",])[,1])))
    relb <- as.vector(summarise_sims(x)["relb",])
    expect_equal(relb, rep(0, length(relb)), tol = tol)
  }
  
}

# right censored models
if (run_sims) {
  set.seed(5050)
  sims_exp <- replicate(n_sims, sim_run(basehaz = "exp"))
  validate_relbias(sims_exp)
}

if (run_sims) {
  set.seed(6060)
  sims_weibull <- replicate(n_sims, sim_run(basehaz = "weibull"))
  validate_relbias(sims_weibull)
}

if (run_sims) {
  set.seed(7070)
  sims_gompertz <- replicate(n_sims, sim_run(basehaz = "gompertz"))
  validate_relbias(sims_gompertz)
}

if (run_sims) {
  set.seed(8080)
  sims_ms <- replicate(n_sims, sim_run(basehaz = "ms"))
  validate_relbias(sims_ms)
}

# delayed entry models
if (run_sims) {
  set.seed(5050)
  sims_exp_d <- replicate(n_sims, sim_run(basehaz = "exp", delay = TRUE))
  validate_relbias(sims_exp_d)
}

if (run_sims) {
  set.seed(6060)
  sims_weibull_d <- replicate(n_sims, sim_run(basehaz = "weibull", delay = TRUE))
  validate_relbias(sims_weibull_d)
}

if (run_sims) {
  set.seed(7070)
  sims_gompertz_d <- replicate(n_sims, sim_run(basehaz = "gompertz", delay = TRUE))
  validate_relbias(sims_gompertz_d)
}

if (run_sims) {
  set.seed(8080)
  sims_ms_d <- replicate(n_sims, sim_run(basehaz = "ms", delay = TRUE))
  validate_relbias(sims_ms_d)
}

# interval censored models
if (run_sims) {
  set.seed(5050)
  sims_exp_i <- replicate(n_sims, sim_run(basehaz = "exp", icens = TRUE))
  validate_relbias(sims_exp_i)
}

if (run_sims) {
  set.seed(6060)
  sims_weibull_i <- replicate(n_sims, sim_run(basehaz = "weibull", icens = TRUE))
  validate_relbias(sims_weibull_i)
}

if (run_sims) {
  set.seed(7070)
  sims_gompertz_i <- replicate(n_sims, sim_run(basehaz = "gompertz", icens = TRUE))
  validate_relbias(sims_gompertz_i)
}

if (run_sims) {
  set.seed(8080)
  sims_ms_i <- replicate(n_sims, sim_run(basehaz = "ms", icens = TRUE))
  validate_relbias(sims_ms_i)
}


# run simulations to check piecewise constant time-varying effects
if (run_sims) {
  
  # number of simulations (for each model specification)
  n_sims <- 250
  
  # define a function to fit the model to one simulated dataset
  sim_run <- function(N = 600, return_relb = FALSE) {         
    
    # simulate data
    covs <- data.frame(id  = 1:N,
                       trt = rbinom(N, 1, 0.5))
    
    dat <- simsurv(dist    = "exponential",
                   x       = covs,
                   lambdas = c(0.15),
                   betas   = c(trt = -0.4),
                   tde     = c(trt = 0.8),
                   tdefun  = function(t) as.numeric(t > 4),
                   maxt    = 15)
    
    dat <- merge(covs, dat)
    
    # define appropriate model formula
    ff <- Surv(eventtime, status) ~ tve(trt, degree = 0, knots = 4)
    
    # fit model
    mod <- stan_surv(formula = ff,
                     data    = dat, 
                     basehaz = "exp", 
                     chains  = 1,
                     refresh = 0,
                     iter    = 1000)
    
    # true parameters (hard coded here)
    true <- c(intercept = log(0.15),
              trt       = -0.4,
              trt_tve   = 0.8)
    
    # extract parameter estimates
    ests <- c(intercept = fixef(mod)[1L],
              trt       = fixef(mod)[2L],
              trt_tve   = fixef(mod)[3L])
    
    # check Rhat
    rhats <- summary(mod)[, "Rhat"]
    rhats <- rhats[!names(rhats) %in% c("lp__", "log-posterior")]
    
    converged <- (all(rhats <= 1.1, na.rm = TRUE))
    
    if (!converged)
      ests <- rep(NA, length(ests)) # set estimates to NA if model didn't converge
    
    if (return_relb)
      return(as.vector((ests - true) / true))
    
    list(true = true,
         ests = ests,
         bias = ests - true,
         relb = (ests - true) / true)
  }
  
  # functions to summarise the simulations and check relative bias
  summarise_sims <- function(x) {
    message("Number of simulations that converged: ", 
            sum(!is.na(do.call(rbind, x["ests",])[,1])))
    rbind(true = colMeans(do.call(rbind, x["true",])),
          ests = colMeans(do.call(rbind, x["ests",])),
          bias = colMeans(do.call(rbind, x["bias",])),
          relb = colMeans(do.call(rbind, x["relb",])))
  }
  
  validate_relbias <- function(x, tol = 0.05) {
    message("Number of simulations that converged: ", 
            sum(!is.na(do.call(rbind, x["ests",])[,1])))
    relb <- as.vector(summarise_sims(x)["relb",])
    expect_equal(relb, rep(0, length(relb)), tol = tol)
  }

}  

# tve models
if (run_sims) {
  set.seed(5050)
  sims_pw <- replicate(n_sims, sim_run())
  validate_relbias(sims_pw)
}
