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
library(rstpm2)
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
source(test_path("helpers", "expect_survfit.R"))
source(test_path("helpers", "expect_ppd.R"))
source(test_path("helpers", "expect_equivalent_loo.R"))
source(test_path("helpers", "SW.R"))
# SW <- function(expr) eval(expr)
source(test_path("helpers", "get_tols_surv.R"))
source(test_path("helpers", "recover_pars_surv.R"))

eo <- function(...) { expect_output (...) }
ee <- function(...) { expect_error  (...) }
ew <- function(...) { expect_warning(...) }
up <- function(...) { update(...) }

#-----------------------------  Models -----------------------------------

#--- Time fixed covariates, time fixed coefficients

cov1 <- data.frame(id = 1:1000,
                   x1 = stats::rbinom(1000, 1, 0.5),
                   x2 = stats::rnorm (1000, -1, 0.5))
dat1 <- simsurv(lambdas = 0.1,
                gammas  = 1.5,
                betas   = c(x1 = -0.5, x2 = -0.3),
                x       = cov1,
                maxt    = 5)
dat1 <- merge(dat1, cov1)
fm1  <- Surv(eventtime, status) ~ x1 + x2
mod1a <- stan_surv(fm1, dat1, chains = 1, refresh = 0L, iter = 1000, basehaz = "ms")
mod1b <- stan_surv(fm1, dat1, chains = 1, refresh = 0L, iter = 1000, basehaz = "bs")
mod1c <- stan_surv(fm1, dat1, chains = 1, refresh = 0L, iter = 1000, basehaz = "exp")
mod1d <- stan_surv(fm1, dat1, chains = 1, refresh = 0L, iter = 1000, basehaz = "weibull")
mod1e <- stan_surv(fm1, dat1, chains = 1, refresh = 0L, iter = 1000, basehaz = "gompertz")


#--------------------------  Arguments -----------------------------------

testmod <- mod1a

test_that("prior_PD argument works", {
  eo(update(testmod, prior_PD = TRUE))
})

test_that("adapt_delta argument works", {
  eo(up(testmod, adapt_delta = NULL))
  eo(up(testmod, adapt_delta = 0.8))
  eo(up(testmod, control = list(adapt_delta = NULL)))
  eo(up(testmod, control = list(adapt_delta = 0.8)))
})

test_that("init argument works", {
  eo(up(testmod, init = "prefit"))
  eo(up(testmod, init = "0"))
  eo(up(testmod, init = 0))
  eo(up(testmod, init = "random"))
})

test_that("qnodes argument works", {
  eo(up(testmod, qnodes = 7))
  eo(up(testmod, qnodes = 11))
  eo(up(testmod, qnodes = 15))
  ee(up(testmod, qnodes = 1),       "must be either 7, 11 or 15")
  ee(up(testmod, qnodes = c(1,2)),  "numeric vector of length 1")
  ee(up(testmod, qnodes = "wrong"), "numeric vector of length 1")
})

test_that("basehaz argument works", {

  eo(up(testmod, basehaz = "exp"))
  eo(up(testmod, basehaz = "weibull"))
  eo(up(testmod, basehaz = "gompertz"))
  eo(up(testmod, basehaz = "ms"))
  eo(up(testmod, basehaz = "bs"))
  eo(up(testmod, basehaz = "piecewise"))

  dfl <- list(df = 5)
  knl <- list(knots = c(1,3,5))
  eo(up(testmod, basehaz = "ms",        basehaz_ops = dfl))
  eo(up(testmod, basehaz = "ms",        basehaz_ops = knl))
  eo(up(testmod, basehaz = "bs",        basehaz_ops = dfl))
  eo(up(testmod, basehaz = "bs",        basehaz_ops = knl))
  eo(up(testmod, basehaz = "piecewise", basehaz_ops = dfl))
  eo(up(testmod, basehaz = "piecewise", basehaz_ops = knl))

  eo(ew(up(testmod, basehaz = "exp",     basehaz_ops = dfl), "'df' will be ignored"))
  eo(ew(up(testmod, basehaz = "exp",     basehaz_ops = knl), "'knots' will be ignored"))
  eo(ew(up(testmod, basehaz = "weibull", basehaz_ops = dfl), "'df' will be ignored"))
  eo(ew(up(testmod, basehaz = "weibull", basehaz_ops = knl), "'knots' will be ignored"))
  eo(ew(up(testmod, basehaz = "gompertz",basehaz_ops = dfl), "'df' will be ignored"))
  eo(ew(up(testmod, basehaz = "gompertz",basehaz_ops = knl), "'knots' will be ignored"))

  ee(up(testmod, basehaz_ops = list(df = 1)), "must be at least 3")
  ee(up(testmod, basehaz_ops = list(knots = -1)), "'knots' must be non-negative")
  ee(up(testmod, basehaz_ops = list(knots = c(1,2,50))), "cannot be greater than the largest event time")

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
                 gammas  = 1.3,
                 betas   = c(X1 = 0.3, X2 = -0.5),
                 x       = covs)
  dat <- merge(dat, covs)

  compare_surv(data = dat, basehaz = "weibull")
  compare_surv(data = dat, basehaz = "ms")

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
  compare_surv(data = dat, basehaz = "ms")


#--------  Check post-estimation functions work

  # fit the models
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

  # new data for predictions
  nd1 <- pbcSurv[pbcSurv$id == 2,]
  nd2 <- pbcSurv[pbcSurv$id %in% c(1,2),]

  # test the models
  for (j in c(1:30)) {

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
