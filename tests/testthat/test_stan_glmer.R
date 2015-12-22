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
stopifnot(require(lme4))
stopifnot(require(gamm4))
SEED <- 123
options(mc.cores = 2L)
if (interactive()) options(mc.cores = parallel::detectCores())
FIXEF_tol <- 0.05
RANEF_tol <- 0.20 


context("stan_lmer")
test_that("stan_lmer returns expected result for slepstudy example", {
  fmla <- Reaction / 10 ~ Days + (Days | Subject)
  fit <- stan_lmer(fmla, data = sleepstudy, 
                   init_r = 0.05, chains = 2, iter = 400, seed = SEED)
  ans <- lmer(fmla, data = sleepstudy)
  expect_equal(fixef(fit), fixef(ans), tol = FIXEF_tol)
  expect_equal(ranef(fit), ranef(ans), tol = RANEF_tol)
  expect_equal(ngrps(fit), ngrps(ans))
})

context("stan_lmer")
test_that("stan_lmer returns expected result for Penicillin example", {
  fmla <- diameter ~ (1|plate) + (1|sample)
  fit <- stan_lmer(fmla, data = Penicillin, chains = 2, iter = 400, seed = SEED)
  ans <- lmer(fmla, data = Penicillin)
  expect_equal(fixef(fit), fixef(ans), tol = FIXEF_tol)
  expect_equal(ranef(fit), ranef(ans), tol = RANEF_tol)
  expect_identical(ngrps(fit), ngrps(ans))
})

context("stan_glmer (binomial)")
test_that("stan_glmer returns expected result for cbpp example", {
  links <- c("logit", "probit", "cauchit", "log", "cloglog")
  # for (i in seq_along(links)) {
  i <- 1L # it seems only logit gives results similar to glmer with same link 
    fmla <- cbind(incidence, size - incidence) ~ period + (1 | herd)
    fit <- stan_glmer(fmla, data = cbpp, family = binomial(links[i]),
                      seed = SEED, chains = 2, iter = 400)
    ans <- glmer(fmla, data = cbpp, family = binomial(links[i]))
    expect_equal(fixef(fit), fixef(ans), tol = FIXEF_tol)
    expect_equal(ranef(fit), ranef(ans), tol = RANEF_tol)
    expect_equal(ngrps(fit), ngrps(ans))
  # }
})

context("stan_gamm4")
test_that("stan_gamm4 returns expected result for sleepstudy example", {
  fit <- stan_gamm4(Reaction / 10 ~ s(Days), data = sleepstudy,
                    random = ~(1|Subject), chains = 2, iter = 400, seed = SEED)
  ans <- gamm4(Reaction / 10 ~ s(Days), data = sleepstudy, 
               random = ~(1|Subject))$mer
  expect_equal(fixef(fit)[-1], fixef(ans)[-1], tol = FIXEF_tol)
  expect_equal(ranef(fit), ranef(ans), tol = RANEF_tol)
  expect_identical(ngrps(fit), ngrps(ans))
})
