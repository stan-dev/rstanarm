# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2015, 2016 Trustees of Columbia University
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

expect_stanreg <- function(x) expect_s3_class(x, "stanreg")
SW <- function(expr) capture.output(suppressWarnings(expr))


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
  SW(fit <- stan_glmer.nb(formula = fmla, data = dd, chains = CHAINS, init_r = 1,
                          iter = ITER, seed = SEED, refresh = REFRESH))
  expect_stanreg(fit)
  
  ans <- glmer.nb(formula = fmla, data = dd)
  # ans is messed up
  # expect_equal(fixef(fit), fixef(ans), tol = FIXEF_tol)
  # expect_equal(ranef(fit), ranef(ans), tol = RANEF_tol)
  expect_equal(ngrps(fit), ngrps(ans))
})

context("stan_lmer")
test_that("stan_lmer returns expected result for slepstudy example", {
  fmla <- Reaction / 10 ~ Days + (Days | Subject)
  SW(fit <- stan_lmer(fmla, data = sleepstudy, refresh = REFRESH,
                      init_r = 0.05, chains = CHAINS, iter = ITER, seed = SEED))
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
expect_gg <- function(x) {
  testthat::expect_s3_class(x, "ggplot")
  invisible(ggplot2::ggplot_build(x))
}
test_that("stan_gamm4 returns stanreg object", {
  SW(fit <- stan_gamm4(Reaction / 10 ~ s(Days), data = sleepstudy, sparse = TRUE,
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

