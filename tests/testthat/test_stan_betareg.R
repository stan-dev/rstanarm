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
library(betareg)
SEED <- 12345
set.seed(SEED)

context("stan_betareg")

SW <- function(expr) capture.output(suppressWarnings(expr))
expect_stanreg <- function(x) expect_s3_class(x, "stanreg")

simple_betareg_data <- function(N, draw_z = FALSE) {
  x <- rnorm(N, 2, 1)
  z <- if (draw_z) rnorm(N, 0, 1) else rep(0, N)
  mu <- binomial(link="logit")$linkinv(1 + 0.2 * x)
  phi <- 20
  y <- rbeta(N, mu * phi, (1 - mu) * phi)
  data.frame(y,x,z)
}

link1 <- c("logit", "probit", "cloglog", "cauchit", "log", "loglog")
link2 <- c("log", "identity", "sqrt")

# sparse currently not used in stan_betareg
test_that("sparse = TRUE errors", {
  dat <- simple_betareg_data(200)
  expect_error(
    stan_betareg(y ~ x, link = "logit", seed = SEED, sparse = TRUE, data = dat),
    "unknown arguments: sparse"
  )
})

# test QR 
test_that("QR errors when number of x and/or z predictors is <= 1", {
  dat <- simple_betareg_data(200, draw_z = TRUE)
  expect_error(
    stan_betareg(y ~ x, link = "logit", seed = SEED, QR = TRUE, data = dat),
    "'QR' can only be specified when there are multiple predictors"
  )
  expect_error(
    stan_betareg(y ~ x | z, link = "logit", seed = SEED, QR = TRUE, data = dat),
    "'QR' can only be specified when there are multiple predictors"
  )
})

test_that("QR works when number of x and/or z predictors is >= 1", {
  dat <- simple_betareg_data(200, draw_z = TRUE)
  SW(fit1 <- stan_betareg(y ~ x + z, link = "logit", seed = SEED, QR = TRUE,
                                   prior = NULL, prior_intercept = NULL,
                                   data = dat, algorithm = "optimizing"))
  expect_stanreg(fit1)
  expect_output(print(prior_summary(fit1)), "Q-space")
  
  SW(fit2 <- stan_betareg(y ~ x + z | z, link = "logit", seed = SEED, QR = TRUE,
                                   prior = NULL, prior_intercept = NULL,
                                   data = dat, algorithm = "optimizing"))
  expect_stanreg(fit2)
})

test_that("stan_betareg returns expected result when modeling x and dispersion", {
  dat <- simple_betareg_data(200)
  for (i in 1:length(link1)) {
    SW(fit <- stan_betareg(y ~ x, link = link1[i], seed = SEED,
                           prior = NULL, prior_intercept = NULL,
                           data = dat, algorithm = "optimizing"))
    expect_stanreg(fit)
    val <- coef(fit)
    ans <- coef(betareg(y ~ x, link = link1[i], data = dat))
    expect_equal(val, ans, tol = 0.1, info = link1[i])
  }
})

test_that("stan_betareg works with QR = TRUE and algorithm = 'optimizing'", {
  dat <- simple_betareg_data(200, draw_z = TRUE)
  SW(fit <- stan_betareg(y ~ x + z, link = "logit", seed = SEED, QR = TRUE,
                         prior = NULL, prior_intercept = NULL,
                         data = dat, algorithm = "optimizing"))
  expect_stanreg(fit)
  val <- coef(fit)
  ans <- coef(betareg(y ~ x + z, link = "logit", data = dat))
  expect_equal(val, ans, tol = 0.1, info = "logit")
})

test_that("stan_betareg works with QR = TRUE and algorithm = 'sampling'", {
  dat <- simple_betareg_data(200, draw_z = TRUE)
  SW(fit <- stan_betareg(y ~ x + z, link = "logit", QR = TRUE,
                         prior = NULL, prior_intercept = NULL, 
                         iter = 100, chains = 2, data = dat))
  expect_stanreg(fit)
  val <- coef(fit)
  ans <- coef(betareg(y ~ x + z, link = "logit", data = dat))
  expect_equal(val, ans, tol = 0.1)
})

test_that("stan_betareg ok when modeling x and z (link.phi = 'log')", {
  N <- 200
  dat <- data.frame(x = rnorm(N, 2, 1), z = rnorm(N, 2, 1))
  mu <- binomial(link="logit")$linkinv(1 + 0.2 * dat$x)
  phi <- poisson(link = link2[1])$linkinv(1.5 + 0.4*dat$z)
  dat$y <- rbeta(N, mu * phi, (1 - mu) * phi)
  
  for (i in 1:length(link1)) {
    SW(fit <- stan_betareg(y ~ x | z, link = link1[i], link.phi = link2[1], 
                           seed = SEED,
                           prior = NULL, prior_intercept = NULL,
                           prior_z = NULL, prior_intercept_z = NULL,
                           data = dat, algorithm = "optimizing"))
    expect_stanreg(fit)
    val <- coef(fit)
    ans <- coef(betareg(y ~ x | z, link = link1[i], link.phi = link2[1], 
                        data = dat))
    expect_equal(val, ans, tol = 0.1, info = c(link1[i], link2[1]))
  }
})

# tests use sampling instead of optimizing (the latter fails)
test_that("stan_betareg ok when modeling x and z (link.phi = 'identity')", {
  N <- 200
  dat <- data.frame(x = rnorm(N, 2, 1), z = rnorm(N, 2, 1))
  mu <- binomial(link = "logit")$linkinv(1 + 0.2*dat$x)
  phi <- dat$z - min(dat$z) + 5.5
  dat$y <- rbeta(N, mu * phi, (1 - mu) * phi)
  for (i in 1:length(link1)) {
    SW(fit <- stan_betareg(y ~ x | z, link = link1[i], link.phi = link2[2],
                           prior = NULL, prior_intercept = NULL,
                           prior_z = NULL, prior_intercept_z = NULL,
                           data = dat, algorithm = "sampling", 
                           chains = 2, iter = 300, seed = SEED))
    expect_stanreg(fit)
    val <- coef(fit)
    ans <- coef(betareg(y ~ x | z, link = link1[i], link.phi = link2[2], data = dat))
    expect_equal(val, ans, tol = 0.15, info = c(link1[i], link2[2]))
  }
})

# sqrt link is unstable so only testing that the model runs. 
test_that("stan_betareg ok when modeling x and z (link.phi = 'sqrt')", {
  for (i in 1:length(link1)) {  # FIXME!
    N <- 300
    dat <- data.frame(x = rnorm(N, 2, 1), z = rep(1, N))
    mu <- binomial(link = "logit")$linkinv(-0.8 + 0.5*dat$x)
    phi <- poisson(link = "sqrt")$linkinv(8 + 2*dat$z)
    dat$y <- rbeta(N, mu * phi, (1 - mu) * phi)

    SW(fit <- stan_betareg(y ~ x | 1, link = link1[i], link.phi = link2[3], 
                           data = dat, algorithm = "sampling",
                           chains = 1, iter = 1)) 
    expect_stanreg(fit)
  }
})

# test weights/offset (make test more comprehensive once the beta_rng() update is in stan math)
test_that("stan_betareg ok when modeling x and dispersion with offset and weights", {
  N <- 200
  weights <- rbeta(N, 2, 2)
  offset <- rep(0.3, N)
  dat <- data.frame(x = rnorm(N, 2, 1))
  mu <- binomial(link="logit")$linkinv(1+0.2*dat$x)
  phi <- 20
  dat$y <- rbeta(N, mu * phi, (1 - mu) * phi)
  SW(fit <- stan_betareg(y ~ x, link = "logit", seed = SEED,
                         prior = NULL, prior_intercept = NULL,
                         data = dat, weights = weights, offset = offset, 
                         algorithm = "optimizing", iter = 2000))
  expect_stanreg(fit)
  val <- coef(fit)
  ans <- coef(betareg(y ~ x, link = "logit", weights = weights, offset = offset, data = dat))
  expect_equal(val, ans, tol = 0.3, info = "logit")
})

test_that("heavy tailed priors work with stan_betareg", {
  dat <- simple_betareg_data(100, TRUE)
  expect_output(stan_betareg(y ~ x | z, data = dat, 
                             prior = product_normal(), prior_z = product_normal(), 
                             chains = 1, iter = 1))
  expect_output(stan_betareg(y ~ x | z, data = dat, 
                             prior = laplace(), prior_z = laplace(), 
                             chains = 1, iter = 1))
  expect_output(stan_betareg(y ~ x | z, data = dat, 
                             prior = lasso(), prior_z = lasso(), 
                             chains = 1, iter = 1))
})
