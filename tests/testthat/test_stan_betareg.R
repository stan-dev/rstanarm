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

SW <- function(expr) capture.output(suppressWarnings(expr))
expect_stanreg <- function(x) expect_s3_class(x, "stanreg")

link1 <- c("logit", "probit", "cloglog", "cauchit", "log", "loglog")
link2 <- c("log", "identity", "sqrt")

# sparse currently not used in stan_betareg
context("stan_betareg sparse error")
test_that("sparse = TRUE errors", {
  dat <- list()
  dat$N <- 200
  dat$x <- rnorm(dat$N, 2, 1)
  dat$mu <- binomial(link="logit")$linkinv(1+0.2*dat$x)
  dat$phi <- 20
  dat$y <- rbeta(dat$N, dat$mu * dat$phi, (1 - dat$mu) * dat$phi)
  dat <- data.frame(dat$y, dat$x)
  colnames(dat) <- c("y", "x")
  expect_error(fit <- stan_betareg(y ~ x, link = "logit", seed = SEED, sparse = TRUE,
                                   prior = NULL, prior_intercept = NULL,
                                   data = dat, algorithm = "optimizing"))
})

# test QR 
context("stan_betareg QR error")
test_that("QR errors when x and/or z predictors are <= 1", {
  dat <- list()
  dat$N <- 200
  dat$x <- rnorm(dat$N, 2, 1)
  dat$z <- rep(0, dat$N)
  dat$mu <- binomial(link="logit")$linkinv(1+0.2*dat$x)
  dat$phi <- 20
  dat$y <- rbeta(dat$N, dat$mu * dat$phi, (1 - dat$mu) * dat$phi)
  dat <- data.frame(dat$y, dat$x, dat$z)
  colnames(dat) <- c("y", "x", "z")
  expect_error(fit <- stan_betareg(y ~ x, link = "logit", seed = SEED, QR = TRUE,
                           prior = NULL, prior_intercept = NULL,
                           data = dat, algorithm = "optimizing"))
  expect_error(fit <- stan_betareg(y ~ x | z, link = "logit", seed = SEED, QR = TRUE,
                                   prior = NULL, prior_intercept = NULL,
                                   data = dat, algorithm = "optimizing"))
})

context("stan_betareg QR tests")
test_that("QR works when x and/or z predictors are >= 1", {
  dat <- list()
  dat$N <- 200
  dat$x <- rnorm(dat$N, 2, 1)
  dat$z <- rep(0, dat$N)
  dat$mu <- binomial(link="logit")$linkinv(1+0.2*dat$x)
  dat$phi <- 20
  dat$y <- rbeta(dat$N, dat$mu * dat$phi, (1 - dat$mu) * dat$phi)
  dat <- data.frame(dat$y, dat$x, dat$z)
  colnames(dat) <- c("y", "x", "z")
  fit1 <- stan_betareg(y ~ x + z, link = "logit", seed = SEED, QR = TRUE,
                                   prior = NULL, prior_intercept = NULL,
                                   data = dat, algorithm = "optimizing")
  fit2 <- stan_betareg(y ~ x | z, link = "logit", seed = SEED, QR = TRUE,
                                   prior = NULL, prior_intercept = NULL,
                                   data = dat, algorithm = "optimizing")
})

context("stan_betareg (dispersion only)")
test_that("stan_betareg returns expected result when modeling x and dispersion", {
  dat <- list()
  dat$N <- 200
  dat$x <- rnorm(dat$N, 2, 1)
  dat$z <- rep(0, dat$N)
  dat$mu <- binomial(link="logit")$linkinv(1+0.2*dat$x)
  dat$phi <- 20
  dat$y <- rbeta(dat$N, dat$mu * dat$phi, (1 - dat$mu) * dat$phi)
  dat <- data.frame(dat$y, dat$x)
  colnames(dat) <- c("y", "x")
  for (i in 1:length(link1)) {
    SW(fit <- stan_betareg(y ~ x, link = link1[i], seed = SEED,
                           prior = NULL, prior_intercept = NULL,
                           data = dat, algorithm = "optimizing"))
    expect_stanreg(fit)
    val <- coef(fit)
    ans <- coef(betareg(y ~ x, link = link1[i], data = dat))
    expect_equal(val, ans, tol = 0.1, info = link1[i])
    # cat("... used link = ", link1[i], "\n")
  }
})

test_that("stan_betareg works with QR = TRUE and algorithm = 'optimizing'", {
  dat <- list()
  dat$N <- 200
  dat$x <- rnorm(dat$N, 2, 1)
  dat$z <- rnorm(dat$N, 0, 1)
  dat$mu <- binomial(link="logit")$linkinv(1+0.2*dat$x + dat$z)
  dat$phi <- 20
  dat$y <- rbeta(dat$N, dat$mu * dat$phi, (1 - dat$mu) * dat$phi)
  dat <- data.frame(dat$y, dat$x, dat$z)
  colnames(dat) <- c("y", "x", "z")
  SW(fit <- stan_betareg(y ~ x + z, link = "logit", seed = SEED, QR = TRUE,
                         prior = NULL, prior_intercept = NULL,
                         data = dat, algorithm = "optimizing"))
  expect_stanreg(fit)
  val <- coef(fit)
  ans <- coef(betareg(y ~ x + z, link = "logit", data = dat))
  expect_equal(val, ans, tol = 0.1, info = "logit")
})

test_that("stan_betareg works with QR = TRUE and algorithm = 'sampling'", {
  dat <- list()
  dat$N <- 200
  dat$x <- rnorm(dat$N, 2, 1)
  dat$z <- rnorm(dat$N, 0, 1)
  dat$mu <- binomial(link="logit")$linkinv(1+0.2*dat$x + dat$z)
  dat$phi <- 20
  dat$y <- rbeta(dat$N, dat$mu * dat$phi, (1 - dat$mu) * dat$phi)
  dat <- data.frame(dat$y, dat$x, dat$z)
  colnames(dat) <- c("y", "x", "z")
    SW(fit <- stan_betareg(y ~ x + z, link = "logit", QR = TRUE,
                           prior = NULL, prior_intercept = NULL, iter = 100, chains = 2,
                           data = dat))
    expect_stanreg(fit)
    val <- coef(fit)
    ans <- coef(betareg(y ~ x + z, link = "logit", data = dat))
    expect_equal(val, ans, tol = 0.1)
})

context("stan_betareg (x and z using link.phi = 'log')")
test_that("stan_betareg returns expected result when modeling x and z", {
  dat <- list()
  dat$N <- 200
  dat$x <- rnorm(dat$N, 2, 1)
  dat$z <- rnorm(dat$N, 2, 1)
  dat$mu <- binomial(link = "logit")$linkinv(0.5 + 0.2*dat$x)
  dat$phi <- poisson(link = link2[1])$linkinv(1.5 + 0.4*dat$z)
  dat$y <- rbeta(dat$N, dat$mu * dat$phi, (1 - dat$mu) * dat$phi)
  dat <- data.frame(dat$y, dat$x, dat$z)
  colnames(dat) <- c("y", "x", "z")
  for (i in 1:length(link1)) {
    # cat("... using link =", link1[i], "and link.phi =", link2[1], "\n")
    SW(fit <- stan_betareg(y ~ x | z, link = link1[i], link.phi = link2[1], 
                           seed = SEED,
                           prior = NULL, prior_intercept = NULL,
                           prior_z = NULL, prior_intercept_z = NULL,
                           data = dat, algorithm = "optimizing"))
    expect_stanreg(fit)
    val <- coef(fit)
    ans <- coef(betareg(y ~ x | z, link = link1[i], link.phi = link2[1], data = dat))
    expect_equal(val, ans, tol = 0.1, info = c(link1[i], link2[1]))
  }
})

context("stan_betareg (x and z using link.phi = 'identity')")
# tests use sampling instead of optimizing (the latter fails)
test_that("stan_betareg returns expected result when modeling x and z", {
  dat <- list()
  dat$N <- 200
  dat$x <- rnorm(dat$N, 2, 1)
  dat$z <- rnorm(dat$N, 2, 1)
  dat$mu <- binomial(link = "logit")$linkinv(1 + 0.2*dat$x)
  dat$phi <- 0.5*dat$z
  dat$phi <- dat$z - min(dat$z) + 5.5
  dat$y <- rbeta(dat$N, dat$mu * dat$phi, (1 - dat$mu) * dat$phi)
  dat <- data.frame(dat$y, dat$x, dat$z)
  colnames(dat) <- c("y", "x", "z")
  for (i in 1:length(link1)) {
    # cat("... using link =", link1[i], "and link.phi =", link2[2], "\n")
    SW(fit <- stan_betareg(y ~ x | z, link = link1[i], link.phi = link2[2], 
                           seed = SEED,
                           prior = NULL, prior_intercept = NULL,
                           prior_z = NULL, prior_intercept_z = NULL,
                           data = dat, algorithm = "sampling", 
                           chains = 2, iter = 300))
    expect_stanreg(fit)
    val <- coef(fit)
    ans <- coef(betareg(y ~ x | z, link = link1[i], link.phi = link2[2], data = dat))
    expect_equal(val, ans, tol = 0.15, info = c(link1[i], link2[2]))
  }
})

# sqrt link is unstable so only testing that the model runs. 
context("stan_betareg (x and z using link.phi = 'sqrt')")
test_that("stan_betareg returns expected result when modeling x and z using link.phi = 'sqrt'", {
  # for (i in 1:length(link1)) {
  for (i in 1:length(link1)) {  # FIXME!
    dat <- list()
    dat$N <- 300
    dat$x <- rnorm(dat$N, 2, 1)
    # dat$z <- rnorm(dat$N, 2, 1)
    dat$z <- rep(1, dat$N)
    dat$mu <- binomial(link = "logit")$linkinv(-0.8 + 0.5*dat$x)
    dat$phi <- poisson(link = "sqrt")$linkinv(8 + 2*dat$z)
    dat$y <- rbeta(dat$N, dat$mu * dat$phi, (1 - dat$mu) * dat$phi)
    dat <- data.frame(dat$y, dat$x, dat$z)
    colnames(dat) <- c("y", "x", "z")

    # cat("... using link =", link1[i], "and link.phi =", link2[3], "\n")
    SW(fit <- stan_betareg(y ~ x | 1, link = link1[i], link.phi = link2[3], 
                        seed = SEED,
                        prior = NULL, prior_intercept = NULL,
                        prior_z = NULL, prior_intercept_z = NULL,
                        data = dat, algorithm = "sampling", 
                        chains = 2, iter = 100))
    expect_stanreg(fit)
  }
})

# test weights/offset (make test more comprehensive once the beta_rng() update is in stan math)
context("stan_betareg (dispersion only) with offset and weights")
test_that("stan_betareg returns expected result when modeling x and dispersion with offset and weights", {
  dat <- list()
  dat$N <- 200
  weights <- rbeta(dat$N, 2, 2)
  offset <- rep(0.3, dat$N)
  dat$x <- rnorm(dat$N, 2, 1)
  dat$mu <- binomial(link="logit")$linkinv(1+0.2*dat$x)
  dat$phi <- 20
  dat$y <- rbeta(dat$N, dat$mu * dat$phi, (1 - dat$mu) * dat$phi)
  dat <- data.frame(dat$y, dat$x)
  colnames(dat) <- c("y", "x")
  SW(fit <- stan_betareg(y ~ x, link = "logit", seed = SEED,
                         prior = NULL, prior_intercept = NULL,
                         data = dat, weights = weights, offset = offset, 
                         algorithm = "optimizing", iter = 2000))
  expect_stanreg(fit)
  val <- coef(fit)
  ans <- coef(betareg(y ~ x, link = "logit", weights = weights, offset = offset, data = dat))
  expect_equal(val, ans, tol = 0.3, info = "logit")
})
