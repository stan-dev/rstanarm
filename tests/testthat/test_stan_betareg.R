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

if (require(betareg)) {
  library(rstanarm)
  SEED <- 12345
  set.seed(SEED)
  ITER <- 10
  CHAINS <- 2
  REFRESH <- 0
  
  context("stan_betareg")
  
  SW <- suppressWarnings
  expect_stanreg <- function(x) expect_s3_class(x, "stanreg")
  
  simple_betareg_data <- function(N, draw_z = FALSE) {
    x <- rnorm(N, 2, 1)
    z <- if (draw_z) rnorm(N, 0, 1) else rep(0, N)
    mu <- binomial(link="logit")$linkinv(1 + 0.2 * x)
    phi <- 20
    y <- rbeta(N, mu * phi, (1 - mu) * phi)
    data.frame(y,x,z)
  }
  dat <- simple_betareg_data(200, draw_z = TRUE)
  
  link1 <- c("logit", "probit", "cloglog", "cauchit", "log", "loglog")
  link2 <- c("log", "identity", "sqrt")
  
  # sparse currently not used in stan_betareg
  test_that("sparse = TRUE errors", {
    expect_error(
      stan_betareg(y ~ x, link = "logit", seed = SEED, sparse = TRUE, data = dat),
      "unknown arguments: sparse"
    )
  })
  
  # test QR 
  test_that("QR errors when number of x and/or z predictors is <= 1", {
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
    SW(fit <- stan_betareg(y ~ x + z, link = "logit", seed = SEED, QR = TRUE,
                           prior = NULL, prior_intercept = NULL,
                           data = dat, algorithm = "optimizing"))
    expect_stanreg(fit)
    val <- coef(fit)
    ans <- coef(betareg(y ~ x + z, link = "logit", data = dat))
    expect_equal(val, ans, tol = 0.1, info = "logit")
  })
  
  test_that("stan_betareg works with QR = TRUE and algorithm = 'sampling'", {
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
                             data = dat, algorithm = "optimizing", 
                             seed = SEED))
      expect_stanreg(fit)
      val <- coef(fit)
      ans <- coef(betareg(y ~ x | z, link = link1[i], link.phi = link2[2], data = dat))
      expect_equal(val, ans, tol = 0.15, info = c(link1[i], link2[2]))
    }
  })
  
  # sqrt link is unstable so only testing that the model runs. 
  test_that("stan_betareg ok when modeling x and z (link.phi = 'sqrt')", {
    for (i in 1:length(link1)) {  # FIXME!
      N <- 1000
      dat <- data.frame(x = rnorm(N, 2, 1), z = rep(1, N))
      mu <- binomial(link = "logit")$linkinv(-0.8 + 0.5*dat$x)
      phi <- poisson(link = "sqrt")$linkinv(8 + 2*dat$z)
      dat$y <- rbeta(N, mu * phi, (1 - mu) * phi)
  
      SW(fit <- stan_betareg(y ~ x | 1, link = link1[i], link.phi = link2[3], 
                             data = dat, algorithm = "sampling", chains = 1, iter = 1)) 
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
  
  expect_equivalent_loo <- function(fit) {
    l <- suppressWarnings(loo(fit))
    w <- suppressWarnings(waic(fit))
    expect_s3_class(l, "loo")
    expect_s3_class(w, "loo")
    expect_s3_class(w, "waic")
    
    att_names <- c("names", "log_lik_dim", "class", "name", "discrete", "yhash")
    expect_named(attributes(l), att_names)
    expect_named(attributes(w), att_names)
    
    discrete <- attr(l, "discrete")
    expect_true(!is.na(discrete) && is.logical(discrete))
    
    expect_equivalent(l, suppressWarnings(loo(log_lik(fit))))
    expect_equivalent(w, suppressWarnings(waic(log_lik(fit))))
  }
  
  test_that("loo/waic for stan_betareg works", {
    data("GasolineYield", package = "betareg")
    SW(fit_logit <- stan_betareg(yield ~ batch + temp | temp, data = GasolineYield,
                                 link = "logit",
                                 chains = CHAINS, iter = ITER,
                                 seed = SEED, refresh = REFRESH))
    expect_equivalent_loo(fit_logit)
    expect_identical(ll_fun(fit_logit), rstanarm:::.ll_beta_i)
  })
  
  N <- 200
  x <- rnorm(N, 2, 1)
  z <- rnorm(N, 2, 1)
  mu <- binomial(link = "logit")$linkinv(1 + 0.2*x)
  phi <- exp(1.5 + 0.4*z)
  y <- rbeta(N, mu * phi, (1 - mu) * phi)
  fake_dat <- data.frame(y, x, z)
  remove(N, x, y, z, mu, phi)
  
  capture.output(
    stan_betareg1 <- stan_betareg(y ~ x | z, data = fake_dat, 
                                  link = "logit", link.phi = "log",
                                  iter = ITER, chains = CHAINS, seed = SEED)
  )
  betareg1 <- betareg(y ~ x | z, data = fake_dat, 
                      link = "logit", link.phi = "log")
  
  test_that("a bunch of methods stuff works properly for stan_betareg", {
    expect_equal(resid(stan_betareg1), stan_betareg1$residuals)
    expect_equal(coef(stan_betareg1), stan_betareg1$coefficients)
    expect_equal(vcov(stan_betareg1), stan_betareg1$covmat)
    expect_equal(fitted(stan_betareg1), stan_betareg1$fitted.values)
    expect_equal(se(stan_betareg1), stan_betareg1$ses)
    
    expect_error(confint(stan_betareg1), regexp = "use posterior_interval")
    expect_silent(ci6 <- posterior_interval(stan_betareg1, prob = 0.5))
    expect_identical(colnames(ci6), c("25%", "75%"))
    expect_equal(log_lik(stan_betareg1), log_lik(stan_betareg1, newdata = fake_dat))
    expect_error(ngrps(stan_betareg1), "stan_glmer and stan_lmer models only")
    expect_equal(nobs(stan_betareg1), nobs(betareg1))
    expect_equal(dimnames(vcov(stan_betareg1)), dimnames(vcov(betareg1)))
    expect_output(print(stan_betareg1, digits = 2), "stan_betareg")
    
    # errors not modeled in stan_betareg
    expect_double <- function(x) expect_type(x, "double")
    expect_double(sig <- rstanarm::sigma(stan_betareg1))
    expect_true(identical(sig, 1))
    expect_error(VarCorr(stan_betareg1), "stan_glmer and stan_lmer models only")
    expect_error(ranef(stan_betareg1), "stan_glmer and stan_lmer models only")
    expect_identical(model.frame(stan_betareg1), model.frame(betareg1))
    expect_identical(terms(stan_betareg1), terms(betareg1))
    expect_identical(formula(stan_betareg1), formula(betareg1))
    expect_output(print(prior_summary(stan_betareg1)),
                  "stan_betareg1")
    expect_named(prior_summary(stan_betareg1),
                 c("prior", "prior_z", "prior_intercept", "prior_intercept_z", "prior_aux"))  
    expect_error(predictive_error(stan_betareg1, draws = 600),
                 "'draws' should be <= posterior sample size")
    expect_error(predictive_interval(stan_betareg1, draws = 600),
                 "'draws' should be <= posterior sample size")
    expect_error(predictive_interval(stan_betareg1, prob = c(0.25, 0.76)),
                 "'prob' should be a single number greater than 0 and less than 1")
    
    # stan_betareg
    expect_warning(s <- summary(stan_betareg1, pars = "varying"),
                   regexp = "No group-specific parameters. 'varying' ignored.")
    expect_silent(s <- summary(stan_betareg1, pars = c("alpha", "beta"), digits = 3))
    expect_s3_class(s, "summary.stanreg")
    expect_output(print(s), "stan_betareg")
    expect_identical(attr(s, "algorithm"), "sampling")
    
    # betareg
    mat <- as.matrix(stan_betareg1)
    df <- as.data.frame(stan_betareg1)
    arr <- as.array(stan_betareg1)
    expect_identical(df, as.data.frame(mat))
    expect_identical(mat[1:2, 1], arr[1:2, 1, 1])
    expect_equal(dim(mat), c(floor(ITER/2) * CHAINS, 4L))
    expect_equal(dim(arr), c(floor(ITER/2), CHAINS, 4L))
    expect_identical(last_dimnames(mat), c("(Intercept)", "x", "(phi)_(Intercept)", "(phi)_z"))
    expect_identical(last_dimnames(arr), last_dimnames(mat))
    
    SW(capture.output(fit3 <- update(stan_betareg1, 
                                     iter = ITER * 2, chains = CHAINS * 2)))
    pss <- rstanarm:::posterior_sample_size
    expect_equal(pss(fit3), 4 * pss(stan_betareg1))
    
    preds <- posterior_predict(stan_betareg1, seed = 123)
    expect_equal(
      predictive_error(stan_betareg1, seed = 123),
      predictive_error(preds, y = stan_betareg1$y)
    )
    preds <- posterior_predict(stan_betareg1, seed = 123)
    expect_equal(
      predictive_interval(stan_betareg1, seed = 123),
      predictive_interval(preds)
    )
    
  })
  
  # These tests just make sure that posterior_predict doesn't throw errors and
  # that result has correct dimensions
  check_for_error <- function(fit, data = NULL, offset = NULL) {
    nsims <- nrow(as.data.frame(fit))
    mf <- if (!is.null(data)) 
      data else model.frame(fit)
    if (identical(deparse(substitute(fit)), "example_model"))
      mf <- lme4::cbpp
    
    expect_silent(yrep1 <- posterior_predict(fit))
    expect_silent(lin1 <- posterior_linpred(fit))
    expect_silent(posterior_linpred(fit, transform = TRUE))
    expect_equal(dim(yrep1), c(nsims, nobs(fit)))
    expect_equal(dim(lin1), c(nsims, nobs(fit)))
    
    expect_silent(yrep2 <- posterior_predict(fit, draws = 1))
    expect_equal(dim(yrep2), c(1, nobs(fit)))
    
    offs <- if (!is.null(offset)) offset[1] else offset
    expect_silent(yrep3 <- posterior_predict(fit, newdata = mf[1,], offset = offs))
    expect_silent(lin3 <- posterior_linpred(fit, newdata = mf[1,], offset = offs))
    expect_equal(dim(yrep3), c(nsims, 1))
    expect_equal(dim(lin3), c(nsims, 1))
    
    expect_silent(yrep4 <- posterior_predict(fit, draws = 2, newdata = mf[1,], offset = offs))
    expect_equal(dim(yrep4), c(2, 1))
    
    offs <- if (!is.null(offset)) offset[1:5] else offset
    expect_silent(yrep5 <- posterior_predict(fit, newdata = mf[1:5,], offset = offs))
    expect_silent(lin5 <- posterior_linpred(fit, newdata = mf[1:5,], offset = offs))
    expect_equal(dim(yrep5), c(nsims, 5))
    expect_equal(dim(lin5), c(nsims, 5))
    
    expect_silent(yrep6 <- posterior_predict(fit, draws = 3, newdata = mf[1:5,], offset = offs))
    expect_equal(dim(yrep6), c(3, 5))
    
    expect_error(posterior_predict(fit, draws = nsims + 1), 
                 regexep = "posterior sample size is only")
  }
  
  expect_linpred_equal <- function(object, tol = 0.1) {
    linpred <- posterior_linpred(object)
    expect_equal(apply(linpred, 2, median), object$linear.predictors, 
                 tolerance = tol, 
                 check.attributes = FALSE)
  }
  
  context("posterior_predict (stan_betareg)")
  test_that("compatible with stan_betareg with z", {
    data("GasolineYield", package = "betareg")
    fit <- SW(stan_betareg(yield ~ pressure + temp | temp, data = GasolineYield,
                           iter = ITER*5, chains = 2*CHAINS, seed = SEED, 
                           refresh = REFRESH))
    check_for_error(fit)
    # expect_linpred_equal(fit)
  })
  
  test_that("compatible with stan_betareg without z", {
    data("GasolineYield", package = "betareg")
    fit <- SW(stan_betareg(yield ~ temp, data = GasolineYield, 
                           iter = ITER, chains = CHAINS, seed = SEED, refresh = REFRESH))
    check_for_error(fit)
    # expect_linpred_equal(fit)
  })
  
  test_that("compatible with betareg with offset", {
    GasolineYield2 <- GasolineYield
    GasolineYield2$offs <- runif(nrow(GasolineYield2))
    fit <- SW(stan_betareg(yield ~ temp, data = GasolineYield2, offset = offs,
                           iter = ITER*5, chains = CHAINS, seed = SEED, refresh = REFRESH))
    fit2 <- SW(stan_betareg(yield ~ temp + offset(offs), data = GasolineYield2,
                            iter = ITER*5, chains = CHAINS, seed = SEED, refresh = REFRESH))
    
    expect_warning(posterior_predict(fit, newdata = GasolineYield), 
                   "offset")
    check_for_error(fit, data = GasolineYield2, offset = GasolineYield2$offs)
    check_for_error(fit2, data = GasolineYield2, offset = GasolineYield2$offs)
    expect_linpred_equal(fit)
    expect_linpred_equal(fit2)
  })
  
  test_that("predict ok for stan_betareg", {
    dat <- list()
    dat$N <- 200
    dat$x <- rnorm(dat$N, 2, 1)
    dat$z <- rnorm(dat$N, 2, 1)
    dat$mu <- binomial(link = "logit")$linkinv(0.5 + 0.2*dat$x)
    dat$phi <- exp(1.5 + 0.4*dat$z)
    dat$y <- rbeta(dat$N, dat$mu * dat$phi, (1 - dat$mu) * dat$phi)
    dat <- data.frame(dat$y, dat$x, dat$z)
    colnames(dat) <- c("y", "x", "z")
    
    betaregfit <- betareg(y ~ x | z, data = dat)
    SW(capture.output(
      stanfit <- stan_betareg(y ~ x | z, data = dat, chains = CHAINS,
                              iter = ITER, seed = SEED, refresh = REFRESH)
    ))
    
    pb <- predict(betaregfit, type = "response")
    ps <- predict(stanfit, type = "response")
    # expect_equal(pb, ps, tol = 0.05)
    expect_error(presp(stanfit))
    
    newd <- data.frame(x = c(300,305))
    pb <- predict(betaregfit, newdata = newd, type = "link")
    ps <- predict(stanfit, newdata = newd, type = "link")
    # expect_equal(pb, ps, tol = 0.05)
  })

}
