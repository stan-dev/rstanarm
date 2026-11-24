# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2015, 2016, 2017 Trustees of Columbia University
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

suppressPackageStartupMessages(library(rstanarm))
stopifnot(require(lme4))
# stopifnot(require(gamm4))
stopifnot(require(HSAUR3))
ITER <- 400
CHAINS <- 2
SEED <- 123
REFRESH <- ITER
set.seed(SEED)
FIXEF_tol <- 0.05
RANEF_tol <- 0.25

if (!exists("example_model")) {
  example_model <- run_example_model()
}

SW(fit <- stan_lmer(Reaction / 10 ~ Days + (Days | Subject), 
                    data = sleepstudy, refresh = 0,
                    init_r = 0.05, chains = CHAINS, iter = ITER, seed = SEED))

test_that("draws from stan_glmer (gaussian) same as from stan_lmer", {
  SW(fit1 <- stan_glmer(mpg ~ wt + (1|cyl), data = mtcars, 
                        iter = 10, chains = 1, seed = SEED, refresh = 0))
  SW(fit2 <- stan_lmer(mpg ~ wt + (1|cyl), data = mtcars, 
                       iter = 10, chains = 1, seed = SEED, refresh = 0))
  expect_identical(as.matrix(fit1), as.matrix(fit2))
})
test_that("stan_glmer returns expected result for binomial cbpp example", {
  links <- c("logit", "probit", "cauchit", "log", "cloglog")
  # for (i in seq_along(links)) {
  i <- 1L # it seems only logit gives results similar to glmer with same link 
    fmla <- cbind(incidence, size - incidence) ~ period + (1 | herd)
    SW(fit <- stan_glmer(fmla, data = cbpp, family = binomial(links[i]),
                      chains = CHAINS, iter = ITER, seed = SEED, refresh = 0))
    expect_stanreg(fit)
    
    ans <- glmer(fmla, data = cbpp, family = binomial(links[i]))
    expect_equal(fixef(fit), fixef(ans), tolerance = FIXEF_tol)
    expect_equal(ranef(fit), ranef(ans), tolerance = RANEF_tol, ignore_attr = TRUE)
    expect_equal(ngrps(fit), ngrps(ans))
  # }
})

test_that("stan_glmer.nb ok", {
  dd <- expand.grid(f1 = factor(1:3),
                    f2 = LETTERS[1:2], g=1:9, rep=1:15,
                    KEEP.OUT.ATTRS=FALSE)
  mu <- 5*(-4 + with(dd, as.integer(f1) + 4*as.numeric(f2)))
  dd$y <- rnbinom(nrow(dd), mu = mu, size = 0.5)
  fmla <- as.formula(y ~ f1*f2 + (1|g))
  SW(fit <- stan_glmer.nb(formula = fmla, data = dd, init_r = 1, refresh = 0,
                          iter = ITER, seed = SEED, algorithm = "meanfield"))
  expect_stanreg(fit)
  
  ans <- glmer.nb(formula = fmla, data = dd)
  # ans is messed up
  # expect_equal(fixef(fit), fixef(ans), tol = FIXEF_tol)
  # expect_equal(ranef(fit), ranef(ans), tol = RANEF_tol)
  expect_equal(ngrps(fit), ngrps(ans))
})

test_that("stan_lmer returns expected result for slepstudy example", {
  fmla <- formula(fit)
  expect_stanreg(fit)
  
  ans <- lmer(fmla, data = sleepstudy)
  expect_equal(fixef(fit), fixef(ans), tolerance = FIXEF_tol)
  # expect_equal(ranef(fit), ranef(ans), tol = RANEF_tol)
  expect_equal(ngrps(fit), ngrps(ans))
})
test_that("stan_lmer returns expected result for Penicillin example", {
  fmla <- as.formula(diameter ~ (1|plate) + (1|sample))
  SW(fit <- stan_lmer(fmla, data = Penicillin, chains = CHAINS, iter = ITER, 
                      seed = SEED, refresh = 0, sparse = TRUE))
  expect_stanreg(fit)
  
  ans <- lmer(fmla, data = Penicillin)
  expect_equal(fixef(fit), fixef(ans), tolerance = FIXEF_tol)
  expect_equal(ranef(fit), ranef(ans), tolerance = RANEF_tol, ignore_attr = TRUE)
  expect_identical(ngrps(fit), ngrps(ans))
})
test_that("stan_lmer ok if global intercept forced to 0", {
  SW(fit <- stan_lmer(mpg ~ 0 + (1|cyl), data = mtcars, iter = 10, 
                      seed = SEED, refresh = 0))
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

test_that("error if y is character", {
  expect_error(
    stan_lmer(as.character(mpg) ~ wt + (1|cyl),  data = mtcars),
    "Outcome variable can't be type 'character'"
  )
  expect_error(
    stan_glmer.nb(as.character(mpg) ~ wt + (1|cyl),  data = mtcars),
    "Outcome variable can't be type 'character'"
  )
})

test_that("stan_gamm4 returns stanreg object", {
  skip_if_not_installed("mgcv")
  sleepstudy$y <- sleepstudy$Reaction / 10
  SW(fit <- stan_gamm4(y ~ s(Days), data = sleepstudy, sparse = TRUE,
                       random = ~(1|Subject), chains = CHAINS, iter = ITER, 
                       seed = SEED, refresh = 0))
  expect_stanreg(fit)
  # ans <- gamm4(Reaction / 10 ~ s(Days), data = sleepstudy, 
  #              random = ~(1|Subject))$mer
  # expect_equal(fixef(fit)[-1], fixef(ans)[-1], tol = FIXEF_tol, ignore_attr = TRUE)
  # expect_equal(ranef(fit), ranef(ans), tol = RANEF_tol)
  # expect_identical(ngrps(fit), ngrps(ans))
  
  p1 <- plot_nonlinear(fit)
  p2 <- plot_nonlinear(fit, smooths = "s(Days)")
  expect_gg(p1)
  expect_gg(p2)
})

test_that("stan_gamm4 doesn't error when bs='cc", {
  # https://github.com/stan-dev/rstanarm/issues/362
  skip_if_not_installed("mgcv")
  N <- 100
  y <- rnorm(N, 0, 1)
  x <- rep(1:(N/2),2)
  x2 <- rnorm(N)
  data <- data.frame(x, x2, y)
  
  # only run a few iter to make sure it doesn't error
  SW(fit1 <- stan_gamm4(y ~ x2 + s(x, bs = "cc"), data=data, iter = 5, chains = 1, init = 0, refresh = 0))
  expect_stanreg(fit1)
  
  # with another smooth term
  SW(fit2 <- stan_gamm4(y ~ s(x2) + s(x, bs = "cc"), data=data, iter = 5, chains = 1, init = 0, refresh = 0))
  expect_stanreg(fit2)
  
  # with another 'cc' smooth term
  SW(fit3 <- stan_gamm4(y ~ s(x2, bs = "cc") + s(x, bs = "cc"), data=data, iter = 5, chains = 1, init = 0, refresh = 0))
  expect_stanreg(fit3)
})

test_that("stan_gamm4 errors if no smooth terms in formula", {
  dat <- data.frame(
    y = rnorm(100), 
    x = rnorm(100), 
    id = gl(5, 20)
  )
  expect_error(
    stan_gamm4(y ~ x, random = ~(1 | id), data = dat),
    "Formula must have at least one smooth term to use stan_gamm4"
  )
})

test_that("stan_gamm4 predictions with t2 terms match the design matrix", {
  # mgcv uses another basis for the prediction of t2 terms than for the fit
  skip_if_not_installed("mgcv")
  N <- 60
  i <- seq_len(N)
  dat <- data.frame(w = (i * 0.2847) %% 1, x = (i * 0.618) %% 1,
                    z = (i * 0.4142) %% 1, g = gl(6, 10))
  dat$y <- cos(2 * pi * dat$w) + sin(2 * pi * dat$x) + cos(3 * dat$z) +
    as.numeric(dat$g) / 3 + 0.3 * sin(7 * i)
  SW(fit1 <- stan_gamm4(y ~ s(w) + t2(x, z), data = dat, chains = 1,
                        iter = 100, seed = SEED, refresh = 0))
  SW(fit2 <- stan_gamm4(y ~ s(w) + t2(x, z), random = ~(1 | g), data = dat,
                        chains = 1, iter = 100, seed = SEED, refresh = 0))
  idx <- c(5:1, 31:35)
  fits <- list(gam = fit1, gamm = fit2)
  for (nm in names(fits)) {
    fit <- fits[[nm]]
    expect_false(is.null(fit$jam$P), info = nm)
    x_all <- as.matrix(fit$x)
    draws <- as.matrix(fit$stanfit)[, colnames(x_all)]
    eta <- tcrossprod(draws, x_all)
    expect_equal(posterior_linpred(fit), eta, tolerance = 1e-8,
                 ignore_attr = TRUE, info = nm)
    expect_equal(posterior_linpred(fit, newdata = dat[idx, ]), eta[, idx],
                 tolerance = 1e-8, ignore_attr = TRUE, info = nm)
    ll <- dnorm(matrix(dat$y, nrow(eta), N, byrow = TRUE), eta,
                as.matrix(fit)[, "sigma"], log = TRUE)
    expect_equal(log_lik(fit, newdata = dat), ll, tolerance = 1e-8,
                 ignore_attr = TRUE, info = nm)
    # the jam object has the posterior mean and covariance of the
    # population-level coefficients
    pop <- colnames(fit$jam$X)
    eta_pop <- tcrossprod(draws[, pop], x_all[idx, pop])
    pred_jam <- predict(fit$jam, newdata = dat[idx, ], se.fit = TRUE)
    expect_equal(pred_jam$fit, colMeans(eta_pop), tolerance = 1e-8,
                 ignore_attr = TRUE, info = nm)
    expect_equal(pred_jam$se.fit, apply(eta_pop, 2, sd), tolerance = 1e-8,
                 ignore_attr = TRUE, info = nm)
  }
  # log_lik without newdata uses fit$x (GAM only: see get_x.gamm4)
  expect_equal(log_lik(fit1, newdata = dat), log_lik(fit1), tolerance = 1e-8)

  # the prediction constraint of t2 centers the term over the data
  expect_lt(max(abs(colSums(predict(fit1$jam, type = "terms")))), 1e-8)

  # plot_nonlinear uses the basis of the prediction matrix, as plot.gam does
  P1 <- fit1$jam$P
  B1 <- as.matrix(fit1)[, colnames(fit1$jam$X)] %*% t(P1)
  X1 <- predict(fit1$jam, type = "lpmatrix")
  sw <- grep("s(w).", colnames(X1), fixed = TRUE)
  p1 <- plot_nonlinear(fit1, smooths = "s(w)")
  ord <- order(dat$w)
  expect_equal(p1$data$middle,
               apply(tcrossprod(B1[, sw], X1[ord, sw]), 2, median),
               tolerance = 1e-8, ignore_attr = TRUE)
  SW(fit3 <- stan_gamm4(y ~ t2(x, z), data = dat, chains = 1, iter = 100,
                        seed = SEED, refresh = 0))
  p3 <- plot_nonlinear(fit3)
  grid <- data.frame(x = p3$data$x, z = p3$data$y)
  X3 <- predict(fit3$jam, newdata = grid, type = "lpmatrix")
  t2c <- grep("t2(x,z).", colnames(X3), fixed = TRUE)
  B3 <- as.matrix(fit3)[, colnames(fit3$jam$X)] %*% t(fit3$jam$P)
  expect_equal(p3$data$z,
               apply(tcrossprod(B3[, t2c], X3[, t2c]), 2, median),
               tolerance = 1e-8, ignore_attr = TRUE)
  # the transformed draws give the linear predictor, and the t2 term of each
  # draw has mean 0 over the data
  expect_equal(tcrossprod(B3, X3), posterior_linpred(fit3, newdata = grid),
               tolerance = 1e-8, ignore_attr = TRUE)
  X3d <- predict(fit3$jam, type = "lpmatrix")
  expect_lt(max(abs(rowMeans(tcrossprod(B3[, t2c], X3d[, t2c])))), 1e-8)
})


test_that("loo/waic for stan_glmer works", {
  ll_fun <- rstanarm:::ll_fun
  
  # gaussian
  expect_equivalent_loo(fit)
  expect_identical(ll_fun(fit), rstanarm:::.ll_gaussian_i)
  
  # binomial
  expect_equivalent_loo(example_model)
  expect_identical(ll_fun(example_model), rstanarm:::.ll_binomial_i)
})

test_that("stan_gamm4 returns expected result for sleepstudy example", {
  skip_if_not_installed("mgcv")
  sleepstudy$y <- sleepstudy$Reaction / 10
  SW(fit <- stan_gamm4(y ~ s(Days), data = sleepstudy,
                       random = ~(1|Subject), chains = CHAINS, iter = ITER, 
                       seed = SEED, refresh = 0))
  expect_silent(yrep1 <- posterior_predict(fit))
  # expect_equal(dim(yrep1), c(nrow(as.data.frame(fit)), nobs(fit)))
  expect_silent(yrep2 <- posterior_predict(fit, draws = 1))
  # expect_equal(dim(yrep2), c(1, nobs(fit)))
  expect_silent(posterior_predict(fit, newdata = sleepstudy))
})

test_that("compatible with stan_lmer", {
  check_for_pp_errors(fit)
  expect_linpred_equal(fit)
})
test_that("compatible with stan_glmer (binomial)", {
  check_for_pp_errors(example_model)
  expect_linpred_equal(example_model)
  predprob <- posterior_linpred(example_model, transform = TRUE)
  expect_true(all(predprob > 0) && all(predprob < 1))
})
test_that("compatible with stan_(g)lmer with transformation in formula", {
  d <- mtcars
  d$cyl <- as.factor(d$cyl)
  args <- list(formula = mpg ~ log1p(wt) + (1|cyl) + (1|gear), data = d, 
               iter = 10, chains = 1,  seed = SEED, refresh = 0)
  SW(fit1 <- do.call("stan_lmer", args))
  SW(fit2 <- do.call("stan_glmer", args))
  nd <- d[6:10, ]
  nd$wt <- runif(5)
  expect_silent(posterior_predict(fit1))
  expect_silent(posterior_predict(fit2))
  expect_silent(posterior_predict(fit1, newdata = nd))
  expect_silent(posterior_predict(fit2, newdata = nd))
  
  expect_silent(posterior_linpred(fit1))
  expect_silent(posterior_linpred(fit2))
  expect_silent(posterior_linpred(fit1, newdata = nd))
  expect_silent(posterior_linpred(fit2, newdata = nd))
})

test_that("compatible with stan_lmer with offset", {
  offs <- rnorm(nrow(mtcars))
  SW(fit <- stan_lmer(mpg ~ wt + (1|cyl) + (1 + wt|gear), data = mtcars, 
                      prior = normal(0,1), iter = 10, chains = 1,
                      seed = SEED, refresh = 0, offset = offs))
  
  expect_warning(posterior_predict(fit, newdata = mtcars[1:2, ], offset = offs),
                 "STATS")
  check_for_pp_errors(fit, offset = offs)
})

test_that("predition with family mgcv::betar doesn't error", {
  test_data <- data.frame(y = c(0.1, 0.3), x = c(TRUE, FALSE))
  SW(fit <- stan_glmer(y ~ (1|x), family=mgcv::betar(link="logit"), 
                       data=test_data,  seed = 101, iter = 10, chains = 1, refresh = 0))
  expect_silent(posterior_linpred(fit, newdata=test_data))
  expect_silent(posterior_predict(fit, newdata=test_data))
})

# compare to lme4 ---------------------------------------------------------

test_that("posterior_predict close to predict.merMod for gaussian", {
  mod1 <- as.formula(mpg ~ wt + (1|cyl) + (1|gear))
  mod2 <- as.formula(mpg ~ log1p(wt) + I(disp/100) + (1|cyl))
  mod3 <- as.formula(mpg ~ wt + (1|cyl) + (1 + wt|gear))
  mod4 <- as.formula(log(mpg) ~ wt + (1 + wt|cyl) + (1 + wt + am|gear))
  
  lfit1 <- lmer(mod1, data = mtcars)
  SW(sfit1 <- stan_glmer(mod1, data = mtcars, iter = 400,
                      chains = CHAINS, seed = SEED, refresh = 0))
  lfit2 <- update(lfit1, formula = mod2)
  SW(sfit2 <- update(sfit1, formula = mod2))
  lfit3 <- update(lfit1, formula = mod3)
  SW(sfit3 <- update(sfit1, formula = mod3))
  lfit4 <- update(lfit1, formula = mod4)
  SW(sfit4 <- update(sfit1, formula = mod4))
  
  nd <- nd2 <- mtcars[1:5, ]
  nd2$cyl[2] <- 5 # add new levels
  nd3 <- nd2
  nd3$gear[2] <- 7
  nd3$gear[5] <- 1
  
  tol <- 0.3
  for (j in 1:4) {
    expect_equal(
      colMeans(posterior_predict(get(paste0("sfit", j)), newdata = nd, seed = SEED)),
      unname(predict(get(paste0("lfit", j)), newdata = nd)),
      tolerance = tol, ignore_attr = TRUE)
    expect_equal(
      colMeans(posterior_predict(get(paste0("sfit", j)), newdata = nd2, seed = SEED,
                                 allow.new.levels = TRUE)),
      unname(predict(get(paste0("lfit", j)), newdata = nd2, allow.new.levels = TRUE)),
      tolerance = tol, ignore_attr = TRUE)
    expect_equal(
      colMeans(posterior_predict(get(paste0("sfit", j)), newdata = nd3, seed = SEED,
                                 allow.new.levels = TRUE)),
      unname(predict(get(paste0("lfit", j)), newdata = nd3, allow.new.levels = TRUE)),
      tolerance = tol, ignore_attr = TRUE)
  }
})

test_that("posterior_predict close to predict.merMod for binomial", {
  d <- nd <- lme4::cbpp
  sfit <- example_model
  lfit <- glmer(formula(example_model), data = d, family = "binomial")
  levels(nd$herd) <- c(levels(nd$herd), "99")
  nd$herd[1:2] <- "99"
  lpred <- simulate(lfit, newdata = nd, re.form = NULL, allow.new.levels = TRUE,
                    nsim = 500, seed = SEED)
  for (j in 1:ncol(lpred)) {
    lpred[, j] <- lpred[, j][, 1] / rowSums(lpred[, j])
  }
  lpred <- t(as.matrix(lpred))
  spred <- posterior_predict(sfit, draws = 500, newdata = nd, 
                             seed = SEED)
  spred <- sweep(spred, 2, rowSums(get_y(sfit)), "/")
  expect_equal(colMeans(spred), unname(colMeans(lpred)),
               tolerance = .125, ignore_attr = TRUE)
})
