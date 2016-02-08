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
SEED <- 12345
set.seed(SEED)
ITER <- 10L
CHAINS <- 2L
REFRESH <- 0

SW <- suppressWarnings

context("helper functions")

test_that("nlist works", {
  nlist <- rstanarm:::nlist
  a <- 1
  b <- 2
  c <- 3
  val <- list(nlist(a, b, c), 
              nlist(a, b, c = "tornado"), 
              nlist(a = -1, b = -2, c))
  ans <- list(list(a = a, b = b, c = c), 
              list(a = a, b = b, c = "tornado"), 
              list(a = -1, b = -2, c = c))
  expect_identical(val, ans)
})

test_that("family checking works", {
  fams <- rstanarm:::nlist("binomial", "gaussian", "poisson", gamma = "Gamma", 
                           ig = "inverse.gaussian", nb = "neg_binomial_2")
  for (j in seq_along(fams)) {
    is.f <- getFromNamespace(paste0("is.", names(fams)[j]), "rstanarm")
    f <- get(fams[[j]])()$family
    expect_true(is.f(f))
    expect_false(is.f("not a family"))
  }
})

test_that("%ORifNULL% works", {
  `%ORifNULL%` <- rstanarm:::`%ORifNULL%`
  a <- list(NULL, NA, NaN, 1, "a", FALSE)
  b <- 1
  ans <- c(b, a[-1])
  for (j in seq_along(a)) {
    expect_identical(a[[j]] %ORifNULL% b, ans[[j]])
  }
})

test_that("%ORifINF% works", {
  `%ORifINF%` <- rstanarm:::`%ORifINF%`
  a <- list(Inf, -Inf, 1, "a", FALSE)
  b <- 0
  ans <- c(b, a[-1])
  for (j in seq_along(a)) {
    expect_identical(a[[j]] %ORifINF% b, ans[[j]])
  }
})

test_that("maybe_broadcast works", {
  maybe_broadcast <- rstanarm:::maybe_broadcast
  n <- 5
  x <- list(numeric(0), NULL, 1, c(1,1))
  ans <- list(rep(0,n), rep(0,n), rep(1,n), c(1,1))
  for (j in seq_along(ans)) {
    expect_equal(maybe_broadcast(x[[j]], n), ans[[j]])  
  }
})

test_that("set_prior_scale works", {
  set_prior_scale <- rstanarm:::set_prior_scale
  expect_error(set_prior_scale("a", "b", "c"))
  expect_error(set_prior_scale(1, 1, 1))
  expect_equal(set_prior_scale(NULL, 1, "a"), 1)
  expect_equal(set_prior_scale(NULL, 1, "probit"), dnorm(0) / dlogis(0))
  expect_equal(set_prior_scale(2, 1, "a"), 2)
  expect_equal(set_prior_scale(2, 1, "probit"), 2 * dnorm(0) / dlogis(0))
})

test_that("validate_parameter_value works", {
  validate_parameter_value <- rstanarm:::validate_parameter_value
  expect_error(validate_parameter_value(-1), "should be positive")
  expect_error(validate_parameter_value(0), "should be positive")
  expect_error(validate_parameter_value("a"), "should be NULL or numeric")
  expect_error(validate_parameter_value(NA), "should be NULL or numeric")
  expect_true(validate_parameter_value(NULL))
  expect_true(validate_parameter_value(.01))
  expect_true(validate_parameter_value(.Machine$double.xmax))
})

test_that("validate_weights works", {
  validate_weights <- rstanarm:::validate_weights
  ff <- function(weights) validate_weights(weights)
  expect_equal(ff(), double(0))
  expect_equal(ff(x <- rexp(10)), x)
  expect_equal(validate_weights(NULL), double(0))
  expect_equal(validate_weights(1:10), 1:10)
  expect_error(validate_weights(LETTERS), regexp = "numeric")
  expect_error(validate_weights(c(-1,2,3)), regexp = "negative", ignore.case = TRUE)
  expect_error(stan_glm(mpg ~ wt, data = mtcars, weights = rep(-1, nrow(mtcars))), 
               regexp = "negative", ignore.case = TRUE)
  expect_is(stan_glm(mpg ~ wt, data = mtcars, algorithm = "optimizing", seed = SEED,
                     weights = rexp(nrow(mtcars))), "stanreg")
})

test_that("validate_offset works", {
  validate_offset <- rstanarm:::validate_offset
  expect_equal(validate_offset(NULL), double(0))
  expect_equal(validate_offset(rep(1, 10), rnorm(10)), rep(1, 10))
  expect_error(validate_offset(rep(1, 10), rnorm(5)))
  expect_error(validate_offset(rep(1, 5), rnorm(10)), 
               regexp = "number of offsets", ignore.case = TRUE)
  fito <- stan_glm(mpg ~ wt, data = mtcars, algorithm = "optimizing", seed = SEED)
  fito2 <- update(fito, offset = rep(5, nrow(mtcars)))
  expect_equal(coef(fito)[1], 5 + coef(fito2)[1], tol = 0.2)
})

test_that("validate_family works", {
  validate_family <- rstanarm:::validate_family
  expect_equal(validate_family("gaussian"), gaussian())
  expect_equal(validate_family(gaussian), gaussian())
  expect_equal(validate_family(gaussian()), gaussian())
  expect_equal(validate_family(gaussian(link = "log")), gaussian(link = "log"))
  expect_equal(validate_family(binomial(link = "probit")), binomial(link = "probit"))
  expect_equal(validate_family(neg_binomial_2()), neg_binomial_2())
  expect_error(validate_family("not a family"))
  expect_error(validate_family(rnorm(10)), "must be a family")
  expect_error(stan_glm(mpg ~ wt, data = mtcars, family = "not a family"))
})

test_that("array1D_check works", {
  array1D_check <- rstanarm:::array1D_check
  y1 <- rnorm(10)
  expect_equal(array1D_check(y1), y1)
  names(y1) <- rep_len(letters, length(y1))
  expect_equal(array1D_check(y1), y1)
  expect_identical(array1D_check(as.array(y1)), y1)
  
  y2 <- cbind(1:10, 11:20)
  expect_equal(array1D_check(y2), y2)
})

test_that("fac2bin works", {
  fac2bin <- rstanarm:::fac2bin
  y <- gl(2, 2, 20, labels = c("lo", "hi"))
  expect_identical(fac2bin(y), rep_len(c(0L, 0L, 1L, 1L), 20))
  y <- gl(2, 8, labels = c("Control", "Treat"))
  expect_identical(fac2bin(y), rep(c(0L, 1L), each = 8))
  expect_identical(fac2bin(factor(c(1,2))), c(0L, 1L))
  expect_error(fac2bin(rnorm(10)))
  expect_error(fac2bin(factor(c(1,2,3))))
  expect_error(fac2bin(factor(mtcars$cyl, labels = c("lo", "mid", "hi"))))
})

test_that("check_constant_vars works", {
  check_constant_vars <- rstanarm:::check_constant_vars
  mf <- model.frame(glm(mpg ~ ., data = mtcars))
  mf2 <- mf
  mf2$wt <- 2
  expect_equal(check_constant_vars(mf), mf)
  expect_error(check_constant_vars(mf2), "wt")
  mf2$gear <- 1
  expect_error(check_constant_vars(mf2), "wt, gear")
  expect_error(stan_glm(mpg ~ ., data = mf2), "wt, gear")
  expect_is(stan_glm(mpg ~ ., data = mf, algorithm = "optimizing", seed = SEED), 
            "stanreg")
  expect_is(stan_glm(mpg ~ ., data = mf, weights = rep(2, nrow(mf)),
                     offset = rep(1, nrow(mf)), algorithm = "optimizing", 
                     seed = SEED), "stanreg")
  
  esoph2 <- esoph
  esoph2$agegp[1:nrow(esoph2)] <- "75+"
  expect_error(stan_polr(tobgp ~ agegp, data = esoph2, iter = 10,
                         prior = R2(0.2, "mean"), init_r = 0.1, seed = SEED, 
                         refresh = 10), 
               regexp = "agegp")
})

test_that("linear_predictor methods work", {
  linpred_vec <- rstanarm:::linear_predictor.default
  linpred_mat <- rstanarm:::linear_predictor.matrix
  x <- cbind(1, 1:4)
  bmat <- matrix(c(-0.5, 0, 0.5, 1), nrow = 2, ncol = 2)
  bvec <- bmat[1, ]
  vec_ans <- seq(0, 1.5, 0.5)
  mat_ans <- rbind(vec_ans, 1:4)
  offset <- rep(2, nrow(x))
  expect_equivalent(linpred_vec(bvec, x), vec_ans)
  expect_equivalent(linpred_vec(bvec, x, offset = NULL), vec_ans)
  expect_equivalent(linpred_vec(bvec, x, offset), vec_ans + offset)
  expect_equivalent(linpred_mat(bmat, x), mat_ans)
  expect_equivalent(linpred_mat(bmat, x, offset = NULL), mat_ans)
  expect_equivalent(linpred_mat(bmat, x, offset), mat_ans + offset)
})

# fits to use in multiple calls to test_that below
fit <- SW(stan_glm(mpg ~ wt, data = mtcars, iter = ITER, 
                   chains = CHAINS, seed = SEED, refresh = REFRESH))
fit2 <- SW(stan_glmer(mpg ~ wt + (1|cyl), data = mtcars, 
                      iter = ITER, chains = CHAINS, seed = SEED, refresh = REFRESH))
fito <- stan_glm(mpg ~ wt, data = mtcars, algorithm = "optimizing", seed = SEED)
fitvb <- update(fito, algorithm = "meanfield", seed = SEED)
fitvb2 <- update(fitvb, algorithm = "fullrank", seed = SEED)

test_that("is.stanreg works", {
  is.stanreg <- rstanarm:::is.stanreg
  expect_true(is.stanreg(fit))
  expect_true(is.stanreg(fit2))
  expect_true(is.stanreg(fito))
  expect_true(is.stanreg(fitvb))
  expect_false(is.stanreg(fit$stanfit))
})

test_that("used.sampling, used.optimizing, and used.variational work", {
  used.sampling <- rstanarm:::used.sampling
  used.optimizing <- rstanarm:::used.optimizing
  used.variational <- rstanarm:::used.variational
  expect_true(used.sampling(fit))
  expect_true(used.sampling(fit2))
  expect_false(used.optimizing(fit))
  expect_false(used.optimizing(fit2))
  expect_false(used.variational(fit))
  expect_false(used.variational(fit2))
  
  expect_true(used.optimizing(fito))
  expect_false(used.sampling(fito))
  expect_false(used.variational(fito))
  
  expect_true(used.variational(fitvb))
  expect_true(used.variational(fitvb2))
  expect_false(used.sampling(fitvb))
  expect_false(used.sampling(fitvb2))
  expect_false(used.optimizing(fitvb))
  expect_false(used.optimizing(fitvb2))
  
  # should return error if passed anything but a stanreg object
  expect_error(used.sampling(fit$stanfit))
  expect_error(used.variational(fitvb$stanfit))
  expect_error(used.optimizing(fito$stanfit))
})

test_that("is.mer works", {
  is.mer <- rstanarm:::is.mer
  bad1 <- bad2 <- example_model
  bad1$glmod <- NULL
  class(bad2) <- "stanreg"
  
  expect_true(is.mer(example_model))
  expect_true(is.mer(fit2))
  expect_false(is.mer(fit))
  expect_false(is.mer(fito))
  expect_false(is.mer(fitvb))
  expect_false(is.mer(fitvb2))
  expect_error(is.mer(bad1), regexp = "Bug found")
  expect_error(is.mer(bad2), regexp = "Bug found")
})

test_that("get_x, get_y, get_z work", {
  x_ans <- cbind("(Intercept)" = 1, wt = mtcars$wt)
  y_ans <- mtcars$mpg
  expect_equivalent(get_x(fit), x_ans)
  expect_equivalent(get_y(fit), y_ans)
  expect_error(get_z(fit), "no applicable method")
  
  z_ans2 <- model.matrix(mpg ~ -1 + factor(cyl), data = mtcars)
  expect_equivalent(get_x(fit2), x_ans)
  expect_equivalent(get_y(fit2), y_ans)
  expect_equivalent(as.matrix(get_z(fit2)), z_ans2)
  
  fit3 <- SW(stan_glmer(mpg ~ wt + (1 + wt|cyl), data = mtcars, 
                        iter = 10, chains = 1, refresh = 5, seed = SEED))
  z_ans3 <- mat.or.vec(nr = nrow(mtcars), nc = 6)
  z_ans3[, c(1, 3, 5)] <- model.matrix(mpg ~ 0 + factor(cyl), data = mtcars)
  z_ans3[, c(2, 4, 6)] <- model.matrix(mpg ~ 0 + wt:factor(cyl), data = mtcars)
  expect_equivalent(get_x(fit3), x_ans)
  expect_equivalent(get_y(fit3), y_ans)
  expect_equivalent(as.matrix(get_z(fit3)), z_ans3)
})

test_that("set_sampling_args works", {
  set_sampling_args <- rstanarm:::set_sampling_args
  
  # user specifies stepsize and also overrides default max_treedepth
  control1 <- list(max_treedepth = 10, stepsize = 0.01)
  # user specifies control but doesn't override max_treedepth
  control2 <- list(stepsize = 0.01)
  # no user 'control' argument 
  no_control <- list()
  
  # normal prior --> adapt_delta = 0.95
  val1 <- set_sampling_args(fit, prior = normal(),
                            user_dots = list(control = control1, iter = 100),  
                            user_adapt_delta = NULL)
  # use fit2 instead of fit to check that it doesn't matter which fit object is used
  val1b <- set_sampling_args(fit2,
                             prior = normal(),
                             user_dots = list(control = control1, iter = 100),  
                             user_adapt_delta = NULL)
  # normal prior --> adapt_delta = 0.95, but user override to 0.9
  val2 <- set_sampling_args(fit, prior = normal(), 
                            user_dots = list(control = control1),  
                            user_adapt_delta = 0.9)
  # cauchy/t_1 prior --> adapt_delta = 0.99
  val3 <- set_sampling_args(fit, prior = student_t(1), 
                            user_dots = list(control = control1),  
                            user_adapt_delta = NULL)
  # cauchy/t_1 prior --> adapt_delta = 0.99, but user override to 0.8
  val4 <- set_sampling_args(fit, prior = cauchy(),
                            user_dots = list(control = control2),  
                            user_adapt_delta = 0.8)
  # hs prior --> adapt_delta = 0.99
  val5 <- set_sampling_args(fit, prior = hs(), 
                            user_dots = no_control,
                            user_adapt_delta = NULL)
  val6 <- set_sampling_args(fit, prior = hs_plus(), 
                            user_dots = no_control,
                            user_adapt_delta = NULL)
  expect_equal(val1$control, c(control1, adapt_delta = 0.95))
  expect_equal(val1$iter, 100)
  expect_equal(val1$control, val1b$control)
  expect_equal(val2$control, c(control1, adapt_delta = 0.9))
  expect_equal(val3$control, c(control1, adapt_delta = 0.99))
  expect_equal(val4$control, c(control2, adapt_delta = 0.8, max_treedepth = 15))
  expect_equal(val5$control, list(adapt_delta = 0.99, max_treedepth = 15))
  expect_equal(val6$control, list(adapt_delta = 0.99, max_treedepth = 15))
})

test_that("linkinv methods work", {
  linkinv.stanreg <- rstanarm:::linkinv.stanreg
  linkinv.character <- rstanarm:::linkinv.character
  linkinv.family <- rstanarm:::linkinv.family
  
  expect_identical(linkinv.family(gaussian()), gaussian()$linkinv)
  expect_identical(linkinv.family(neg_binomial_2()), neg_binomial_2()$linkinv)
  expect_identical(linkinv.family(binomial(link = "probit")), 
                   binomial(link = "probit")$linkinv)
  
  fit_polr <- SW(stan_polr(tobgp ~ agegp, data = esoph, method = "loglog",
                           prior = R2(0.2, "mean"), init_r = 0.1, 
                           chains = CHAINS, iter = ITER, seed = SEED, 
                           refresh = REFRESH))
  expect_identical(linkinv.stanreg(fit_polr), rstanarm:::pgumbel)
  expect_identical(linkinv.character(fit_polr$family), rstanarm:::pgumbel)
  expect_identical(linkinv.stanreg(example_model), binomial()$linkinv)
  expect_identical(linkinv.stanreg(fit), gaussian()$linkinv)
  expect_error(rstanarm:::polr_linkinv(example_model), 
               regexp = "should be a stanreg object created by stan_polr")
})

test_that("collect_pars and grep_for_pars work", {
  fit <- example_model
  collect_pars <- rstanarm:::collect_pars
  grep_for_pars <- rstanarm:::grep_for_pars
  
  all_period <- paste0("period", 2:4)
  all_varying <- rstanarm:::b_names(rownames(fit$stan_summary), value = TRUE)
  
  expect_identical(grep_for_pars(fit, "period"), all_period)
  expect_identical(grep_for_pars(fit, c("period", "size")), c(all_period, "size"))
  expect_identical(grep_for_pars(fit, "period|size"), c("size", all_period))
  expect_identical(grep_for_pars(fit, "(2|3)$"), all_period[1:2])
  expect_identical(grep_for_pars(fit, "herd"), all_varying)
  expect_identical(grep_for_pars(fit, "b\\["), all_varying)
  expect_identical(grep_for_pars(fit, "Intercept"), c("(Intercept)", all_varying))
  expect_identical(grep_for_pars(fit, "herd:[3,5]"), all_varying[c(3,5)])
  expect_identical(grep_for_pars(fit, "herd:[3-5]"), all_varying[3:5])
  expect_error(grep_for_pars(fit, "NOT A PARAMETER"), regexp = "No matches")
  expect_error(grep_for_pars(fit, "b["))
  
  expect_identical(collect_pars(fit, regex_pars = "period"), all_period)
  expect_identical(collect_pars(fit, pars = "size", regex_pars = "period"), 
                   c("size", all_period))
  expect_identical(collect_pars(fit, pars = c("(Intercept)", "size")), 
                   c("(Intercept)", "size"))
  expect_identical(collect_pars(fit, pars = "period2", regex_pars = "herd:[[1]]"), 
                   c("period2", all_varying[1]))
  expect_identical(collect_pars(fit, pars = "size", regex_pars = "size"), "size")
  expect_identical(collect_pars(fit, regex_pars = c("period", "herd")), 
                   c(all_period, all_varying))
})

test_that("posterior_sample_size works", {
  pss <- rstanarm:::posterior_sample_size
  expect_equal(pss(example_model), 500)
  expect_equal(pss(fit), nrow(as.matrix(fit)))
  expect_equal(pss(fit2), ITER * CHAINS / 2)
  expect_equal(pss(fitvb), 1000)
  expect_equal(pss(fitvb2), 1000)
  expect_null(pss(fito))
  
  fit3 <- suppressWarnings(stan_glm(mpg ~ wt, data = mtcars, iter = 20, 
                                    chains = 1, thin = 2))
  expect_equal(pss(fit3), nrow(as.matrix(fit3)))
})
