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

library(rstanarm)
library(lme4)
library(MASS)
SEED <- 12345
set.seed(SEED)
ITER <- 10
CHAINS <- 2
REFRESH <- 0

SW <- suppressWarnings

N <- 200
x <- rnorm(N, 2, 1)
z <- rnorm(N, 2, 1)
mu <- binomial(link = "logit")$linkinv(1 + 0.2*x)
phi <- exp(1.5 + 0.4*z)
y <- rbeta(N, mu * phi, (1 - mu) * phi)
fake_dat <- data.frame(y, x, z)
remove(N, x, y, z, mu, phi)

capture.output(
  stan_glm1 <- SW(stan_glm(mpg ~ wt + cyl, data = mtcars, iter = ITER,
                           chains = CHAINS, seed = SEED, refresh = 0)),
  stan_glm_opt1 <- stan_glm(mpg ~ wt + cyl, data = mtcars, algorithm = "optimizing",
                            seed = SEED, refresh = 0),
  stan_glm_vb1 <- update(stan_glm_opt1, algorithm = "meanfield", QR = TRUE, iter = 10000),
  glm1 <- glm(mpg ~ wt + cyl, data = mtcars),
  
  lmer1 <- lmer(diameter ~ (1|plate) + (1|sample), data = Penicillin),
  stan_lmer1 <- SW(stan_lmer(diameter ~ (1|plate) + (1|sample), data = Penicillin,
                             prior_intercept = normal(0, 50, autoscale = FALSE),
                             prior_aux = normal(0, 10),
                             iter = ITER, chains = CHAINS, seed = SEED, refresh = 0)),
  lmer2 <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy),
  stan_lmer2 <- SW(stan_lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy,
                             iter = ITER, chains = CHAINS, seed = SEED,
                             refresh = 0)),
  
  stan_polr1 <- SW(stan_polr(tobgp ~ agegp, data = esoph, prior = R2(0.2, "mean"),
                             init_r = 0.1, iter = ITER, chains = CHAINS,
                             seed = SEED, refresh = 0)),
  polr1 <- polr(tobgp ~ agegp, data = esoph, Hess = TRUE),
  
  stan_gamm41 <- SW(stan_gamm4(mpg ~ s(wt) + cyl, data = mtcars, iter = ITER,
                               chains = CHAINS, seed = SEED, refresh = 0)),

  stan_betareg1 <- SW(stan_betareg(y ~ x | z, data = fake_dat, 
                                   link = "logit", link.phi = "log", refresh = 0,
                                   iter = ITER, chains = CHAINS, seed = SEED)),
  betareg1 <- betareg::betareg(y ~ x | z, data = fake_dat, link = "logit", link.phi = "log")
)

att_names <- function(object) {
  nms <- names(object)
  att_nms <- names(attributes(object))
  att_nms2 <- lapply(object, function(x) sort(names(attributes(x))))
  c(nms, att_nms, att_nms2)
}
check_att_names <- function(x,y) {
  expect_identical(att_names(x), att_names(y))
}
check_sizes <- function(x,y) {
  expect_equal(length(x), length(y))
  expect_equal(lapply(x, dim), lapply(y, dim))
}


context("methods for stanreg objects")


# extractors --------------------------------------------------------------
test_that("stanreg extractor methods work properly", {
  expect_equal(resid(stan_glm1), stan_glm1$residuals)
  expect_equal(coef(stan_glm1), stan_glm1$coefficients)
  expect_equal(vcov(stan_glm1), stan_glm1$covmat)
  expect_equal(fitted(stan_glm1), stan_glm1$fitted.values)
  expect_equal(se(stan_glm1), stan_glm1$ses)

  expect_equal(resid(stan_polr1), stan_polr1$residuals)
  expect_equal(coef(stan_polr1), stan_polr1$coefficients)
  expect_equal(vcov(stan_polr1), stan_polr1$covmat)
  expect_equal(fitted(stan_polr1), stan_polr1$fitted.values)
  expect_equal(se(stan_polr1), stan_polr1$ses)

  expect_equal(vcov(stan_glm_opt1), stan_glm_opt1$covmat)
  expect_equal(vcov(stan_glm_opt1, correlation = TRUE),
               cov2cor(stan_glm_opt1$covmat))
  expect_equal(resid(stan_glm_opt1), stan_glm_opt1$residuals)
  expect_equal(coef(stan_glm_opt1), stan_glm_opt1$coefficients)
  expect_equal(fitted(stan_glm_opt1), stan_glm_opt1$fitted.values)
  expect_equal(se(stan_glm_opt1), stan_glm_opt1$ses)

  expect_equal(resid(stan_lmer1), stan_lmer1$residuals)
  expect_equal(fitted(stan_lmer1), stan_lmer1$fitted.values)
  expect_equal(se(stan_lmer1), stan_lmer1$ses)
  expect_equal(resid(example_model), example_model$residuals)
  expect_equal(fitted(example_model), example_model$fitted.values)
  expect_equal(se(example_model), example_model$ses)
  # coef and vcov are different for stan_(g)lmer models and are tested
  # separately later in this file
  
  expect_equal(resid(stan_betareg1), stan_betareg1$residuals)
  expect_equal(coef(stan_betareg1), stan_betareg1$coefficients)
  expect_equal(vcov(stan_betareg1), stan_betareg1$covmat)
  expect_equal(fitted(stan_betareg1), stan_betareg1$fitted.values)
  expect_equal(se(stan_betareg1), stan_betareg1$ses)
  
})


# confint -----------------------------------------------------------------
test_that("confint method returns correct structure", {
  expect_silent(ci <- confint(stan_glm_opt1))
  expect_silent(ci2 <- confint(stan_glm_opt1, parm = "wt", level = 0.9))
  expect_equal(rownames(ci), c("(Intercept)", "wt", "cyl"))
  expect_equal(colnames(ci), c("2.5 %", "97.5 %"))
  expect_equal(rownames(ci2), c("wt"))
  expect_equal(colnames(ci2), c("5 %", "95 %"))

  expect_error(confint(stan_glm1), regexp = "use posterior_interval")
  expect_error(confint(stan_glm_vb1), regexp = "use posterior_interval")
  expect_error(confint(stan_polr1), regexp = "use posterior_interval")
  expect_error(confint(stan_lmer1), regexp = "use posterior_interval")
  expect_error(confint(stan_lmer2), regexp = "use posterior_interval")
  expect_error(confint(stan_betareg1), regexp = "use posterior_interval")
})


# posterior_interval -----------------------------------------------------
test_that("posterior_interval returns correct structure", {
  expect_silent(ci <- posterior_interval(stan_glm1, prob = 0.5))
  expect_silent(ci2 <- posterior_interval(stan_glm_vb1, pars = "wt", prob = 0.95))
  expect_silent(ci3 <- posterior_interval(example_model, prob = 0.95, regex_pars = "herd"))
  expect_silent(ci4 <- posterior_interval(example_model, prob = 0.8, pars = "(Intercept)",
                               regex_pars = "period"))
  expect_silent(ci5 <- posterior_interval(stan_polr1, prob = 0.9))
  
  expect_identical(rownames(ci), c("(Intercept)", "wt", "cyl", "sigma"))
  expect_identical(rownames(ci2), "wt")
  expect_identical(rownames(ci3), c(paste0("b[(Intercept) herd:", 1:15, "]"), 
                                    "Sigma[herd:(Intercept),(Intercept)]"))
  expect_identical(rownames(ci4), c("(Intercept)", paste0("period", 2:4)))
  expect_identical(colnames(ci), c("25%", "75%"))
  expect_identical(colnames(ci2), c("2.5%", "97.5%"))
  expect_identical(colnames(ci3), c("2.5%", "97.5%"))
  expect_identical(colnames(ci4), c("10%", "90%"))
  expect_identical(colnames(ci5), c("5%", "95%"))

  expect_silent(ci6 <- posterior_interval(stan_betareg1, prob = 0.5))
  expect_identical(colnames(ci6), c("25%", "75%"))
  
  expect_error(posterior_interval(stan_glm1, type = "HPD"),
               regexp = "only option for 'type' is 'central'")
  expect_identical(colnames(posterior_interval(stan_glm_opt1)), c("5%", "95%"))
  expect_error(posterior_interval(lm(mpg ~ wt, data = mtcars)),
               regexp = "should be a matrix")

  prob_msg <- "'prob' should be a single number greater than 0 and less than 1."
  expect_error(posterior_interval(stan_glm1, prob = c(0.25, 0.75)), regexp = prob_msg)
  expect_error(posterior_interval(stan_glm1, prob = 0), regexp = prob_msg)
  expect_error(posterior_interval(stan_glm1, prob = 1), regexp = prob_msg)
  expect_error(posterior_interval(stan_glm1, prob = 2), regexp = prob_msg)
})



# log_lik -----------------------------------------------------------------
test_that("log_lik method works", {
  expect_silent(log_lik(stan_glm_opt1))
  expect_silent(log_lik(stan_glm_vb1))
  expect_silent(log_lik(stan_glm1))

  expect_silent(log_lik(stan_polr1))
  expect_silent(log_lik(stan_gamm41))
  expect_equal(dim(log_lik(stan_polr1)), c(ITER, nobs(stan_polr1)))
  expect_equal(dim(log_lik(stan_lmer1)), c(ITER, nobs(stan_lmer1)))
  expect_equal(log_lik(stan_betareg1), log_lik(stan_betareg1, newdata = fake_dat))
  
  # Compute log-lik matrix using different method than log_lik.stanreg
  # and compare
  samp <- as.matrix(stan_glm1)
  y <- get_y(stan_glm1)
  y_new <- y[1:10] + rnorm(10)
  x <- get_x(stan_glm1)
  x_new <- cbind(1, x[1:10, 2:3] + rnorm(10))
  sigma <- samp[, 4]
  eta <- tcrossprod(x, samp[, 1:3])
  eta_new <- tcrossprod(x_new, samp[, 1:3])
  llmat <- matrix(NA, nrow = nrow(samp), ncol = nrow(eta))
  llmat_new <- matrix(NA, nrow = nrow(samp), ncol = nrow(eta_new))
  for (i in 1:nrow(llmat)) {
    llmat[i, ] <- dnorm(y, mean = eta[, i], sd = sigma[i], log = TRUE)
    llmat_new[i, ] <- dnorm(y_new, mean = eta_new[, i], sd = sigma[i], log = TRUE)
  }
  expect_equal(log_lik(stan_glm1), llmat, check.attributes = FALSE)
  nd <- data.frame(mpg = y_new, wt = x_new[, 2], cyl = x_new[, 3])
  expect_equal(log_lik(stan_glm1, newdata = nd), llmat_new, check.attributes = FALSE)


  # make sure log_lik with newdata equals log_lik if newdata is the same as the
  # data used to fit the model
  expect_equal(log_lik(example_model), log_lik(example_model, newdata = cbpp))
  expect_equal(log_lik(stan_lmer2), log_lik(stan_lmer2, newdata = sleepstudy))
  expect_equal(log_lik(stan_glm1), log_lik(stan_glm1, newdata = mtcars))
  expect_equal(log_lik(stan_polr1), log_lik(stan_polr1, newdata = esoph))
  expect_equal(log_lik(stan_gamm41), log_lik(stan_gamm41, newdata = mtcars))
})


# ngrps, nobs -------------------------------------------------------------
test_that("ngrps is right", {
  expect_equal(ngrps(lmer1), ngrps(stan_lmer1))
  expect_equal(ngrps(lmer2), ngrps(stan_lmer2))
  expect_error(ngrps(stan_glm1), "stan_glmer and stan_lmer models only")
  expect_error(ngrps(stan_betareg1), "stan_glmer and stan_lmer models only")
  expect_equal(nobs(stan_betareg1), nobs(betareg1))
})

test_that("nobs is right", {
  expect_equal(nobs(lmer1), nobs(stan_lmer1))
  expect_equal(nobs(lmer2), nobs(stan_lmer2))
  expect_equal(nobs(glm1), nobs(stan_glm_opt1))
  expect_equal(nobs(glm1), nobs(stan_glm1))
  expect_equal(nobs(polr1), nobs(stan_polr1))
})


# vcov --------------------------------------------------------------
test_that("vcov returns correct structure", {
  expect_equal(dimnames(vcov(stan_glm1)), dimnames(vcov(glm1)))
  expect_equal(dimnames(vcov(stan_polr1)), dimnames(vcov(polr1)))
  expect_equal(dimnames(vcov(stan_lmer1)), dimnames(vcov(lmer1)))
  expect_equal(dimnames(vcov(stan_lmer2)), dimnames(vcov(lmer2)))
  expect_equal(dimnames(vcov(stan_betareg1)), dimnames(vcov(betareg1)))
})

# sigma --------------------------------------------------------------
test_that("sigma method works", {
  # need to use :: because sigma is masked by lme4's sigma
  rsigma <- rstanarm::sigma
  expect_identical(rsigma(stan_polr1), 1)
  expect_identical(rsigma(example_model), 1)

  expect_double <- function(x) expect_type(x, "double")

  expect_double(sig <- rsigma(stan_lmer1))
  expect_false(identical(sig, 1))
  expect_double(sig <- rsigma(stan_lmer2))
  expect_false(identical(sig, 1))
  expect_double(sig <- rsigma(stan_glm1))
  expect_false(identical(sig, 1))
  expect_double(sig <- rsigma(stan_glm_vb1))
  expect_false(identical(sig, 1))
  expect_double(sig <- rsigma(stan_glm_opt1))
  expect_false(identical(sig, 1))
  
  expect_double(sig <- rsigma(stan_betareg1))
  expect_true(identical(sig, 1))
})


# VarCorr -----------------------------------------------------------------
test_that("VarCorr returns correct structure", {
  vc_lmer1 <- VarCorr(lmer1); vc_stan1 <- VarCorr(stan_lmer1)
  vc_lmer2 <- VarCorr(lmer2); vc_stan2 <- VarCorr(stan_lmer2)
  expect_s3_class(vc_stan1, class(vc_lmer1))
  expect_s3_class(vc_stan2, class(vc_lmer2))
  check_att_names(vc_stan1, vc_lmer1)
  check_att_names(vc_stan2, vc_lmer2)
  v <- sapply(vc_stan1, "[[", 1)
  expect_true(length(unique(v)) == length(v))
  expect_error(VarCorr(stan_glm1), "stan_glmer and stan_lmer models only")
  expect_error(VarCorr(stan_betareg1), "stan_glmer and stan_lmer models only")
})


# ranef,fixef,coef -----------------------------------------------------------
test_that("ranef returns correct structure", {
  re_stan1 <- ranef(stan_lmer1); re_lmer1 <- ranef(lmer1)
  re_stan2 <- ranef(stan_lmer1); re_lmer2 <- ranef(lmer1)
  expect_s3_class(re_stan1, class(re_lmer1))
  expect_s3_class(re_stan2, class(re_lmer2))
  check_att_names(re_stan1, re_lmer1)
  check_att_names(re_stan2, re_lmer2)
  check_sizes(re_stan1, re_lmer1)
  check_sizes(re_stan2, re_lmer2)
  expect_error(ranef(stan_glm1), "stan_glmer and stan_lmer models only")
  expect_error(ranef(stan_betareg1), "stan_glmer and stan_lmer models only")
})
test_that("fixef returns the right coefs", {
  expect_identical(names(fixef(stan_lmer1)), names(fixef(lmer1)))
  expect_identical(names(fixef(stan_lmer2)), names(fixef(lmer2)))
})
test_that("coef returns the right structure", {
  coef_stan1 <- coef(stan_lmer1); coef_lmer1 <- coef(lmer1)
  coef_stan2 <- coef(stan_lmer1); coef_lmer2 <- coef(lmer1)
  check_att_names(coef_stan1, coef_lmer1)
  check_att_names(coef_stan2, coef_lmer2)
  check_sizes(coef_stan1, coef_lmer1)
  check_sizes(coef_stan2, coef_lmer2)
})
test_that("coef ok if any 'ranef' missing from 'fixef'", {
  SW(capture.output(
    stan_lmer3 <- update(stan_lmer2, formula = . ~ (Days | Subject))
  ))
  lmer3 <- update(lmer2, formula = . ~ (Days | Subject))
  coef_stan3 <- coef(stan_lmer3); coef_lmer3 <- coef(lmer3)
  check_att_names(coef_stan3, coef_lmer3)
  check_sizes(coef_stan3, coef_lmer3)
})



# as.matrix,as.data.frame,as.array ----------------------------------------

test_that("as.matrix, as.data.frame, as.array methods work for MCMC", {
  last_dimnames <- rstanarm:::last_dimnames
  # glm
  mat <- as.matrix(stan_glm1)
  df <- as.data.frame(stan_glm1)
  arr <- as.array(stan_glm1)
  expect_identical(df, as.data.frame(mat))
  expect_identical(mat[1:2, 1], arr[1:2, 1, 1])
  expect_equal(dim(mat), c(floor(ITER/2) * CHAINS, 4L))
  expect_equal(dim(arr), c(floor(ITER/2), CHAINS, 4L))
  expect_identical(last_dimnames(mat), c("(Intercept)", "wt", "cyl", "sigma"))
  expect_identical(last_dimnames(arr), last_dimnames(mat))

  # selecting only 1 parameter
  mat <- as.matrix(stan_glm1, pars = "wt")
  df <- as.data.frame(stan_glm1, pars = "wt")
  arr <- as.array(stan_glm1, pars = "wt")
  expect_identical(df, as.data.frame(mat))
  expect_identical(mat[1:2, 1], arr[1:2, 1, 1])
  expect_equal(dim(mat), c(floor(ITER/2) * CHAINS, 1L))
  expect_equal(dim(arr), c(floor(ITER/2), CHAINS, 1L))
  expect_identical(last_dimnames(mat), "wt")
  expect_identical(last_dimnames(arr), last_dimnames(mat))

  # glmer
  mat <- as.matrix(example_model)
  df <- as.data.frame(example_model)
  arr <- as.array(example_model)
  expect_identical(df, as.data.frame(mat))
  expect_identical(mat[1:2, 1], arr[1:2, 1, 1])
  nc <- length(c(fixef(example_model), unlist(ranef(example_model)))) + 1L
  nr <- rstanarm:::posterior_sample_size(example_model)
  nms <- rownames(summary(example_model))[seq_len(nc)]
  expect_equal(dim(mat), c(nr, nc))
  expect_equal(dim(arr), c(nr / 2, 2, nc))
  expect_identical(last_dimnames(mat), nms)
  expect_identical(last_dimnames(mat), last_dimnames(arr))

  # pars & regex_pars
  mat <- as.matrix(example_model, pars = "mean_PPD", regex_pars = "period")
  df <- as.data.frame(example_model, pars = "mean_PPD", regex_pars = "period")
  arr <- as.array(example_model, pars = "mean_PPD", regex_pars = "period")
  expect_identical(df, as.data.frame(mat))
  expect_identical(mat[1:2, 1], arr[1:2, 1, 1])
  expect_equal(dim(mat), c(nr, 4L))
  expect_equal(dim(arr), c(nr/2, 2, 4L))
  expect_identical(last_dimnames(mat), c("mean_PPD", paste0("period", 2:4)))
  expect_identical(last_dimnames(mat), last_dimnames(arr))

  # lmer
  mat <- as.matrix(stan_lmer2)
  df <- as.data.frame(stan_lmer2)
  arr <- as.array(stan_lmer2)
  expect_identical(df, as.data.frame(mat))
  expect_identical(mat[1:2, 1], arr[1:2, 1, 1])
  # +1 for "sigma" and +3 for "Sigma" 
  nc <- length(c(fixef(stan_lmer2), unlist(ranef(stan_lmer2)))) + 4
  nms <- rownames(summary(stan_lmer2))[seq_len(nc)]
  expect_equal(dim(mat), c(floor(ITER/2) * CHAINS, nc))
  expect_equal(dim(arr), c(floor(ITER/2), CHAINS, nc))
  expect_identical(last_dimnames(mat), nms)
  expect_identical(last_dimnames(mat), last_dimnames(arr))
  mat <- as.matrix(stan_lmer2, pars = "(Intercept)", regex_pars = "b\\[Days Subject")
  df <- as.data.frame(stan_lmer2, pars = "(Intercept)", regex_pars = "b\\[Days Subject")
  expect_identical(df, as.data.frame(mat))
  s <- summary(stan_lmer2, pars = "(Intercept)", regex_pars = "b\\[Days Subject")
  expect_equal(dim(mat), c(floor(ITER/2) * CHAINS, nrow(s)))
  expect_identical(colnames(mat), rownames(s))

  # polr
  mat <- as.matrix(stan_polr1)
  df <- as.data.frame(stan_polr1)
  arr <- as.array(stan_polr1)
  expect_identical(df, as.data.frame(mat))
  expect_identical(mat[1:2, 1], arr[1:2, 1, 1])
  nms <- names(c(stan_polr1$coefficients, stan_polr1$zeta))
  expect_equal(dim(mat), c(floor(ITER/2) * CHAINS, length(nms)))
  expect_equal(dim(arr), c(floor(ITER/2), CHAINS, length(nms)))
  expect_identical(last_dimnames(mat), nms)
  expect_identical(last_dimnames(mat), last_dimnames(arr))
  mat <- as.matrix(stan_polr1, regex_pars = "agegp")
  df <- as.data.frame(stan_polr1, regex_pars = "agegp")
  expect_identical(df, as.data.frame(mat))

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
  
})

test_that("as.matrix and as.data.frame work for optimization and vb", {
  # optimization
  mat <- as.matrix(stan_glm_opt1)
  df <- as.data.frame(stan_glm_opt1)
  expect_identical(df, as.data.frame(mat))
  expect_equal(dim(mat), c(1000L, 4L))
  expect_identical(colnames(mat), c("(Intercept)", "wt", "cyl", "sigma"))
  mat <- as.matrix(stan_glm_opt1, pars = "sigma")
  df <- as.data.frame(stan_glm_opt1, pars = "sigma")
  expect_identical(df, as.data.frame(mat))
  expect_equal(dim(mat), c(1000, 1L))
  expect_identical(colnames(mat), "sigma")

  # vb
  mat <- as.matrix(stan_glm_vb1)
  df <- as.data.frame(stan_glm_vb1)
  expect_identical(df, as.data.frame(mat))
  expect_equal(dim(mat), c(1000L, 4L))
  expect_identical(colnames(mat), c("(Intercept)", "wt", "cyl", "sigma"))
  mat <- as.matrix(stan_glm_vb1, pars = c("(Intercept)", "sigma"))
  df <- as.data.frame(stan_glm_vb1, pars = c("(Intercept)", "sigma"))
  expect_identical(df, as.data.frame(mat))
  expect_equal(dim(mat), c(1000, 2L))
  expect_identical(colnames(mat), c("(Intercept)", "sigma"))
})
test_that("as.matrix and as.array errors & warnings", {
  # optimization and vb errors
  expect_error(as.array(stan_glm_opt1),
               regexp = "use 'as.matrix' instead")
  expect_error(as.array(stan_glm_vb1),
               regexp = "use 'as.matrix' instead")

  # pars and regex_pars errors
  expect_error(as.matrix(stan_glm1, pars = c("bad1", "sigma")),
               regexp = "No parameter(s) bad1", fixed = TRUE)
  expect_error(as.matrix(stan_glm1, regex_pars = "not a parameter"),
               regexp = "No matches for 'regex_pars'")
  expect_warning(as.matrix(stan_glm_opt1, regex_pars = "wt"),
                 regexp = "'regex_pars' ignored")
})




# terms, formula, model.frame, model.matrix, update methods -----------------
context("model.frame methods")
test_that("model.frame works properly", {
  expect_identical(model.frame(stan_glm1), model.frame(glm1))
  expect_identical(model.frame(stan_glm_opt1), model.frame(glm1))
  expect_identical(model.frame(stan_glm_vb1), model.frame(glm1))
  expect_identical(model.frame(stan_polr1), model.frame(polr1))
  expect_identical(model.frame(stan_lmer1), model.frame(lmer1))
  expect_identical(model.frame(stan_lmer2), model.frame(lmer2))
  # lme4 is doing something different with the names
  # expect_identical(model.frame(stan_lmer1, fixed.only = TRUE),
  #                  model.frame(lmer1, fixed.only = TRUE))
  # expect_identical(model.frame(stan_lmer2, fixed.only = TRUE),
  #                  model.frame(lmer2, fixed.only = TRUE))
  expect_identical(model.frame(stan_betareg1), model.frame(betareg1))
})

context("terms methods")
test_that("terms works properly", {
  expect_identical(terms(stan_glm1), terms(glm1))
  expect_identical(terms(stan_glm_opt1), terms(glm1))
  expect_identical(terms(stan_glm_vb1), terms(glm1))
  expect_identical(terms(stan_polr1), terms(polr1))
  expect_identical(terms(stan_lmer1), terms(lmer1))
  expect_identical(terms(stan_lmer2), terms(lmer2))
  expect_identical(terms(stan_lmer1, fixed.only = TRUE),
                   terms(lmer1, fixed.only = TRUE))
  expect_identical(terms(stan_lmer2, fixed.only = TRUE),
                   terms(lmer2, fixed.only = TRUE))
  expect_equal(terms(stan_lmer1, random.only = TRUE),
                   terms(lmer1, random.only = TRUE))
  expect_equal(terms(stan_lmer2, random.only = TRUE),
               terms(lmer2, random.only = TRUE))
  expect_error(terms(stan_lmer1, fixed.only = TRUE, random.only = TRUE),
               regexp = "can't both be TRUE")
  expect_identical(terms(stan_betareg1), terms(betareg1))
})

context("formula methods")
test_that("formula works properly", {
  expect_identical(formula(stan_glm1), formula(glm1))
  expect_identical(formula(stan_glm_opt1), formula(glm1))
  expect_identical(formula(stan_glm_vb1), formula(glm1))
  expect_identical(formula(stan_betareg1), formula(betareg1))
  expect_equal(terms(stan_polr1), formula(polr1))
  expect_identical(formula(stan_lmer1), formula(lmer1))
  expect_identical(formula(stan_lmer2), formula(lmer2))
  expect_identical(formula(stan_lmer1, fixed.only = TRUE),
                   formula(lmer1, fixed.only = TRUE))
  expect_identical(formula(stan_lmer2, fixed.only = TRUE),
                   formula(lmer2, fixed.only = TRUE))
  expect_equal(formula(stan_lmer1, random.only = TRUE),
               formula(lmer1, random.only = TRUE))
  expect_equal(formula(stan_lmer2, random.only = TRUE),
               formula(lmer2, random.only = TRUE))
  expect_error(formula(stan_lmer1, fixed.only = TRUE, random.only = TRUE),
               regexp = "can't both be TRUE")
  

  tmp <- stan_lmer1
  tmp$formula <- NULL
  attr(tmp$glmod$fr, "formula") <- NULL
  expect_equal(formula(tmp), formula(lmer1))
  tmp$call <- NULL
  expect_error(formula(tmp), regexp = "can't find formula", ignore.case = TRUE)
})

context("update methods")
test_that("update works properly", {
  pss <- rstanarm:::posterior_sample_size

  SW(fit1 <- update(stan_lmer2, iter = ITER * 2, chains = 2 * CHAINS))
  SW(fit2 <- update(stan_glm1, iter = ITER * 2, chains = 2 * CHAINS))
  SW(fit3 <- update(stan_betareg1, iter = ITER * 2, chains = CHAINS * 2))
  expect_equal(pss(fit1), 4 * pss(stan_lmer2))
  expect_equal(pss(fit2), 4 * pss(stan_glm1))
  expect_equal(pss(fit3), 4 * pss(stan_betareg1))
  
  call_only <- update(fit1, evaluate = FALSE)
  expect_is(call_only, "call")
  expect_identical(call_only, getCall(fit1))

  # expect_error(fit2 <- update(fit2, algorithm = "optimizing"),
  #              regexp = "unknown arguments: chains")
  expect_identical(fit2$algorithm, "sampling")

  fit2$call <- NULL
  expect_error(update(fit2), regexp = "does not contain a 'call' component")
})


# print and summary -------------------------------------------------------
context("print and summary methods")
test_that("print and summary methods ok for mcmc and vb", {
  expect_output(print(example_model, digits = 2), "stan_glmer")
  expect_output(print(example_model, digits = 2), "Error terms")
  expect_output(print(stan_lmer1, digits = 2), "stan_lmer")
  expect_output(print(stan_lmer2), "stan_lmer")
  expect_output(print(stan_polr1), "stan_polr")
  expect_output(print(stan_polr1), "Cutpoints")
  expect_output(print(stan_glm_opt1, digits = 5), "stan_glm")
  expect_output(print(stan_glm_vb1, digits = 5), "stan_glm")
  expect_output(print(stan_betareg1, digits = 2), "stan_betareg")
  
  expect_silent(s <- summary(stan_lmer1, pars = "varying", regex_pars = "Sigma"))
  expect_silent(s_alt <- summary(stan_lmer1, regex_pars = c("plate", "sample")))
  expect_identical(s, s_alt)
  expect_silent(s <- summary(stan_lmer1))
  expect_silent(d <- as.data.frame(s))
  expect_s3_class(s, "summary.stanreg")
  expect_output(print(s), "stan_lmer")
  expect_identical(attr(s, "algorithm"), "sampling")
  expect_identical(colnames(s), colnames(d))
  expect_identical(rownames(s), rownames(d))

  expect_silent(s <- summary(example_model, pars = "beta", regex_pars = "herd"))
  expect_silent(s_alt <- summary(example_model, pars = c("beta", "varying"), regex_pars = "Sigma"))
  expect_identical(s, s_alt)
  expect_silent(d <- as.data.frame(s))
  expect_s3_class(s, "summary.stanreg")
  expect_output(print(s), "stan_glmer")
  expect_output(
    print(s), 
    paste(rstanarm:::posterior_sample_size(example_model)), "(posterior sample size)"
  )
  expect_identical(attr(s, "algorithm"), "sampling")
  expect_identical(colnames(s), colnames(d))
  expect_identical(rownames(s), rownames(d))

  expect_silent(s <- summary(stan_polr1, pars = "beta", probs = c(0.25, 0.75)))
  expect_silent(d <- as.data.frame(s))
  expect_identical(colnames(s), c("mean", "mcse", "sd", "25%", "75%", "n_eff", "Rhat"))
  expect_identical(colnames(s), colnames(d))
  expect_identical(rownames(s), rownames(d))
  expect_s3_class(s, "summary.stanreg")
  expect_output(print(s), "stan_polr")

  expect_warning(s <- summary(stan_glm1, pars = "varying"),
                 regexp = "No group-specific parameters. 'varying' ignored.")
  expect_silent(s <- summary(stan_glm1, pars = c("alpha", "beta"), digits = 3))
  expect_s3_class(s, "summary.stanreg")
  expect_output(print(s), "stan_glm")
  expect_identical(attr(s, "algorithm"), "sampling")
  
  expect_silent(s <- summary(stan_glm_vb1, pars = c("alpha", "beta")))
  expect_silent(d <- as.data.frame(s))
  expect_s3_class(s, "summary.stanreg")
  expect_output(print(s), "stan_glm")
  expect_identical(attr(s, "algorithm"), "meanfield")
  
  expect_warning(s <- summary(stan_betareg1, pars = "varying"),
                 regexp = "No group-specific parameters. 'varying' ignored.")
  expect_silent(s <- summary(stan_betareg1, pars = c("alpha", "beta"), digits = 3))
  expect_s3_class(s, "summary.stanreg")
  expect_output(print(s), "stan_betareg")
  expect_identical(attr(s, "algorithm"), "sampling")
  
})

test_that("print and summary methods ok for optimization", {
  expect_silent(s <- summary(stan_glm_opt1))
  expect_silent(s <- summary(stan_glm_opt1, pars = c("wt", "sigma"), digits = 8))
  expect_warning(s <- summary(stan_glm_opt1, regex_pars = c("wt", "sigma")),
                 regexp = "'regex_pars' ignored")
  expect_silent(d <- as.data.frame(s))
  expect_s3_class(s, "summary.stanreg")
  expect_output(print(s), "stan_glm")
  expect_identical(attr(s, "algorithm"), "optimizing")
  expect_identical(colnames(s), colnames(d))
  expect_identical(rownames(s), rownames(d))

  counts <- c(18,17,15,20,10,20,25,13,12)
  outcome <- gl(3,1,9)
  treatment <- gl(3,3)  
  capture.output(
    fit <- stan_glm.nb(counts ~ outcome + treatment, algorithm = "optimizing",
                       seed = SEED, refresh = 0)
  )
  expect_output(print(fit), "reciprocal_dispersion")

  clotting <- data.frame(log_u = log(c(5,10,15,20,30,40,60,80,100)),
                         lot1 = c(118,58,42,35,27,25,21,19,18),
                         lot2 = c(69,35,26,21,18,16,13,12,12))
  capture.output(
    fit2 <- stan_glm(lot1 ~ log_u, data = clotting, family = Gamma(link="log"),
                     algorithm = "optimizing", seed = SEED, refresh = 0),
    fit3 <- update(fit2, family = inverse.gaussian(link = "log"))
  )
  expect_output(print(fit2), "shape")
  expect_output(print(fit3), "lambda")
})

# prior_summary -----------------------------------------------------------
test_that("prior_summary errors if info not found", {
  tmp <- example_model
  tmp$prior.info <- NULL
  expect_message(s <- prior_summary(tmp), "Priors not found in stanreg object")
  expect_null(s)
})
test_that("prior_summary doesn't error", {
  expect_output(print(prior_summary(example_model, digits = 2)),
                "Priors for model 'example_model'")
  expect_output(print(prior_summary(stan_lmer1, digits = 2)),
                "stan_lmer1")
  expect_output(print(prior_summary(stan_lmer2)),
                "stan_lmer2")
  expect_output(print(prior_summary(stan_polr1)),
                "stan_polr1")
  expect_output(print(prior_summary(stan_glm_opt1)),
                "stan_glm_opt1")
  expect_output(print(prior_summary(stan_glm_vb1)),
                "stan_glm_vb1")
  expect_output(print(prior_summary(stan_betareg1)),
                "stan_betareg1")
  
})
test_that("prior_summary returns correctly named list", {
  expect_named(prior_summary(example_model),
               c("prior", "prior_intercept", "prior_covariance"))
  expect_named(prior_summary(stan_lmer1),
               c("prior", "prior_intercept", "prior_covariance", "prior_aux"))
  expect_named(prior_summary(stan_lmer2),
               c("prior", "prior_intercept", "prior_covariance", "prior_aux"))
  expect_named(prior_summary(stan_polr1),
               c("prior", "prior_counts"))
  expect_named(prior_summary(stan_glm_opt1),
               c("prior", "prior_intercept", "prior_aux"))
  expect_named(prior_summary(stan_glm_vb1),
               c("prior", "prior_intercept", "prior_aux"))
  expect_named(prior_summary(stan_betareg1),
               c("prior", "prior_z", "prior_intercept", "prior_intercept_z", "prior_aux"))  
})


# predictive_error,predictive_interval ------------------------------------
context("predictive error and interval methods")
test_that("predictive_error works", {
  expect_error(predictive_error(stan_glm1, draws = 100),
               "'draws' should be <= posterior sample size")
  expect_error(predictive_error(stan_polr1),
               "not currently available for stan_polr")
  expect_error(predictive_error(stan_betareg1, draws = 600),
               "'draws' should be <= posterior sample size")
  
  mods <- c("stan_glm1", "stan_glm_vb1", "stan_lmer1",
            "stan_lmer2", "example_model")
  for (m in seq_along(mods)) {
    mod <- get(mods[m])
    err <- predictive_error(mod, draws = 5)
    expect_equal(dim(err), c(5, nobs(mod)), info = mods[m])
  }

  err2 <- predictive_error(stan_glm1, newdata = model.frame(stan_glm1)[1:10, ],
                           draws = 7)
  expect_equal(dim(err2), c(7, 10))

  err3 <- predictive_error(example_model, draws = 5,
                           newdata = data.frame(
                             size = c(10, 20),
                             incidence = c(5, 10),
                             period = factor(c(1,2)),
                             herd = c(1, 15)
                            ))
  expect_equal(dim(err3), c(5, 2))
})
test_that("predictive_interval works", {
  expect_error(predictive_interval(stan_glm1, draws = 100),
               "'draws' should be <= posterior sample size")
  expect_error(predictive_interval(stan_glm1, prob = c(0.25, 0.76)),
               "'prob' should be a single number greater than 0 and less than 1")
  expect_error(predictive_interval(stan_polr1),
               "not currently available for stan_polr")
  expect_error(predictive_interval(stan_betareg1, draws = 600),
               "'draws' should be <= posterior sample size")
  expect_error(predictive_interval(stan_betareg1, prob = c(0.25, 0.76)),
               "'prob' should be a single number greater than 0 and less than 1")
  
  mods <- c("stan_glm1", "stan_glm_vb1", "stan_lmer1",
            "stan_lmer2", "example_model")
  for (m in seq_along(mods)) {
    mod <- get(mods[m])
    pint1 <- predictive_interval(mod, draws = 5)
    expect_equal(dim(pint1), c(nobs(mod), 2), info = mods[m])
    expect_identical(colnames(pint1), c("5%", "95%"), info = mods[m])
  }

  pint2 <- predictive_interval(stan_glm1, prob = 0.5, newdata = model.frame(stan_glm1)[1:2, ])
  expect_equal(dim(pint2), c(2, 2))
  expect_identical(colnames(pint2), c("25%", "75%"))

  pint3 <- predictive_interval(example_model, prob = 0.8, newdata = lme4::cbpp[1:10, ])
  expect_equal(dim(pint3), c(10, 2))
  expect_identical(colnames(pint3), c("10%", "90%"))
})

test_that("predictive_error stanreg and ppd methods return the same thing", {
  preds <- posterior_predict(stan_glm1, seed = 123)
  expect_equal(
    predictive_error(stan_glm1, seed = 123),
    predictive_error(preds, y = stan_glm1$y)
  )
  preds <- posterior_predict(stan_betareg1, seed = 123)
  expect_equal(
    predictive_error(stan_betareg1, seed = 123),
    predictive_error(preds, y = stan_betareg1$y)
  )
})
test_that("predictive_interval stanreg and ppd methods return the same thing", {
  preds <- posterior_predict(stan_glm1, seed = 123)
  expect_equal(
    predictive_interval(stan_glm1, seed = 123),
    predictive_interval(preds)
  )
  preds <- posterior_predict(stan_betareg1, seed = 123)
  expect_equal(
    predictive_interval(stan_betareg1, seed = 123),
    predictive_interval(preds)
  )
})



# stanreg lists -----------------------------------------------------------
test_that("stan*_list functions throw proper errors", {
  expect_error(stanreg_list(), ">= 1 is not TRUE")
  expect_error(stanreg_list(stan_glm1, glm1), "For stanreg_list")
  expect_error(stanmvreg_list(stan_glm1, glm1), "For stanmvreg_list")
  expect_error(stanjm_list(stan_glm1, glm1), "For stanjm_list")
})

test_that("stanreg_list works", {
  list1 <- stanreg_list(stan_lmer1, stan_lmer2)
  expect_named(list1, c("stan_lmer1", "stan_lmer2"))
  expect_equivalent(attr(list1, "families"), c("gaussian", "gaussian"))
  expect_identical(list1$stan_lmer1, stan_lmer1)
  expect_identical(list1$stan_lmer2, stan_lmer2)
})
