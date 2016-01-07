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

library(rstanarm)
library(lme4)
library(MASS)
SEED <- 12345
set.seed(SEED)
ITER <- 10
CHAINS <- 2
CORES <- 1
REFRESH <- ITER

stan_glm1 <- suppressWarnings(stan_glm(mpg ~ wt, data = mtcars, iter = ITER, 
                                 chains = CHAINS, cores = CORES, seed = SEED, 
                                 refresh = REFRESH))
stan_glm_opt1 <- stan_glm(mpg ~ wt, data = mtcars, algorithm = "optimizing")
stan_glm_vb1 <- update(stan_glm_opt1, algorithm = "meanfield", iter = 10000, 
                       seed = SEED)
glm1 <- glm(mpg ~ wt, data = mtcars)

lmer1 <- lmer(diameter ~ (1|plate) + (1|sample), data = Penicillin)
stan_lmer1 <- suppressWarnings(stan_lmer(diameter ~ (1|plate) + (1|sample), 
                                        data = Penicillin, iter = ITER, 
                                        chains = CHAINS, cores = CORES,
                                        seed = SEED, refresh = REFRESH))
lmer2 <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy)
stan_lmer2 <- suppressWarnings(stan_lmer(Reaction ~ Days + (Days | Subject), 
                                        data = sleepstudy, iter = ITER, 
                                        chains = CHAINS, cores = CORES,
                                        seed = SEED, refresh = REFRESH))

stan_polr1 <- suppressWarnings(stan_polr(tobgp ~ agegp, data = esoph,
                                       prior = R2(0.2, "mean"), init_r = 0.1, 
                                       iter = ITER, chains = CHAINS, 
                                       cores = CORES, seed = SEED, refresh = REFRESH))
polr1 <- polr(tobgp ~ agegp, data = esoph, Hess = TRUE)

att_names <- function(object) {
  nms <- names(object)
  att_nms <- names(attributes(object))
  att_nms2 <- lapply(object, function(x) names(attributes(x)))
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
  expect_equal(vcov(stan_glm_opt1, correlation = TRUE), cov2cor(stan_glm_opt1$covmat))
  expect_equal(resid(stan_glm_opt1), stan_glm_opt1$residuals)
  expect_equal(coef(stan_glm_opt1), stan_glm_opt1$coefficients)
  expect_equal(fitted(stan_glm_opt1), stan_glm_opt1$fitted.values)
  expect_equal(se(stan_glm_opt1), stan_glm_opt1$ses)
})

test_that("confint method returns correct structure", {
  expect_silent(ci <- confint(stan_glm_opt1))
  expect_silent(ci2 <- confint(stan_glm_opt1, parm = "wt", level = 0.9))
  expect_equal(rownames(ci), c("(Intercept)", "wt"))
  expect_equal(colnames(ci), c("2.5 %", "97.5 %"))
  expect_equal(rownames(ci2), "wt")
  expect_equal(colnames(ci2), c("5 %", "95 %"))
  
  expect_error(confint(stan_glm1), regexp = "use prob_int")
  expect_error(confint(stan_glm_vb1), regexp = "use prob_int")
  expect_error(confint(stan_polr1), regexp = "use prob_int")
  expect_error(confint(stan_lmer1), regexp = "use prob_int")
  expect_error(confint(stan_lmer2), regexp = "use prob_int")
})

test_that("prob_int returns correct structure", {
  # NOTE: prob_int is not an S3 method but is tested here because it's the 
  # analog of confint but for MCMC and VB
  expect_silent(ci <- prob_int(stan_glm1, prob = 0.5))
  expect_silent(ci2 <- prob_int(stan_glm_vb1, pars = "wt", prob = 0.95))
  expect_silent(ci3 <- prob_int(example_model, prob = 0.95, 
                                regex_pars = "herd", type = "spin"))
  expect_silent(ci4 <- prob_int(example_model, prob = 0.75, pars = "(Intercept)", 
                               regex_pars = "period", type = "spin"))
  expect_silent(ci5 <- prob_int(stan_polr1, prob = 0.8))
  expect_equal(rownames(ci), c("(Intercept)", "wt", "sigma"))
  expect_equal(rownames(ci2), "wt")
  expect_equal(rownames(ci3), rstanarm:::.bnames(rownames(example_model$stan_summary), value = TRUE)[1:15])
  expect_equal(rownames(ci4), c("(Intercept)", paste0("period", 2:4)))
  expect_equal(colnames(ci), c("25%", "75%"))
  expect_equal(colnames(ci2), c("2.5%", "97.5%"))
  expect_equal(colnames(ci2), c("2.5%", "97.5%"))
  expect_equal(colnames(ci3), c("lower95", "upper95"))
  expect_equal(colnames(ci4), c("lower75", "upper75"))
  
  expect_error(prob_int(stan_glm_opt1), regexp = "not available")
  expect_error(prob_int(lm(mpg ~ wt, data = mtcars)), 
               regexp = "not a stanreg object")
  
  prob_msg <- "'prob' should be a single number greater than 0 and less than 1."
  expect_error(prob_int(stan_glm1, prob = c(0.25, 0.75)), regexp = prob_msg)
  expect_error(prob_int(stan_glm1, prob = 0), regexp = prob_msg)
  expect_error(prob_int(stan_glm1, prob = 1), regexp = prob_msg)
  expect_error(prob_int(stan_glm1, prob = 2), regexp = prob_msg)
})

test_that("log_lik method works", {
  expect_error(log_lik(stan_glm_opt1))
  expect_error(log_lik(stan_glm_vb1))
  expect_silent(log_lik(stan_glm1))
  
  expect_silent(log_lik(stan_polr1))
  expect_equal(dim(log_lik(stan_polr1)), c(ITER, nobs(stan_polr1)))
  expect_equal(dim(log_lik(stan_lmer1)), c(ITER, nobs(stan_lmer1)))
  
  # Compute log-lik matrix using different method than log_lik.stanreg
  # and compare
  samp <- as.matrix(stan_glm1)
  y <- get_y(stan_glm1)
  eta <- tcrossprod(get_x(stan_glm1), samp[, 1:2])
  sigma <- samp[, 3]
  llmat <- matrix(NA, nrow = nrow(samp), ncol = nrow(eta))
  for (i in 1:nrow(llmat)) {
    llmat[i, ] <- dnorm(y, mean = eta[, i], sd = sigma[i], log = TRUE)
  }
  expect_equal(llmat, log_lik(stan_glm1))
})

test_that("ngrps is right", {
  expect_equal(ngrps(lmer1), ngrps(stan_lmer1))
  expect_equal(ngrps(lmer2), ngrps(stan_lmer2))
  expect_error(ngrps(stan_glm1), "stan_glmer and stan_lmer models only")
})

test_that("nobs is right", {
  expect_equal(nobs(lmer1), nobs(stan_lmer1))
  expect_equal(nobs(lmer2), nobs(stan_lmer2))
  expect_equal(nobs(glm1), nobs(stan_glm_opt1))
  expect_equal(nobs(glm1), nobs(stan_glm1))
  expect_equal(nobs(polr1), nobs(stan_polr1))
})

test_that("vcov returns correct structure", {
  expect_equal(dimnames(vcov(stan_glm1)), dimnames(vcov(glm1)))
  expect_equal(dimnames(vcov(stan_polr1)), dimnames(vcov(polr1)))
  expect_equal(dimnames(vcov(stan_lmer1)), dimnames(vcov(lmer1)))
  expect_equal(dimnames(vcov(stan_lmer2)), dimnames(vcov(lmer2)))
})

test_that("VarCorr returns correct structure", {
  vc_lmer1 <- VarCorr(lmer1); vc_stan1 <- VarCorr(stan_lmer1)
  vc_lmer2 <- VarCorr(lmer2); vc_stan2 <- VarCorr(stan_lmer2)
  expect_is(vc_stan1, class(vc_lmer1))
  expect_is(vc_stan2, class(vc_lmer2))
  check_att_names(vc_stan1, vc_lmer1)
  check_att_names(vc_stan2, vc_lmer2)
  expect_error(VarCorr(stan_glm1), "stan_glmer and stan_lmer models only")
})

test_that("ranef returns correct structure", {
  re_stan1 <- ranef(stan_lmer1); re_lmer1 <- ranef(lmer1)
  re_stan2 <- ranef(stan_lmer1); re_lmer2 <- ranef(lmer1)
  expect_is(re_stan1, class(re_lmer1))
  expect_is(re_stan2, class(re_lmer2))
  check_att_names(re_stan1, re_lmer1)
  check_att_names(re_stan2, re_lmer2)
  check_sizes(re_stan1, re_lmer1)
  check_sizes(re_stan2, re_lmer2)
  expect_error(ranef(stan_glm1), "stan_glmer and stan_lmer models only")
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

test_that("print and summary methods don't throw errors", {
  expect_output(print(stan_lmer1, digits = 2), "stan_lmer")
  expect_output(print(stan_lmer2), "stan_lmer")
  expect_output(print(stan_polr1), "stan_polr")
  expect_output(print(stan_glm_opt1, digits = 5), "stan_glm")
  expect_output(print(stan_glm_vb1, digits = 5), "stan_glm")
  
  expect_silent(summary(stan_lmer1, pars = "varying"))
  expect_silent(s <- summary(stan_lmer1))
  expect_silent(d <- as.data.frame(s))
  expect_is(s, "summary.stanreg")
  expect_output(print(s), "stan_lmer")
  expect_equal(attr(s, "algorithm"), "sampling")
  expect_equal(colnames(s), colnames(d))
  expect_equal(rownames(s), rownames(d))
  
  expect_silent(s <- summary(stan_glm_opt1))
  expect_silent(s <- summary(stan_glm_opt1, pars = c("wt", "sigma"), digits = 8))
  expect_silent(d <- as.data.frame(s))
  expect_is(s, "summary.stanreg")
  expect_output(print(s), "stan_glm")
  expect_equal(attr(s, "algorithm"), "optimizing")
  expect_equal(colnames(s), colnames(d))
  expect_equal(rownames(s), rownames(d))
  
  expect_silent(s <- summary(stan_polr1, pars = "beta", probs = c(0.25, 0.75)))
  expect_silent(d <- as.data.frame(s))
  expect_equal(colnames(s), c("mean", "mcse", "sd", "25%", "75%", "n_eff", "Rhat"))
  expect_equal(rownames(s), c("agegp.L", "agegp.Q", "agegp.C", "agegp^4", "agegp^5"))
  expect_equal(colnames(s), colnames(d))
  expect_equal(rownames(s), rownames(d))
  expect_is(s, "summary.stanreg")
  expect_output(print(s), "stan_polr")
  
  expect_silent(s <- summary(stan_glm1, pars = c("alpha", "beta"), digits = 3))
  expect_is(s, "summary.stanreg")
  expect_output(print(s), "stan_glm")
  expect_equal(attr(s, "algorithm"), "sampling")
  
  expect_silent(s <- summary(stan_glm_vb1, pars = c("alpha", "beta")))
  expect_silent(d <- as.data.frame(s))
  expect_is(s, "summary.stanreg")
  expect_output(print(s), "stan_glm")
  expect_equal(attr(s, "algorithm"), "meanfield")
})


context("terms, formula, model.frame, and model.matrix methods")

test_that("model.frame works properly", {
  expect_identical(model.frame(stan_glm1), model.frame(glm1))
  expect_identical(model.frame(stan_glm_opt1), model.frame(glm1))
  expect_identical(model.frame(stan_glm_vb1), model.frame(glm1))
  expect_identical(model.frame(stan_polr1), model.frame(polr1))
  expect_identical(model.frame(stan_lmer1), model.frame(lmer1))
  expect_identical(model.frame(stan_lmer2), model.frame(lmer2))
  expect_identical(model.frame(stan_lmer1, fixed.only = TRUE), 
                   model.frame(lmer1, fixed.only = TRUE))
  expect_identical(model.frame(stan_lmer2, fixed.only = TRUE), 
                   model.frame(lmer2, fixed.only = TRUE))
})

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
})

test_that("formula works properly", {
  expect_identical(formula(stan_glm1), formula(glm1))
  expect_identical(formula(stan_glm_opt1), formula(glm1))
  expect_identical(formula(stan_glm_vb1), formula(glm1))
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
