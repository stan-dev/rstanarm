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

# tests can be run using devtools::test() or manually by loading testthat 
# package and then running the code below possibly with options(mc.cores = 4).

library(rstanarm)
SEED <- 12345
set.seed(SEED)
CHAINS <- 2
ITER <- 40 # small iter for speed but large enough for psis
REFRESH <- 0

source(file.path("helpers", "expect_stanreg.R"))
source(file.path("helpers", "SW.R"))

SW(
  fit_gaus <- stan_glm(mpg ~ wt, data = mtcars, 
                       chains = CHAINS, iter = ITER,
                       seed = SEED, refresh = REFRESH)
)
dat <- data.frame(ldose = rep(0:5, 2),
                  sex = factor(rep(c("M", "F"), c(6, 6))))
numdead <- c(1, 4, 9, 13, 18, 20, 0, 2, 6, 10, 12, 16)
SF <- cbind(numdead, numalive = 20-numdead)
SW(
  fit_binom <- stan_glm(SF ~ sex*ldose, data = dat, family = binomial,
                        chains = CHAINS, iter = ITER, seed = SEED,
                        refresh = REFRESH)
)
dead <- rbinom(length(numdead), 1, prob = 0.5)
SW(fit_binom2 <- update(fit_binom, formula = factor(dead) ~ .))

d.AD <- data.frame(treatment = gl(3,3), outcome =  gl(3,1,9),
                   counts = c(18,17,15,20,10,20,25,13,12))
SW(fit_pois <- stan_glm(counts ~ outcome + treatment, data = d.AD,
                        family = poisson, chains = CHAINS, iter = ITER,
                        seed = SEED, refresh = REFRESH))
SW(fit_negbin <- update(fit_pois, family = neg_binomial_2))

clotting <- data.frame(log_u = log(c(5,10,15,20,30,40,60,80,100)),
                       lot1 = c(118,58,42,35,27,25,21,19,18),
                       lot2 = c(69,35,26,21,18,16,13,12,12))
SW(fit_gamma <- stan_glm(lot1 ~ log_u, data = clotting, family = Gamma,
                         chains = CHAINS, iter = ITER, seed = SEED,
                         refresh = REFRESH))
SW(fit_igaus <- update(fit_gamma, family = inverse.gaussian))

test_that("loo/waic for stan_glm works", {
  ll_fun <- rstanarm:::ll_fun
  source(file.path("helpers", "expect_equivalent_loo.R"))
  
  # gaussian
  expect_equivalent_loo(fit_gaus)
  expect_identical(ll_fun(fit_gaus), rstanarm:::.ll_gaussian_i)
  
  # binomial
  expect_equivalent_loo(fit_binom)
  expect_equivalent_loo(fit_binom2)
  expect_identical(ll_fun(fit_binom), rstanarm:::.ll_binomial_i)
  expect_identical(ll_fun(fit_binom2), rstanarm:::.ll_binomial_i)
  
  # poisson
  expect_equivalent_loo(fit_pois)
  expect_identical(ll_fun(fit_pois), rstanarm:::.ll_poisson_i)
  
  # negative binomial
  expect_equivalent_loo(fit_negbin)
  expect_identical(ll_fun(fit_negbin), rstanarm:::.ll_neg_binomial_2_i)
  
  # gamma
  expect_equivalent_loo(fit_gamma)
  expect_identical(ll_fun(fit_gamma), rstanarm:::.ll_Gamma_i)
  
  # inverse gaussian
  expect_equivalent_loo(fit_igaus)
  expect_identical(ll_fun(fit_igaus), rstanarm:::.ll_inverse.gaussian_i)
})

context("stan_glm (errors, warnings, messages)")
test_that("stan_glm throws appropriate errors, warnings, and messages", {
  counts <- c(18,17,15,20,10,20,25,13,12)
  outcome <- gl(3,1,9)
  treatment <- gl(3,3)
  f <- as.formula(counts ~ outcome + treatment)
  
  # error: glmer syntax
  expect_error(stan_glm(counts ~ treatment + (1|outcome)), 
               regexp = "model formula not allowed")
  
  # error: empty model
  expect_error(stan_glm(counts ~ 0), 
               regexp = "No intercept or predictors specified")
  
  # error: stan_glm.nb with family argument
  expect_error(stan_glm.nb(f, family = "neg_binomial_2"), 
               regexp = "'family' should not be specified.")
  
  # error: prior and prior_intercept not lists
  expect_error(stan_glm(f, family = "poisson", prior = normal), 
               regexp = "should be a named list")
  expect_error(stan_glm(f, family = "poisson", prior_intercept = normal), 
               regexp = "should be a named list")
  
  # error: QR only with more than 1 predictor
  expect_error(stan_glm(counts ~ 1, family = "poisson", QR = TRUE), 
               regexp = "'QR' can only be specified when there are multiple predictors")
  
  # error: QR and sparse
  expect_error(stan_glm(f, family = "poisson", QR = TRUE, sparse = TRUE), 
               regexp = "'QR' and 'sparse' cannot both be TRUE")
  
  # require intercept for certain family and link combinations
  expect_error(stan_glm(counts ~ -1 + outcome + treatment, 
                        family = poisson(link="identity"), seed = SEED), 
               regexp = "model must have an intercept")
  expect_error(stan_glm(I(counts > 20) ~ -1 + outcome + treatment, 
                        family = binomial(link="log"), seed = SEED), 
               regexp = "model must have an intercept")
  
  # support of outcome variable
  expect_error(stan_glm(cbind(1:10, runif(10)) ~ 1, family = "binomial"), 
               "outcome values must be counts")
  expect_error(stan_glm(c(1,2,1,2) ~ 1, family = "binomial"), 
               "outcome values must be 0 or 1")
  expect_error(stan_glm((-1):3 ~ 1, family = "poisson"), 
               "outcome values must be counts")
  expect_error(stan_glm.nb(runif(3) ~ 1), 
               "outcome values must be counts")
  expect_error(stan_glm(0:3 ~ 1, family = "Gamma"), 
               "outcome values must be positive")
  expect_error(stan_glm(runif(3, -2, -1) ~ 1, family = "inverse.gaussian"), 
               "outcome values must be positive")
  expect_error(stan_glm(cbind(1:10, 1:10) ~ 1, family = "gaussian"), 
               "should not have multiple columns")
  
  # prior_aux can't be NULL if prior_PD is TRUE
  expect_error(stan_glm(mpg ~ wt, data = mtcars, prior_aux = NULL, prior_PD = TRUE),
               "'prior_aux' cannot be NULL if 'prior_PD' is TRUE")
})

context("stan_glm (gaussian)")
test_that("gaussian returns expected result for trees example", {
  # example using trees dataset
  links <- c("identity", "log", "inverse")
  for (i in 1:length(links)) {
    if (links[i] == "inverse") next # unreliable
    fit <- stan_glm(Volume ~ log(Girth) + log(Height), data = trees, 
                    family = gaussian(link = links[i]), algorithm = "optimizing",
                    prior = NULL, prior_intercept = NULL,
                    QR = TRUE, tol_rel_grad = 1e-16, seed = SEED)
    expect_stanreg(fit)
    
    ans <- glm(Volume ~ log(Girth) + log(Height),data = trees, 
               family = gaussian(link = links[i]))
    expect_equal(coef(fit), coef(ans), tol = 0.021)
  }
  
  expect_error(update(fit, prior = dnorm), 
               regexp = "should be a named list")
  expect_error(update(fit, prior_intercept = dnorm), 
               regexp = "should be a named list")
  expect_error(update(fit, prior = R2(0.5)), 
               regexp = "should be one of")
  expect_error(update(fit, prior_intercept = R2(0.5)), 
               regexp = "should be one of")
})

context("stan_glm (poisson)")
links <- c("log", "identity", "sqrt")
test_that("stan_glm returns expected result for glm poisson example", {
  # example from help("glm")
  for (i in 1:length(links)) {
    fit <- stan_glm(counts ~ outcome + treatment, data = d.AD,
                    family = poisson(links[i]), 
                    prior = NULL, prior_intercept = NULL, QR = TRUE,
                    algorithm = "optimizing", tol_rel_grad = 1e-16, seed = SEED)
    expect_stanreg(fit)
    
    ans <- glm(counts ~ outcome + treatment, data = d.AD,
               family = poisson(links[i]), start = coef(fit))
    if (links[i] == "log") expect_equal(coef(fit), coef(ans), tol = 0.01)
    # if (links[i] == "identity") expect_equal(coef(fit)[-1], coef(ans)[-1], tol = 0.03)
    if (links[i] == "sqrt") { # this is weird
      if (coef(ans)[1] > 0)
        expect_equal(coef(fit)[-1], coef(ans)[-1], tol = 0.03)
      else
        expect_equal(-coef(fit)[-1], coef(ans)[-1], tol = 0.03)
    }
  }
})

context("stan_glm (negative binomial)")
if (require(MASS)) 
  test_that("stan_glm returns something for glm negative binomial example", {
  # example from MASS::glm.nb
  
  for (i in 1:length(links)) {
    fit1 <- stan_glm(Days ~ Sex/(Age + Eth*Lrn), data = quine, 
                     family = neg_binomial_2(links[i]), 
                     seed = SEED, chains = 1, iter = 100,
                     prior_PD = TRUE, QR = TRUE, refresh = 100)
    fit2 <- stan_glm.nb(Days ~ Sex/(Age + Eth*Lrn), data = quine, 
                        link = links[i],
                        seed = SEED, chains = 1, iter = 100,
                        prior_PD = TRUE, QR = TRUE, refresh = 100)
    expect_stanreg(fit1)
    expect_stanreg(fit2)
    expect_equal(as.matrix(fit1), as.matrix(fit2))
  }
  # testing results against MASS::glm.nb is unreliable
})

context("stan_glm (gaussian)")
test_that("stan_glm returns expected result for cars example", {
  # example using cars dataset
  fit <- stan_glm(log(dist) ~ log(speed), data = cars, sparse = TRUE,
                  family = gaussian(link = "identity"), seed  = SEED,
                  prior = NULL, prior_intercept = NULL,
                  tol_rel_obj = .Machine$double.eps, algorithm = "optimizing")
  expect_stanreg(fit)
  
  ans <- glm(log(dist) ~ log(speed), data = cars, family = gaussian(link = "identity"))
  expect_equal(coef(fit), coef(ans), tol = 0.1)
})
test_that("stan_glm returns expected result with no intercept for mtcars example", {
  f <- as.formula(mpg ~ -1 + wt + cyl + disp + am + carb)
  fit <- stan_glm(f, data = mtcars,
                  prior = NULL, prior_intercept = NULL,
                  tol_rel_obj = .Machine$double.eps, algorithm = "optimizing",
                  seed  = SEED, sparse = TRUE)
  expect_stanreg(fit)
  
  ans <- glm(f, data = mtcars, family = gaussian(link = "identity"))
  expect_equal(coef(fit), coef(ans), tol = 0.04)
})

context("stan_glm (bernoulli)")
links <- c("logit", "probit", "cauchit", "log", "cloglog")
test_that("stan_glm returns expected result for bernoulli", {
  # bernoulli example
  sd1 <- 1; sd2 <- 0.5; corr_12 <- -0.4
  Sigma <- matrix(c(sd1^2, rep(prod(corr_12, sd1, sd2), 2), sd2^2), 2, 2)
  x <- t(t(chol(Sigma)) %*% matrix(rnorm(50), 2, 250))
  b <- c(2, 1) / 10
  for (i in 1:length(links)) {
    fam <- binomial(links[i])
    theta <- fam$linkinv(-1 + x %*% b)
    y <- rbinom(length(theta), size = 1, prob = theta)
  
    capture.output(
      fit <- stan_glm(y ~ x, family = fam, seed  = SEED, QR = TRUE,
                    prior = NULL, prior_intercept = NULL,
                    tol_rel_obj = .Machine$double.eps, algorithm = "optimizing")
    )
    expect_stanreg(fit)
    
    val <- coef(fit)
    ans <- coef(glm(y ~ x, family = fam, start = val))
    if (links[i] != "log") expect_equal(val, ans, 0.03, info = links[i])
    else expect_equal(val[-1], ans[-1], 0.06, info = links[i])
  }
})

context("stan_glm (binomial)")
test_that("stan_glm returns expected result for binomial example", {
  # example using simulated data
  N <- 200
  trials <- rpois(N, lambda = 30)
  trials <<- trials
  X <- cbind(1, matrix(rnorm(N * 3, sd = 0.5), N, 3))
  X <<- X
  for (i in 1:length(links)) {
    fam <- binomial(links[i])
    if (i == 4) {
      b <- c(0, 0.5, 0.1, -1.0)
      eta <- X %*% b
      b[1] <- -max(eta) - 0.05
    }  
    else b <- c(0, 0.5, 0.1, -1.0)
    yes <- rbinom(N, size = trials, prob = fam$linkinv(X %*% b))
    y <- cbind(yes, trials - yes)
    capture.output(
      fit <- stan_glm(y ~ X[,-1], family = fam, seed  = SEED, QR = TRUE,
                    prior = NULL, prior_intercept = NULL,
                    tol_rel_obj = .Machine$double.eps, algorithm = "optimizing")
    )
    expect_stanreg(fit)
    
    val <- coef(fit)
    ans <- coef(glm(y ~ X[,-1], family = fam, start = val))
    if (links[i] != "log") expect_equal(val, ans, 0.017, info = links[i])
    else expect_equal(val[-1], ans[-1], 0.008, info = links[i])

    prop <- yes / trials
    capture.output(
      fit2 <- stan_glm(prop ~ X[,-1], weights = trials, family = fam, seed  = SEED,
                     prior = NULL, prior_intercept = NULL,
                     tol_rel_obj = .Machine$double.eps, algorithm = "optimizing")
    )
    expect_stanreg(fit2)
    
    val2 <- coef(fit2)
    if (links[i] != "log") expect_equal(val2, ans, 0.018, info = links[i])
    else expect_equal(val2[-1], ans[-1], 0.01, info = links[i])
  }
})

context("stan_glm (other tests)")
test_that("model with hs prior doesn't error", {
  expect_output(fit <- stan_glm(mpg ~ ., data = mtcars, prior = hs(4, 2, .5), 
                         seed = SEED, algorithm = "meanfield", QR = TRUE), 
                regexp = "Begin stochastic gradient ascent")
  expect_output(print(prior_summary(fit)), "~ hs(df = ", fixed = TRUE)
})

context("stan_glm (other tests)")
# test_that("model with hs_plus prior doesn't error", { # this works except on 32bit Windows 
#   expect_output(fit <- stan_glm(mpg ~ ., data = mtcars, prior = hs_plus(4, 1, 2, .5), 
#                                 seed = SEED, algorithm = "meanfield", QR = TRUE), 
#                 regexp = "Begin stochastic gradient ascent")
#   expect_output(print(prior_summary(fit)), "~ hs_plus(df1 = ", fixed = TRUE)
# })

test_that("model with laplace prior doesn't error", {
  expect_output(fit <- stan_glm(mpg ~ ., data = mtcars, prior = laplace(), 
                         seed = SEED, algorithm = "meanfield", QR = FALSE), 
                regexp = "Begin stochastic gradient ascent")
  expect_output(print(prior_summary(fit)), 
                "~ laplace(", fixed = TRUE)
})

test_that("model with lasso prior doesn't error", {
  expect_output(fit <- stan_glm(mpg ~ ., data = mtcars, prior = lasso(), 
                         seed = SEED, algorithm = "meanfield", QR = FALSE), 
                regexp = "Begin stochastic gradient ascent")
  expect_output(print(prior_summary(fit)), 
                "~ lasso(", fixed = TRUE)
}) 

test_that("model with product_normal prior doesn't error", {
  expect_output(fit <- stan_glm(mpg ~ ., data = mtcars, 
                                prior = product_normal(df = 3, scale = 0.5), 
                                seed = SEED, algorithm = "meanfield", QR = FALSE), 
                regexp = "Begin stochastic gradient ascent")
  expect_output(print(prior_summary(fit)), "~ product_normal(df = ", fixed = TRUE)
})

test_that("prior_aux argument is detected properly", {
  fit <- stan_glm(mpg ~ wt, data = mtcars, iter = 10, chains = 1, seed = SEED, 
                  refresh = -1, prior_aux = exponential(5), 
                  prior = normal(autoscale=FALSE), 
                  prior_intercept = normal(autoscale=FALSE))
  expect_identical(
    fit$prior.info$prior_aux, 
    list(dist = "exponential", 
         location = NULL, scale = NULL, 
         adjusted_scale = 1/5 * sd(mtcars$mpg),
         df = NULL, rate = 5, 
         aux_name = "sigma")
  )
  expect_output(print(prior_summary(fit)), 
                "~ exponential(rate = ", fixed = TRUE)
  expect_output(print(prior_summary(fit)), 
                "**adjusted scale", fixed = TRUE)
})

test_that("prior_aux can be NULL", {
  fit <- stan_glm(mpg ~ wt, data = mtcars, iter = 10, chains = 1, seed = SEED, 
                  refresh = -1, prior_aux = NULL)
  expect_output(print(prior_summary(fit)), 
                "~ flat", fixed = TRUE)
})

test_that("autoscale works (insofar as it's reported by prior_summary)", {
  suppressWarnings(capture.output(
    fit <- stan_glm(mpg ~ wt, data = mtcars, iter = 5, 
                    prior = normal(autoscale=FALSE), 
                    prior_intercept = normal(autoscale=FALSE), 
                    prior_aux = cauchy(autoscale=FALSE)), 
    fit2 <- update(fit, prior = normal())
  ))
  
  out <- capture.output(print(prior_summary(fit)))
  expect_false(any(grepl("adjusted", out)))
  
  expect_output(
    print(prior_summary(fit2)), 
    "**adjusted scale", 
    fixed = TRUE
  )
})
test_that("prior_options is deprecated", {
  expect_warning(
    ops <- prior_options(scaled = FALSE, prior_scale_for_dispersion = 3), 
    "deprecated and will be removed"
  )
  expect_warning(
    capture.output(fit <- stan_glm(mpg ~ wt, data = mtcars, iter = 5, prior_ops = ops)),
    "Setting prior scale for aux to value specified in 'prior_options'"
  )
  expect_output(
    print(prior_summary(fit)), 
    "~ half-cauchy(location = 0, scale = 3)", 
    fixed = TRUE
  )
})

test_that("empty interaction levels dropped", {
  x1 <- gl(3, 5, 100)
  x2 <- gl(4, 6, 100)
  x1[x2 == 1] <- 1
  x1[x2 == 2] <- 1
  y <- rnorm(100)
  expect_warning(stan_glm(y ~ x1*x2, chains = 2, iter = 20, refresh = 0), 
                 regexp = "Dropped empty interaction levels")
})

context("posterior_predict (stan_glm)")

test_that("posterior_predict compatible with glms", {
  source(file.path("helpers", "check_for_error.R"))
  source(file.path("helpers", "expect_linpred_equal.R"))
  SW <- suppressWarnings
  
  check_for_error(fit_gaus)
  expect_linpred_equal(fit_gaus)
  
  mtcars2 <- mtcars
  mtcars2$offs <- runif(nrow(mtcars))
  fit2 <- SW(stan_glm(mpg ~ wt + offset(offs), data = mtcars2,
                      iter = ITER, chains = CHAINS, seed = SEED, refresh = REFRESH))
  expect_warning(posterior_predict(fit2, newdata = mtcars2[1:5, ]), 
                 "offset")
  check_for_error(fit_gaus, data = mtcars2, offset = mtcars2$offs)
  check_for_error(fit2, data = mtcars2, offset = mtcars2$offs)
  expect_linpred_equal(fit_gaus)
  expect_linpred_equal(fit2)
  
  check_for_error(fit_pois)
  check_for_error(fit_negbin)
  expect_linpred_equal(fit_pois)
  expect_linpred_equal(fit_negbin)

  check_for_error(fit_gamma)
  check_for_error(fit_igaus)
  expect_linpred_equal(fit_gamma)
  expect_linpred_equal(fit_igaus)
  
})
