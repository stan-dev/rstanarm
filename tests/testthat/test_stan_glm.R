# tests can be run using devtools::test() or manually by loading testthat 
# package and then running the code below possibly with options(mc.cores = 4).

library(rstanarm)
SEED <- 12345
set.seed(SEED)

context("stan_glm (gaussian)")
test_that("gaussian returns expected result for trees example", {
  # example using trees dataset
  links <- c("identity", "log", "inverse")
  for (i in 1:length(links)) {
    if (links[i] == "inverse") next # unreliable
    fit <- stan_glm(Volume ~ log(Girth) + log(Height), data = trees, 
                    family = gaussian(link = links[i]), algorithm = "optimizing",
                    prior = NULL, prior_intercept = NULL, prior_ops = NULL,
                    QR = TRUE, tol_rel_grad = 1e-16, seed  = SEED)
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
  counts <- c(18,17,15,20,10,20,25,13,12)
  outcome <- gl(3,1,9)
  treatment <- gl(3,3)
  for (i in 1:length(links)) {
    fit <- stan_glm(counts ~ outcome + treatment, family = poisson(links[i]), 
                    prior = NULL, prior_intercept = NULL, prior_ops = NULL, QR = TRUE,
                    algorithm = "optimizing", tol_rel_grad = 1e-16, seed  = SEED)
    ans <- glm(counts ~ outcome + treatment, family = poisson(links[i]), start = coef(fit))
    if (links[i] == "log") expect_equal(coef(fit), coef(ans), tol = 0.01)
    if (links[i] == "identity") expect_equal(coef(fit)[-1], coef(ans)[-1], tol = 0.03)
    if (links[i] == "sqrt") { # this is weird
      if (coef(ans)[1] > 0)
        expect_equal(coef(fit)[-1], coef(ans)[-1], tol = 0.03)
      else
        expect_equal(-coef(fit)[-1], coef(ans)[-1], tol = 0.03)
    }
  }
})

context("stan_glm (negative binomial)")
test_that("stan_glm returns something for glm negative binomial example", {
  # example from MASS::glm.nb
  require(MASS)
  for (i in 1:length(links)) {
    fit1 <- stan_glm(Days ~ Sex/(Age + Eth*Lrn), data = quine, 
                     family = neg_binomial_2(links[i]), 
                     seed  = SEED, chains = 1, iter = 500,
                     prior_PD = TRUE, QR = TRUE)
    fit2 <- stan_glm.nb(Days ~ Sex/(Age + Eth*Lrn), data = quine, 
                        link = links[i],
                        seed  = SEED, chains = 1, iter = 500,
                        prior_PD = TRUE, QR = TRUE)
    expect_is(fit1, "stanreg")
    expect_is(fit2, "stanreg")
    expect_equal(as.matrix(fit1), as.matrix(fit2))
  }
  # testing results against MASS::glm.nb is unreliable
})

context("stan_glm (gaussian)")
test_that("stan_glm returns expected result for cars example", {
  # example using cars dataset
  fit <- stan_glm(log(dist) ~ log(speed), data = cars, 
                  family = gaussian(link = "identity"), seed  = SEED,
                  prior = NULL, prior_intercept = NULL, prior_ops = NULL,
                  tol_rel_obj = .Machine$double.eps, algorithm = "optimizing")
  ans <- glm(log(dist) ~ log(speed), data = cars, family = gaussian(link = "identity"))
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
  
    fit <- stan_glm(y ~ x, family = fam, seed  = SEED, QR = TRUE,
                    prior = NULL, prior_intercept = NULL, prior_ops = NULL,
                    tol_rel_obj = .Machine$double.eps, algorithm = "optimizing")
    val <- coef(fit)
    ans <- coef(glm(y ~ x, family = fam, start = val))
    if (links[i] != "log") expect_equal(val, ans, 0.03)
    else expect_equal(val[-1], ans[-1], 0.06)
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
    fit <- stan_glm(y ~ X[,-1], family = fam, seed  = SEED, QR = TRUE,
                    prior = NULL, prior_intercept = NULL, prior_ops = NULL,
                    tol_rel_obj = .Machine$double.eps, algorithm = "optimizing")
    val <- coef(fit)
    ans <- coef(glm(y ~ X[,-1], family = fam, start = val))
    if (links[i] != "log") expect_equal(val, ans, 0.017, info = links[i])
    else expect_equal(val[-1], ans[-1], 0.008, info = links[i])

    prop <- yes / trials
    fit2 <- stan_glm(prop ~ X[,-1], weights = trials, family = fam, seed  = SEED,
                     prior = NULL, prior_intercept = NULL, prior_ops = NULL,
                     tol_rel_obj = .Machine$double.eps, algorithm = "optimizing")
    val2 <- coef(fit2)
    if (links[i] != "log") expect_equal(val2, ans, 0.018, info = links[i])
    else expect_equal(val2[-1], ans[-1], 0.005, info = links[i])
  }
})
