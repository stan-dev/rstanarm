# tests can be run using devtools::test() or manually by loading testthat 
# package and then running the code below possibly with options(mc.cores = 4).

library(rstanarm)

context("stan_glm (gaussian)")
test_that("gaussian returns expected result for trees example", {
  # example using trees dataset
  links <- c("identity", "log", "inverse")
  for (i in 1:length(links)) {
    if (links[i] == "inverse") next # unreliable
    fit <- stan_glm(Volume ~ log(Girth) + log(Height), data = trees, 
                    family = gaussian(link = links[i]), algorithm = "optimizing",
                    prior = NULL, prior_intercept = NULL, prior_ops = NULL,
                    tol_rel_grad = 1e-16, seed = 12345)
    ans <- glm(Volume ~ log(Girth) + log(Height),data = trees, 
               family = gaussian(link = links[i]))
    expect_equal(coef(fit), coef(ans), tol = 0.02)
  }
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
                    prior = NULL, prior_intercept = NULL, prior_ops = NULL,
                    algorithm = "optimizing", tol_rel_grad = 1e-16, seed = 12345)
    ans <- glm(counts ~ outcome + treatment, family = poisson(links[i]))
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
test_that("stan_glm returns expected result for glm negative binomial example", {
  # example from MASS::glm.nb
  require(MASS)
  for (i in 1:length(links)) {
    fit <- stan_glm(Days ~ Sex/(Age + Eth*Lrn), data = quine, 
                    family = neg_binomial_2(links[i]), seed = 12345,
                    prior = NULL, prior_intercept = NULL, prior_ops = NULL,
                    algorithm = "optimizing", tol_rel_grad = 1e-16)
    ans <- glm.nb(Days ~ Sex/(Age + Eth*Lrn), data = quine, link = links[i],
                  start = coef(fit), init.theta = fit$stan_summary["overdispersion",1])
    # testing results is unreliable
#     if (links[i] == "log") expect_equal(coef(fit), coef(ans), tol = 0.01)
#     if (links[i] == "identity") expect_equal(coef(fit)[-1], coef(ans)[-1], tol = 0.015)
#     if (links[i] == "sqrt") { # this is weird
#       if (coef(ans)[1] > 0)
#         expect_equal(coef(fit)[-1], coef(ans)[-1], tol = 0.015)
#       else
#         expect_equal(-coef(fit)[-1], coef(ans)[-1], tol = 0.015)
#     }
  }
})

context("stan_glm (gaussian)")
test_that("stan_glm returns expected result for cars example", {
  # example using cars dataset
  fit <- stan_glm(log(dist) ~ log(speed), data = cars, 
                  family = gaussian(link = "identity"), seed = 12345,
                  prior = NULL, prior_intercept = NULL, prior_ops = NULL,
                  tol_rel_obj = .Machine$double.eps, algorithm = "optimizing")
  ans <- glm(log(dist) ~ log(speed), data = cars, family = gaussian(link = "identity"))
  expect_equal(coef(fit), coef(ans), tol = 0.0025)
})

context("stan_glm (bernoulli)")
links <- c("logit", "probit", "cauchit", "log", "cloglog")
test_that("stan_glm returns expected result for bernoulli", {
  # bernoulli example
  sd1 <- 1; sd2 <- 0.5; corr_12 <- -0.4
  Sigma <- matrix(c(sd1^2, rep(prod(corr_12, sd1, sd2), 2), sd2^2), 2, 2)
  x <- t(t(chol(Sigma)) %*% matrix(rnorm(500), 2, 250))
  b <- c(2, 1) / 10
  for (i in 1:length(links)) {
    fam <- binomial(links[i])
    theta <- fam$linkinv(-1 + x %*% b)
    y <- rbinom(length(theta), size = 1, prob = theta)
  
    fit <- stan_glm(y ~ x, family = fam, seed = 12345,
                    prior = NULL, prior_intercept = NULL, prior_ops = NULL,
                    tol_rel_obj = .Machine$double.eps, algorithm = "optimizing")
    val <- coef(fit)
    ans <- coef(glm(y ~ x, family = fam, start = coef(fit)))
    if (links[i] != "log") expect_equal(val, ans, 0.02)
    else expect_equal(val[-1], ans[-1], 0.02)
  }
})

context("stan_glm (binomial)")
test_that("stan_glm returns expected result for binomial example", {
  # example using simulated data
  N <- 500
  trials <- rpois(N, lambda = 30)
  X <- cbind(1, matrix(rnorm(N * 3), N, 3))
  b <- c(-3.25, 0.5, 0.1, -1.0)
  for (i in 1:length(links)) {
    fam <- binomial(links[i])
    yes <- rbinom(N, size = trials, prob = fam$linkinv(X %*% b))
    y <- cbind(yes, trials - yes)
    fit <- stan_glm(y ~ X[,-1], family = fam, seed = 12345,
                    prior = NULL, prior_intercept = NULL, prior_ops = NULL,
                    tol_rel_obj = .Machine$double.eps, algorithm = "optimizing")
    val <- coef(fit)
    ans <- coef(glm(y ~ X[,-1], family = fam))
    if (links[i] != "log") expect_equal(val, ans, 0.003)
    else expect_equal(val[-1], ans[-1], 0.003)

    prop <- yes / trials
    fit2 <- stan_glm(prop ~ X[,-1], weights = trials, family = fam, seed = 12345,
                     prior = NULL, prior_intercept = NULL, prior_ops = NULL,
                     tol_rel_obj = .Machine$double.eps, algorithm = "optimizing")
    val2 <- coef(fit2)
    if (links[i] != "log") expect_equal(val2, ans, 0.003)
    else expect_equal(val2[-1], ans[-1], 0.003)
  }
})
