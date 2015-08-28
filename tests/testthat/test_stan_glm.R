# tests can be run using devtools::test() or manually by loading testthat 
# package and then running the code below possibly with options(mc.cores = 4).

library(rstanarm)

context("stan_glm (gaussian)")
test_that("gaussian returns expected result for trees example", {
  # example using trees dataset
  links <- c("identity", "log", "inverse")
  for (i in 1:length(links)) {
    fit <- stan_glm(Volume ~ log(Girth) + log(Height), data = trees, 
                    family = gaussian(link = links[i]), algorithm = "optimizing",
                    prior = NULL, prior.for.intercept = NULL, prior.options = NULL,
                    tol_rel_grad = 1e-16)
    ans <- glm(Volume ~ log(Girth) + log(Height),data = trees, 
               family = gaussian(link = links[i]))
    expect_equal(coef(fit), coef(ans), tol = 2e-5)
  }
})


context("stan_glm (poisson)")
test_that("stan_glm returns expected result for glm poisson example", {
  # example from help("glm")
  counts <- c(18,17,15,20,10,20,25,13,12)
  outcome <- gl(3,1,9)
  treatment <- gl(3,3)
  links <- c("log", "identity", "sqrt")
  for (i in 1:length(links)) {
    fit <- stan_glm(counts ~ outcome + treatment, family = poisson(links[i]), 
                    prior = NULL, prior.for.intercept = NULL, prior.options = NULL,
                    algorithm = "optimizing", tol_rel_grad = 1e-16)
    ans <- glm(counts ~ outcome + treatment, family = poisson(links[i]))
    expect_equal(coef(fit), coef(ans), tol = 2e-4)
  }
})

context("stan_glm (negative binomial)")
test_that("stan_glm returns expected result for glm negative binomial example", {
  # example from MASS::glm.nb
  links <- c("log", "identity", "sqrt")
  require(MASS)
  for (i in 1:length(links)) {
    fit <- stan_glm(Days ~ Sex/(Age + Eth*Lrn), data = quine, 
                    family = neg_binomial_2(links[i]),
                    prior = NULL, prior.for.intercept = NULL, prior.options = NULL,
                    algorithm = "optimizing", tol_rel_grad = 1e-16)
    ans <- glm.nb(Days ~ Sex/(Age + Eth*Lrn), data = quine, link = links[i])
    expect_equal(coef(fit), coef(ans), tol = 0.0011)
  }
})

context("stan_glm (gaussian)")
test_that("stan_glm returns expected result for cars example", {
  # example using cars dataset
  fit <- stan_glm(log(dist) ~ log(speed), data = cars, 
                  family = gaussian(link = "identity"), 
                  prior = NULL, prior.for.intercept = NULL, prior.options = NULL,
                  tol_rel_obj = .Machine$double.eps, algorithm = "optimizing")
  ans <- glm(log(dist) ~ log(speed), data = cars, family = gaussian(link = "identity"))
  expect_equal(coef(fit), coef(ans), tol = 1e-6)
})

context("stan_glm (bernoulli)")
test_that("stan_glm returns expected result for bernoulli (logit and probit)", {
  # bernoulli example
  sd1 <- 1; sd2 <- 3; corr_12 <- -0.4
  Sigma <- matrix(c(sd1^2, rep(prod(corr_12, sd1, sd2), 2), sd2^2), 2, 2)
  x <- t(t(chol(Sigma)) %*% matrix(rnorm(500), 2, 250))
  b <- c(-0.5, 1)
  theta <- 1/(1 + exp(-x %*% b))
  y <- rbinom(length(theta), size = 1, prob = theta)
  
  fit <- stan_glm(y ~ x, family = "binomial", iter = 400, seed = 123)
  fit2 <- stan_glm(y ~ x, family = binomial(link = "probit"), iter = 400, 
                   seed = 123)
  val <- f1(fit) 
  val2 <- f1(fit2)
  ans <- f2(glm(y ~ x, family = "binomial"))
  ans2 <- f2(glm(y ~ x, family = binomial(link = "probit"))) 
  diff <- abs(val - ans)
  diff2 <- abs(val2 - ans2)
  testthat::expect_true(all(diff < threshold))
  testthat::expect_true(all(diff2 < threshold))
})

context("stan_glm (binomial)")
test_that("stan_glm returns expected result for binomial example", {
  # example using simulated data
  N <- 50
  trials <- rpois(N, lambda = 30)
  X <- cbind(1, matrix(rnorm(N * 3), N, 3))
  b <- c(-0.5, 0.5, 0.1, -0.75)
  yes <- rbinom(N, size = trials, prob = 1 / (1 + exp(- X %*% b)))
  y <- cbind(yes, trials - yes)
  X <- X[,-1]
  fit <- stan_glm(y ~ X, family = binomial, iter = 400, seed = 123)
  val <- f1(fit)
  ans <- f2(glm(y ~ X, family = binomial))
  diff <- abs(val - ans)
  expect_true(all(diff < threshold))
  
  prop <- yes / trials
  fit2 <- stan_glm(prop ~ X, weights = trials, family = binomial, 
                   iter = 400, seed = 123)
  val2 <- f1(fit2)
  diff2 <- abs(val2 - ans)
  expect_true(all(diff2 < threshold))
})
