# The tests in this file check that for very simple models with default priors 
# the point estimates and standard error estimates are similar to the results 
# from the corresponding R function (lm or glm).

# Tests can be run using devtools::test() or manually by loading testthat 
# package and then running the code below.

library(rstanarm)
set.seed(123)

threshold <- 0.1

f1 <- function(x) cbind(coef(x), se(x))
f2 <- function(x) summary(x)$coefficients[,1:2]

context("stan_lm")
test_that("stan_lm returns expected result for simulated example", {
  # example using fake data
  N <- 100
  X <- cbind(rnorm(N), rnorm(N))
  b <- c( -1, .1)
  a <- .5
  y <- a + X %*% b + rnorm(N)
  fit <- stan_lm(y ~ X, iter = 400, seed = 123)
  val <- f1(fit)
  ans <- f2(lm(y ~ X))
  diff <- abs(val - ans)
  expect_true(all(diff < threshold))
})

context("stan_lm and stan_glm")
test_that("gaussian(link = 'log') returns expected result for trees example", {
  # example using trees dataset
  fit1 <- stan_glm(Volume ~ log(Girth) + log(Height), data = trees, 
                 family = gaussian(link = "log"), iter = 400, seed = 123)
  fit2 <- stan_lm(log(Volume) ~ log(Girth) + log(Height), data = trees, 
                  iter = 400, seed = 123)
  val1 <- f1(fit1); val2 <- f1(fit2)
  ans <- f2(lm(log(Volume) ~ log(Girth) + log(Height),data = trees))
  diff1 <- abs(val1 - ans); diff2 <- abs(val2 - ans)
  expect_true(all(diff1 < threshold & diff2 < threshold))
})


context("stan_glm (poisson)")
test_that("stan_glm returns expected result for glm poisson example", {
  # example from help("glm")
  counts <- c(18,17,15,20,10,20,25,13,12)
  outcome <- gl(3,1,9)
  treatment <- gl(3,3)
  fit <- stan_glm(counts ~ outcome + treatment, family = poisson(), 
                  iter = 400, seed = 123)
  val <- f1(fit)
  ans <- f2(glm(counts ~ outcome + treatment, family = poisson()))
  diff <- abs(val - ans)
  expect_true(all(diff < threshold))
})

context("stan_glm (gaussian)")
test_that("stan_glm returns expected result for cars example", {
  # example using cars dataset
  fit <- stan_glm(log(dist) ~ log(speed), data = cars, 
                  family = gaussian(link = "identity"), 
                  iter = 400, seed = 123)
  val <- f1(fit)
  ans <- f2(glm(log(dist) ~ log(speed), data = cars, family = gaussian(link = "identity")))
  diff <- abs(val - ans)
  testthat::expect_true(all(diff < threshold))
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
