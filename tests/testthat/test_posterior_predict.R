# tests can be run using devtools::test() or manually by loading testthat 
# package and then running the code below possibly with options(mc.cores = 4).

library(rstanarm)
SEED <- 123
ITER <- 10
CHAINS <- 2
CORES <- 1

# These tests just make sure that posterior_predict doesn't throw errors and
# that result has correct dimensions
check_for_error <- function(fit) {
  nsims <- nrow(as.matrix(fit))
  
  expect_silent(yrep1 <- posterior_predict(fit))
  expect_equal(dim(yrep1), c(nsims, nobs(fit)))

  expect_silent(yrep2 <- posterior_predict(fit, draws = 1))
  expect_equal(dim(yrep2), c(1, nobs(fit)))
  
  expect_silent(yrep3 <- posterior_predict(fit, newdata = model.frame(fit)[1,]))
  expect_equal(dim(yrep3), c(nsims, 1))
  
  expect_silent(yrep4 <- posterior_predict(fit, draws = 2, newdata = model.frame(fit)[1,]))
  expect_equal(dim(yrep4), c(2, 1))
  
  expect_silent(yrep5 <- posterior_predict(fit, newdata = model.frame(fit)[1:5,]))
  expect_equal(dim(yrep5), c(nsims, 5))
  
  expect_silent(yrep5 <- posterior_predict(fit, draws = 3, newdata = model.frame(fit)[1:5,]))
  expect_equal(dim(yrep5), c(3, 5))
}

context("posterior_predict (stan_lm)")
test_that("posterior_predict compatible with stan_lm", {
  fit <- stan_lm(mpg ~ wt + cyl + am, data = mtcars, prior = R2(0.5), 
                 iter = ITER, chains = CHAINS, cores = CORES, seed = SEED)
  check_for_error(fit)
})

context("posterior_predict (stan_glm)")
test_that("posterior_predict compatible with gaussian glm", {
  fit <- stan_glm(mpg ~ wt, data = mtcars, 
                  iter = ITER, chains = CHAINS, cores = CORES, seed = SEED)
  check_for_error(fit)
})
test_that("posterior_predict compatible with poisson & negbin glm", {
  counts <- c(18,17,15,20,10,20,25,13,12)
  outcome <- gl(3,1,9)
  treatment <- gl(3,3)
  fit <- stan_glm(counts ~ outcome + treatment, family = poisson(), 
                  iter = ITER, chains = CHAINS, cores = CORES, seed = SEED)
  fitnb <- update(fit, family = neg_binomial_2)
  check_for_error(fit)
  check_for_error(fitnb)
})

context("posterior_predict (stan_polr)")
test_that("posterior_predict compatible with stan_polr", {
  fit <- stan_polr(tobgp ~ agegp + alcgp, data = esoph, 
                   prior = R2(location = 0.4),
                   iter = ITER, chains = CHAINS, cores = CORES, seed = SEED)
  check_for_error(fit)
})

context("posterior_predict (stan_(g)lmer)")
test_that("posterior_predict compatible with stan_lmer", {
  fit <- stan_lmer(mpg ~ wt + (1|cyl) + (1 + wt|gear), data = mtcars, 
                   prior = normal(0,1), 
                   iter = ITER, chains = CHAINS, cores = CORES, seed = SEED)
  check_for_error(fit)
})
test_that("posterior_predict compatible with stan_glmer (binomial)", {
  check_for_error(example_model)
})

