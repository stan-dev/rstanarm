# tests can be run using devtools::test() or manually by loading testthat 
# package and then running the code below possibly with options(mc.cores = 4).

library(rstanarm)
library(lme4)
SEED <- 123
set.seed(SEED)
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
test_that("compatible with gaussian glm", {
  fit <- stan_glm(mpg ~ wt, data = mtcars, 
                  iter = ITER, chains = CHAINS, cores = CORES, seed = SEED)
  check_for_error(fit)
})
test_that("compatible with poisson & negbin glm", {
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
test_that("compatible with stan_polr", {
  fit <- stan_polr(tobgp ~ agegp + alcgp, data = esoph, 
                   prior = R2(location = 0.4),
                   iter = ITER, chains = CHAINS, cores = CORES, seed = SEED)
  check_for_error(fit)
})

context("posterior_predict (stan_(g)lmer)")
test_that("compatible with stan_lmer", {
  fit <- stan_lmer(mpg ~ wt + (1|cyl) + (1 + wt|gear), data = mtcars, 
                   prior = normal(0,1), 
                   iter = ITER, chains = CHAINS, cores = CORES, seed = SEED)
  check_for_error(fit)
})
test_that("compatible with stan_glmer (binomial)", {
  check_for_error(example_model)
})
test_that("compatible with stan_glmer with transformation in formula", {
  fit <- stan_lmer(mpg ~ log1p(wt) + (1|cyl) + (1|gear), data = mtcars, 
                   iter = ITER, chains = CHAINS, cores = CORES, seed = SEED)
  expect_silent(yrep1 <- posterior_predict(fit))
  expect_silent(yrep1 <- posterior_predict(fit, newdata = mtcars[1:5, ]))
})


context("posterior_predict (optimizing and vb)")
test_that("errors for optimizing and silent for vb", {
  fit <- stan_glm(mpg ~ wt + cyl + am, data = mtcars, algorithm = "optimizing")
  fit2 <- update(fit, algorithm = "meanfield")
  fit3 <- update(fit, algorithm = "fullrank")
  expect_error(posterior_predict(fit), regexp = "optimizing")
  expect_silent(posterior_predict(fit2))
  expect_silent(posterior_predict(fit3))
})


context("posterior_predict (compare to lme4)")
test_that("posterior_predict close to predict.merMod", {
  mod1 <- as.formula(mpg ~ wt + (1|cyl) + (1|gear))
  mod2 <- as.formula(mpg ~ log1p(wt) + (1|cyl))
  mod3 <- as.formula(mpg ~ wt + (1|cyl) + (1 + wt|gear))

  lfit1 <- lmer(mod1, data = mtcars)
  sfit1 <- stan_glmer(mod1, data = mtcars, cores = CORES, chains = CHAINS, 
                      iter = 400, seed = SEED)
  lfit2 <- update(lfit1, formula = mod2)
  sfit2 <- update(sfit1, formula = mod2)
  lfit3 <- update(lfit1, formula = mod3)
  sfit3 <- update(sfit1, formula = mod3)
  
  nd <- nd2 <- mtcars[1:5, ]
  nd2$cyl[2] <- 5 # new level
  nd3 <- nd2
  nd3$gear[2] <- 7
  nd3$gear[5] <- 1
  
  for (j in 1:3) {
    expect_equal(
      predict(get(paste0("sfit", j))),
      unname(predict(get(paste0("lfit", j)))),
      tol = 0.5)
    expect_equal(
      predict(get(paste0("sfit", j)), newdata = nd),
      predict(get(paste0("lfit", j)), newdata = nd),
      tol = 0.5)
    expect_equal(
      colMeans(posterior_predict(get(paste0("sfit", j)), newdata = nd)),
      unname(predict(get(paste0("lfit", j)), newdata = nd)),
      tol = 0.5)
    expect_equal(
      colMeans(posterior_predict(get(paste0("sfit", j)), newdata = nd2, 
                                 allow.new.levels = TRUE)),
      unname(predict(lfit1, newdata = nd2, allow.new.levels = TRUE)),
      tol = 0.5)
    expect_equal(
      colMeans(posterior_predict(get(paste0("sfit", j)), newdata = nd3, 
                                 allow.new.levels = TRUE)),
      unname(predict(lfit1, newdata = nd3, allow.new.levels = TRUE)),
      tol = 0.5)
    
    expect_error(
      posterior_predict(get(paste0("sfit", j)), newdata = nd3, 
                        allow.new.levels = FALSE),
      regexp = "new levels", ignore.case = TRUE)
  }
})

