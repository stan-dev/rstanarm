# tests can be run using devtools::test() or manually by loading testthat 
# package and then running the code below possibly with options(mc.cores = 4).

library(rstanarm)
SEED <- 123

threshold <- 0.1

plink <- function(fit, nd = NULL) 
  predict(fit, newdata = nd, type = "link", se.fit = TRUE)
presp <- function(fit, nd = NULL) 
  predict(fit, newdata = nd, type = "response", se.fit = TRUE)
get_diffs <- function(glm_preds, stan_preds) {
  fit_diff <- glm_preds$fit - stan_preds$fit
  se_diff <- glm_preds$se.fit - stan_preds$se.fit
  abs(cbind(fit_diff, se_diff))
}


context("predict")
test_that("predict ok for binomial", {
  # example from help(predict.glm)
  ldose <- rep(0:5, 2)
  numdead <- c(1, 4, 9, 13, 18, 20, 0, 2, 6, 10, 12, 16)
  sex <- factor(rep(c("M", "F"), c(6, 6)))
  SF <- cbind(numdead, numalive = 20-numdead)
  
  glmfit <- glm(SF ~ sex*ldose, family = binomial)
  stanfit <- stan_glm(SF ~ sex*ldose, family = binomial, 
                      chains = 2, iter = 50, seed = SEED)
  pg <- plink(glmfit)
  ps <- plink(stanfit)
  expect_equal(pg$fit, ps$fit, tol = 0.2)
  expect_equal(pg$se.fit, ps$se.fit, tol = 0.2)
  expect_error(presp(stanfit))
  
  ld <- seq(0, 5, 0.1)
  newd <- data.frame(ldose = ld, sex = factor(rep("M", length(ld)), levels = levels(sex)))
  pg <- plink(glmfit, newd)
  ps <- plink(stanfit, newd)
  expect_equal(pg$fit, ps$fit, tol = 0.2)
  expect_equal(pg$se.fit, ps$se.fit, tol = 0.15)
})

test_that("predict ok for gaussian", {
  glmfit <- glm(mpg ~ wt, data = mtcars)
  stanfit <- stan_glm(mpg ~ wt, data = mtcars, 
                      chains = 2, iter = 400, seed = SEED)
  pg <- plink(glmfit)
  ps <- plink(stanfit)
  expect_equal(pg$fit, ps$fit, tol = 0.05)
  expect_equal(pg$se.fit, ps$se.fit, tol = 0.1)
  expect_error(presp(stanfit))

  newd <- data.frame(wt = c(1,5))
  pg <- plink(glmfit, newd)
  ps <- plink(stanfit, newd)
  expect_equal(pg$fit, ps$fit, tol = 0.1)
  expect_equal(pg$se.fit, ps$se.fit, tol = 0.1)
})

test_that("predict ok for Poisson", {
  counts <- c(18,17,15,20,10,20,25,13,12)
  outcome <- gl(3,1,9)
  treatment <- gl(3,3)

  glmfit <- glm(counts ~ outcome + treatment, family = poisson())
  stanfit <- stan_glm(counts ~ outcome + treatment, family = poisson(), 
                      chains = 2, iter = 50, seed = SEED)
  pg <- plink(glmfit)
  ps <- plink(stanfit)
  expect_equal(pg$fit, ps$fit, tol = 0.05)
  expect_equal(pg$se.fit, ps$se.fit, tol = 0.2)
  expect_error(presp(stanfit))
})

