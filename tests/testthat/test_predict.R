# tests can be run using devtools::test() or manually by loading testthat 
# package and then running the code below possibly with options(mc.cores = 4).

library(rstanarm)
set.seed(123)

threshold <- 0.3

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
  stanfit <- stan_glm(SF ~ sex*ldose, family = binomial, iter = 400, seed = 123)
  diffs_link <- get_diffs(plink(glmfit), plink(stanfit))
  diffs_resp <- get_diffs(presp(glmfit), presp(stanfit))
  expect_true(all(diffs_link < threshold))
  expect_true(all(diffs_resp < threshold))
  
  ld <- seq(0, 5, 0.1)
  newd <- data.frame(ldose = ld, sex = factor(rep("M", length(ld)), levels = levels(sex)))
  diffs_link <- get_diffs(plink(glmfit, newd), plink(stanfit, newd))
  diffs_resp <- get_diffs(presp(glmfit, newd), presp(stanfit, newd))
  expect_true(all(diffs_link < threshold))
  expect_true(all(diffs_resp < threshold))
})
