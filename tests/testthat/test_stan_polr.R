# tests can be run using devtools::test() or manually by loading testthat 
# package and then running the code below possibly with options(mc.cores = 4).

library(rstanarm)
SEED <- 123

threshold <- 0.03

context("stan_polr")
test_that("stan_polr returns expected result for esoph example", {
  library(MASS)
  f <- tobgp ~ agegp + alcgp
  fit <- stan_polr(f, data = esoph, prior = R2(location = 0.4),
                   chains = 2, iter = 400, seed = SEED)
  # fit <- stan_polr(f, data = esoph, prior = NULL, 
  #                  algorithm = "fullrank", seed = SEED)
  # check <- polr(f, data = esoph)
  # expect_equal(coef(fit), coef(check), threshold)
})
