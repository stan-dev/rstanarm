library(rstanarm)
library(loo)
options(loo.cores = 1)
SEED <- 1234

# These tests just check that the loo.stanreg method (which calls loo.function
# method) results are identical to the loo.matrix results. Since for these tests 
# the log-likelihood matrix is computed using the log-likelihood function, the 
# only thing these tests really do is make sure that loo.stanreg and all the 
# log-likelihood functions don't return any errors and whatnot (it does not check
# that the results returned by loo are actually correct). 

loo_with_fn <- function(fit) loo(fit)
loo_with_mat <- function(fit) loo(log_lik(fit))
expect_identical_loo <- function(fit) {
  expect_identical(loo_with_fn(fit), loo_with_mat(fit))
}


context("loo")

test_that("loo for gaussian works", {
  fit_gaus <- stan_glm(mpg ~ wt, data = mtcars, chains = 2, iter = 50, seed = SEED)
  expect_identical_loo(fit_gaus)
})
test_that("loo for binomial works", {
  dat <- data.frame(ldose = rep(0:5, 2), sex = factor(rep(c("M", "F"), c(6, 6))))
  numdead <- c(1, 4, 9, 13, 18, 20, 0, 2, 6, 10, 12, 16)
  SF <- cbind(numdead, numalive = 20-numdead)
  fit_binom <- stan_glm(SF ~ sex*ldose, data = dat, family = binomial, 
                        chains = 2, iter = 50, seed = SEED)
  expect_identical_loo(fit_binom)
})
test_that("loo for poisson works", {
  d.AD <- data.frame(treatment = gl(3,3), outcome =  gl(3,1,9), 
                     counts = c(18,17,15,20,10,20,25,13,12))
  fit_pois <- stan_glm(counts ~ outcome + treatment, data = d.AD, 
                  family = poisson, chains = 2, iter = 50, seed = SEED)
  expect_identical_loo(fit_pois)
})
test_that("loo for gamma works", {
  clotting <- data.frame(u = c(5,10,15,20,30,40,60,80,100),
                         lot1 = c(118,58,42,35,27,25,21,19,18),
                         lot2 = c(69,35,26,21,18,16,13,12,12))
  fit_gamma <- stan_glm(lot1 ~ log(u), data = clotting, family = Gamma, 
                        chains = 2, iter = 50, seed = SEED)
  expect_identical_loo(fit_gamma)
})
test_that("loo for stan_polr (logistic) works", {
  fit_polr <- stan_polr(tobgp ~ agegp, data = esoph, prior = R2(0.2, "mean"),  
                   chains = 2, iter = 50, init_r = 0.1, seed = SEED)
  expect_identical_loo(fit_polr)
})
test_that("loo for stan_lm works", {
  fit_lm <- stan_lm(mpg ~ ., data = mtcars, prior = R2(0.75), 
                    chains = 2, iter = 50, seed = SEED)
  expect_identical_loo(fit_lm)
})
