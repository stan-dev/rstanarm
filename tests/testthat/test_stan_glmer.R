# tests can be run using devtools::test() or manually by loading testthat 
# package and then running the code below possibly with options(mc.cores = 4).

library(rstanarm)
require(lme4)
options(mc.cores = 2L)
if (interactive()) options(mc.cores = parallel::detectCores())
set.seed(123)

threshold <- 0.3

f1 <- function(x) cbind(coef(x), se(x))
f2 <- function(x) summary(x)$coefficients[,1:2]

context("stan_glmer")
test_that("stan_glmer returns expected result for slepstudy example", {
  # example using sleepstudy dataset
  fit <- stan_lmer(Reaction / 10 ~ Days + (Days | Subject), sleepstudy, init_r = 0.05)
  fm1 <- lmer(Reaction / 10 ~ Days + (Days | Subject), sleepstudy)
})
