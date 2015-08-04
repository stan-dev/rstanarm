# tests can be run using devtools::test() or manually by loading testthat 
# package and then running the code below possibly with options(mc.cores = 4).

library(rstanarm)
set.seed(123)

threshold <- 0.3

f1 <- function(x) cbind(coef(x), se(x))
f2 <- function(x) summary(x)$coefficients[,1:2]

context("stan_polr")
test_that("stan_polr returns expected result for mtcars example", {
  # example using mtcars dataset
  f <- as.ordered(cyl) ~ vs + am + carb
  fit <- stan_polr(f, data = mtcars, prior = LKJ(location = 0.5))
  fit <- stan_polr(f, data = mtcars, prior = NULL, algorithm = "opt")
  check <- MASS::polr(f, data = mtcars)
  expect_true(all(abs(coef(fit) - coef(check)) < threshold))
})
