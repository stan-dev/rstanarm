# tests can be run using devtools::test() or manually by loading testthat 
# package and then running the code below possibly with options(mc.cores = 4).

library(rstanarm)
SEED <- 123
ITER <- 10
CHAINS <- 2
CORES <- 1

fit <- example_model

context("plot.stanreg")
test_that("plot method doesn't throw bad errors and creates ggplot objects", {
  expect_default_message <- function(fit, message = "default", ...) {
    expect_message(p <- plot(fit, ...), regexp = message)
  }
  
  plotters1 <- paste0("stan_", c("plot", "trace", "hist", "dens", "ac")) # scat
  plotters2 <- paste0("stan_", c("rhat", "ess", "mcse")) # diag, par
  
  for (j in seq_along(plotters1)) {
    expect_message(p <- plot(fit, plotfun = plotters1[j]), regexp = "default")
    expect_is(p, "ggplot")
    expect_is(p <- plot(fit, plotfun = plotters1[j], pars = "beta"), "ggplot")
  }
  for (j in seq_along(plotters2)) {
    expect_silent(p <- plot(fit, plotfun = plotters2[j]))
    expect_is(p, "ggplot")
  }
})
