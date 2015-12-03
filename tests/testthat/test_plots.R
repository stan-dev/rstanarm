# tests can be run using devtools::test() or manually by loading testthat 
# package and then running the code below possibly with options(mc.cores = 4).

library(rstanarm)
SEED <- 123
ITER <- 10
CHAINS <- 2
CORES <- 1

fit <- example_model
fito <- stan_glm(mpg ~ ., data = mtcars, algorithm = "optimizing")

context("plot.stanreg helpers")
test_that(".grep_for_pars works", {
  .grep_for_pars <- rstanarm:::.grep_for_pars
  .bnames <- rstanarm:::.bnames
  
  all_period <- paste0("period", 2:4)
  all_varying <- .bnames(rownames(fit$stan_summary), value = TRUE)
  
  expect_equal(.grep_for_pars(fit, "period"), all_period)
  expect_equal(.grep_for_pars(fit, c("period", "size")), c(all_period, "size"))
  expect_equal(.grep_for_pars(fit, "period|size"), c("size", all_period))
  expect_equal(.grep_for_pars(fit, "(2|3)$"), all_period[1:2])
  expect_equal(.grep_for_pars(fit, "herd"), all_varying)
  expect_equal(.grep_for_pars(fit, "b\\["), all_varying)
  expect_equal(.grep_for_pars(fit, "Intercept"), c("(Intercept)", all_varying))
  expect_equal(.grep_for_pars(fit, "herd:[3,5]"), all_varying[c(3,5)])
  expect_equal(.grep_for_pars(fit, "herd:[3-5]"), all_varying[3:5])
  expect_error(.grep_for_pars(fit, "NOT A PARAMETER"), regexp = "No matches")
  expect_error(.grep_for_pars(fit, "b["))
})

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
    expect_is(p <- plot(fit, plotfun = plotters1[j], pars = "alpha", 
                        regex_pars = "period"), "ggplot")
    expect_is(p <- plot(fit, plotfun = plotters1[j], regex_pars = "period"), 
              "ggplot")
  }
  for (j in seq_along(plotters2)) {
    expect_silent(p <- plot(fit, plotfun = plotters2[j]))
    expect_is(p, "ggplot")
  }
})

test_that("plot.stanreg ok for optimization", {
  expect_silent(plot(fito))
  expect_silent(plot(fito, pars = "alpha"))
  expect_silent(plot(fito, pars = "beta"))
  expect_silent(plot(fito, pars = c("wt", "gear")))
  expect_warning(plot(fito, regex_pars = "wt"), regexp = "'regex_pars' ignored")
  expect_warning(plot(fito, "trace"), regexp = "'plotfun' ignored")
})

context("pairs.stanreg")
test_that("pairs method ok", {
  requireNamespace("rstan")
  requireNamespace("KernSmooth")
  expect_silent(pairs(fit, pars = c("period2", "log-posterior")))
  expect_error(pairs(fito), regexp = "only available for models fit using MCMC")
})
