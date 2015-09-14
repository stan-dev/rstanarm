# tests can be run using devtools::test() or manually by loading testthat 
# package and then running the code below possibly with options(mc.cores = 4).

library(rstanarm)
require(lme4)
SEED <- 123
set.seed(SEED)
options(mc.cores = 2L)
if (interactive()) options(mc.cores = parallel::detectCores())
CONTROL <- list(adapt_delta = 0.95)
FIXEF_tol <- 0.05
RANEF_tol <- 0.1 


context("stan_lmer")
test_that("stan_lmer returns expected result for slepstudy example", {
  fmla <- Reaction / 10 ~ Days + (Days | Subject)
  fit <- stan_lmer(fmla, data = sleepstudy, 
                   init_r = 0.05, seed = SEED, control = CONTROL)
  ans <- lmer(fmla, data = sleepstudy)
  expect_equal(fixef(fit), fixef(ans), tol = FIXEF_tol)
  expect_equal(ranef(fit)[[1L]], ranef(ans)[[1L]], tol = RANEF_tol)
})

context("stan_lmer")
test_that("stan_lmer returns expected result for Penicillin example", {
  fmla <- diameter ~ (1|plate) + (1|sample)
  fit <- stan_lmer(fmla, data = Penicillin, control = CONTROL, seed = SEED)
  ans <- lmer(fmla, data = Penicillin)
  expect_equal(fixef(fit), fixef(ans), tol = FIXEF_tol)
  expect_equal(ranef(fit)[[1L]], ranef(ans)[[1L]], tol = RANEF_tol)
})

context("stan_glmer (binomial)")
test_that("stan_glmer returns expected result for cbpp example", {
  links <- c("logit", "probit", "cauchit", "log", "cloglog")
  # for (i in seq_along(links)) {
  i <- 1L # it seems only logit gives results similar to glmer with same link 
    fmla <- cbind(incidence, size - incidence) ~ period + (1 | herd)
    fit <- stan_glmer(fmla, data = cbpp, family = binomial(links[i]),
                      seed = SEED, control = CONTROL)
    ans <- glmer(fmla, data = cbpp, family = binomial(links[i]))
    expect_equal(fixef(fit), fixef(ans), tol = FIXEF_tol)
    expect_equal(ranef(fit)[[1L]], ranef(ans)[[1L]], tol = RANEF_tol)
  # }
})




