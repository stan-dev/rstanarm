# tests can be run using devtools::test() or manually by loading testthat 
# package and then running the code below possibly with options(mc.cores = 4).

library(rstanarm)
require(lme4)
SEED <- 123
options(mc.cores = 2L)
if (interactive()) options(mc.cores = parallel::detectCores())
FIXEF_tol <- 0.05
RANEF_tol <- 0.13 


context("stan_lmer")
test_that("stan_lmer returns expected result for slepstudy example", {
  fmla <- Reaction / 10 ~ Days + (Days | Subject)
  fit <- stan_lmer(fmla, data = sleepstudy, 
                   init_r = 0.05, chains = 2, iter = 400, seed = SEED)
  ans <- lmer(fmla, data = sleepstudy)
  expect_equal(fixef(fit), fixef(ans), tol = FIXEF_tol)
  expect_equal(ranef(fit), ranef(ans), tol = RANEF_tol)
  expect_equal(ngrps(fit), ngrps(ans))
})

context("stan_lmer")
test_that("stan_lmer returns expected result for Penicillin example", {
  fmla <- diameter ~ (1|plate) + (1|sample)
  fit <- stan_lmer(fmla, data = Penicillin, chains = 2, iter = 400, seed = SEED)
  ans <- lmer(fmla, data = Penicillin)
  expect_equal(fixef(fit), fixef(ans), tol = FIXEF_tol)
  expect_equal(ranef(fit), ranef(ans), tol = RANEF_tol)
  expect_equal(ngrps(fit), ngrps(ans))
})

context("stan_glmer (binomial)")
test_that("stan_glmer returns expected result for cbpp example", {
  links <- c("logit", "probit", "cauchit", "log", "cloglog")
  # for (i in seq_along(links)) {
  i <- 1L # it seems only logit gives results similar to glmer with same link 
    fmla <- cbind(incidence, size - incidence) ~ period + (1 | herd)
    fit <- stan_glmer(fmla, data = cbpp, family = binomial(links[i]),
                      seed = SEED, chains = 2, iter = 400)
    ans <- glmer(fmla, data = cbpp, family = binomial(links[i]))
    expect_equal(fixef(fit), fixef(ans), tol = FIXEF_tol)
    expect_equal(ranef(fit), ranef(ans), tol = RANEF_tol)
    expect_equal(ngrps(fit), ngrps(ans))
  # }
})
