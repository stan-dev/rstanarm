library(rstanarm)
library(lme4)

set.seed(123)

context("methods for stan_lmer models")

fit_lmer1 <- lmer(diameter ~ (1|plate) + (1|sample), data = Penicillin)
fit_stan1 <- stan_lmer(diameter ~ (1|plate) + (1|sample), data = Penicillin, 
                      iter = 100, chains = 1)

fit_lmer2 <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy)
fit_stan2 <- stan_lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy, 
                       iter = 100, chains = 1)

check_all_names <- function(object) {
  nms <- names(object)
  att_nms <- names(attributes(object))
  att_nms2 <- lapply(object, function(x) names(attributes(x)))
  c(nms, att_nms, att_nms2)
}

test_that("ngrps is right", {
  expect_equal(ngrps(fit_lmer1), ngrps(fit_stan1))
  expect_equal(ngrps(fit_lmer2), ngrps(fit_stan2))
})
test_that("VarCorr returns with correct structure", {
  vc_lmer1 <- VarCorr(fit_lmer1); vc_stan1 <- VarCorr(fit_stan1)
  vc_lmer2 <- VarCorr(fit_lmer2); vc_stan2 <- VarCorr(fit_stan2)
  expect_is(vc_stan1, class(vc_lmer1))
  expect_is(vc_stan2, class(vc_lmer2))
  expect_identical(check_att_names(vc_stan1), check_att_names(vc_lmer1))
  expect_identical(check_att_names(vc_stan2), check_att_names(vc_lmer2))
})
test_that("ranef returns with correct structure", {
  re_stan1 <- ranef(fit_stan1); re_lmer1 <- ranef(fit_lmer1)
  re_stan2 <- ranef(fit_stan1); re_lmer2 <- ranef(fit_lmer1)
  expect_is(re_stan1, class(re_lmer1))
  expect_is(re_stan2, class(re_lmer2))
  expect_identical(check_att_names(re_stan1), check_att_names(re_lmer1))
  expect_identical(check_att_names(re_stan2), check_att_names(re_lmer2))
})
test_that("fixef returns the right coefs", {
  expect_identical(names(fixef(fit_stan1)), names(fixef(fit_lmer1)))
  expect_identical(names(fixef(fit_stan2)), names(fixef(fit_lmer2)))
})

