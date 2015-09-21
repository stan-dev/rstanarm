# loads packages, creates ROOT, SEED, and DATA_ENV
demo("SETUP", package = "rstanarm", verbose = FALSE, echo = FALSE, ask = FALSE)

source(paste0(ROOT, "ARM/Ch.14/election88.data.R"), local = DATA_ENV, verbose = FALSE) 
election88 <- with(DATA_ENV, data.frame(y, black, v.prev.full = v_prev_full,
                                        region.full = region_full, age.edu = age_edu,
                                        age, edu, female, state))

t_prior <- student_t(df = 7)
M1 <- stan_glmer(y ~ black + female + (1|state), data = election88, 
                 family=binomial(link="logit"),
                 prior = t_prior, prior_intercept = t_prior)
fixef(M1)
ranef(M1)
VarCorr(M1)

fmla <- y ~ black + female + black:female + v.prev.full + 
  (1 | age) + (1 | edu) + (1 | age.edu) + (1 | state) + (1 | region.full)
M2 <- stan_glmer(fmla, data = election88, family = binomial(link = "logit"),
                 prior = t_prior, prior_intercept = t_prior)


ANSWER <- tolower(readline("Do you want to remove the objects this demo created? (y/n) "))
if (ANSWER != "n") {
  rm(election88)
  # removes stanreg and loo objects, plus what was created by STARTUP
  demo("CLEANUP", package = "rstanarm", verbose = FALSE, echo = FALSE, ask = FALSE)
}
