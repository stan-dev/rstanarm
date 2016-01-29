# loads packages, creates ROOT, SEED, and DATA_ENV
demo("SETUP", package = "rstanarm", verbose = FALSE, echo = FALSE, ask = FALSE)

source(paste0(ROOT, "ARM/Ch.14/election88.data.R"), local = DATA_ENV, verbose = FALSE) 
election88 <- with(DATA_ENV, data.frame(y, black, v.prev.full = v_prev_full,
                                        region.full = region_full, age.edu = age_edu,
                                        age, edu, female, state))

t_prior <- student_t(df = 7)
fmla1 <- y ~ black + female + (1 | state)
M1 <- stan_glmer(fmla1, data = election88, family = binomial(link="logit"),
                 prior = t_prior, prior_intercept = t_prior, 
                 seed = SEED, iter = 250, refresh = 125) # this model is a bit slow to run
print(M1, digits = 2) # can also do fixef(M1), ranef(M1), VarCorr(M1), etc. 

fmla2 <- y ~ black + female + black:female + v.prev.full + 
  (1 | age) + (1 | edu) + (1 | age.edu) + (1 | state) + (1 | region.full)
M2 <- update(M1, formula = fmla2)
print(M2)

ANSWER <- tolower(readline("Do you want to remove the objects this demo created? (y/n) "))
if (ANSWER != "n") {
  rm(election88, t_prior, fmla1, fmla2, ANSWER)
  # removes stanreg and loo objects, plus what was created by STARTUP
  demo("CLEANUP", package = "rstanarm", verbose = FALSE, echo = FALSE, ask = FALSE)
}
