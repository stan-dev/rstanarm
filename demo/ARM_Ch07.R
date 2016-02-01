# loads packages, creates ROOT, SEED, and DATA_ENV
demo("SETUP", package = "rstanarm", verbose = FALSE, echo = FALSE, ask = FALSE)

source(paste0(ROOT, "ARM/Ch.7/congress.data.R"), local = DATA_ENV, verbose = FALSE)

# The stuff in sections 7.0 -- 7.2 is not very relevant 

post1 <- stan_lm(vote_88 ~ vote_86 + incumbency_88, data = DATA_ENV, 
                 prior = R2(0.9, what = "mean"), 
                 seed = SEED, refresh = REFRESH)
post1 # badly underfitting
y_tilde <- posterior_predict(post1) # incumbency_90 is not available
summary(rowSums(y_tilde > 0.5))

source(paste0(ROOT, "ARM/Ch.6/wells.data.R"), local = DATA_ENV, verbose = FALSE)
post2 <- stan_glm(switch ~ I(dist / 100), data = DATA_ENV, family = "binomial", 
                  seed = SEED, refresh = REFRESH)
prop.table(table(c(ppd)))

# the compound model is not good because it assumes the two errors are 
# independent. rstanarm will eventually support Heckman models, which
# would be a better choice here.

ANSWER <- tolower(readline("Do you want to remove the objects this demo created? (y/n) "))
if (ANSWER != "n") {
  rm(y_tilde, ppd, ANSWER)
  # removes stanreg and loo objects, plus what was created by STARTUP
  demo("CLEANUP", package = "rstanarm", verbose = FALSE, echo = FALSE, ask = FALSE)
}
