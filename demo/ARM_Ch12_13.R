# loads packages, creates ROOT, SEED, and DATA_ENV
demo("SETUP", package = "rstanarm", verbose = FALSE, echo = FALSE, ask = FALSE)

### Radon data
source(paste0(ROOT, "ARM/Ch.12/radon.data.R"), local = DATA_ENV, verbose = FALSE) 
radon <- with(DATA_ENV, data.frame(y, x, u, radon, county))

# complete pooling
(pool <- stan_glm(y ~ x, data = radon, seed = SEED, refresh = REFRESH))
# no pooling 
(no_pool <- update(pool, formula = y ~ x + factor(county) - 1))

# varying intercept with no predictors
M0 <- stan_lmer(y ~ 1 + (1 | county), data = radon, 
                seed = SEED, refresh = REFRESH)
# varying intercept with individual-level predictor
M1 <- update(M0, formula = y ~ x + (1 | county))
# include group-level predictor
M2 <- update(M0, formula = y ~ x + u + (1 | county))
# varying intercepts and slopes
M3 <- update(M0, formula = y ~ x + (1 + x | county))
# varying intercepts and slopes with group-level predictor
M4 <- update(M0, formula = y ~ x + u + x:u + (1 + x | county))

# Can use VarCorr, coef, fixef, ranef just like after using lmer, e.g.
VarCorr(M2)
coef(M2)
fixef(M2)
ranef(M2)


### Pilots data
source(paste0(ROOT, "ARM/Ch.13/pilots.data.R"), local = DATA_ENV, verbose = FALSE) 
pilots <- with(DATA_ENV, data.frame(y, scenario_id, group_id))
M5 <- stan_lmer(y ~ 1 + (1 | group_id) + (1 | scenario_id), data = pilots, 
                seed = SEED, refresh = REFRESH)
VarCorr(M5)

### Earnings data
# regressions of earnings on ethnicity categories, age categories, and height
source(paste0(ROOT, "ARM/Ch.13/earnings.data.R"), local = DATA_ENV, verbose = FALSE)
earnings <- with(DATA_ENV, data.frame(earn = earn / 1e4, 
                                      height = scale(height), 
                                      eth, age))

f1 <- log(earn) ~ 1 + (1 + height | eth)
f2 <- log(earn) ~ 1 + (1 + height | eth) + (1 + height | age) + (1 + height | eth:age)
(fit1 <- stan_lmer(f1, data = earnings, seed = SEED, refresh = REFRESH))
(fit2 <- update(fit1, formula = f2))

ANSWER <- tolower(readline("Do you want to remove the objects this demo created? (y/n) "))
if (ANSWER != "n") {
  rm(radon, pilots, earnings, f1, f2, ANSWER)
  # removes stanreg and loo objects, plus what was created by STARTUP
  demo("CLEANUP", package = "rstanarm", verbose = FALSE, echo = FALSE, ask = FALSE)
}
