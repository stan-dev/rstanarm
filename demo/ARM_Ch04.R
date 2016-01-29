# loads packages, creates ROOT, SEED, and DATA_ENV
demo("SETUP", package = "rstanarm", verbose = FALSE, echo = FALSE, ask = FALSE)

source(paste0(ROOT, "ARM/Ch.4/earnings.data.R"), local = DATA_ENV, verbose = FALSE)

# The stuff in sections 4.0 -- 4.3 is not very relevant 
# Moreover, centering predictors is NOT recommended in the rstanarm package
# Just look at the posterior predictive distribution 
# over a range of values to interpret the effect of a predictor

# These two models are essentially equivalent in the likelihood
# But the "same" priors affect the posterior differently
post1 <- stan_glm(log(earn) ~ height, data = DATA_ENV, 
                  family = gaussian(link = "identity"), 
                  seed = SEED, refresh = REFRESH)
# post2 <- stan_glm(earn ~ height, data = DATA_ENV, 
#                   family = gaussian(link = "log"), 
#                   seed = SEED, refresh = REFRESH)
# and this does not even converge

# These models add terms to the right-hand side
post3 <- stan_lm(log(earn) ~ height + male, data = DATA_ENV,
                 prior = R2(location = 0.3, what = "mean"), 
                 seed = SEED, refresh = REFRESH)
post4 <- update(post3, formula = log(earn) ~ height * male)

# Compare them with loo
loo1 <- loo(post1)
# post2 is not comparable to the others anyway
loo3 <- loo(post3)
loo4 <- loo(post4)
compare(loo1, loo3, loo4) # loo1 is dominated

# Generate predictions to interpret
WOMEN_SEQ <- seq(from = 58, to = 75, by = 1)
MEN_SEQ <- seq(from = 60, to = 77, by = 1)
YLIM <- c(500, 100000)
y_women <- posterior_predict(post4, fun = exp,
                             newdata = data.frame(male = 0, height = WOMEN_SEQ))
y_men <- posterior_predict(post4, fun = exp,
                           newdata = data.frame(male = 1, height = MEN_SEQ))
par(mfrow = c(1:2), mar = c(5,4,2,1) + .1)
boxplot(y_women, axes = FALSE, outline = FALSE, log = "y", ylim = YLIM,
        xlab = "Height in Inches", ylab = "", main = "Predicted Earnings of Women")
axis(1, at = 1:ncol(y_women), labels = WOMEN_SEQ, las = 3)
axis(2, las = 1)
boxplot(y_men, outline = FALSE, col = "red", axes = FALSE, log = "y", ylim = YLIM,
        xlab = "Height in Inches", ylab = "", main = "Predicted Earnings of Men")
axis(1, at = 1:ncol(y_men), labels = MEN_SEQ, las = 3)

# Prediction of the weight of mesquite trees 
source(paste0(ROOT, "ARM/Ch.4/mesquite.data.R"), local = DATA_ENV, verbose = FALSE)
CONTINUE1 <- tolower(readline(
  paste("A heads up: the next part of the demo (Predicting weight of mesquite trees )",
        "prints many lines \nto the console as it runs many models and compares the results", 
        "Proceed? (y/n)")
))
if (CONTINUE1 != "n") {
  post5 <- stan_lm(weight ~ diam1 + diam2 + canopy_height + total_height +
                     density + group, data = DATA_ENV,
                   prior = R2(0.9), seed = SEED, refresh = REFRESH)
  post6 <- update(post5, formula = log(weight) ~ log(diam1) + log(diam2) + 
                    log(canopy_height) + log(total_height) + log(density) + group)
  post7 <- update(post5, formula = log(weight) ~ log(diam1 * diam2 * canopy_height), 
                  prior = R2(0.75, what = "mean"))
  post8 <- update(post5, formula = log(weight) ~ log(diam1 * diam2 * canopy_height) + 
                    log(diam1 * diam2) + group, prior = R2(0.8))
  post9 <- update(post5, formula = log(weight) ~ log(diam1 * diam2 * canopy_height) + 
                    log(diam1 * diam2) + log(diam1 / diam2) + group, 
                  prior = R2(0.85))
  
  # Compare them with loo
  compare(loo(post5), loo(post6), loo(post7), loo(post8), loo(post9))
}

# Predicting "continuous" party ID over time without multilevel stuff
CONTINUE2 <- tolower(readline(
  paste("A heads up: the next part of the demo (Predicting party ID over time)",
        "prints many lines \nto the console as it runs many models and compares the results", 
        "Proceed? (y/n)")
))
if (CONTINUE2 != "n") {
  YEARS <- as.character(seq(from = 1972, to = 2000, by = 4))
  round(digits = 2, x = sapply(YEARS, FUN = function(YEAR) {
    source(paste0(ROOT, "ARM/Ch.4/nes", YEAR, ".data.R"), local = DATA_ENV, verbose = FALSE)
    coef(stan_lm(partyid7 ~ real_ideo + I(race_adj == 1) + as.factor(age_discrete) + 
                   educ1 + gender + income, data = DATA_ENV, prior = R2(0.5),
                 seed = SEED, refresh = 0))
  }))
}

ANSWER <- tolower(readline("Do you want to remove the objects this demo created? (y/n) "))
if (ANSWER != "n") {
  if (CONTINUE2 != "n") rm(YEARS)
  rm(WOMEN_SEQ, MEN_SEQ, y_women, y_men, YLIM, ANSWER, CONTINUE1, CONTINUE2)
  # removes stanreg and loo objects, plus what was created by STARTUP
  demo("CLEANUP", package = "rstanarm", verbose = FALSE, echo = FALSE, ask = FALSE)
}
