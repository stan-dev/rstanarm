# loads packages, creates ROOT, SEED, and DATA_ENV
demo("SETUP", package = "rstanarm", verbose = FALSE, echo = FALSE, ask = FALSE)

source(paste0(ROOT, "ARM/Ch.5/nes1992_vote.data.R"), local = DATA_ENV, verbose = FALSE)
nes1992 <- with(DATA_ENV, data.frame(vote, income))

invlogit <- plogis

# We'll use a Student t distribution with 7 degrees of freedom as our default 
# weakly informative prior for logistic regression coefficients. This prior 
# reflects the belief that the coefficients  are probably close to zero, are as
# likely to be positive as they are to be negative, but do have a small chance
# of being quite far from zero. Using a normal distribution instead of the t 
# distribution would be a more informative prior, as the tails of the normal
# distribution as less heavy. The t is therefore a bit more robust. 
t_prior <- student_t(df = 7, location = 0, scale = 2.5)

# Logistic regression with one predictor
vote_fit <- stan_glm(vote ~ income, data = nes1992, family=binomial(link="logit"),
                  prior = t_prior, prior_intercept = t_prior, 
                  seed = SEED, refresh = REFRESH)
print(vote_fit, digits = 2)
b <- coef(vote_fit)
plot(vote_fit, "hist", pars = names(b))

# Probability of Bush vote at various values of income
pr_bush <- function(x, coefs) invlogit(coefs[[1]] + coefs[[2]] * x)
income_vals <- with(nes1992, c(min(income), median(income), max(income)))
pr_bush(income_vals, b)
#  How the probability differs with a unit difference in x near the central value
pr_bush(3, b) - pr_bush(2, b)


# Wells in Bangladesh
source(paste0(ROOT, "ARM/Ch.5/wells.data.R"), local = DATA_ENV, verbose = FALSE)
wells <- with(DATA_ENV, data.frame(switch = switched, dist100 = dist/100, arsenic))

# Only use distance (in 100m) as predictor
post1 <- stan_glm(switch ~ dist100, data = wells, family = "binomial", 
                  prior = t_prior, prior_intercept = t_prior, 
                  seed = SEED, refresh = REFRESH)

# Add arsenic as predictor
post2 <- update(post1, formula = switch ~ dist100 + arsenic)

# Add interaction of dist100 and arsenic
post3 <- update(post2, formula = .~. + dist100:arsenic)
plot(post3, "areas", prob = 0.9, prob_outer = 1)

# Compare them with loo
loo1 <- loo(post1)
loo2 <- loo(post2)
loo3 <- loo(post3)
compare(loo1, loo2, loo3) # loo1 is dominated


# Graphing the fitted models
op <- par('mfrow')
par(mfrow = c(1,2))

jitter.binary <- function(a, jitt=.05){
  ifelse(a==0, runif(length(a), 0, jitt), runif (length(a), 1-jitt, 1))
}
b2 <- coef(post2)
b3 <- coef(post3)

# As function of dist100 
with(wells, plot(dist100, jitter.binary(switch), 
                 xlim=c(0, max(dist100)), ylab = "Prob"))
# Model with two predictors in red
curve(invlogit(cbind(1, x, .5) %*% b2), add = TRUE, col = "red", lty = 2) 
curve(invlogit(cbind(1, x, 1) %*% b2), add = TRUE, col = "red", lty = 2)
# Model with interaction in blue
curve(invlogit(cbind(1, x, .5, .5 * x) %*% b3), add = TRUE, col = "blue") 
curve(invlogit(cbind(1, x, 1, 1 * x) %*% b3), add = TRUE, col = "blue")

# As function of arsenic 
with(wells, plot(arsenic, jitter.binary(switch), 
                 xlim=c(0, max(arsenic)), ylab = "Prob"))
curve(invlogit(cbind (1, 0, x) %*% b2), add = TRUE, col = "red", lty = 2) 
curve(invlogit(cbind (1,.5, x) %*% b2), add = TRUE, col = "red", lty = 2)
curve(invlogit(cbind(1, 0, x, 0 * x) %*% b3), add = TRUE, col = "blue") 
curve(invlogit(cbind(1, .5, x, .5 * x) %*% b3), add = TRUE, col = "blue")
par(mfrow = op)



ANSWER <- tolower(readline("Do you want to remove the objects this demo created? (y/n) "))
if (ANSWER != "n") {
  rm(nes1992, invlogit, t_prior, b, pr_bush, income_vals, 
     wells, jitter.binary, b2, b3, op, ANSWER)
  # removes stanreg and loo objects, plus what was created by STARTUP
  demo("CLEANUP", package = "rstanarm", verbose = FALSE, echo = FALSE, ask = FALSE)
}
