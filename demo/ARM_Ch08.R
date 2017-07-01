# loads packages, creates ROOT, SEED, and DATA_ENV
demo("SETUP", package = "rstanarm", verbose = FALSE, echo = FALSE, ask = FALSE)

source(paste0(ROOT, "ARM/Ch.8/lightspeed.data.R"), local = DATA_ENV, verbose = FALSE)
light_dat <- with(DATA_ENV, data.frame(y))
# The stuff in sections 8.0 -- 8.2 is not very relevant 

(post1 <- stan_glm(y ~ 1, data = light_dat, seed = SEED, refresh = REFRESH))
y_rep <- posterior_predict(post1)

pp_check(post1, plotfun = "stat", stat = "min") + 
  ggtitle("Minimum Predicted Measurement Error")

# make similar plot manually
hist(apply(y_rep, 1, min), prob = TRUE, main = "", las = 1,
     xlab = "Minimum Predicted Measurement Error", xlim = c(-45,20))
abline(v = min(DATA_ENV$y), col = "red")

# Compare observed y to several replicated y (y_rep) from posterior predictive
# distribution
ttl <- paste("Measurement Error for the Speed of Light", 
             "\nvs Predicted Measurement Error")
pp_check(post1, plotfun = "hist") + ggtitle(ttl)

# Make similar plot manually but combine all y_rep
op <- par('mfrow')
par(mfrow = 1:2, mar = c(5,4,1,1) + .1)
hist(light_dat$y, prob = TRUE, main = "", las = 1,
     xlab = "Measurement Error for the Speed of Light")
hist(y_rep, prob = TRUE, main = "", las = 1,
     xlab = "Predicted Measurement Error")
par(mfrow = op)

# Roaches example
data(roaches, package = "rstanarm")
post2 <- stan_glm(y ~ roach1 + treatment + senior, data = roaches, 
                  family = poisson(link = "log"), seed = SEED, refresh = REFRESH)
y_rep <- posterior_predict(post2)

# Compare observed proportion of zeros to predicted proportion of zeros
mean(y_rep == 0)
mean(roaches$y == 0)
summary(apply(y_rep == 0, 1, mean))
prop0 <- function(x) mean(x == 0)
pp_check(post2, plotfun = "stat", stat = "prop0") # model doesn't predict enough zeros

# Negative binomial model does a much better job handling the zeros
post3 <- update(post2, family = neg_binomial_2())
pp_check(post3, plotfun = "stat", stat = "prop0")


# rstanarm does not yet support time-series models

ANSWER <- tolower(readline("Do you want to remove the objects this demo created? (y/n) "))
if (ANSWER != "n") {
  rm(y_rep, prop0, ttl, op, ANSWER)
  # removes stanreg and loo objects, plus what was created by STARTUP
  demo("CLEANUP", package = "rstanarm", verbose = FALSE, echo = FALSE, ask = FALSE)
}
