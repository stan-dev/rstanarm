# loads packages, creates ROOT, SEED, and DATA_ENV
demo("SETUP", package = "rstanarm", verbose = FALSE, echo = FALSE, ask = FALSE)

source(paste0(ROOT, "ARM/Ch.8/lightspeed.data.R"), local = DATA_ENV, verbose = FALSE)

# The stuff in sections 8.0 -- 8.2 is not very relevant 

post1 <- stan_glm(y ~ 1, data = DATA_ENV, seed = SEED, family = gaussian())
post1
y_rep <- posterior_predict(post1)

hist(apply(y_rep, 1, min), prob = TRUE, main = "", las = 1,
     xlab = "Minimum Predicted Measurement Error", xlim = c(-45,20))
abline(v = min(DATA_ENV$y), col = "red")
par(mfrow = 1:2, mar = c(5,4,1,1) + .1)
hist(DATA_ENV$y, prob = TRUE, main = "", las = 1,
     xlab = "Measurement Error for the Speed of Light")
hist(y_rep, prob = TRUE, main = "", las = 1,
     xlab = "Predicted Measurement Error")

source(paste0(ROOT, "ARM/Ch.8/roaches.data.R"), local = DATA_ENV, verbose = FALSE)
post2 <- stan_glm(y ~ roach1 + treatment + senior, data = DATA_ENV, 
                  family = poisson(link = "log"), seed = SEED)
y_rep <- posterior_predict(post2)
mean(y_rep == 0)
mean(DATA_ENV$y == 0)
summary(apply(y_rep == 0, 1, mean))

# rstanarm does not currently support time-series models

ANSWER <- tolower(readline("Do you want to remove the objects this demo created? (y/n) "))
if (ANSWER != "n") {
  rm(y_rep, ANSWER)
  # removes stanreg and loo objects, plus what was created by STARTUP
  demo("CLEANUP", package = "rstanarm", verbose = FALSE, echo = FALSE, ask = FALSE)
}
