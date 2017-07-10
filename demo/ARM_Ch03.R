# loads packages, creates ROOT, SEED, and DATA_ENV
demo("SETUP", package = "rstanarm", verbose = FALSE, echo = FALSE, ask = FALSE)

source(paste0(ROOT, "ARM/Ch.3/kidiq.data.R"), local = DATA_ENV, verbose = FALSE)
dat <- with(DATA_ENV, data.frame(kid_score, mom_hs, mom_iq))

# Estimate four contending models
post1 <- stan_glm(kid_score ~ mom_hs, data = dat, 
                  family = gaussian(link = "identity"), 
                  seed = SEED, refresh = REFRESH)
post2 <- update(post1, formula = kid_score ~ mom_iq)
post3 <- stan_lm(kid_score ~ mom_hs + mom_iq, data = dat,
                 prior = R2(location = 0.25, what = "mean"), 
                 seed = SEED, refresh = REFRESH)
post4 <- update(post3, formula = kid_score ~ mom_hs * mom_iq,
                 prior = R2(location = 0.30, what = "mean"))

# Compare them with loo
loo1 <- loo(post1)
loo2 <- loo(post2)
loo3 <- loo(post3)
loo4 <- loo(post4)
par(mfrow = c(2,2), mar = c(4,4,2,1) + .1)
plot(loo1, label_points = TRUE); title(main = "Model 1")
plot(loo2, label_points = TRUE); title(main = "Model 2")
plot(loo3, label_points = TRUE); title(main = "Model 3")
plot(loo4, label_points = TRUE); title(main = "Model 4")
compare(loo1, loo2, loo3, loo4) # fourth model dominates

# Generate predictions
IQ_SEQ <- seq(from = 75, to = 135, by = 5)
y_nohs <- posterior_predict(post4, newdata = data.frame(mom_hs = 0, mom_iq = IQ_SEQ))
y_hs <- posterior_predict(post4, newdata = data.frame(mom_hs = 1, mom_iq = IQ_SEQ))
par(mfrow = c(1:2), mar = c(5,4,2,1))
boxplot(y_hs, axes = FALSE, outline = FALSE, ylim = c(30,160),
        xlab = "Mom IQ", ylab = "Predicted Kid IQ", main = "Mom HS")
axis(1, at = 1:ncol(y_hs), labels = IQ_SEQ, las = 3)
axis(2, las = 1)
boxplot(y_nohs, outline = FALSE, col = "red", axes = FALSE, ylim = c(30,160),
        xlab = "Mom IQ", ylab = "Predicted Kid IQ", main = "Mom No HS")
axis(1, at = 1:ncol(y_hs), labels = IQ_SEQ, las = 3)
axis(2, las = 1)

# External Validation
source(paste0(ROOT, "ARM/Ch.3/kids_before1987.data.R"), 
       local = DATA_ENV, verbose = FALSE)
source(paste0(ROOT, "ARM/Ch.3/kids_after1987.data.R"), 
       local = DATA_ENV, verbose = FALSE)

fit_data <- with(DATA_ENV, data.frame(ppvt, hs, afqt))
pred_data <- with(DATA_ENV, data.frame(ppvt_ev, hs_ev, afqt_ev))

post5 <- stan_lm(ppvt ~ hs + afqt, data = fit_data,
                 prior = R2(location = 0.25, what = "mean"), 
                 seed = SEED, refresh = REFRESH)
y_ev <- posterior_predict(
  post5, 
  newdata = with(pred_data, data.frame(hs = hs_ev, afqt = afqt_ev))
)
par(mfrow = c(1,1))
hist(-sweep(y_ev, 2, STATS = pred_data$ppvt_ev, FUN = "-"), prob = TRUE,
     xlab = "Predictive Errors in ppvt", main = "", las = 2)


ANSWER <- tolower(readline("Do you want to remove the objects this demo created? (y/n) "))
if (ANSWER != "n") {
  rm(IQ_SEQ, y_nohs, y_hs, y_ev, ANSWER)
  # removes stanreg and loo objects, plus what was created by STARTUP
  demo("CLEANUP", package = "rstanarm", verbose = FALSE, echo = FALSE, ask = FALSE)
}
