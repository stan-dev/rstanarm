# loads packages, creates ROOT, SEED, and DATA_ENV
demo("SETUP", package = "rstanarm", verbose = FALSE, echo = FALSE, ask = FALSE)
# read data into DATA_ENV environment
source(paste0(ROOT, "ARM/Ch.9/electric_grade4.data.R"), local = DATA_ENV, 
       verbose = FALSE)

post1 <- stan_lm(post_test ~ treatment * pre_test, data = DATA_ENV, 
                 prior = R2(0.75), seed = SEED, refresh = REFRESH)
post1 # underfitting but ok because it is an experiment
plot(post1)

y_0 <- posterior_predict(post1, data.frame(treatment = 0, pre_test = DATA_ENV$pre_test))
y_1 <- posterior_predict(post1, data.frame(treatment = 1, pre_test = DATA_ENV$pre_test))
diff <- y_1 - y_0
mean(diff)
sd(diff) # much larger than in ARM
hist(diff, prob = TRUE, main = "", xlab = "Estimated Average Treatment Effect", las = 1)


stopifnot(require(gridExtra))
plots <- sapply(1:4, simplify = FALSE, FUN = function(k) {
  source(paste0(ROOT, "ARM/Ch.9/electric_grade", k, "_supp.data.R"), 
         local = DATA_ENV, verbose = FALSE)
  out <- plot(stan_lm(post_test ~ supp + pre_test, data = DATA_ENV, 
                    seed = SEED, refresh = REFRESH,
                    prior = R2(0.75, what = "mean")))
  out + ggtitle(paste("Grade =", k))
})
marrangeGrob(plots, nrow = 2, ncol = 2)

ANSWER <- tolower(readline("Do you want to remove the objects this demo created? (y/n) "))
if (ANSWER != "n") {
  rm(y_0, y_1, diff, plots, ANSWER)
  # removes stanreg and loo objects, plus what was created by STARTUP
  demo("CLEANUP", package = "rstanarm", verbose = FALSE, echo = FALSE, ask = FALSE)
}
