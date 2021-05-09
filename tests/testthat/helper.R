run_example_model <- function() {
  o <- capture.output(suppressWarnings(
    fit <- stan_glmer(cbind(incidence, size - incidence) ~ size + period + (1|herd),
             data = lme4::cbpp, family = binomial, QR = TRUE,
             # this next line is only to keep the example small in size!
             chains = 2, cores = 1, seed = 12345, iter = 1000, refresh = 0)
  ))
  fit
}