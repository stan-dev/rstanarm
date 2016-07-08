stan_betareg.fit <- function (x, y, z = NULL, weights = NULL, offset = NULL, link = "logit", 
          link.phi = "log", ...
          prior = normal(), prior_intercept = normal(),
          prior_ops = prior_options(), prior_PD = FALSE, 
          algorithm = c("sampling", "optimizing", 
                        "meanfield", "fullrank"),
          adapt_delta = NULL, QR = FALSE, sparse = FALSE) {
  
  # lots of tedious but simple stuff including standata which is a big list to pass to data {}
  # process the prior information like stan_glm.fit() does
  stanfit <- stanmodels$continuous
  stanfit <- sampling(stanfit, data = standata, ...)
  return(stanfit)
}
