stan_betareg.fit <- function (x, y, z = NULL, weights = NULL, offset = NULL,
                              link = "logit", link.phi = "log", ...,
                              prior = normal(), prior_intercept = normal(),
                              prior_ops = prior_options(), prior_PD = FALSE, 
                              algorithm = c("sampling", "optimizing", "meanfield", "fullrank"),
                              adapt_delta = NULL, QR = FALSE, sparse = FALSE) {
  
  # lots of tedious but simple stuff including standata which is a big list to pass to data {}
  # process the prior information like stan_glm.fit() does
  
  algorithm <- match.arg(algorithm)
  
  # no family argument
  
  supported_links <- c("logit", "probit","cloglog", "cauchit", "log", "loglog")
  
  link <- which(supported_links == link)
  if (!length(link)) 
    stop("'link' must be one of ", paste(supported_links, collapse = ", "))
  
  # useless assignments to pass R CMD check (need to include?)
  
  # prior distributions (handle_glm_prior() from data_block.R)
  prior_stuff <- handle_glm_prior(prior, nvars, family$link, default_scale = 2.5)
  for (i in names(prior_stuff)) # prior_{dist, mean, scale, df}
    assign(i, prior_stuff[[i]])
  prior_intercept_stuff <- handle_glm_prior(prior_intercept, nvars = 1, default_scale = 10,
                                            family$link, ok_dists = 
                                              nlist("normal", student_t = "t", "cauchy"))
  names(prior_intercept_stuff) <- paste0(names(prior_intercept_stuff), "_for_intercept")
  for (i in names(prior_intercept_stuff)) # prior_{dist, mean, scale, df}_for_intercept
    assign(i, prior_intercept_stuff[[i]])
  
  # create entries in the data block of the .stan file
  standata <- nlist(
    N = nrow(xtemp), K = ncol(xtemp), xbar = as.array(xbar), dense_X = !sparse,
    link, has_weights = length(weights) > 0, has_offset = length(offset) > 0,
    prior_dist, prior_mean, prior_scale, prior_df,
    family = 4L)
  
  # call stan() to draw from posterior distribution
  standata$prior_scale_for_dispersion <- 
    prior_scale_for_dispersion %ORifINF% 0
  
  stanfit <- stanmodels$continuous
  stanfit <- sampling(stanfit, data = standata, ...)
  return(stanfit)
}
