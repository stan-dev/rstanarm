#' Draw from posterior predictive distribution
#' 
#' @export
#' @param object a model fit with \code{\link{stan_lm}} or 
#'   \code{\link{stan_glm}}.
#' @param newdata optionally, a data frame in which to look for variables with 
#'   which to predict. If omitted, the model matrix is used.
#' @param draws the number of draws to return. The default and maximum number of
#'   draws is the size of the posterior sample.
#' @param fun optional function to apply to the results. See examples below. 
#' 
#' @return a matrix of draws from the posterior predictive distribution.
#' 
#' @examples 
#' fit <- stan_glm(mpg ~ wt, data = mtcars)
#' ppd <- posterior_predict(fit) 
#' mean_wt <- mean(mtcars$wt)
#' hist(posterior_predict(fit, newdata = data.frame(wt = mean_wt)))
#'  
#' fit <- stan_glm(I(log(mpg)) ~ wt, data = mtcars)
#' ppd <- posterior_predict(fit, fun = exp)
#' 
posterior_predict <- function(object, newdata = NULL, draws = NULL, fun) {
  if (inherits(object, "stanreg-mle")) # replace with whatever name we end up using
    stop("posterior_predict only available for MCMC")
  family <- object$family
  famname <- family$family
  ppfun <- paste0(".pp_", famname)
  dat <- .pp_data(object, newdata)
  stanmat <- as.matrix(object$stanfit)
  S <- nrow(stanmat)
  if (is.null(draws)) 
    draws <- S
  else {
    if (draws > S)
      stop(paste("draws =", draws, "but only", S, 
                 "posterior draws found."), call. = FALSE)
  } 
  beta <- stanmat[, 1:ncol(dat$x)]
  eta <- linear_predictor(beta, dat$x, dat$offset)
  if (draws < S)
    eta <- eta[sample(S, draws), , drop = FALSE]
  ppargs <- list(mu = family$linkinv(eta))
  if (famname == "gaussian")
    ppargs$sigma <- stanmat[, "sigma"]
  if (famname == "binomial") {
    y <- if (!is.null(object$y)) 
      object$y else model.response(model.frame(object))
    ppargs$trials <- if (NCOL(y) == 2L) rowSums(y) else rep(1, NROW(y))
  }
  ytilde <- do.call(ppfun, ppargs)
  if (missing(fun)) ytilde
  else do.call(fun, list(ytilde))
}

.pp_gaussian <- function(mu, sigma) {
  t(sapply(1:nrow(mu), function(s) {
    rnorm(ncol(mu), mu[s,], sigma[s])
  }))
}
.pp_poisson <- function(mu) {
  t(sapply(1:nrow(mu), function(s) {
    rpois(ncol(mu), mu[s,])
  }))
}
.pp_binomial <- function(mu, trials) {
  t(sapply(1:nrow(mu), function(s) {
    rbinom(ncol(mu), size = trials, prob = mu[s,])
  }))
}
