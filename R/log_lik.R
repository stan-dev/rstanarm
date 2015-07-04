#' Pointwise log-likelihood
#' 
#' Compute the \eqn{S} by \eqn{N} pointwise log-likelihood matrix, where \eqn{S}
#' is the size of the posterior sample (the number of simulations) and \eqn{N}
#' is the number of data points.
#' 
#' @export 
#' @param object a model fit with \code{\link{stan_lm}} or
#'   \code{\link{stan_glm}}.
#' @return a matrix.
#' 
log_lik <- function(object) {
  fam <- object$family
  famname <- fam$family
  llfun <- paste0(".ll_", famname)
  y <- model.response(object$model)
  X <- model.matrix(object)
  beta <- as.matrix(rstan::extract(object$stanfit, pars = "beta")$beta)
  theta <- fam$linkinv(X %*% t(beta))
  llargs <- nlist(y, theta)
  if (famname == "gaussian") {
    llargs$sigma <- rstan::extract(object$stanfit, pars = "sigma")$sigma  
  }
  ll <- do.call(llfun, llargs)
  if (all(object$weights == 1)) 
    ll 
  else 
    sweep(ll, MARGIN = 2, object$weights,`*`)
}

.ll_gaussian <- function(y, theta, sigma) {
  t(sapply(1:ncol(theta), function(s) {
    dnorm(y, mean = theta[,s], sd = sigma[s], log = TRUE)
  }))
}
.ll_poisson <- function(y, theta) {
  t(sapply(1:ncol(theta), function(s) {
    dpois(y, lambda = theta[,s], log = TRUE)
  }))
}
.ll_binomial <- function(y, theta) {
  if (NCOL(y) == 2) {
    trials <- y[,1] + y[,2]
    y <- y[,1]
  } else {
    trials <- 1
    if (is.factor(y)) 
      y <- y != levels(y)[1L]
    if (!all(y %in% c(0L, 1L)))
      stop("Not yet supported")
  }
  t(sapply(1:ncol(theta), function(s) {
    dbinom(y, size = trials, prob = theta[,s], log = TRUE)
  }))
}