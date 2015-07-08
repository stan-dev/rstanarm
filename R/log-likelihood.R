# Pointwise log-likelihood
# 
# Internal functions to compute the S by N pointwise log-likelihood matrix,
# where S is the size of the posterior sample (the number of simulations) and N
# is the number of data points.
# 
# @param family,x,y,weights,offset,beta,sigma (note: sigma should be NULL if not
#   gaussian model)
#
# @return a matrix.
#
pw_log_lik <- function(family, x, y, weights, offset, beta, sigma = NULL) {
  f <- family
  llfun <- paste0(".ll_", f$family)
  eta <- x %*% t(beta)
  if (any(offset != 0))
    eta <- sweep(eta, MARGIN = 1L, offset, `+`)
  theta <- f$linkinv(eta)
  args <- nlist(y, theta)
  if (!is.null(sigma))
    args$sigma <- sigma
  ll <- do.call(llfun, args)
  if (all(weights == 1)) ll 
  else sweep(ll, MARGIN = 2L, weights, `*`)
}

# pw_log_lik calls one of the functions below depending on the model
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
  if (NCOL(y) == 2L) {
    trials <- rowSums(y)
    y <- y[, 1L]
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