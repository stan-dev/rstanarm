# Pointwise log-likelihood
# 
# Internal functions to compute the S by N pointwise log-likelihood matrix,
# where S is the size of the posterior sample (the number of simulations) and N
# is the number of data points.
# 
# @param llargs list with components family, x, y, weights, offset, beta, and
# sigma (note: sigma should be NULL if not gaussian model)
#
# @return a matrix.
#
pw_log_lik <- function(llargs) {
  fam <- llargs$family
  famname <- fam$family
  llfun <- paste0(".ll_", famname)
  eta <- llargs$x %*% t(llargs$beta)
  if (any(llargs$offset != 0))
    eta <- sweep(eta, MARGIN = 1L, llargs$offset, `+`)
  theta <- fam$linkinv(eta)
  args <- list(y = llargs$y, theta = theta)
  if (!is.null(llargs$sigma))
    args$sigma <- llargs$sigma
  ll <- do.call(llfun, args)
  if (all(llargs$weights == 1)) ll 
  else sweep(ll, MARGIN = 2L, llargs$weights,`*`)
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