# Pointwise log-likelihood
# 
# Internal functions to compute the S by N pointwise log-likelihood matrix,
# where S is the size of the posterior sample (the number of simulations) and N
# is the number of data points.
# 
# @param family,x,y,weights,offset,beta,sigma
#
# @return a matrix.
#
pw_log_lik <- function(family, x, y, weights, beta,
                       # remaining arguments might be applicable
                       sigma = NULL, zeta = NULL, offset = NULL) {
  eta <- linear_predictor(beta, x, offset)
  f <- family
  if (is(f, "family")) {
    llfun <- paste0(".ll_", f$family)
    linkinv <- f$linkinv
    mu <- linkinv(eta)
    args <- nlist(y, mu)
  }
  else if (is.character(f)) {
    llfun <- ".ll_polr"
    args <- nlist(y, eta, f)
  }
  else stop("'family' must be a family or a character string")
  if (!is.null(sigma)) args$sigma <- sigma
  if (!is.null(zeta))  args$zeta  <- zeta
  ll <- do.call(llfun, args)
  if (all(weights == 1)) ll 
  else sweep(ll, MARGIN = 2L, weights, `*`)
}

# pw_log_lik calls one of the functions below depending on the model
.ll_gaussian <- function(y, mu, sigma) {
  t(sapply(1:nrow(mu), function(s) {
    dnorm(y, mean = mu[s,], sd = sigma[s], log = TRUE)
  }))
}
.ll_poisson <- function(y, mu) {
  t(sapply(1:nrow(mu), function(s) {
    dpois(y, lambda = mu[s,], log = TRUE)
  }))
}
.ll_nb <- function(y, rho, theta) {
  stop(".ll_nb not working")
  t(sapply(1:nrow(rho), function(s) {
    alpha <- rho[s,]
    beta <- theta[s]
    lgamma(y + alpha - 1) - lgamma(alpha - 1) - lgamma(y) +
      alpha * log(beta / (beta + 1)) - y * log(beta + 1)
  }))
}
.ll_binomial <- function(y, mu) {
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
  t(sapply(1:nrow(mu), function(s) {
    dbinom(y, size = trials, prob = mu[s,], log = TRUE)
  }))
}
.ll_polr <- function(y, eta, f, zeta) {
  if (f == "logistic")    linkinv <- make.link("logit")$linkinv
  else if (f == "loglog") linkinv <- pgumbel
  else                    linkinv <- make.link(f)$linkinv
  y <- as.integer(y)
  N <- NROW(y)
  J <- max(y)
  ll <- matrix(NA_real_, nrow = nrow(zeta), ncol = N)
  for (i in 1:N) {
    y_i <- y[i]
    if      (y_i == 1) ll[,i] <- log(linkinv(zeta[,1] - eta[,i]))
    else if (y_i == J) ll[,i] <- log1p(-linkinv(zeta[,J-1] - eta[,i]))
    else ll[,i] <- log(linkinv(zeta[,y_i] - eta[,i]) - 
                       linkinv(zeta[,y_i - 1L] - eta[,i]))
  }
  return(ll)
}
