stanreg <- function(object) {
  stanfit <- object$stanfit
  weights <- object$weights
  offset <- object$offset
  family <- object$family
  y <- object$y
  x <- object$x
  nvars <- ncol(x)
  nobs <- NROW(y)
  ynames <- if (is.matrix(y)) rownames(y) else names(y)
  rank <- qr(x, tol = .Machine$double.eps, LAPACK = TRUE)$rank 
  
  # rstan::summary
  levs <- c(0.5, 0.8, 0.95, 0.99)
  qq <- (1 - levs) / 2
  probs <- sort(c(0.5, qq, 1 - qq))
  if (is.list(stanfit)) { # used optimization
    L <- t(chol(stanfit$cov.scaled))
    k <- nrow(L)
    unconstrained <- stanfit$par[1:k] + L %*% matrix(rnorm(4000 * k), k)
    stanmat <- t(apply(unconstrained, 2, FUN = function(u)
      unlist(constrain_pars(stanfit$stanfit, u))))
    stan_summary <- cbind(Estimate = stanfit$par, 
                          "Std. Error" = apply(stanmat, 2, sd),
                          t(apply(stanmat, 2, quantile, 
                                  probs = c(0.025, .975))))
    covmat <- cov(stanmat)
    coefs <- stanfit$par[grep("^gamma|^sigma|^overdispersion|^mean_PPD",
                              names(stanfit$par), invert = TRUE)]
    if ("(Intercept)" %in% names(coefs)) coefs <- c(tail(coefs, 1), head(coefs, -1))
  }
  else {
    stan_summary <- rstan::summary(stanfit, probs = probs, digits = 10)$summary
    coefs <- stan_summary[1:nvars, 1]
  }    

  eta <- linear_predictor(coefs, x, offset)
  mu <- family$linkinv(eta)
  
  # residuals (of type 'response', unlike glm which does type 'deviance' by
  # default)
  residuals <- if (NCOL(y) == 2L)
    y[, 1] / rowSums(y) - mu else y - mu
  df.residual <- nobs - sum(weights == 0) - rank
  
  if (!is.list(stanfit)) {
    stanmat <- as.matrix(stanfit)
    covmat <- cov(stanmat[,1:nvars])
    rownames(covmat) <- colnames(covmat) <- rownames(stan_summary)[1:nvars]
  }
  
  # pointwise log-likelihood
  mark <- 1:nvars
  if (is.list(stanfit)) {
    mark <- c(grep("^alpha", colnames(stanmat)), 
              grep("^beta",  colnames(stanmat)))
  }
  llargs <- nlist(family, x, y, weights, offset, 
                  beta = stanmat[,mark])
  if (family$family == "gaussian") llargs$sigma <- stanmat[, "sigma"]
  if (pmatch("Negative Binomial", family$family, nomatch = 0L) == 1) {
    llargs$theta <- stanmat[,"overdispersion"]
    family$family <- "nb"
  }
  else log_lik <- do.call("pw_log_lik", llargs)
  
  names(eta) <- names(mu) <- names(residuals) <- ynames
  offset <- if (any(offset != 0)) offset else NULL
  out <- list(
    coefficients = coefs, fitted.values = mu, linear.predictors = eta,
    residuals = residuals, df.residual = df.residual, covmat = covmat,
    y = y, x = x, model = object$model, data = object$data, rank = rank,
    offset = offset, weights = weights, prior.weights = weights, 
    family = family, contrasts = object$contrasts, na.action = object$na.action,
    call = object$call, formula = object$formula, terms = object$terms,
    prior.info = object$prior.info, log_lik = log_lik,
    stan_summary = stan_summary,  
    stanfit = if (is.list(stanfit)) stanfit$stanfit else stanfit
  )
  class(out) <- c("stanreg", "glm", "lm")
  out
}
