stanreg <- function(object) {
  stanfit <- object$stanfit
  weights <- object$weights
  offset <- if (is.null(object$offset)) 0 else object$offset
  family <- object$family
  y <- object$y
  x <- object$x
  nvars <- ncol(x)
  nobs <- NROW(y)
  ynames <- if (is.matrix(y)) rownames(y) else names(y)
  rank <- qr(x, tol = .Machine$double.eps, LAPACK = TRUE)$rank 
  
  # rstan::summary
  levs <- c(0.5, 0.8, 0.95, 0.99)
  qq <- (1 - levs)/2
  probs <- sort(c(0.5, c(qq, 1 - qq)))
  stan_summary <- rstan::summary(stanfit, probs = probs, digits = 10)$summary
  
  # linear predictors and fitted values
  mu <- stan_summary[1:nvars, "mean"]
  eta <- if (NCOL(x) == 1L) x * mu else x %*% mu
  eta <- as.vector(eta) + offset
  mu <- family$linkinv(eta)
  
  # residuals (of type 'response', unlike glm which does type 'deviance' by
  # default)
  residuals <- if (NCOL(y) == 2L)
    y[, 1] / rowSums(y) - mu else y - mu
  df.residual <- nobs - sum(weights == 0) - rank
  
  # covariance matrix
  beta <- rstan::extract(stanfit, pars = "beta")$beta
  covmat <- cov(beta)
  
  # pointwise log-likelihood
  llargs <- nlist(family, x, y, weights, offset, beta, sigma = NULL)
  if (family$family == "gaussian") {
    llargs$sigma <- rstan::extract(stanfit, pars = "sigma")$sigma  
  } 
  log_lik <- do.call("pw_log_lik", llargs)
  
  names(eta) <- names(mu) <- names(residuals) <- ynames
  rownames(covmat) <- colnames(covmat) <- rownames(stan_summary)[1:nvars]
  
  out <- list(
    coefficients = stan_summary[1:nvars, "mean"],
    fitted.values = mu, linear.predictors = eta,
    residuals = residuals, df.residual = df.residual, covmat = covmat,
    y = y, x = x, model = object$model, data = object$data, rank = rank,
    offset = object$offset, weights = weights, prior.weights = weights, 
    family = family, contrasts = object$contrasts, na.action = object$na.action,
    call = object$call, formula = object$formula, terms = object$terms,
    prior.info = object$prior.info, log_lik = log_lik,
    stan_summary = stan_summary, stanfit = stanfit
  )
  class(out) <- c("stanreg", "glm", "lm")
  out
}
