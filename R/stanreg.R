# Create a stanreg object
#
# @param object A list provided by one of the \code{stan_*} modeling functions.
# @return A stanreg object.
#
stanreg <- function(object) {
  stanfit <- object$stanfit
  family <- object$family
  y <- object$y
  x <- object$x
  nvars <- ncol(x)
  nobs <- NROW(y)
  ynames <- if (is.matrix(y)) rownames(y) else names(y)
  rank <- qr(x, tol = .Machine$double.eps, LAPACK = TRUE)$rank 
  
  opt <- object$algorithm == "optimizing"
  mer <- !is.null(object$glmod) # used stan_(g)lmer
  
  levs <- c(0.5, 0.8, 0.95, 0.99)
  qq <- (1 - levs) / 2
  probs <- sort(c(0.5, qq, 1 - qq))
  if (opt) {
    stanmat <- stanfit$theta_tilde
    stan_summary <- cbind(Median = apply(stanmat, 2L, median), 
                          MAD_SD = apply(stanmat, 2L, mad),
                          t(apply(stanmat, 2L, quantile, 
                                  probs = c(0.025, .975))))
    covmat <- cov(stanmat)
    coefs <- apply(stanmat[, colnames(x), drop = FALSE], 2L, median)
    ses <- apply(stanmat[, colnames(x), drop = FALSE], 2L, mad)
  } else {
    stan_summary <- rstan::summary(stanfit, probs = probs, digits = 10)$summary
    coefs <- stan_summary[1:nvars, .select_median(object$algorithm)]
    if (length(coefs) == 1L) # ensures that if only a single coef it still gets a name
      names(coefs) <- rownames(stan_summary)[1L]
    
    if (any(stan_summary[,"Rhat"] > 1.1, na.rm = TRUE)) 
      warning("Markov chains did not converge! Do not analyze results!", 
              call. = FALSE, noBreaks. = TRUE)
    
    stanmat <- as.matrix(stanfit)[, 1:nvars, drop = FALSE]
    covmat <- cov(stanmat)
    rownames(covmat) <- colnames(covmat) <- rownames(stan_summary)[1:nvars]
    ses <- apply(stanmat, 2L, mad)
  }
  
  # linear predictor, fitted values, and residuals (of type 'response', unlike
  # glm which does 'deviance' residuals by default)
  eta <- linear_predictor(coefs, x, object$offset)
  mu <- family$linkinv(eta)

  if (NCOL(y) == 2L) {
    residuals <- y[, 1L] / rowSums(y) - mu 
  } else {
    ytmp <- if (is.factor(y)) as.integer(y != levels(y)[1L]) else y
    residuals <- ytmp - mu
  }
  df.residual <- nobs - sum(object$weights == 0) - rank
  names(eta) <- names(mu) <- names(residuals) <- ynames
  
  out <- nlist(
    coefficients = coefs, 
    ses,
    fitted.values = mu, 
    linear.predictors = eta,
    residuals, 
    df.residual, 
    covmat,
    y, 
    x, 
    model = object$model, 
    data = object$data, 
    family, 
    rank,
    offset = if (any(object$offset != 0)) object$offset else NULL,
    weights = object$weights, 
    prior.weights = object$weights, 
    contrasts = object$contrasts, 
    na.action = object$na.action,
    call = object$call, 
    formula = object$formula, 
    terms = object$terms,
    xlevels = object$xlevels,
    prior.info = object$prior.info,
    algorithm = object$algorithm,
    stan_summary,  
    stanfit = if (opt) stanfit$stanfit else stanfit,
    asymptotic_sampling_dist = if (opt) stanmat else NULL,
    glmod = if (mer) object$glmod else NULL
  )
  structure(out, class = c("stanreg", "glm", "lm"))
}
