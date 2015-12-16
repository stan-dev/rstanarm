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
  
  opt <- object$algorithm == "optimizing" # used optimization
  mer <- !is.null(object$glmod) # used stan_(g)lmer
  
  # rstan::summary
  levs <- c(0.5, 0.8, 0.95, 0.99)
  qq <- (1 - levs) / 2
  probs <- sort(c(0.5, qq, 1 - qq))
  if (opt) {
    stanmat <- stanfit$theta_tilde
    stan_summary <- cbind(Median = apply(stanmat, 2, median), 
                          MAD_SD = apply(stanmat, 2, mad),
                          t(apply(stanmat, 2, quantile, 
                                  probs = c(0.025, .975))))
    covmat <- cov(stanmat)
    coefs <- apply(stanmat[,colnames(x),drop=FALSE], 2, median)
    ses <- apply(stanmat[,colnames(x),drop=FALSE], 2, mad)
  }
  else {
    stan_summary <- rstan::summary(stanfit, probs = probs, digits = 10)$summary
    coefs <- stan_summary[1:nvars, .select_median(object$algorithm)]
    if (length(coefs) == 1L) # ensures that if only a single coef it still gets a name
      names(coefs) <- rownames(stan_summary)[1L]
    if (any(stan_summary[,"Rhat"] > 1.1, na.rm = TRUE)) 
      warning("Markov chains did not converge! Do not analyze results!", 
              call. = FALSE, noBreaks. = TRUE)
  }    
  eta <- linear_predictor(coefs, x, offset)
  mu <- family$linkinv(eta)
  
  # residuals (of type 'response', unlike glm which does type 'deviance' by
  # default)
  residuals <- if (NCOL(y) == 2L)
    y[, 1] / rowSums(y) - mu else y - mu
  df.residual <- nobs - sum(weights == 0) - rank
  
  if (!opt) {
    stanmat <- as.matrix(stanfit)[, 1:nvars, drop = FALSE]
    covmat <- cov(stanmat)
    rownames(covmat) <- colnames(covmat) <- rownames(stan_summary)[1:nvars]
    ses <- apply(stanmat, 2, mad)
  }
  
  names(eta) <- names(mu) <- names(residuals) <- ynames
  offset <- if (any(offset != 0)) offset else NULL
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
    offset, 
    weights, 
    prior.weights = weights, 
    contrasts = object$contrasts, 
    na.action = object$na.action,
    call = object$call, 
    formula = object$formula, 
    terms = object$terms,
    prior.info = object$prior.info,
    algorithm = object$algorithm,
    stan_summary,  
    stanfit = if (opt) stanfit$stanfit else stanfit
  )
  if (opt) out$asymptotic_sampling_dist <- stanmat
  if (mer) out$glmod <- object$glmod
  structure(out, class = c("stanreg", "glm", "lm"))
}
