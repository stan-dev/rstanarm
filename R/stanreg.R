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
  fit_summary <- rstan::summary(stanfit, probs = probs, digits = 5)$summary
  
  # linear predictors and fitted values
  mu <- fit_summary[1:nvars, "mean"]
  eta <- if (NCOL(x) == 1L) x * mu else x %*% mu
  eta <- as.vector(eta) + offset
  mu <- family$linkinv(eta)
  
  # residuals (default to response type for linear models and deviance residuals
  # otherwise. this mimics lm and glm behavior but maybe we don't want to do
  # that)
  if (family$family == "gaussian" && family$link == "identity") {
    residuals <- y - mu
    attr(residuals, "type") <- "response"
  } else {
    d.res <- sqrt(pmax((object$family$dev.resids)(y, mu, weights), 0))
    residuals <- ifelse(y > mu, d.res, -d.res)
    attr(residuals, "type") <- "deviance"
  }
  df.residual <- nobs - sum(weights == 0) - rank
  
  # covariance matrix
  beta <- rstan::extract(stanfit, pars = "beta")$beta
  covmat <- cov(beta)
  
  names(eta) <- names(mu) <- names(residuals) <- ynames
  rownames(covmat) <- colnames(covmat) <- rownames(fit_summary)[1:nvars]
  
  out <- list(
    coefficients = fit_summary[1:nvars, "mean"],
    fitted.values = mu, linear.predictors = eta,
    residuals = residuals, df.residual = df.residual, covmat = covmat,
    y = y, x = x, model = object$model, data = object$data,
    offset = object$offset, prior.weights = weights, rank = rank,
    family = family, contrasts = object$contrasts, na.action = object$na.action,
    call = object$call, formula = object$formula, terms = object$terms,
    prior.info = object$prior.info, stanfit = stanfit
  )
  class(out) <- c("stanreg", "glm", "lm")
}
