# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2015, 2016 Trustees of Columbia University
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

# Create a stanreg object
#
# @param object A list provided by one of the \code{stan_*} modeling functions.
# @return A stanreg object.
#
stanreg <- function(object) {
  opt <- object$algorithm == "optimizing"
  mer <- !is.null(object$glmod) # used stan_(g)lmer
  stanfit <- object$stanfit
  family <- object$family
  y <- object$y
  x <- object$x
  nvars <- ncol(x)
  nobs <- NROW(y)
  ynames <- if (is.matrix(y)) rownames(y) else names(y)
  if (opt) {
    stanmat <- stanfit$theta_tilde
    probs <- c(0.025, .975)
    stan_summary <- cbind(Median = apply(stanmat, 2L, median), 
                          MAD_SD = apply(stanmat, 2L, mad),
                          t(apply(stanmat, 2L, quantile, probs)))
    xnms <- colnames(x)
    covmat <- cov(stanmat)[xnms, xnms]
    coefs <- apply(stanmat[, xnms, drop = FALSE], 2L, median)
    ses <- apply(stanmat[, xnms, drop = FALSE], 2L, mad)
    rank <- qr(x, tol = .Machine$double.eps, LAPACK = TRUE)$rank
    df.residual <- nobs - sum(object$weights == 0) - rank
  } else {
    levs <- c(0.5, 0.8, 0.95, 0.99)
    qq <- (1 - levs) / 2
    probs <- sort(c(0.5, qq, 1 - qq))
    stan_summary <- rstan::summary(stanfit, probs = probs, digits = 10)$summary
    coefs <- stan_summary[1:nvars, select_median(object$algorithm)]
    if (length(coefs) == 1L) # ensures that if only a single coef it still gets a name
      names(coefs) <- rownames(stan_summary)[1L]
    
    stanmat <- as.matrix(stanfit)[, 1:nvars, drop = FALSE]
    covmat <- cov(stanmat)
    rownames(covmat) <- colnames(covmat) <- rownames(stan_summary)[1:nvars]
    ses <- apply(stanmat, 2L, mad)
    if (object$algorithm == "sampling") 
      check_rhats(stan_summary[, "Rhat"])
  }
  
  # linear predictor, fitted values
  eta <- linear_predictor(coefs, x, object$offset)
  mu <- family$linkinv(eta)

  if (NCOL(y) == 2L) {
    # residuals of type 'response', (glm which does 'deviance' residuals by default)
    residuals <- y[, 1L] / rowSums(y) - mu 
  } else {
    ytmp <- if (is.factor(y)) fac2bin(y) else y
    residuals <- ytmp - mu
  }
  names(eta) <- names(mu) <- names(residuals) <- ynames
  
  out <- nlist(
    coefficients = coefs, 
    ses,
    fitted.values = mu,
    linear.predictors = eta,
    residuals, 
    df.residual = if (opt) df.residual else NA_integer_, 
    covmat,
    y, 
    x, 
    model = object$model, 
    data = object$data, 
    family, 
    offset = if (any(object$offset != 0)) object$offset else NULL,
    weights = object$weights, 
    prior.weights = object$weights, 
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
  if (opt) 
    out$asymptotic_sampling_dist <- stanmat
  if (mer) 
    out$glmod <- object$glmod
  
  structure(out, class = c("stanreg", "glm", "lm"))
}
