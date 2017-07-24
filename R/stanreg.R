# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2015, 2016, 2017 Trustees of Columbia University
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
  
  is_betareg <- is.beta(family$family)
  if (is_betareg) { 
    family_phi <- object$family_phi  # pull out phi family/link
    if (is.null(family_phi)) {
      family_phi <- beta_fam("log")
      z <- matrix(1, nrow = nobs, ncol = 1, dimnames = list(NULL, "(Intercept)"))
    }
    else z <- object$z   # pull out betareg z vars so that they can be used in posterior_predict/loo
    nvars_z <- NCOL(z)   # used so that all coefficients are printed with coef()
  }
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
    if (is_betareg) {
      if (length(colnames(z)) == 1)
        coefs_z <- apply(stanmat[, grepl("(phi)", colnames(stanmat), fixed = TRUE), drop = FALSE], 2L, median)
      else
        coefs_z <- apply(stanmat[, paste0("(phi)_",colnames(z)), drop = FALSE], 2L, median)
    }
  } else {
    stan_summary <- make_stan_summary(stanfit)
    coefs <- stan_summary[1:nvars, select_median(object$algorithm)]
    if (is_betareg) {
      coefs_z <- stan_summary[(nvars + 1):(nvars + nvars_z), select_median(object$algorithm)]
      if (length(coefs_z) == 1L)
        names(coefs_z) <- rownames(stan_summary)[nvars + 1]
    }
    if (length(coefs) == 1L) # ensures that if only a single coef it still gets a name
      names(coefs) <- rownames(stan_summary)[1L]

    if (is_betareg) {
      stanmat <- as.matrix(stanfit)[,c(names(coefs),names(coefs_z)), drop = FALSE]
      colnames(stanmat) <- c(names(coefs),names(coefs_z))
    } else {
      stanmat <- as.matrix(stanfit)[, 1:nvars, drop = FALSE]
      colnames(stanmat) <- colnames(x)
    }
    ses <- apply(stanmat, 2L, mad)
    if (mer) {
      mark <- sum(sapply(object$stanfit@par_dims[c("alpha", "beta")], prod))
      stanmat <- stanmat[,1:mark, drop = FALSE]
    }
    covmat <- cov(stanmat)
    # rownames(covmat) <- colnames(covmat) <- rownames(stan_summary)[1:nrow(covmat)]
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
  if (is_betareg) {
    eta_z <- linear_predictor(coefs_z, z, object$offset)
    phi <- family_phi$linkinv(eta_z)
  }
  
  out <- nlist(
    coefficients = unpad_reTrms(coefs), 
    ses = unpad_reTrms(ses),
    fitted.values = mu,
    linear.predictors = eta,
    residuals, 
    df.residual = if (opt) df.residual else NA_integer_, 
    # covmat = unpad_reTrms(unpad_reTrms(covmat, col = TRUE), col = FALSE),
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
    formula = object$formula, 
    terms = object$terms,
    prior.info = attr(stanfit, "prior.info"),
    algorithm = object$algorithm,
    stan_summary,  
    stanfit = if (opt) stanfit$stanfit else stanfit,
    rstan_version = utils::packageVersion("rstan"),
    call = object$call, 
    # sometimes 'call' is no good (e.g. if using do.call(stan_glm, args)) so
    # also include the name of the modeling function (for use when printing,
    # etc.)
    stan_function = object$stan_function
  )

  if (opt) 
    out$asymptotic_sampling_dist <- stanmat
  if (mer) 
    out$glmod <- object$glmod
  if (is_betareg) {
    out$coefficients <- unpad_reTrms(c(coefs, coefs_z))
    out$z <- z
    out$family_phi <- family_phi
    out$eta_z <- eta_z
    out$phi <- phi
  }
  
  structure(out, class = c("stanreg", "glm", "lm"))
}
