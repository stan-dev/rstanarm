# This file is part of rstanarm.
# Copyright 2013 Stan Development Team
# rstanarm is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# rstanarm is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with rstanarm.  If not, see <http://www.gnu.org/licenses/>.

#' Fitting regularized linear models via Stan
#'
#' Full Bayesian inference for linear modeling with regularizing priors on the
#' model parameters that are driven by prior beliefs about \eqn{R^2}, the
#' proportion of variance in the outcome attributable to the predictors.
#' @export
#'
#' @param formula,data,subset Same as in \code{\link[stats]{lm}}
#' @param weights,na.action,method,model,x,y,qr,singular.ok,contrasts,offset 
#'   Also the same as in \code{\link[stats]{lm}} but rarely specified
#' @param eta Either a positive scalar or \code{NULL} (the default)
#' @param prior.R2 A numeric scalar; see the Details section
#' @param prior.what A character vector of length one; see the Details section
#' @param ... Further arguments passed to \code{\link[rstan]{stan}}
#'
#'
#' @details The \code{stan_glm} function is similar in syntax to the 
#'   \code{\link[stats]{lm}} function but rather than choosing the parameters
#'   to minimize the sum of squared residuals, samples from the posterior
#'   distribution are drawn using Markov Chain Monte Carlo. 
#'   
#'   The prior on all the parameters hinges on the prior beliefs about 
#'   \eqn{R^2}, the proportion of variance in the outcome attributable to the
#'   predictors, which is given a \code{\link[stats]{Beta}} prior with first shape 
#'   hyperparameter equal to \eqn{K / 2} --- where \eqn{K} is the number of 
#'   predictors excluding any constant term --- and second shape hyperparameter
#'   \code{eta}. Any supplied characteristic of this \code{\link[stats]{Beta}}
#'   distribution implies the value of \code{eta}. Typically, the user would
#'   specify a value of \code{R2}, in which case the \code{prior.what}
#'   indicates whether \code{R2} is indended to be the prior mode, mean,
#'   median, or natural logarithm of \eqn{R^2}. For example, if 
#'   \code{R2 = 0.5}, then the mode, mean, and median of the \code{\link[stats]{Beta}} 
#'   distribution are all the same and by implication, \code{eta = K / 2}. If
#'   \code{prior.what = "log"}, then \code{R2} must be a negative number,
#'   and \code{eta} is derived such that the expected natural logarithm of
#'   \eqn{R^2} is \code{R2}. Otherwise, \code{R2} should be positive and
#'   at most one.
#'   
#'   The smaller is \eqn{R^2}, the larger is the value of \code{eta}, and the 
#'   more concentrated near zero is the prior density for the regression 
#'   coefficients. Hence, the prior on the coefficients is regularizing and
#'   should yield a posterior distribution with good out-of-sample predictions
#'   \emph{if} the value of \code{eta} implied by \code{R2} is reasonable. In 
#'   addition to estimating \code{sigma} --- the standard deviation of the 
#'   normally-distributed errors --- this model estimates a positive parameter
#'   called \code{fit_ratio}. If it is greater than \eqn{1}, the posterior 
#'   predictive variance of the outcome will exceed the sample variance of the
#'   outcome by a multiplicative factor equal to the square of 
#'   \code{fit_ratio}. Conversely if \code{fit_ratio} is less than \eqn{1}, 
#'   then the model underfits.
#'   
#' @return An object of class \code{"stanreg"}, which is a list containing the 
#'   components
#' 
#' \describe{
#'   \item{coefficients}{named vector of coefficients (posterior means)}
#'   \item{residuals}{the residuals. For linear models \code{residuals}
#'    contains the response minus fitted values. Otherwise \code{residuals}
#'    contains the deviance residuals. See also \code{\link{residuals.stanreg}}}.
#'   \item{fitted.values}{the fitted mean values (for glms 
#'   the linear predictors are transformed by the invserse link function).}
#'   \item{linear.predictors}{the linear fit on the link scale (for linear models
#'   this is the same as \code{fitted.values}).}
#'   \item{covmat}{variance-covariance matrix for the coefficients (estimated
#'   from the posterior draws.)}
#'   \item{y}{if requested, the \code{y} vector used.}
#'   \item{x}{if requested, the model matrix.}
#'   \item{model}{if requested, the model frame.}
#'   \item{family}{the \code{\link[stats]{family}} object used.}
#'   \item{prior.weights}{any weights supplied by the user.}
#'   \item{df.residual}{the residual degrees of freedom}
#'   \item{call}{the matched call.}
#'   \item{formula}{the formula supplied.}
#'   \item{data}{the \code{data} argument.}
#'   \item{prior.info}{a list with information about the prior distributions
#'   used.}
#'   \item{stanfit}{the stanfit object returned by \code{\link[rstan]{stan}}}
#' } 
#' The accessor functions \code{coef}, \code{fitted}, and \code{resid}
#' can be used with objects of class \code{"stanreg"}. There are also  
#' \code{vcov}, \code{confint} and \code{\link{se}} methods.
#'
#' @seealso \code{\link[stats]{lm}}, \code{\link[rstan]{stan}}, and
#'   \code{\link{stan_glm}}, which --- if 
#'   \code{family = gaussian(link = "identity")} --- also estimates
#'   a linear model with normally-distributed errors but specifies
#'   different priors
#' @examples 
#' \dontrun{
#' stan_lm(mpg ~ ., data = mtcars, prior.R2 = 0.75)
#' }


stan_lm <- function(formula, data, subset, weights, na.action, method = "qr",
                    model = TRUE, x = FALSE, y = FALSE, qr = TRUE, 
                    singular.ok = TRUE, contrasts = NULL, offset, 
                    eta = NULL, prior.R2 = NULL, 
                    prior.what = c("mode", "mean", "median", "log"), ...) {
  
  call <- match.call()
  mf <- match.call(expand.dots = FALSE)
  mf[[1L]] <- as.name("lm")
  mf$x <- mf$y <- mf$singular.ok <- TRUE
  mf$qr <- FALSE
  mf$eta <- mf$prior.R2 <- mf$prior.what <- NULL
  
  ols <- suppressWarnings(eval(mf, parent.frame()))
  mt <- attr(ols, "terms")
  
  X <- ols$x
  if (colnames(X)[1] == "(Intercept)") {
    has_intercept <- 1L
    X <- X[,-1]
  }
  else has_intercept <- 0L
  
  J <- 1L
  N <- array(nrow(X), c(J))
  K <- ncol(X)
  if (K == 0) stop("'stan_lm' is not suitable for estimating a mean ",
                   "use 'stan_glm' with 'family = gaussian()' instead")
  sqrt_weights <- model.weights(mf)
  if (is.null(sqrt_weights)) sqrt_weights <- rep(1,N) 
  else {
    sqrt_weights <- sqrt(weights)
    X <- sqrt_weights * X
  }
  
  b <- coef(ols)
  b[is.na(b)] <- 0.0
  b <- array(b[-1], c(J,K))
  SSR <- array(crossprod(sqrt_weights * residuals(ols))[1], J)
  
  s_X <- array(apply(X, 2, sd), c(J,K))
  xbar <- array(colMeans(X), c(J,K))
  X <- sweep(X, 2, xbar, FUN = "-")
  XtX <- crossprod(X)
  dim(XtX) <- c(J, K, K)
  s_Y <- array(sd(ols$y), J)
  ybar <- array(mean(ols$y), J)
  
  if (is.null(eta)) {
    if (is.null(prior.R2)) 
      stop("the 'prior.R2' argument must be specified if 'eta' is unspecified")
    stopifnot(is.numeric(prior.R2))
    stopifnot(is.numeric(K), K > 0, K == as.integer(K))
    prior.what <- match.arg(prior.what)
    half_K <- K / 2
    if (prior.what == "mode") {
      stopifnot(prior.R2 > 0, prior.R2 <= 1)
      if (K <= 2) stop("mode of beta distribution does not exist when K <= 2 ",
                       "specify 'prior.what' as 'mean' or 'log' instead")
      eta <- (half_K - 1  - prior.R2 * half_K + prior.R2 * 2) / prior.R2
    }
    else if (prior.what == "mean") {
      stopifnot(prior.R2 > 0, prior.R2 <= 1)
      eta <- (half_K - prior.R2 * half_K) / prior.R2
    }
    else if (prior.what == "median") {
      stopifnot(prior.R2 > 0, prior.R2 <= 1)
      FUN <- function(eta) qbeta(0.5, half_K, qexp(eta)) - prior.R2
      eta <- qexp(uniroot(FUN, interval = 0:1)$root)
    }
    else { # prior.what == "log"
      stopifnot(prior.R2 < 0)
      FUN <- function(eta) digamma(half_K) - digamma(half_K + qexp(eta)) - prior.R2
      eta <- qexp(uniroot(FUN, interval = 0:1, 
                          f.lower = -prior.R2, f.upper = -.Machine$double.xmax)$root)
    }
  }
  else if (!is.numeric(eta) || length(eta) != 1L || eta <= 0) {
    stop("'eta' must be a positive scalar")
  }
  
  stanfit <- get("stanfit_lm")
  standata <- nlist(K, has_intercept, J, N, xbar, s_X, XtX, ybar, s_Y, b, SSR, eta)
  pars <- c(if (has_intercept) "alpha", "beta", "sigma", "omega")
  stanfit <- rstan::stan(fit = stanfit, data = standata, pars = pars, ...)
  parameters <- dimnames(stanfit)$parameters
  new_names <- c(if (has_intercept) "(Intercept)", colnames(X), 
                 "sigma", "fit_ratio", "log-posterior")
  stanfit@sim$fnames_oi <- new_names

  fit <- nlist(stanfit, family = gaussian(), formula, offset = model.offset(mf), 
               weights = ols$weights, x = ols$x, y = ols$y, 
               data, prior.info = eta, call = call, terms = mt, 
               model = if (model) mf else NULL,
               na.action = attr(ols, "na.action"), 
               contrasts = attr(X, "contrasts"))
  fit <- stanreg(fit)
  if (!x) fit$x <- NULL
  if (!y) fit$y <- NULL
  fit
}
