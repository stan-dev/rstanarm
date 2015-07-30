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
#' @param weights,na.action,method,model,qr,singular.ok,contrasts,offset 
#'   Also the same as in \code{\link[stats]{lm}} but rarely specified
#' @param x,y In \code{stan_lm}, logical scalars indicating whether to
#'   return the design matrix and response vector. In \code{stan_lm.fit},
#'   a design matrix and response vector.
#' @param w Same as in \code{\link[stats]{lm.wfit}} but rarely specified   
#' @param eta Either a positive scalar or \code{NULL} (the default)
#' @param prior.R2 A numeric scalar; see the Details section
#' @param prior.what A character vector of length one; see the Details section
#' @param prior_PD A logical scalar (defaulting to \code{FALSE}) indicating
#'   whether to draw from the prior predictive distribution instead of
#'   conditioning on the outcome. Note that if \code{TRUE}, the draws are
#'   merely proportional to the actual distribution because of an improper
#'   prior on a scale parameter
#' @param ... Further arguments passed to \code{\link[rstan]{stan}}
#'
#'
#' @details The \code{stan_lm} function is similar in syntax to the 
#'   \code{\link[stats]{lm}} function but rather than choosing the parameters
#'   to minimize the sum of squared residuals, samples from the posterior
#'   distribution are drawn using Markov Chain Monte Carlo. The 
#'   \code{stan_lm} function has a formula-based interface and would usually
#'   be called by users but the \code{stan_lm.fit} and \code{stan_lm.wfit}
#'   functions might be called by other functions that parse the data 
#'   themselves and are analagous to \code{\link[stats]{lm.fit}} and
#'   \code{\link[stats]{lm.wfit}} respectively.
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
#'   \emph{if} the value of \code{eta} implied by \code{R2} is reasonable. 
#'   
#'   In addition to estimating \code{sigma} --- the standard deviation of the
#'   normally-distributed errors --- this model estimates a positive parameter
#'   called \code{log-fit_ratio}. If it is positive, the marginal posterior 
#'   variance of the outcome will exceed the sample variance of the outcome
#'   by a multiplicative factor equal to the square of \code{fit_ratio}.
#'   Conversely if \code{log-fit_ratio} is negative, then the model underfits.
#'   Given the regularizing nature of the priors, a slight underfit is good.
#'   However, even if the \emph{marginal} posterior variance is off, the 
#'   \emph{conditional} variance may still be reasonable.
#'   
#'   Finally, the posterior predictive distribution is generated with the
#'   predictors fixed at their sample means. This quantity is useful for
#'   checking convergence because it is reasonably normally distributed
#'   and a function of all the parameters in the model.
#'   
#' @return The \code{stan_lm.fit} and \code{stan_lm.wfit} functions return an 
#'   object of class \code{\link[rstan]{stanfit-class}}. The more typically
#'   used \code{stan_lm} function returns an object of class \code{"stanreg"}, 
#'   which is a list containing the components
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
#'   \item{stanfit}{an object of \code{\link[rstan]{stanfit-class}}}
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
  
  modelframe <- suppressWarnings(eval(mf, parent.frame()))
  mt <- attr(modelframe, "terms")
  
  Y <- modelframe$y
  X <- modelframe$x
  w <- modelframe$weights
  offset <- model.offset(mf)
  stanfit <- stan_lm.wfit(y = Y, x = X, w, offset, method = "qr", singular.ok = TRUE,
                          eta = eta, prior.R2 = prior.R2, 
                          prior.what = prior.what, ...)
  
  K <- ncol(X) - (colnames(X)[1] == "(Intercept)")
  if (is.null(eta)) eta <- make_eta(prior.R2, prior.what, K)
  
  fit <- nlist(stanfit, family = gaussian(), formula, offset, 
               weights = w, x = X, y = Y, data,
               prior.info = nlist(dist = "beta", shape1 = K / 2, shape2 = eta), 
               call = call, terms = mt,
               model = if (model) mf else NULL,
               na.action = attr(modelframe, "na.action"),
               contrasts = attr(X, "contrasts"))
  fit <- stanreg(fit)
  if (!x) fit$x <- NULL
  if (!y) fit$y <- NULL
  fit
}
