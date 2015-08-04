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
#' proportion of variance in the outcome attributable to the predictors. See
#' \code{\link{priors}} for an explanation of this critical point.
#' @export
#'
#' @param formula,data,subset Same as in \code{\link[stats]{lm}}
#' @param weights,na.action,method,model,qr,singular.ok,contrasts,offset 
#'   Also the same as in \code{\link[stats]{lm}} but rarely specified
#' @param x,y In \code{stan_lm}, logical scalars indicating whether to
#'   return the design matrix and response vector. In \code{stan_lm.fit},
#'   a design matrix and response vector.
#' @param w Same as in \code{\link[stats]{lm.wfit}} but rarely specified   
#' @param prior Must be a call to \code{\link{LKJ}} with its 
#'   \code{location} argument specified
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
#' stan_lm(mpg ~ ., data = mtcars, prior = LKJ(location = 0.75))
#' }

stan_lm <- function(formula, data, subset, weights, na.action, method = "qr",
                    model = TRUE, x = FALSE, y = FALSE, qr = TRUE, 
                    singular.ok = TRUE, contrasts = NULL, offset, 
                    prior = LKJ(), prior_PD = FALSE, ...) {
  
  call <- match.call()
  mf <- match.call(expand.dots = FALSE)
  mf[[1L]] <- as.name("lm")
  mf$x <- mf$y <- mf$singular.ok <- TRUE
  mf$qr <- FALSE
  mf$prior <- NULL
  
  modelframe <- suppressWarnings(eval(mf, parent.frame()))
  mt <- attr(modelframe, "terms")
  
  Y <- modelframe$y
  X <- modelframe$x
  w <- modelframe$weights
  offset <- model.offset(mf)

  stanfit <- stan_lm.wfit(y = Y, x = X, w, offset, method = "qr", singular.ok = TRUE,
                          prior = prior,  prior_PD = prior_PD, ...)
  
  fit <- nlist(stanfit, family = gaussian(), formula, offset, 
               weights = w, x = X, y = Y, data,
               prior.info = prior, 
               call = call, terms = mt,
               model = if (model) mf else NULL,
               na.action = attr(modelframe, "na.action"),
               contrasts = attr(X, "contrasts"))
  fit <- stanreg(fit)
  if (!x) fit$x <- NULL
  if (!y) fit$y <- NULL
  fit
}
