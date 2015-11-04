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

#' Regularized linear models via Stan
#'
#' Bayesian inference for linear modeling with regularizing priors on the 
#' model parameters that are driven by prior beliefs about \eqn{R^2}, the 
#' proportion of variance in the outcome attributable to the predictors. See 
#' \code{\link{priors}} for an explanation of this critical point.
#' 
#' 
#' @export
#' @templateVar fun stan_lm, stan_aov
#' @templateVar fitfun stan_lm.fit or stan_lm.wfit
#' @templateVar pkg stats
#' @templateVar pkgfun lm
#' @templateVar rareargs model,offset,weights
#' @templateVar rareargs2 na.action,singular.ok,contrasts
#' @template return-stanreg-object
#' @template return-stanfit-object
#' @template args-formula-data-subset
#' @template args-same-as-rarely
#' @template args-same-as-rarely-2
#' @template args-x-y
#' @template args-dots
#' @template args-algorithm
#' @template args-adapt_delta
#'
#' @param w Same as in \code{\link[stats]{lm.wfit}} but rarely specified.
#' @param prior Must be a call to \code{\link{R2}} with its 
#'   \code{location} argument specified.
#' @param prior_PD A logical scalar (defaulting to \code{FALSE}) indicating
#'   whether to draw from the prior predictive distribution instead of
#'   conditioning on the outcome. Note that if \code{TRUE}, the draws are
#'   merely proportional to the actual distribution because of an improper
#'   prior on a scale parameter.
#'
#'
#' @details The \code{stan_lm} function is similar in syntax to the 
#'   \code{\link[stats]{lm}} function but rather than choosing the parameters
#'   to minimize the sum of squared residuals, samples from the posterior
#'   distribution (if \code{algorithm = "sampling"}) are drawn using MCMC. The
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
#'   \emph{conditional} variance (\code{sigma}) may still be reasonable.
#'   
#'   Finally, the posterior predictive distribution is generated with the
#'   predictors fixed at their sample means. This quantity is useful for
#'   checking convergence because it is reasonably normally distributed
#'   and a function of all the parameters in the model.
#'   
#'   The \code{stan_aov} function is similar to \code{\link[stats]{aov}} and
#'   has a somewhat customized \code{\link{print}} method but basically just 
#'   calls \code{stan_lm} with dummy variables to do a Bayesian analysis of
#'   variance.
#'   
#' 
#' @seealso 
#' The vignettes for \code{stan_lm} and \code{stan_aov}, which have more
#' thorough descriptions and examples.
#' 
#' Also see \code{\link{stan_glm}}, which --- if \code{family =
#' gaussian(link = "identity")} --- also estimates a linear model with
#' normally-distributed errors but specifies different priors.
#'   
#'   
#' @examples
#' \dontrun{
#' (fit <- stan_lm(mpg ~ wt + qsec + am, data = mtcars, prior = R2(0.75), seed = 12345))
#' plot(fit)
#' ppcheck(fit, check = "dist", nreps = 25)
#' }
#' 
stan_lm <- function(formula, data, subset, weights, na.action,
                    model = TRUE, x = FALSE, y = FALSE, 
                    singular.ok = TRUE, contrasts = NULL, offset, ...,
                    prior = R2(stop("'location' must be specified")), 
                    prior_PD = FALSE, algorithm = c("sampling", "optimizing"), 
                    adapt_delta = 0.95) {
  
  call <- match.call()
  mf <- match.call(expand.dots = FALSE)
  mf[[1L]] <- as.name("lm")
  mf$x <- mf$y <- mf$singular.ok <- TRUE
  mf$qr <- FALSE
  mf$prior <- NULL
  
  modelframe <- suppressWarnings(eval(mf, parent.frame()))
  mt <- modelframe$terms
  Y <- modelframe$y
  X <- modelframe$x
  if (!singular.ok) X <- X[,!is.na(modelframe$coefficients),drop = FALSE]
  w <- modelframe$weights
  offset <- model.offset(mf)
  stanfit <- stan_lm.wfit(y = Y, x = X, w, offset, singular.ok = TRUE,
                          prior = prior,  prior_PD = prior_PD, 
                          algorithm = algorithm, adapt_delta = adapt_delta, ...)

  fit <- nlist(stanfit, family = gaussian(), formula, offset, 
               weights = w, x = X, y = Y, data,
               prior.info = prior, algorithm = "sampling",
               call = call, terms = mt,
               model = if (model) model.frame(modelframe) else NULL,
               na.action = attr(modelframe, "na.action"),
               contrasts = attr(X, "contrasts"))
  fit <- stanreg(fit)
  if (!x) fit$x <- NULL
  if (!y) fit$y <- NULL
  fit
}
