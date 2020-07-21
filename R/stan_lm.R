# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2013, 2014, 2015, 2016, 2017 Trustees of Columbia University
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

#' Bayesian regularized linear models via Stan
#' 
#' \if{html}{\figure{stanlogo.png}{options: width="25px" alt="http://mc-stan.org/about/logo/"}}
#' Bayesian inference for linear modeling with regularizing priors on the model
#' parameters that are driven by prior beliefs about \eqn{R^2}, the proportion
#' of variance in the outcome attributable to the predictors. See
#' \code{\link{priors}} for an explanation of this critical point.
#' \code{\link{stan_glm}} with \code{family="gaussian"} also estimates a linear
#' model with normally-distributed errors and allows for various other priors on
#' the coefficients.
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
#' @template args-prior_PD
#' @template args-algorithm
#' @template args-adapt_delta
#'
#' @param w Same as in \code{lm.wfit} but rarely specified.
#' @param prior Must be a call to \code{\link{R2}} with its 
#'   \code{location} argument specified or \code{NULL}, which would
#'   indicate a standard uniform prior for the \eqn{R^2}.
#' @param prior_intercept Either \code{NULL} (the default) or a call to
#'   \code{\link{normal}}. If a \code{\link{normal}} prior is specified
#'   without a \code{scale}, then the standard deviation is taken to be
#'   the marginal standard deviation of the outcome divided by the square
#'   root of the sample size, which is legitimate because the marginal
#'   standard deviation of the outcome is a primitive parameter being
#'   estimated.
#'   
#'   \strong{Note:} If using a dense representation of the design matrix
#'   ---i.e., if the \code{sparse} argument is left at its default value of
#'   \code{FALSE}--- then the prior distribution for the intercept is set so it
#'   applies to the value \emph{when all predictors are centered}. If you prefer
#'   to specify a prior on the intercept without the predictors being
#'   auto-centered, then you have to omit the intercept from the
#'   \code{\link[stats]{formula}} and include a column of ones as a predictor,
#'   in which case some element of \code{prior} specifies the prior on it,
#'   rather than \code{prior_intercept}. Regardless of how
#'   \code{prior_intercept} is specified, the reported \emph{estimates} of the
#'   intercept always correspond to a parameterization without centered
#'   predictors (i.e., same as in \code{glm}).
#'
#'
#' @details The \code{stan_lm} function is similar in syntax to the 
#'   \code{\link[stats]{lm}} function but rather than choosing the parameters to
#'   minimize the sum of squared residuals, samples from the posterior 
#'   distribution are drawn using MCMC (if \code{algorithm} is
#'   \code{"sampling"}). The \code{stan_lm} function has a formula-based
#'   interface and would usually be called by users but the \code{stan_lm.fit}
#'   and \code{stan_lm.wfit} functions might be called by other functions that
#'   parse the data themselves and are analogous to \code{lm.fit}
#'   and \code{lm.wfit} respectively.
#'      
#'   In addition to estimating \code{sigma} --- the standard deviation of the
#'   normally-distributed errors --- this model estimates a positive parameter
#'   called \code{log-fit_ratio}. If it is positive, the marginal posterior 
#'   variance of the outcome will exceed the sample variance of the outcome
#'   by a multiplicative factor equal to the square of \code{fit_ratio}.
#'   Conversely if \code{log-fit_ratio} is negative, then the model underfits.
#'   Given the regularizing nature of the priors, a slight underfit is good.
#'   
#'   Finally, the posterior predictive distribution is generated with the
#'   predictors fixed at their sample means. This quantity is useful for
#'   checking convergence because it is reasonably normally distributed
#'   and a function of all the parameters in the model.
#'   
#'   The \code{stan_aov} function is similar to \code{\link[stats]{aov}}, but
#'   does a Bayesian analysis of variance that is basically equivalent to
#'   \code{stan_lm} with dummy variables. \code{stan_aov} has a somewhat
#'   customized \code{\link{print}} method that prints an ANOVA-like table in
#'   addition to the output printed for \code{stan_lm} models.
#'   
#'   
#' @references 
#' Lewandowski, D., Kurowicka D., and Joe, H. (2009). Generating random
#' correlation matrices based on vines and extended onion method. 
#' \emph{Journal of Multivariate Analysis}. \strong{100}(9), 1989--2001.
#' 
#' @seealso 
#' The vignettes for \code{stan_lm} and \code{stan_aov}, which have more
#' thorough descriptions and examples.
#' \url{http://mc-stan.org/rstanarm/articles/}
#' 
#' Also see \code{\link{stan_glm}}, which --- if \code{family =
#' gaussian(link="identity")} --- also estimates a linear model with
#' normally-distributed errors but specifies different priors.
#'   
#'   
#' @examples
#' (fit <- stan_lm(mpg ~ wt + qsec + am, data = mtcars, prior = R2(0.75), 
#'                 # the next line is only to make the example go fast enough
#'                 chains = 1, iter = 300, seed = 12345, refresh = 0))
#' plot(fit, "hist", pars = c("wt", "am", "qsec", "sigma"), 
#'      transformations = list(sigma = "log"))
#' 
stan_lm <- function(formula, data, subset, weights, na.action,
                    model = TRUE, x = FALSE, y = FALSE, 
                    singular.ok = TRUE, contrasts = NULL, offset, ...,
                    prior = R2(stop("'location' must be specified")), 
                    prior_intercept = NULL,
                    prior_PD = FALSE, 
                    algorithm = c("sampling", "meanfield", "fullrank"), 
                    adapt_delta = NULL) {
  
  algorithm <- match.arg(algorithm)
  validate_glm_formula(formula)
  data <- validate_data(data, if_missing = environment(formula))
  
  call <- match.call(expand.dots = TRUE)
  mf <- match.call(expand.dots = FALSE)
  mf[[1L]] <- as.name("lm")
  mf$data <- data
  mf$x <- mf$y <- mf$singular.ok <- TRUE
  mf$qr <- FALSE
  mf$prior <- mf$prior_intercept <- mf$prior_PD <- mf$algorithm <- 
    mf$adapt_delta <- NULL
  mf$method <- "model.frame"
  modelframe <- suppressWarnings(eval(mf, parent.frame()))
  mt <- attr(modelframe, "terms")
  Y <- model.response(modelframe, "numeric")
  X <- model.matrix(mt, modelframe, contrasts)
  w <- as.vector(model.weights(modelframe))
  offset <- as.vector(model.offset(modelframe))
  stanfit <- stan_lm.wfit(y = Y, x = X, w, offset, singular.ok = singular.ok,
                          prior = prior, prior_intercept = prior_intercept, 
                          prior_PD = prior_PD, 
                          algorithm = algorithm, adapt_delta = adapt_delta, 
                          ...)
  if (algorithm != "optimizing" && !is(stanfit, "stanfit")) return(stanfit)
  fit <- nlist(stanfit, family = gaussian(), formula, offset, weights = w,
               x = X[,intersect(colnames(X), dimnames(stanfit)[[3]]), drop = FALSE], 
               y = Y, 
               data = data,
               prior.info = prior, 
               algorithm, call, terms = mt,
               model = if (model) modelframe else NULL,
               na.action = attr(modelframe, "na.action"),
               contrasts = attr(X, "contrasts"), 
               stan_function = "stan_lm")
  out <- stanreg(fit)
  out$xlevels <- .getXlevels(mt, modelframe)
  if (!x) 
    out$x <- NULL
  if (!y) 
    out$y <- NULL
  
  return(out)
}
