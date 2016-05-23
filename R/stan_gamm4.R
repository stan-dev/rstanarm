# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2015, 2016 Trustees of Columbia University
# Copyright (C) 2009, 2010, 2011, 2012, 2013 Simon N. Wood
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

#' Bayesian generalized linear additive models with group-specific terms via
#' Stan
#' 
#' Bayesian inference for GAMMs with flexible priors.
#' 
#' @export
#' @templateVar fun stan_gamm4
#' @templateVar pkg gamm4
#' @templateVar pkgfun gamm4
#' @template return-stanreg-object
#' @template see-also
#' @template args-priors
#' @template args-prior_PD
#' @template args-algorithm
#' @template args-adapt_delta
#' @template args-QR
#' 
#' @param formula,random,family,data,knots,drop.unused.levels Same as for 
#'   \code{\link[gamm4]{gamm4}}.
#' @param subset,weights,na.action Same as \code{\link[stats]{glm}}, 
#'   but rarely specified.
#' @param ... Further arguments passed to \code{\link[rstan]{sampling}} (e.g. 
#'   \code{iter}, \code{chains}, \code{cores}, etc.) or to
#'   \code{\link[rstan]{vb}} (if \code{algorithm} is \code{"meanfield"} or
#'   \code{"fullrank"}).
#' @param prior_covariance Cannot be \code{NULL}; see \code{\link{decov}} for
#'   more information about the default arguments.
#'
#' @details The \code{stan_gamm4} function is similar in syntax to 
#'   \code{\link[gamm4]{gamm4}}, which accepts a syntax that is similar to (but 
#'   not quite as extensive as) that for \code{\link[mgcv]{gamm}} and converts 
#'   it internally into the syntax accepted by \code{\link[lme4]{glmer}}. But 
#'   rather than performing (restricted) maximum likelihood estimation, the 
#'   \code{stan_gamm4} function utilizes MCMC to perform Bayesian estimation. 
#'   The Bayesian model adds independent priors on the common regression 
#'   coefficients (in the same way as \code{\link{stan_glm}}) and priors on the 
#'   terms of a decomposition of the covariance matrices of the group-specific 
#'   parameters, including the smooths. Estimating these models via MCMC avoids
#'   the optimization issues that often crop up with GAMMs and provides better
#'   estimates for the uncertainty in the parameter estimates. 
#'   
#'   See \code{\link[gamm4]{gamm4}} for more information about the model
#'   specicification and \code{\link{priors}} for more information about the
#'   priors.
#' @examples
#' # see example(gamm4, package = "gamm4") but prefix gamm4() calls with stan_
#' 
stan_gamm4 <- function(formula, random = NULL, family = gaussian(), data = list(), 
                       weights = NULL, subset = NULL, na.action, knots = NULL, 
                       drop.unused.levels = TRUE, ..., 
                       prior = normal(), prior_intercept = normal(),
                       prior_ops = prior_options(),
                       prior_covariance = decov(), prior_PD = FALSE, 
                       algorithm = c("sampling", "meanfield", "fullrank"), 
                       adapt_delta = NULL, QR = FALSE) {

  mc <- match.call(expand.dots = FALSE)
  glmod <- gamm4_to_glmer(formula, random, family, data, weights, subset, 
                          na.action, knots, drop.unused.levels)

  X <- glmod$X
  y <- glmod$fr[, as.character(glmod$formula[2L])]
  if (is.matrix(y) && ncol(y) == 1L) y <- as.vector(y)
  offset <- model.offset(glmod$fr) %ORifNULL% double(0)
  weights <- validate_weights(glmod$fr$weights)
  
  if (is.null(prior))
    prior <- list()
  if (is.null(prior_intercept)) 
    prior_intercept <- list()
  if (!length(prior_ops)) 
    prior_ops <- list(scaled = FALSE, prior_scale_for_dispersion = Inf)
  

  group <- glmod$reTrms
  group$decov <- prior_covariance
  algorithm <- match.arg(algorithm)
  
  stanfit <- stan_glm.fit(x = X, y = y, weights = weights,
                          offset = offset, family = family,
                          prior = prior, prior_intercept = prior_intercept,
                          prior_ops = prior_ops, prior_PD = prior_PD, 
                          algorithm = algorithm, adapt_delta = adapt_delta,
                          group = group, QR = QR, ...)
  
  Z <- pad_reTrms(Z = t(group$Zt), cnms = group$cnms, 
                  flist = group$flist)$Z
  colnames(Z) <- b_names(names(stanfit), value = TRUE)
  fit <- nlist(stanfit, family, formula, offset, weights, 
               x = if (getRversion() < "3.2.0") cBind(X, Z) else cbind2(X, Z), 
               y = y, data, call = mc, terms = NULL, model = NULL, 
               prior.info = get_prior_info(call, formals()),
               algorithm, glmod = glmod)
  out <- stanreg(fit)
  class(out) <- c(class(out), "lmerMod", "gamm4")
  return(out)
  # TODO: maybe convert back to gam parameterization?
}

