# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2013, 2014, 2015, 2016, 2017 Trustees of Columbia University
# Copyright (C) 2016, 2017 Sam Brilleman
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

#' Bayesian multivariate generalized linear models with correlated 
#' group-specific terms via Stan
#' 
#' \if{html}{\figure{stanlogo.png}{options: width="25px" alt="http://mc-stan.org/about/logo/"}}
#' Bayesian inference for multivariate GLMs with group-specific coefficients 
#' that are assumed to be correlated across the GLM submodels.
#' 
#' @export
#' @template args-dots
#' @template args-prior_PD
#' @template args-algorithm
#' @template args-adapt_delta
#' @template args-max_treedepth
#' @template args-QR
#' @template args-sparse
#' 
#' @param formula A two-sided linear formula object describing both the 
#'   fixed-effects and random-effects parts of the longitudinal submodel  
#'   similar in vein to formula specification in the \strong{lme4} package
#'   (see \code{\link[lme4]{glmer}} or the \strong{lme4} vignette for details). 
#'   Note however that the double bar (\code{||}) notation is not allowed 
#'   when specifying the random-effects parts of the formula, and neither
#'   are nested grouping factors (e.g. \code{(1 | g1/g2))} or 
#'   \code{(1 | g1:g2)}, where \code{g1}, \code{g2} are grouping factors. 
#'   For a multivariate GLM this should be a list of such formula objects, 
#'   with each element of the list providing the formula for one of the 
#'   GLM submodels.
#' @param data A data frame containing the variables specified in
#'   \code{formula}. For a multivariate GLM, this can
#'   be either a single data frame which contains the data for all 
#'   GLM submodels, or it can be a list of data frames where each
#'   element of the list provides the data for one of the GLM submodels.
#' @param family The family (and possibly also the link function) for the 
#'   GLM submodel(s). See \code{\link[lme4]{glmer}} for details. 
#'   If fitting a multivariate GLM, then this can optionally be a
#'   list of families, in which case each element of the list specifies the
#'   family for one of the GLM submodels. In other words, a different family
#'   can be specified for each GLM submodel. 
#' @param weights Same as in \code{\link[stats]{glm}},
#'   except that when fitting a multivariate GLM and a list of data frames 
#'   is provided in \code{data} then a corresponding list of weights 
#'   must be provided. If weights are 
#'   provided for one of the GLM submodels, then they must be provided for 
#'   all GLM submodels.
#' @param prior,prior_intercept,prior_aux Same as in \code{\link{stan_glmer}}
#'   except that for a multivariate GLM a list of priors can be provided for 
#'   any of \code{prior}, \code{prior_intercept} or \code{prior_aux} arguments. 
#'   That is, different priors can optionally be specified for each of the GLM  
#'   submodels. If a list is not provided, then the same prior distributions are 
#'   used for each GLM submodel. Note that the \code{"product_normal"} prior is
#'   not allowed for \code{stan_mvmer}.
#' @param prior_covariance Cannot be \code{NULL}; see \code{\link{priors}} for
#'   more information about the prior distributions on covariance matrices.
#'   Note however that the default prior for covariance matrices in 
#'   \code{stan_mvmer} is slightly different to that in \code{\link{stan_glmer}} 
#'   (the details of which are described on the \code{\link{priors}} page).
#' @param init The method for generating initial values. See
#'   \code{\link[rstan]{stan}}.
#'   
#' @details The \code{stan_mvmer} function can be used to fit a multivariate
#'   generalized linear model (GLM) with group-specific terms. The model consists
#'   of distinct GLM submodels, each which contains group-specific terms; within
#'   a grouping factor (for example, patient ID) the grouping-specific terms are
#'   assumed to be correlated across the different GLM submodels. It is 
#'   possible to specify a different outcome type (for example a different
#'   family and/or link function) for each of the GLM submodels. \cr
#'   \cr
#'   Bayesian estimation of the model is performed via MCMC, in the same way as 
#'   for \code{\link{stan_glmer}}. Also, similar to \code{\link{stan_glmer}},
#'   an unstructured covariance matrix is used for the group-specific terms 
#'   within a given grouping factor, with priors on the terms of a decomposition
#'   of the covariance matrix.See \code{\link{priors}} for more information about 
#'   the priors distributions that are available for the covariance matrices, 
#'   the regression coefficients and the intercept and auxiliary parameters.
#'
#' @return A \link[=stanreg-objects]{stanmvreg} object is returned.
#' 
#' @seealso \code{\link{stan_glmer}}, \code{\link{stan_jm}},
#'   \code{\link{stanreg-objects}}, \code{\link{stanmvreg-methods}}, 
#'   \code{\link{print.stanmvreg}}, \code{\link{summary.stanmvreg}},
#'   \code{\link{posterior_predict}}, \code{\link{posterior_interval}}.
#'    
#' @examples
#' \donttest{
#' #####
#' # A multivariate GLM with two submodels. For the grouping factor 'id', the 
#' # group-specific intercept from the first submodel (logBili) is assumed to
#' # be correlated with the group-specific intercept and linear slope in the 
#' # second submodel (albumin)
#' f1 <- stan_mvmer(
#'         formula = list(
#'           logBili ~ year + (1 | id), 
#'           albumin ~ sex + year + (year | id)),
#'         data = pbcLong, 
#'         # this next line is only to keep the example small in size!
#'         chains = 1, cores = 1, seed = 12345, iter = 1000)
#' summary(f1) 
#' 
#' #####
#' # A multivariate GLM with one bernoulli outcome and one
#' # gaussian outcome. We will artificially create the bernoulli
#' # outcome by dichotomising log serum bilirubin
#' pbcLong$ybern <- as.integer(pbcLong$logBili >= mean(pbcLong$logBili))
#' f2 <- stan_mvmer(
#'         formula = list(
#'           ybern ~ year + (1 | id), 
#'           albumin ~ sex + year + (year | id)),
#'         data = pbcLong,
#'         family = list(binomial, gaussian),
#'         chains = 1, cores = 1, seed = 12345, iter = 1000)
#' }
#' 
stan_mvmer <- function(formula, data, family = gaussian, weights,				          
                       prior = normal(autoscale=TRUE), prior_intercept = normal(autoscale=TRUE), 
                       prior_aux = cauchy(0, 5, autoscale=TRUE),
                       prior_covariance = lkj(autoscale=TRUE), prior_PD = FALSE, 
                       algorithm = c("sampling", "meanfield", "fullrank"), 
                       adapt_delta = NULL, max_treedepth = 10L, 
                       init = "random", QR = FALSE, sparse = FALSE, ...) {
  
  #-----------------------------
  # Pre-processing of arguments
  #-----------------------------  
  
  algorithm <- match.arg(algorithm)
  
  if (missing(weights)) weights <- NULL
  
  if (!is.null(weights)) 
    stop("'weights' are not yet implemented.")
  if (QR)               
    stop("'QR' decomposition is not yet implemented.")
  if (sparse)
    stop("'sparse' option is not yet implemented.")
  
  # Formula
  formula <- validate_arg(formula, "formula"); M <- length(formula)
	if (M > 3L)
	  stop("'stan_mvmer' is currently limited to a maximum of 3 outcomes.")
  
  # Data
  data <- validate_arg(data, "data.frame", validate_length = M)  
  data <- xapply(formula, data, FUN = get_all_vars) # drop additional vars
  
  # Family
  ok_classes <- c("function", "family", "character")
  ok_families <- c("binomial", "gaussian", "Gamma", 
                   "inverse.gaussian", "poisson", "neg_binomial_2")
  family <- validate_arg(family, ok_classes, validate_length = M)
  family <- lapply(family, validate_famlink, ok_families)

  # Observation weights
  if (!is.null(weights)) {
    if (!is(weights, "list")) 
      weights <- rep(list(weights), M)
    weights <- lapply(weights, validate_weights)
  }
  
  # Is prior* already a list?
  prior <- broadcast_prior(prior, M)
  prior_intercept <- broadcast_prior(prior_intercept, M)
  prior_aux <- broadcast_prior(prior_aux, M)
  
  #-----------
  # Fit model
  #----------- 
  
  stanfit <- stan_jm.fit(formulaLong = formula, dataLong = data, family = family,
                         weights = weights, priorLong = prior, 
                         priorLong_intercept = prior_intercept, priorLong_aux = prior_aux, 
                         prior_covariance = prior_covariance, prior_PD = prior_PD, 
                         algorithm = algorithm, adapt_delta = adapt_delta, 
                         max_treedepth = max_treedepth, init = init, 
                         QR = QR, sparse = sparse, ...)
  if (algorithm != "optimizing" && !is(stanfit, "stanfit")) return(stanfit)
  y_mod <- attr(stanfit, "y_mod")
  cnms  <- attr(stanfit, "cnms")
  flevels <- attr(stanfit, "flevels")
  prior_info <- attr(stanfit, "prior_info")
  stanfit <- drop_attributes(stanfit, "y_mod", "cnms", "flevels", "prior_info")
  
  terms <- fetch(y_mod, "terms")
  n_yobs <- fetch_(y_mod, "x", "N")
  n_grps <- sapply(flevels, n_distinct)
  
  fit <- nlist(stanfit, formula, family, weights, M, cnms, flevels, n_grps, n_yobs, 
               algorithm, terms, glmod = y_mod, data, prior.info = prior_info, 
               stan_function = "stan_mvmer", call = match.call(expand.dots = TRUE))
  
  out <- stanmvreg(fit)
  return(out)
}

