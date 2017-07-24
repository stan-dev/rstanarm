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

#' Bayesian generalized linear models with group-specific terms via Stan
#' 
#' Bayesian inference for GLMs with group-specific coefficients that have 
#' unknown covariance matrices with flexible priors.
#' 
#' @export
#' @templateVar armRef (Ch. 11-15)
#' @templateVar fun stan_glmer, stan_lmer, stan_glmer.nb
#' @templateVar pkg lme4
#' @templateVar pkgfun glmer
#' @template return-stanreg-object
#' @template see-also
#' @template args-priors
#' @template args-prior_aux
#' @template args-prior_covariance
#' @template args-prior_PD
#' @template args-algorithm
#' @template args-adapt_delta
#' @template args-QR
#' @template args-sparse
#' @template reference-gelman-hill
#' 
#' @param formula,data,family Same as for \code{\link[lme4]{glmer}}. \emph{We
#'   strongly advise against omitting the \code{data} argument}. Unless 
#'   \code{data} is specified (and is a data frame) many post-estimation 
#'   functions (including \code{update}, \code{loo}, \code{kfold}) are not 
#'   guaranteed to work properly.
#' @param subset,weights,offset Same as \code{\link[stats]{glm}}.
#' @param na.action,contrasts Same as \code{\link[stats]{glm}}, but rarely 
#'   specified.
#' @param ... For \code{stan_glmer}, further arguments passed to 
#'   \code{\link[rstan]{sampling}} (e.g. \code{iter}, \code{chains}, 
#'   \code{cores}, etc.) or to \code{\link[rstan]{vb}} (if \code{algorithm} is 
#'   \code{"meanfield"} or \code{"fullrank"}). For \code{stan_lmer} and 
#'   \code{stan_glmer.nb}, \code{...} should also contain all relevant arguments
#'   to pass to \code{stan_glmer} (except \code{family}).
#'
#' @details The \code{stan_glmer} function is similar in syntax to 
#'   \code{\link[lme4]{glmer}} but rather than performing (restricted) maximum 
#'   likelihood estimation of generalized linear models, Bayesian estimation is 
#'   performed via MCMC. The Bayesian model adds priors on the 
#'   regression coefficients (in the same way as \code{\link{stan_glm}}) and
#'   priors on the terms of a decomposition of the covariance matrices of the
#'   group-specific parameters. See \code{\link{priors}} for more information
#'   about the priors.
#'   
#'   The \code{stan_lmer} function is equivalent to \code{stan_glmer} with 
#'   \code{family = gaussian(link = "identity")}. 
#'   
#'   The \code{stan_glmer.nb} function, which takes the extra argument 
#'   \code{link}, is a wrapper for \code{stan_glmer} with \code{family = 
#'   \link{neg_binomial_2}(link)}.
#'   
#'   
#' @seealso The vignette for \code{stan_glmer} and the \emph{Hierarchical 
#'   Partial Pooling} vignette.
#'    
#' @examples
#' # see help(example_model) for details on the model below
#' if (!exists("example_model")) example(example_model) 
#' print(example_model, digits = 1)
#' 
#' @importFrom lme4 glFormula
#' @importFrom Matrix Matrix t cBind
stan_glmer <- function(formula, data = NULL, family = gaussian, 
                       subset, weights, 
                       na.action = getOption("na.action", "na.omit"), 
                       offset, contrasts = NULL, ...,
                       prior = normal(), prior_intercept = normal(),
                       prior_aux = cauchy(0, 5),
                       prior_covariance = decov(), prior_PD = FALSE, 
                       algorithm = c("sampling", "meanfield", "fullrank"), 
                       adapt_delta = NULL, QR = FALSE, sparse = FALSE) {
  
  call <- match.call(expand.dots = TRUE)
  mc <- match.call(expand.dots = FALSE)
  data <- validate_data(data) #, if_missing = environment(formula))
  family <- validate_family(family)
  mc[[1]] <- quote(lme4::glFormula)
  mc$control <- make_glmerControl()
  mc$prior <- mc$prior_intercept <- mc$prior_covariance <- mc$prior_aux <-
    mc$prior_PD <- mc$algorithm <- mc$scale <- mc$concentration <- mc$shape <-
    mc$adapt_delta <- mc$... <- mc$QR <- mc$sparse <- NULL
  glmod <- eval(mc, parent.frame())
  X <- glmod$X
  y <- glmod$fr[, as.character(glmod$formula[2L])]
  if (is.matrix(y) && ncol(y) == 1L)
    y <- as.vector(y)

  offset <- model.offset(glmod$fr) %ORifNULL% double(0)
  weights <- validate_weights(weights)
  if (is.null(prior)) 
    prior <- list()
  if (is.null(prior_intercept)) 
    prior_intercept <- list()
  if (is.null(prior_aux)) 
    prior_aux <- list()
  if (is.null(prior_covariance))
    stop("'prior_covariance' can't be NULL.", call. = FALSE)
  group <- glmod$reTrms
  group$decov <- prior_covariance
  algorithm <- match.arg(algorithm)
  stanfit <- stan_glm.fit(x = X, y = y, weights = weights,
                          offset = offset, family = family,
                          prior = prior, prior_intercept = prior_intercept,
                          prior_aux = prior_aux, prior_PD = prior_PD, 
                          algorithm = algorithm, adapt_delta = adapt_delta,
                          group = group, QR = QR, sparse = sparse, ...)
  if (family$family == "Beta regression") family$family <- "beta"
  Z <- pad_reTrms(Ztlist = group$Ztlist, cnms = group$cnms, 
                  flist = group$flist)$Z
  colnames(Z) <- b_names(names(stanfit), value = TRUE)
  
  fit <- nlist(stanfit, family, formula, offset, weights, 
               x = if (getRversion() < "3.2.0") cBind(X, Z) else cbind2(X, Z), 
               y = y, data, call, terms = NULL, model = NULL, 
               na.action = attr(glmod$fr, "na.action"), contrasts, algorithm, glmod, 
               stan_function = "stan_glmer")
  out <- stanreg(fit)
  class(out) <- c(class(out), "lmerMod")
  
  return(out)
}


#' @rdname stan_glmer
#' @export
stan_lmer <- function(formula,
                      data = NULL,
                      subset,
                      weights,
                      na.action = getOption("na.action", "na.omit"),
                      offset,
                      contrasts = NULL,
                      ...,
                      prior = normal(),
                      prior_intercept = normal(),
                      prior_aux = cauchy(0, 5),
                      prior_covariance = decov(),
                      prior_PD = FALSE,
                      algorithm = c("sampling", "meanfield", "fullrank"),
                      adapt_delta = NULL,
                      QR = FALSE) {
  if ("family" %in% names(list(...)))
    stop("'family' should not be specified.")
  mc <- call <- match.call(expand.dots = TRUE)
  if (!"formula" %in% names(call))
    names(call)[2L] <- "formula"
  mc[[1L]] <- quote(stan_glmer)
  mc$REML <- NULL
  mc$family <- "gaussian"
  out <- eval(mc, parent.frame())
  out$call <- call
  out$stan_function <- "stan_lmer"
  return(out)
}


#' @rdname stan_glmer
#' @export
#' @param link For \code{stan_glmer.nb} only, the link function to use. See 
#'   \code{\link{neg_binomial_2}}.
#' 
stan_glmer.nb <- function(formula,
                          data = NULL,
                          subset,
                          weights,
                          na.action = getOption("na.action", "na.omit"),
                          offset,
                          contrasts = NULL,
                          link = "log",
                          ...,
                          prior = normal(),
                          prior_intercept = normal(),
                          prior_aux = cauchy(0, 5),
                          prior_covariance = decov(),
                          prior_PD = FALSE,
                          algorithm = c("sampling", "meanfield", "fullrank"),
                          adapt_delta = NULL,
                          QR = FALSE) {
  if ("family" %in% names(list(...)))
    stop("'family' should not be specified.")
  mc <- call <- match.call(expand.dots = TRUE)
  if (!"formula" %in% names(call))
    names(call)[2L] <- "formula"
  mc[[1]] <- quote(stan_glmer)
  mc$REML <- mc$link <- NULL
  mc$family <- neg_binomial_2(link = link)
  out <- eval(mc, parent.frame())
  out$call <- call
  out$stan_function <- "stan_glmer.nb"
  return(out)
}
