# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2016 Trustees of Columbia University
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

#' Bayesian regularized linear but stratified models via Stan
#' 
#' This model uses the same likelihood and syntax as the  \code{\link[lme4]{lmList}} 
#' function in the \pkg{lme4} package, which allows the parameters of a linear
#' regression model to vary according to a mutually exclusive and exhaustive grouping
#' factor, but differs fundamentally in that group-level parameters have hyperpriors
#' that shrink --- to some estimated extent --- the group-level parameters toward
#' common parameters. Thus, \code{stan_lmList} is more like a \code{\link[lme4]{lmer}}
#' model with a sinking grouping variable that pertains too all the regression
#' coefficients and the intercept (and, unlike \code{\link[lme4]{lmer}}, the
#' standard deviation of the errors).
#' 
#' @export
#' @param formula,data,family Same as in \code{\link[lme4]{lmList}} but 
#'   \code{family} must be omitted or \code{\link[stats]{gaussian}}
#' @param subset,weights,na.action,offset Same as in \code{\link[lme4]{lmList}}
#'   but rarely specified
#' @template args-dots
#' @param prior,prior_intercept Same as in \code{\link{stan_lm}}
#' @template args-kappa_mean
#' @template args-prior_PD
#' @template args-algorithm
#' @template args-adapt_delta
#' 
#' @examples 
#' \dontrun{
#' post <- stan_lmList(mpg ~ disp + wt + hp | cyl, data = mtcars, 
#'                     prior = R2(0.25))
#' coef(post) # posterior medians stratified by number of cylinders
#' }
stan_lmList <- function(formula, data, family, subset, weights, na.action, offset, 
                        ..., prior = R2(stop("'location' must be specified")), 
                        prior_intercept = NULL, kappa_mean = 1, prior_PD = FALSE, 
                        algorithm = c("sampling", "meanfield", "fullrank"),
                        adapt_delta = NULL) {
  f <- as.character(formula)
  f <- as.formula(paste(f[2], "~", "-1 + (", f[3], ")"))
  isGLM <- !missing(family)
  if (isGLM) stop("GLMs are not currently supported so omit the 'family' argument")
  mCall <- match.call(expand.dots = FALSE)
  mCall$prior <- mCall$prior_intercept <- mCall$prior_PD <- mCall$algorithm <- 
    mCall$adapt_delta <- mCall$`...` <- NULL
  mCall$pool <- FALSE
  mCall$y <- TRUE
  mCall[[1L]] <- quote(lme4::lmList)
  lms <- eval.parent(mCall)
  group_names <- levels(lms@groups)
  has_intercept <- identical(attr(lms@.Data[[1L]]$terms, "intercept"), 1L)
  if (has_intercept) {
    b <- lapply(lms@.Data, FUN = function(ols) coef(ols)[-1L])
    X <- lapply(lms@.Data, FUN = function(ols) model.matrix(ols)[,-1L, drop = FALSE])
  }
  else {
    b <- lapply(lms@.Data, FUN = coef)
    X <- lapply(lms@.Data, FUN = model.matrix)
  }
  names(b) <- names(X) <- group_names
  R <- sapply(X, simplify = FALSE, FUN = function(x) {
    qr.R(qr(sweep(x, MARGIN = 2, STATS = colMeans(x))))
  })
  SSR <- sapply(lms@.Data, FUN = function(ols) crossprod(residuals(ols))[1])
  names(SSR) <- group_names
  N <- sapply(X, FUN = NROW)
  xbar <- sapply(X, simplify = FALSE, FUN = colMeans)
  ybar <- sapply(lms@.Data, FUN = function(ols) mean(ols$y))
  s_y  <- sapply(lms@.Data, FUN = function(ols)   sd(ols$y))
  names(ybar) <- names(s_y) <- group_names
  algorithm <- match.arg(algorithm)
  stanfit <- stan_biglm.fit(b, R, SSR, N, xbar, ybar, s_y, has_intercept, ..., 
                            prior = prior, prior_intercept = prior_intercept,
                            kappa_mean = kappa_mean, prior_PD = prior_PD, 
                            algorithm = algorithm, adapt_delta = adapt_delta)
  g <- glFormula(f, data)
  modelframe <- g$fr
  Z <- t(g$reTrms$Zt)
  K <- NCOL(X[[1]]) + has_intercept
  colnames(Z) <- c(sapply(group_names, FUN = function(g) paste0(1:K, ":", g)))
  Z <- Z[,sort(colnames(Z))]
  fit <- nlist(stanfit, family = gaussian(), formula, offset = NULL, weights = NULL,
               x = Z, y = modelframe[,1], data = if (missing("data")) environment(formula) else data,
               prior.info = prior, 
               algorithm, call = match.call(), terms = attr(g$fr, "terms"),
               model = NULL,
               na.action = attr(modelframe, "na.action"),
               contrasts = attr(modelframe, "contrasts"))
  out <- stanreg(fit)
  out$groups <- lms@groups
  class(out) <- c("stanreg", "lmList")
  return(out)
}
