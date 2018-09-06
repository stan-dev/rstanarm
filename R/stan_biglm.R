# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2016, 2017 Trustees of Columbia University
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

#' Bayesian regularized linear but big models via Stan
#' 
#' \if{html}{\figure{stanlogo.png}{options: width="25px" alt="http://mc-stan.org/about/logo/"}}
#' This is the same model as with \code{\link{stan_lm}} but it utilizes the
#' output from \code{\link[biglm]{biglm}} in the \pkg{biglm} package in order to
#' proceed when the data is too large to fit in memory.
#' 
#' @export
#' @param biglm The list output by \code{\link[biglm]{biglm}} in the \pkg{biglm}
#'   package.
#' @param xbar A numeric vector of column means in the implicit design matrix 
#'   excluding the intercept for the observations included in the model.
#' @param ybar A numeric scalar indicating the mean of the outcome for the
#'   observations included in the model.
#' @param s_y A numeric scalar indicating the unbiased sample standard deviation
#'   of the outcome for the observations included in the model.
#' @template args-dots
#' @param prior Must be a call to \code{\link{R2}} with its \code{location}
#'   argument specified or \code{NULL}, which would indicate a standard uniform
#'   prior for the \eqn{R^2}.
#' @inheritParams stan_lm
#' @template args-prior_PD
#' @template args-algorithm
#' @template args-adapt_delta
#' 
#' @details The \code{stan_biglm} function is intended to be used in the same 
#'   circumstances as the \code{\link[biglm]{biglm}} function in the \pkg{biglm}
#'   package but with an informative prior on the \eqn{R^2} of the regression. 
#'   Like \code{\link[biglm]{biglm}}, the memory required to estimate the model 
#'   depends largely on the number of predictors rather than the number of 
#'   observations. However, \code{stan_biglm} and \code{stan_biglm.fit} have 
#'   additional required arguments that are not necessary in 
#'   \code{\link[biglm]{biglm}}, namely \code{xbar}, \code{ybar}, and \code{s_y}.
#'   If any observations have any missing values on any of the predictors or the 
#'   outcome, such observations do not contribute to these statistics.
#'   
#' @return The output of both \code{stan_biglm} and \code{stan_biglm.fit} is an
#'   object of \code{\link[rstan]{stanfit-class}} rather than
#'   \code{\link{stanreg-objects}}, which is more limited and less convenient
#'   but necessitated by the fact that \code{stan_biglm} does not bring the full
#'   design matrix into memory. Without the full design matrix,some of the
#'   elements of a \code{\link{stanreg-objects}} object cannot be calculated,
#'   such as residuals. Thus, the functions in the \pkg{rstanarm} package that
#'   input \code{\link{stanreg-objects}}, such as 
#'   \code{\link{posterior_predict}} cannot be used.
#'   
stan_biglm <- function(biglm, xbar, ybar, s_y, ...,
                       prior = R2(stop("'location' must be specified")), 
                       prior_intercept = NULL, prior_PD = FALSE, 
                       algorithm = c("sampling", "meanfield", "fullrank"),
                       adapt_delta = NULL) {
  if (!inherits(biglm, "biglm")    || is.null(biglm$qr) ||
      !inherits(biglm$qr, "bigqr") || is.null(biglm$terms))
    stop("'biglm' must be of S3 class biglm as defined by the biglm package.")

  b <- coef(biglm)
  R <- diag(length(b))
  R[upper.tri(R)] <- biglm$qr$rbar
  R <- sqrt(biglm$qr$D) * R
  if (identical(attr(biglm$terms, "intercept"), 1L)) {
    b <- b[-1]
    R <- R[-1,-1]
    has_intercept <- TRUE
  }
  else has_intercept <- FALSE

  return(stan_biglm.fit(b, R, SSR = biglm$qr$ss, N = biglm$n, xbar, ybar, s_y, 
                        has_intercept, ...,
                        prior = prior, prior_intercept = prior_intercept,
                        prior_PD = prior_PD, algorithm = algorithm, 
                        adapt_delta = adapt_delta))
}

