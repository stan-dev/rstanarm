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

#' Bayesian regularized linear but big models via Stan
#' 
#' This is the same model as with \code{\link{stan_lm}} but it utilizes
#' the output from \code{\link[biglm]{biglm}} in the \pkg{biglm} package
#' in order to proceed when the data is too large to fit in memory.
#' 
#' @export
#' @param biglm The list output by \code{\link[biglm]{biglm}} in the \pkg{biglm}
#'   package. The original call to \code{\link[biglm]{biglm}} must not have an
#'   intercept and must utilize \emph{centered} but not \emph{standardized}
#'   predictors. See the Details section or the Example.
#' @param xbar A numeric vector of means in the implicit design matrix for
#'   the observations included in the model
#' @param ybar A numeric scalar indicating the same mean of the outcome for
#'   the observations included in the model
#' @param s_y A numeric scalar indicating the unbiased sample standard deviation
#'   of the outcome for the observations included in the model
#' @param has_intercept A logical scalar indicating whether to add an intercept
#'   to the model when estimating it
#' @template args-dots
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
#' @template args-prior_PD
#' @template args-algorithm
#' @template args-adapt_delta
#' 
#' @details The \code{stan_biglm} function is intended to be used in the same
#'   circumstances as the \code{\link[biglm]{biglm}} function in the \pkg{biglm}
#'   package but with an informative prior on the \eqn{R^2} of the regression.
#'   Like \code{\link[biglm]{biglm}}, the memory required to estimate the model
#'   depends largely on the number of predictors rather than the number of 
#'   observations. However, the original call to \code{\link[biglm]{biglm}} must
#'   be a little unconventional. The original \code{\link[stats]{formula}} must
#'   not include an intercept and all the columns of the implicit design matrix
#'   must be expressed as deviations from the sample mean. If the design matrix
#'   is on the hard disk, the column sums must be accumulated, divided by the
#'   sample size to produce the column means, and then the column means must be
#'   swept from the design matrix on disk. If any observations have any missing
#'   values on any of the predictors or the outcome, such observations do not
#'   contribute to the column means, which must be passed as the \code{xbar}
#'   argument. If the outcome is also expressed as the deviation from its 
#'   sample mean, then the coefficients produced by \code{\link[biglm]{biglm}}
#'   are the same as if the raw data were used and an intercept were included.
#'   The sample mean and sample standard deviation of the outcome must also
#'   be passed.
#'   
#' @return The output of both \code{stan_biglm} and \code{stan_biglm.fit} is an object of
#'   \code{\link[rstan]{stanfit-class}} rather than \code{\link{stanreg-objects}}, 
#'   which is more limited and less convenient but necessitated by the fact that 
#'   \code{stan_biglm} does not bring the full design matrix into memory. Without the 
#'   full design matrix,some of the elements of a \code{\link{stanreg-objects}} object 
#'   cannot be calculated, such as residuals. Thus, the functions in the \pkg{rstanarm}
#'   package that input \code{\link{stanreg-objects}}, such as
#'   \code{\link{posterior_predict}} cannot be used.
stan_biglm <- function(biglm, xbar, ybar, s_y, has_intercept = TRUE, ...,
                       prior = R2(stop("'location' must be specified")), 
                       prior_intercept = NULL, prior_PD = FALSE, 
                       algorithm = c("sampling", "meanfield", "fullrank"),
                       adapt_delta = NULL) {
  if (!inherits(biglm, "biglm")    || is.null(biglm$qr) ||
      !inherits(biglm$qr, "bigqr") || is.null(biglm$terms))
    stop("'biglm' must be of S3 class biglm as defined by the biglm package")
  if (identical(attr(biglm$terms, "intercept"), 1L))
    stop("The original biglm must not include an intercept", 
         "Rather, the predictors should be centered so that the intercept is the sample mean of the outcome")
  b <- coef(biglm)
  R <- diag(length(b))
  pos <- 1L
  for (i in 1:(ncol(R) - 1L)) for(j in (i+1):ncol(R)) {
    R[i,j] <- biglm$qr$rbar[pos]
    pos <- pos + 1L
  }
  R <- sqrt(biglm$qr$D) * R
  return(stan_biglm.fit(b, R, SSR = biglm$qr$ss, N = biglm$n, xbar, ybar, s_y, has_intercept, 
                        ...,
                        prior = prior, prior_intercept = prior_intercept,
                        prior_PD = prior_PD, algorithm = algorithm, 
                        adapt_delta = adapt_delta))
}

