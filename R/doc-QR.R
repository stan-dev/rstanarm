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

#' The \code{QR} argument 
#' 
#' Details about the \code{QR} argument to \pkg{rstanarm}'s modeling
#' functions.
#' 
#' @name QR-argument
#' @template reference-stan-manual
#'   
#' @details The \code{QR} argument is a logical scalar defaulting to
#'   \code{FALSE}, but if \code{TRUE} applies a scaled \code{\link{qr}}
#'   decomposition to the design matrix, \eqn{X = Q^\ast R^\ast}{X = Q* R*}. 
#'   If \code{autoscale = TRUE} (the default) 
#'   in the call to the function passed to the \code{prior} argument, then
#'   \eqn{Q^\ast = Q \sqrt{n-1}}{Q* = Q (n-1)^0.5} and 
#'   \eqn{R^\ast = \frac{1}{\sqrt{n-1}} R}{R* = (n-1)^(-0.5) R}. When 
#'   \code{autoscale = FALSE}, \eqn{R} is scaled such that the lower-right
#'   element of \eqn{R^\ast}{R*} is \eqn{1}.
#'   
#'   The coefficients relative to \eqn{Q^\ast}{Q*} are obtained and then 
#'   premultiplied by the inverse of \eqn{R^{\ast}}{R*} to obtain coefficients 
#'   relative to the original predictors, \eqn{X}. Thus, when 
#'   \code{autoscale = FALSE}, the coefficient on the last column of \eqn{X} 
#'   is the same as the coefficient on the last column of \eqn{Q^\ast}{Q*}.
#'   
#'   These transformations do not change the likelihood of the data but are 
#'   recommended for computational reasons when there are multiple predictors. 
#'   Importantly, while the columns of \eqn{X} are almost generally correlated, 
#'   the columns of \eqn{Q^\ast}{Q*} are uncorrelated by design, which often makes 
#'   sampling from the posterior easier. However, because when \code{QR} is 
#'   \code{TRUE} the \code{prior} argument applies to the coefficients relative to 
#'   \eqn{Q^\ast}{Q*} (and those are not very interpretable), setting \code{QR=TRUE} 
#'   is only recommended if you do not have an informative prior for the regression
#'   coefficients or if the only informative prior is on the last regression
#'   coefficient (in which case you should set \code{autoscale = FALSE} when
#'   specifying such priors). 
#'   
#'   For more details see the Stan case study 
#'   \emph{The QR Decomposition For Regression Models} at 
#'   \url{http://mc-stan.org/users/documentation/case-studies/qr_regression.html}.
#'   
NULL
