# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2015, 2016 Trustees of Columbia University
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

#' Family function for negative binomial GLMs
#' 
#' Specifies the information required to fit a Negative Binomial GLM in a 
#' similar way to \code{\link[MASS]{negative.binomial}}. However, here the 
#' overdispersion parameter \code{theta} is not specified by the user and always
#' estimated. A call to this function can be passed to the \code{family}
#' argument of \code{\link{stan_glm}} or \code{\link{stan_glmer}} to estimate a
#' Negative Binomial model. Alternatively, the \code{\link{stan_glm.nb}} and 
#' \code{\link{stan_glmer.nb}} wrapper functions may be used, which call 
#' \code{neg_binomial_2} internally.
#' 
#' @export
#' @param link The same as for \code{\link{poisson}}, typically a character
#'   vector of length one among \code{"log"}, \code{"identity"}, and
#'   \code{"sqrt"}.
#' @return An object of class \code{\link[stats]{family}} very similar to
#'   that of \code{\link[stats]{poisson}} but with a different family name.
#' @examples
#' if (!grepl("^sparc",  R.version$platform))
#' stan_glm(Days ~ Sex/(Age + Eth*Lrn), data = MASS::quine, seed = 123,
#'          family = neg_binomial_2, QR = TRUE, algorithm = "fullrank") 
#'                 
#' # or, equivalently, call stan_glm.nb() without specifying the family
#'
neg_binomial_2 <- function(link = "log") {
  out <- poisson(link)
  out$family <- "neg_binomial_2"
  out$variance <- function(mu, theta) mu + mu^2 / theta
  out$dev.resids <- function(y, mu, wt)
    stop("'dev.resids' function should not be called")
  out$aic <- function(y, n, mu, wt, dev)
    stop("'aic' function should not have been called")
  out$simulate <- function(object, nsim)
    stop("'simulate' function should not have been called")
  return(out)
}
