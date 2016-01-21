# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2015 Trustees of Columbia University
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
#' @name neg_binomial_2
#' @export
#' @param link The same as for \code{\link{poisson}}, typically a character
#'   vector of length one among \code{"log"}, \code{"identity"}, and
#'   \code{"sqrt"}.
#' @return An object of class \code{\link[stats]{family}} very similar to
#'   that of \code{\link[stats]{poisson}} but with a different family name.
#' @examples
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


#' Family function for Student t GLMs
#' 
#' Estimates of regression coefficients are less sensitive to outliers if the
#' Student t distribution is used in place of the normal distribution in
#' settings where some errors may be large. These models are sometimes referred
#' to as \emph{robust regression models}. To estimate such a model, a call to
#' \code{t_family} can be passed to the \code{family} argument of 
#' \code{\link{stan_glm}} or \code{\link{stan_glmer}} instead of a call to 
#' \code{\link{gaussian}}.
#' 
#' @name t_family
#' @export
#' @templateVar armRef (Ch. 6)
#' @template reference-gelman-hill
#' @param link The same as for \code{\link[stats]{gaussian}}.
#' @return An object of class \code{\link[stats]{family}} very similar to
#'   that of \code{\link[stats]{gaussian}} but with a different family name.
#' @examples 
#' \dontrun{
#' SEED <- 1234
#' set.seed(SEED)
#' x <- rnorm(1000)
#' alpha <- 2; beta <- 0.5; df <- 4
#' y <- alpha + beta * x + rt(1000, df)
#' (fit <- stan_glm(y ~ x, family = t_family(), seed = SEED, cores = 4))
#' }
#' 
t_family <- function(link = "identity") {
  out <- gaussian(link)
  out$family <- "t_family"
  out$variance <- NULL
  out$dev.resids <- function(y, mu, wt)
    stop("'dev.resids' function should not be called")
  out$aic <- function(y, n, mu, wt, dev)
    stop("'aic' function should not have been called")
  out$simulate <- function(object, nsim)
    stop("'simulate' function should not have been called")
  return(out)
}
