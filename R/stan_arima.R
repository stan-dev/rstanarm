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

#' Bayesian ARIMA modeling of time series
#'
#' @export
#' @param x,order,seasonal,xreg,include.mean Same as \code{\link[stats]{arima}}.
#' @templateVar pkg stats
#' @templateVar pkgfun arima
#' @templateVar fun stan_arima
#' @templateVar fitfun stan_arima.fit
#' @template return-stanreg-object
#' @template return-stanfit-object
#' @template see-also
#' @template args-dots
# #' @template args-priors
# #' @template args-prior_PD
#' @template args-algorithm
#' @template args-adapt_delta
#' 
#' @details The \code{stan_arima} function is similar in syntax to 
#'   \code{\link[stats]{arima}} but rather than performing maximum 
#'   likelihood estimation, full Bayesian estimation is performed (if 
#'   \code{algorithm} is \code{"sampling"}) via MCMC.
#'   
#' @seealso The vignette for \code{stan_arima}.
#' 
#' 
#' @examples 
#' # add example of stan_arima
#'
stan_arima <-
  function(x, 
           order = c(0L, 0L, 0L),
           seasonal = list(order = c(0L, 0L, 0L), period = NA),
           xreg = NULL, 
           include.mean = TRUE,
           ...,
           prior = normal(),
           prior_intercept = normal(),
           # prior_PD = FALSE,
           algorithm = c("sampling", "optimizing", "meanfield", "fullrank"),
           adapt_delta = NULL
           #QR = FALSE
           ) {
    
    mc <- match.call(expand.dots = FALSE)
    algorithm <- match.arg(algorithm)

    stanfit <-
      stan_arima.fit(
        y = x,
        x = xreg,
        ...,
        has_intercept = include.mean,
        prior = prior,
        prior_intercept = prior_intercept,
        # prior_PD = prior_PD,
        algorithm = algorithm,
        adapt_delta = adapt_delta
        # QR = QR
      )
    fit <-
      nlist(
        stanfit,
        algorithm,
        y = x,
        x = xreg,
        family = gaussian(),
        call = match.call(),
        stan_function = "stan_arima"
      )
    out <- stanreg(fit)
    structure(out, class = c("stanreg", "Arima"))
  }

