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

#' Predictive intervals
#' 
#' For models fit using MCMC (\code{algorithm="sampling"}) or one of the 
#' variational approximations (\code{"meanfield"} or \code{"fullrank"}), the 
#' \code{predictive_interval} function computes Bayesian predictive intervals. 
#' The method for stanreg objects calls \code{\link{posterior_predict}}
#' internally, whereas the method for objects of class \code{"ppd"} accepts the
#' matrix returned by \code{posterior_predict} as input and can be used to avoid
#' multiple calls to \code{posterior_predict}.
#' 
#' @export
#' @aliases predictive_interval
#' 
#' @param object Either a fitted model object returned by one of the 
#'   \pkg{rstanarm} modeling functions (a \link[=stanreg-objects]{stanreg 
#'   object}) or, for the \code{"ppd"} method, a matrix of draws from the 
#'   posterior predictive distribution returned by 
#'   \code{\link{posterior_predict}}.
#' @template args-dots-ignored
#' @inheritParams posterior_interval.stanreg
#' @param newdata,draws,fun,offset,re.form,seed Passed to 
#'   \code{\link[=posterior_predict]{posterior_predict}}.
#' 
#' @return A matrix with two columns and as many rows as are in \code{newdata}. 
#'   If \code{newdata} is not provided then the matrix will have as many rows as
#'   the data used to fit the model. For a given value of \code{prob}, \eqn{p},
#'   the columns correspond to the lower and upper \eqn{100p}\% central interval
#'   limits and have the names \eqn{100\alpha/2}\% and \eqn{100(1 -
#'   \alpha/2)}\%, where \eqn{\alpha = 1-p}. For example, if \code{prob=0.9} is
#'   specified (a \eqn{90}\% interval), then the column names will be
#'   \code{"5\%"} and \code{"95\%"}, respectively.
#'   
#' @seealso \code{\link{predictive_error}}, \code{\link{posterior_predict}}, 
#'   \code{\link{posterior_interval}}
#' 
#' @examples 
#' fit <- stan_glm(mpg ~ wt, data = mtcars, iter = 300)
#' predictive_interval(fit)
#' predictive_interval(fit, newdata = data.frame(wt = range(mtcars$wt)), 
#'                     prob = 0.5)
#' 
#' # stanreg vs ppd methods
#' preds <- posterior_predict(fit, seed = 123)
#' all.equal(
#'   predictive_interval(fit, seed = 123),
#'   predictive_interval(preds)
#' )
#' 
predictive_interval.stanreg <-
  function(object,
           prob = 0.9,
           newdata = NULL,
           draws = NULL,
           re.form = NULL,
           fun = NULL,
           seed = NULL,
           offset = NULL,
           ...) {
    if (used.optimizing(object))
      STOP_not_optimizing("posterior_interval")
    if (inherits(object, "polr"))
      stop("'predictive_interval' is not currently available for stan_polr.")
    
    ytilde <- posterior_predict(
      object,
      newdata = newdata,
      draws = draws,
      seed = seed,
      re.form = re.form,
      offset = offset,
      fun = fun
    )
    predictive_interval.ppd(ytilde, prob = prob)
  }

#' @rdname predictive_interval.stanreg
#' @export
predictive_interval.ppd <- function(object, prob = 0.9, ...) {
  ytilde <- unclass(object)
  rstantools::predictive_interval(ytilde, prob = prob)
}
