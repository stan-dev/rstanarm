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

#' Predictive intervals
#' 
#' For models fit using MCMC (\code{algorithm="sampling"}) or one of the 
#' variational approximations (\code{"meanfield"} or \code{"fullrank"}), the 
#' \code{predictive_interval} function computes Bayesian predictive intervals. 
#' 
#' @export
#' @templateVar stanregArg object
#' @template args-stanreg-object
#' @inheritParams posterior_interval
#' @param newdata,draws,fun,offset,re.form,seed Passed to
#'   \code{\link[=posterior_predict.stanreg]{posterior_predict}}.
#' @param ... Currently ignored by the method for stanreg objects. The S3
#'   generic uses \code{...} to pass arguments to any defined methods.
#' 
#' @return A matrix with two columns and as many rows as are in \code{newdata}. 
#'   If \code{newdata} is not provided then the matrix will as many rows as in
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
#' predictive_interval(fit, newdata = data.frame(wt = range(mtcars$wt)))
#' 
predictive_interval <- function(object, ...) {
  UseMethod("predictive_interval")
}

#' @rdname predictive_interval
#' @export
predictive_interval.stanreg <- function(object,
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
    fun = fun,
    ...
  )
  central_intervals(ytilde, prob)
}
