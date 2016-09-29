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

#' Pointwise log-likelihood matrix
#' 
#' For models fit using MCMC only, the \code{log_lik} function returns the 
#' \eqn{S} by \eqn{N} pointwise log-likelihood matrix, where \eqn{S} is the size
#' of the posterior sample and \eqn{N} is the number of data points. We define a
#' new function \code{log_lik} rather than a \code{\link[stats]{logLik}} method
#' because (in addition to the conceptual difference) the documentation for
#' \code{logLik} states that the return value will be a single number, whereas
#' \code{log_lik} returns a matrix.
#' 
#' @export
#' @keywords internal
#' @param object Fitted model object.
#' @param newdata An optional data frame of new data (e.g. holdout data) to use
#'   when evaluating the log-likelihood. See the description of \code{newdata}
#'   for \code{\link{posterior_predict}}.
#' @param offset A vector of offsets. Only required if \code{newdata} is
#'   specified and an \code{offset} was specified when fitting the model.
#' @param ... Currently ignored by the method for stanreg objects. The S3
#'   generic uses \code{...} to pass arguments to any defined methods.
#'   
#' @return An \eqn{S} by \eqn{N} matrix, where \eqn{S} is the size of the
#'   posterior sample and \eqn{N} is the number of data points.
#' 
log_lik <- function(object, ...) {
  UseMethod("log_lik")
}

#' @rdname log_lik
#' @export
log_lik.stanreg <- function(object, newdata = NULL, offset = NULL, ...) {
  if (!used.sampling(object)) 
    STOP_sampling_only("Pointwise log-likelihood matrix")
  if (!is.null(newdata)) {
    if ("gam" %in% names(object))
      stop("'log_lik' with 'newdata' not yet supported ", 
           "for models estimated via 'stan_gamm4'.")
    newdata <- as.data.frame(newdata)
  }
  fun <- ll_fun(object)
  args <- ll_args(object, newdata = newdata, offset = offset)
  sapply(seq_len(args$N), function(i) {
    as.vector(fun(i = i, data = args$data[i, , drop = FALSE], 
                  draws = args$draws))
  })
}