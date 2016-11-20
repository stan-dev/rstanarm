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
#' For models fit using MCMC only, the \code{log_lik} method returns the
#' \eqn{S} by \eqn{N} pointwise log-likelihood matrix, where \eqn{S} is the size
#' of the posterior sample and \eqn{N} is the number of data points.
#'
#' @aliases log_lik
#' @export
#'
#' @templateVar stanregArg object
#' @template args-stanreg-object
#' @template args-dots-ignored
#' @param newdata An optional data frame of new data (e.g. holdout data) to use
#'   when evaluating the log-likelihood. See the description of \code{newdata}
#'   for \code{\link{posterior_predict}}.
#' @param offset A vector of offsets. Only required if \code{newdata} is
#'   specified and an \code{offset} was specified when fitting the model.
#'
#' @return An \eqn{S} by \eqn{N} matrix, where \eqn{S} is the size of the
#'   posterior sample and \eqn{N} is the number of data points.
#'   
#'   
#' @examples 
#' \donttest{
#'  roaches$roach100 <- roaches$roach1 / 100
#'  fit <- stan_glm(
#'     y ~ roach100 + treatment + senior,
#'     offset = log(exposure2),
#'     data = roaches,
#'     family = poisson(link = "log"),
#'     prior = normal(0, 2.5),
#'     prior_intercept = normal(0, 10),
#'     iter = 500 # to speed up example
#'  )
#'  ll <- log_lik(fit)
#'  dim(ll)
#'  all.equal(ncol(ll), nobs(fit))
#'
#'  # using newdata argument
#'  nd <- roaches[1:2, ]
#'  nd$treatment[1:2] <- c(0, 1)
#'  ll2 <- log_lik(fit, newdata = nd, offset = c(0, 0))
#'  head(ll2)
#'  dim(ll2)
#'  all.equal(ncol(ll2), nrow(nd))
#' }
#'
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
