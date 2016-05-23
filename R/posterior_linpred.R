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

#' Posterior distribution of the linear predictor
#' 
#' Extract the posterior draws of the linear predictor, possibly transformed by 
#' the inverse-link function. This function is occasionally convenient, but it 
#' should be used sparingly. Inference and model checking should generally be 
#' carried out using the posterior predictive distribution (see
#' \code{\link{posterior_predict}}).
#' 
#' @export
#' @templateVar stanregArg object
#' @template args-stanreg-object
#' @param transform Should the linear predictor be transformed using the
#'   inverse-link function? The default is \code{FALSE}, in which case the 
#'   untransformed linear predictor is returned.
#' @param newdata,re.form Same as for \code{\link{posterior_predict}}.
#' @param ... Currently unused.
#' 
#' @return A \code{draws} by \code{nrow(newdata)} matrix of simulations from the
#'   posterior distribution of the (possibly transformed) linear predictor.
#' 
#' @seealso \code{\link{posterior_predict}} to draw from the posterior 
#'   predictive distribution of the outcome, which is almost always preferable.
#' 
#' @examples
#' if (!exists("example_model")) example(example_model)
#' 
#' # linear predictor on log-odds scale
#' linpred <- posterior_linpred(example_model)
#' # probabilities
#' probs <- posterior_linpred(example_model, transform = TRUE)
#' 
#' # not conditioning on any group-level parameters
#' probs2 <- posterior_linpred(example_model, transform = TRUE, re.form = NA)
#' 
posterior_linpred <- function(object, transform = FALSE, newdata = NULL, 
                              re.form = NULL, ...) {
  validate_stanreg_object(object)
  if (used.optimizing(object))
    STOP_not_optimizing("posterior_linpred")
  if (!is.null(newdata)) {
    if ("gam" %in% names(object))
      stop("'posterior_linpred' with 'newdata' not yet supported ", 
           "for models estimated via 'stan_gamm4'.")
    newdata <- as.data.frame(newdata)
    if (any(is.na(newdata))) 
      stop("Currently NAs are not allowed in 'newdata'.")
  }
  dat <- pp_data(object, newdata = newdata, re.form = re.form, ...)
  eta <- pp_eta(object, data = dat, draws = NULL)[["eta"]]
  if (!transform)
    return(eta)
  linkinv(object)(eta)
}
