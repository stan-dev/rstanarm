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
#' the inverse-link function. This function is occasionally useful, but it 
#' should be used sparingly. Inference and model checking should generally be 
#' carried out using the posterior predictive distribution (i.e., using 
#' \code{\link{posterior_predict}}).
#'
#' @export
#' @templateVar stanregArg object
#' @template args-stanreg-object
#' @template args-dots-ignored
#' @param transform Should the linear predictor be transformed using the
#'   inverse-link function? The default is \code{FALSE}, in which case the
#'   untransformed linear predictor is returned.
#' @param newdata,re.form,offset Same as for \code{\link{posterior_predict}}.
#' @param XZ If \code{TRUE} then instead of computing the linear predictor the 
#'   design matrix \code{X} (or \code{cbind(X,Z)} for models with group-specific
#'   terms) constructed from \code{newdata} is returned. The default is 
#'   \code{FALSE}.
#'   
#' @return The default is to return a \code{draws} by \code{nrow(newdata)} 
#'   matrix of simulations from the posterior distribution of the (possibly 
#'   transformed) linear predictor. The exception is if the argument \code{XZ} 
#'   is set to \code{TRUE} (see the \code{XZ} argument description above).
#'   
#' @seealso \code{\link{posterior_predict}} to draw from the posterior 
#'   predictive distribution of the outcome, which is typically preferable.
#'
#' @examples
#' if (!exists("example_model")) example(example_model)
#' print(family(example_model))
#' 
#' # linear predictor on log-odds scale
#' linpred <- posterior_linpred(example_model)
#' # probabilities
#' probs <- posterior_linpred(example_model, transform = TRUE)
#'
#' # not conditioning on any group-level parameters
#' probs2 <- posterior_linpred(example_model, transform = TRUE, re.form = NA)
#'
posterior_linpred <- function(object, ...) {
  UseMethod("posterior_linpred")
}

#' @rdname posterior_linpred
#' @export 
posterior_linpred.stanreg <-
  function(object,
           transform = FALSE,
           newdata = NULL,
           re.form = NULL,
           offset = NULL,
           XZ = FALSE,
           ...) {
    if (used.optimizing(object))
      STOP_not_optimizing("posterior_linpred")
    
    newdata <- validate_newdata(newdata)
    dat <- pp_data(object,
                   newdata = newdata,
                   re.form = re.form,
                   offset = offset)
    if (XZ) {
      XZ <- dat[["x"]]
      if (is.mer(object))
        XZ <- cbind(XZ, t(dat[["Zt"]]))
      return(XZ)
    }
    eta <- pp_eta(object, data = dat, draws = NULL)[["eta"]]
    if (!transform)
      return(eta)
    
    linkinv(object)(eta)
  }
