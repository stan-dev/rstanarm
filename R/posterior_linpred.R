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

#' Posterior distribution of the linear predictor
#'
#' Extract the posterior draws of the linear predictor, possibly transformed by 
#' the inverse-link function. This function is occasionally useful, but it 
#' should be used sparingly. Inference and model checking should generally be 
#' carried out using the posterior predictive distribution (i.e., using 
#' \code{\link{posterior_predict}}).
#'
#' @aliases posterior_linpred
#' @export
#' 
#' @templateVar stanregArg object
#' @template args-stanreg-object
#' @param transform Should the linear predictor be transformed using the
#'   inverse-link function? The default is \code{FALSE}, in which case the
#'   untransformed linear predictor is returned.
#' @param newdata,draws,re.form,offset Same as for \code{\link{posterior_predict}}.
#' @param XZ If \code{TRUE} then instead of computing the linear predictor the 
#'   design matrix \code{X} (or \code{cbind(X,Z)} for models with group-specific
#'   terms) constructed from \code{newdata} is returned. The default is 
#'   \code{FALSE}.
#' @param ... Currently ignored.   
#'   
#' @return The default is to return a \code{draws} by \code{nrow(newdata)}
#'   matrix of simulations from the posterior distribution of the (possibly
#'   transformed) linear predictor. The exception is if the argument \code{XZ}
#'   is set to \code{TRUE} (see the \code{XZ} argument description above).
#'   
#' @note For models estimated with \code{\link{stan_clogit}}, the number of 
#'   successes per stratum is ostensibly fixed by the research design. Thus,
#'   when calling \code{posterior_linpred} with new data and \code{transform =
#'   TRUE}, the \code{data.frame} passed to the \code{newdata} argument must
#'   contain an outcome variable and a stratifying factor, both with the same
#'   name as in the original \code{data.frame}. Then, the probabilities will
#'   condition on this outcome in the new data.
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
#' colMeans(linpred)
#' 
#' # probabilities
#' probs <- posterior_linpred(example_model, transform = TRUE)
#' colMeans(probs)
#' 
#' # not conditioning on any group-level parameters
#' probs2 <- posterior_linpred(example_model, transform = TRUE, re.form = NA)
#' apply(probs2, 2, median)
#' 
posterior_linpred.stanreg <-
  function(object,
           transform = FALSE,
           newdata = NULL,
           draws = NULL,
           re.form = NULL,
           offset = NULL,
           XZ = FALSE,
           ...) {

    if (is.stanmvreg(object)) {
      STOP_if_stanmvreg("'posterior_linpred'")
    }
    
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
    
    eta <- pp_eta(object, data = dat, draws = draws)[["eta"]]
    if (is.null(newdata)) {
      colnames(eta) <- rownames(model.frame(object))
    } else {
      colnames(eta) <- rownames(newdata)
    }
    
    if (!transform || is.nlmer(object)) {
      return(eta)
    }
    
    if (is_clogit(object)) {
      return(clogit_linpred_transform(object, newdata = newdata, eta = eta))
    }
    
    g <- linkinv(object)
    return(g(eta))
  }



# internal ----------------------------------------------------------------
clogit_linpred_transform <- function(object, newdata = NULL, eta = NULL) {
  g <- linkinv(object)
  if (!is.null(newdata)) {
    y <- eval(formula(object)[[2L]], newdata)
    strata <- as.factor(eval(object$call$strata, newdata))
    formals(g)$g <- strata
    formals(g)$successes <- aggregate(y, by = list(strata), FUN = sum)$x
  }
  return(t(apply(eta, 1, FUN = g)))
}
