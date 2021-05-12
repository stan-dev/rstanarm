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

#' Predict method for stanreg objects
#' 
#' This method is primarily intended to be used only for models fit using 
#' optimization. For models fit using MCMC or one of the variational
#' approximations, see \code{\link{posterior_predict}}.
#' 
#' @export
#' @templateVar stanregArg object
#' @template args-stanreg-object
#' @param ... Ignored.
#' @param newdata Optionally, a data frame in which to look for variables with 
#'   which to predict. If omitted, the model matrix is used.
#' @param type The type of prediction. The default \code{'link'} is on the scale
#'   of the linear predictors; the alternative \code{'response'} is on the scale
#'   of the response variable.
#' @param se.fit A logical scalar indicating if standard errors should be 
#'   returned. The default is \code{FALSE}.
#'   
#' @return A vector if \code{se.fit} is \code{FALSE} and a list if \code{se.fit}
#'   is \code{TRUE}.
#'
#' @seealso \code{\link{posterior_predict}}
#' 
predict.stanreg <- function(object,
                            ...,
                            newdata = NULL,
                            type = c("link", "response"),
                            se.fit = FALSE) {
  if (is.mer(object)) {
    stop(
      "'predict' is not available for models fit with ",
      object$stan_function,
      ". Please use the 'posterior_predict' function instead.",
      call. = FALSE
    )
  }
  
  type <- match.arg(type)
  if (!se.fit && is.null(newdata)) {
    preds <- if (type == "link") 
      object$linear.predictors else object$fitted.values
    return(preds)
  }

  if (isTRUE(object$stan_function == "stan_betareg") && 
      !is.null(newdata)) {
    # avoid false positive warnings about missing z variables in newdata
    zvars <- all.vars(object$terms$precision)
    for (var in zvars) {
      if (!var %in% colnames(newdata)) newdata[[var]] <- NA
    }
  }
  
  dat <- pp_data(object, newdata)
  stanmat <- as.matrix.stanreg(object)
  beta <- stanmat[, seq_len(ncol(dat$x))]
  eta <- linear_predictor(beta, dat$x, dat$offset)
  if (type == "response") {
    inverse_link <- linkinv(object)
    eta <- inverse_link(eta)
    if (is(object, "polr") && ("alpha" %in% colnames(stanmat)))
      eta <- apply(eta, 1L, FUN = `^`, e2 = stanmat[, "alpha"])
  }
  fit <- colMeans(eta)
  if (!se.fit)
    return(fit)
  
  se.fit <- apply(eta, 2L, sd)
  nlist(fit, se.fit)
}
