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

#' In-sample or out-of-sample predictive errors
#' 
#' This is a convenience function for computing \eqn{y - y^{rep}}{y - yrep} 
#' (in-sample, for observed \eqn{y}) or \eqn{y - \tilde{y}}{y - ytilde} 
#' (out-of-sample, for new or held-out \eqn{y}). The method for stanreg objects 
#' calls \code{\link{posterior_predict}} internally, whereas the method for
#' objects with class \code{"ppd"} accepts the matrix returned by 
#' \code{posterior_predict} as input and can be used to avoid multiple calls to 
#' \code{posterior_predict}.
#' 
#' @aliases predictive_error
#' @export
#' 
#' @param object Either a fitted model object returned by one of the 
#'   \pkg{rstanarm} modeling functions (a \link[=stanreg-objects]{stanreg 
#'   object}) or, for the \code{"ppd"} method, a matrix of draws from the 
#'   posterior predictive distribution returned by 
#'   \code{\link{posterior_predict}}.
#' @param newdata,draws,seed,offset,re.form Optional arguments passed to 
#'   \code{\link{posterior_predict}}. For binomial models, please see the
#'   \strong{Note} section below if \code{newdata} will be specified.
#' @template args-dots-ignored
#' 
#' @return A \code{draws} by \code{nrow(newdata)} matrix. If \code{newdata} is 
#'   not specified then it will be \code{draws} by \code{nobs(object)}.
#'   
#' @note The \strong{Note} section in \code{\link{posterior_predict}} about 
#'   \code{newdata} for binomial models also applies for
#'   \code{predictive_error}, with one important difference. For
#'   \code{posterior_predict} if the left-hand side of the model formula is 
#'   \code{cbind(successes, failures)} then the particular values of 
#'   \code{successes} and \code{failures} in \code{newdata} don't matter, only 
#'   that they add to the desired number of trials. \strong{This is not the case
#'   for} \code{predictive_error}. For \code{predictive_error} the particular
#'   value of \code{successes} matters because it is used as \eqn{y} when
#'   computing the error.
#' 
#' @seealso \code{\link[=posterior_predict.stanreg]{posterior_predict}} to draw
#'   from the posterior predictive distribution without computing predictive
#'   errors.
#'   
#' @examples
#' if (!exists("example_model")) example(example_model)
#' err1 <- predictive_error(example_model, draws = 50)
#' hist(err1)
#' 
#' # Using newdata with a binomial model
#' formula(example_model)
#' nd <- data.frame(
#'  size = c(10, 20), 
#'  incidence = c(5, 10), 
#'  period = factor(c(1,2)), 
#'  herd = c(1, 15)
#' )
#' err2 <- predictive_error(example_model, newdata = nd, draws = 10, seed = 1234)
#' 
#' # stanreg vs ppd methods
#' fit <- stan_glm(mpg ~ wt, data = mtcars, iter = 300)
#' preds <- posterior_predict(fit, seed = 123)
#' all.equal(
#'   predictive_error(fit, seed = 123),
#'   predictive_error(preds, y = fit$y)
#' )
#' 
predictive_error.stanreg <-
  function(object,
           newdata = NULL,
           draws = NULL,
           re.form = NULL,
           seed = NULL,
           offset = NULL,
           ...) {
    if (used.optimizing(object))
      STOP_not_optimizing("predictive_error")
    if (inherits(object, "polr"))
      stop("'predictive_error' is not currently available for stan_polr.")
    if ("y" %in% names(list(...)))
      stop("Argument 'y' should not be specified if 'object' is a stanreg object.")
    
    y <- if (is.null(newdata))
      get_y(object) else eval(formula(object)[[2L]], newdata)
    
    fam <- family(object)$family
    if (is.binomial(fam) && NCOL(y) == 2)
      y <- y[, 1]
    
    ytilde <- posterior_predict(
      object,
      newdata = newdata,
      draws = draws,
      offset = offset,
      seed = seed,
      re.form = re.form
    )
    predictive_error.ppd(ytilde, y = y)
  }

#' @rdname predictive_error.stanreg
#' @export
#' @param y For the \code{"ppd"} method only, a vector of \eqn{y} values the 
#'   same length as the number of columns in the matrix used as \code{object}. 
#'   The method for stanreg objects takes \code{y} directly from the fitted 
#'   model object.
#'   
predictive_error.ppd <- function(object, y, ...) {
  ytilde <- unclass(object)
  rstantools::predictive_error(ytilde, y = y)
}

# @rdname predictive_error.stanreg
# @export
# @param m For \code{stanmvreg} models, the submodel for which to calculate
#   the prediction error. Can be an integer, or for \code{\link{stan_mvmer}}
#   models it can be \code{"y1"}, \code{"y2"}, etc, or for \code{\link{stan_jm}}
#   models it can be \code{"Event"}, \code{"Long1"}, \code{"Long2"}, etc.
# @param t,u Only relevant for \code{\link{stan_jm}} models and when \code{m = "Event"}. 
#   The argument \code{t} specifies the time up to which individuals must have survived
#   as well as being the time up to which the longitudinal data in \code{newdata}
#   is available. The argument \code{u} specifies the time at which the 
#   prediction error should be calculated (i.e. the time horizon).
#   
predictive_error.stanmvreg <-
  function(object,
           newdataLong = NULL,
           newdataEvent = NULL,
           m = "Event",
           draws = NULL,
           re.form = NULL,
           seed = NULL,
           offset = NULL,
           t, u,
           lossfn = "square",
           ...) {
    stop("This function is not yet implemented for stanmvreg objects.")
    if ("y" %in% names(list(...)))
      stop("Argument 'y' should not be specified if 'object' is a stanmvreg object.")
    if (!is.jm(object))
      stop("This function is currently only implemented for stan_jm models.")
    if (missing(t)) 
      t <- NULL
    if (missing(u))
      u <- NULL
    M <- get_M(object)
    
    if (m == "Event") { # prediction error for event submodel
            
      if (!is.surv(object))
        stop("No event submodel was found in the fitted object.")
      if (is.null(t) || is.null(u))
        stop("'t' and 'u' must be specified when calculating the ",
             "prediction error for the event submodel.")
      if (u <= t)
        stop("'u' must be greater than 't'.")
      
      # Construct prediction data
      # ndL: dataLong to be used in predictions
      # ndE: dataEvent to be used in predictions
      if (!identical(is.null(newdataLong), is.null(newdataEvent)))
        stop("Both newdataLong and newdataEvent must be supplied together.")
      if (is.null(newdataLong)) { # user did not specify newdata
        dats <- get_model_data(object)
        ndL <- dats[1:M]
        ndE <- dats[["Event"]]
      } else { # user specified newdata
        newdatas <- validate_newdatas(object, newdataLong, newdataEvent)
        ndL <- newdatas[1:M]
        ndE <- newdatas[["Event"]]   
      }
      
      # Subset prediction data to only include
      # observations prior to time t
      fm_LHS <- formula(object, m = "Event")[[2L]]
      event_tvar <- as.character(fm_LHS[[length(fm_LHS) - 1L]])
      sel <- which(ndE[[event_tvar]] > t)
      ndE <- ndE[sel, , drop = FALSE]
      ndL <- lapply(ndL, function(x) {
        sel <- which(x[[object$time_var]] > t)
        x[sel, , drop = FALSE]
      })      
      id_var <- object$id_var
      ids <- ndE[[id_var]]
      for (i in 1:length(ndL))
        ids <- intersect(ndL[[i]][[id_var]], ids)
      if (!length(ids))
        stop("No individuals still at risk at time 't' and ",
             "with longitudinal measurements prior to 't'.")
      ndE <- ndE[ndE[[id_var]] %in% ids, , drop = FALSE]
      ndL <- lapply(ndL, function(x) {
        x[x[[id_var]] %in% ids, , drop = FALSE]
      })
      
      # Observed y: event status at time u
      event_dvar <- as.character(fm_LHS[[length(fm_LHS)]])
      y <- ndE[, c(id_var, event_tvar, event_dvar), drop = FALSE]

      # Predicted y: conditional survival probability at time u
      ytilde <- posterior_survfit(
        object, 
        newdataLong = ndL, 
        newdataEvent = ndE,
        times = u,
        last_time = t,
        last_time2 = event_tvar,
        condition = TRUE,
        extrapolate = FALSE,
        draws = draws,
        seed = seed)
      ytilde <- ytilde[, c(id_var, "survpred", "survpred_eventtime"), drop = FALSE]

      y <- merge(y, ytilde, by = id_var)

      loss <- switch(lossfn,
                     square = function(x) {x*x},
                     absolute = function(x) {abs(x)})
      
      y$dummy <- as.integer(y[[event_tvar]] > u)
      y$status <- as.integer(y[[event_dvar]])
      y$res <- 
        y$dummy * loss(1 - y$survpred) +
        y$status * (1 - y$dummy) * loss(0 - y$survpred) +
        (1 - y$status) * (1 - y$dummy) * (
          y$survpred_eventtime * loss(1- y$survpred) + 
            (1 - y$survpred_eventtime) + loss(0- y$survpred)
        )
      return(list(PE = mean(y$res), N = nrow(y)))
      
    } else { # prediction error for longitudinal submodel
      
      y <- if (is.null(newdataLong))
        get_y(object, m = m) else 
          eval(formula(object, m = m)[[2L]], newdataLong)
      
      fam <- family(object, m = m)$family
      if (is.binomial(fam) && NCOL(y) == 2)
        y <- y[, 1]      
      
      ytilde <- posterior_predict(
        object,
        m = m,
        newdata = newdataLong,
        draws = draws,
        offset = offset,
        seed = seed,
        re.form = re.form
      )
      
      return(predictive_error.ppd(ytilde, y = y))
    }
  }
