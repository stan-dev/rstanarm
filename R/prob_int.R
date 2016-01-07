# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2015 Trustees of Columbia University
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

#' Bayesian uncertainty intervals
#' 
#' For models fit using MCMC (\code{algorithm="sampling"}) or one of the 
#' variational approximations (\code{"meanfield"} or \code{"fullrank"}), the
#' \code{prob_int} function computes Bayesian uncertainty intervals.
#' 
#' @export
#' @templateVar stanregArg object
#' @template args-stanreg-object
#' @template args-regex-pars
#' @param type The type of interval to compute. If \code{type="central"} (the 
#'   default), the intervals are computed from posterior quantiles. The central 
#'   \eqn{100p}\% interval is defined by the \eqn{\alpha} and \eqn{1 - \alpha} 
#'   quantiles, where \eqn{\alpha = (1 - p) / 2}. If \code{type="spin"}, 
#'   shortest probability intervals are computed via \code{\link[SPIn]{SPIn}}.
#'   A shortest probability interval is a particular type of highest probability
#'   density (HPD) interval, i.e., the shortest interval with probability 
#'   coverage \eqn{p}. See Liu, Gelman, and Zheng (2015) for a detailed
#'   discussion.
#'   
#' @param pars An optional character vector of parameter names.
#' @param prob A number between 0 and 1 indicating the desired posterior
#'   probability mass \eqn{p} to include in the interval. The default is 0.5,
#'   yielding a 50\% interval.
#' @param ... Currently ignored.
#' 
#' @return A two-column matrix.
#' 
#' @references 
#' Liu, Y., Gelman, A., and Zheng, T. (2015). Simulation-efficient shortest
#' probability intervals. \emph{Statistics and Computing}. 25(4), 809--819.
#'   
#' @seealso \code{\link{confint.stanreg}}, which, for models fit using 
#'   optimization, can be used to compute traditional confidence intervals.
#' 
#' @examples 
#' prob_int(example_model)
#' prob_int(example_model, regex_pars = "herd")
#' prob_int(example_model, type = "hpd", pars = "period2", prob = 0.5)
#' 
prob_int <- function(object, type = c("central", "spin"), prob = 0.5, 
                    pars = NULL, regex_pars = NULL, ...) {
  if (!is.stanreg(object))
    stop(deparse(substitute(object)), " is not a stanreg object.")
  if (used.optimizing(object))
    STOP_not_optimizing("credint")
  if (!identical(length(prob), 1L) || prob <= 0 || prob >= 1)
    stop("'prob' should be a single number greater than 0 and less than 1.", 
         call. = FALSE)
  
  type <- match.arg(type)
  mat <- as.matrix.stanreg(object)
  pars <- .collect_pars(object, pars, regex_pars)
  if (!is.null(pars)) {
    mat <- mat[, colnames(mat) %in% pars, drop = FALSE]
  } else {
    pars <- grep("mean_PPD|log-posterior", colnames(mat), invert = TRUE, 
                 value = TRUE)
    mat <- mat[, pars, drop = FALSE]
  }
  if (type == "central") {
    alpha <- (1 - prob) / 2
    probs <- c(alpha, 1 - alpha)
    labs <- paste0(100 * probs, "%")
    ci <- t(apply(mat, 2L, quantile, probs = probs))
  } 
  else {
    stopifnot(type == "spin") 
    if (!requireNamespace("SPIn", quietly = TRUE)) {
      stop("Please install the 'SPIn' package.", call. = FALSE)
    }
    spin <- function(x, prob) SPIn::SPIn(x, conf = prob)[["spin"]]
    ci <- t(apply(mat, 2L, spin, prob = prob))
    labs <- paste0(c("lower", "upper"), 100 * prob)
  }
  return(structure(ci, dimnames = list(colnames(mat), labs)))
}
