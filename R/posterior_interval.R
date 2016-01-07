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
#' \code{posterior_interval} function computes Bayesian uncertainty intervals. These
#' intervals are often referred to as \emph{credible} intervals, but we use the
#' term \emph{uncertainty} intervals to highlight the fact that wider intervals
#' correspond to greater uncertainty.
#' 
#' @export
#' @templateVar stanregArg object
#' @template args-stanreg-object
#' @template args-regex-pars
#' @param prob A number between 0 and 1 indicating the desired posterior 
#'   probability mass \eqn{p} to include in the intervals. The default is to 
#'   report 50\% intervals (\code{prob=0.5}) rather than the traditionally used
#'   95\%, although the latter can be computed by specifying \code{prob=0.95}. 
#'   We use this default for several reasons:
#'   \itemize{
#'   \item Computational stability: 50\% intervals are more stable than 95\%
#'   intervals (for which each end relies on only 2.5\% of the posterior draws).
#'   \item More intuitive evaluation: (roughly) half of the 50\% intervals
#'   should contain the true value.
#'   \item In aplications, a good first step is to get a sense of where the 
#'   parameters and predicted values will be, not to attempt an unrealistic 
#'   near-certainty.
#'   }
#' @param type The type of interval to compute. Currently the only option is 
#'   \code{"central"}, although other possibilities may be available in future
#'   releases. A central \eqn{100p}\% interval is defined by the \eqn{\alpha}
#'   and \eqn{1 - \alpha} quantiles, where \eqn{\alpha = (1 - p)/2}.
#' @param pars An optional character vector of parameter names.
#' @param ... Currently ignored.
#' 
#' @return A two-column matrix.
#'   
#' @seealso \code{\link{confint.stanreg}}, which, for models fit using 
#'   optimization, can be used to compute traditional confidence intervals.
#'   
#' @references 
#' Morey, R. D., Hoekstra, R., Rouder, J., Lee, M. D., and Wagenmakers, E.
#' (2015). The fallacy of placing confidence in confidence intervals.
#' \emph{Psychonomic Bulletin & Review}. 1 -- 21.
#' 
#' @examples 
#' posterior_interval(example_model)
#' posterior_interval(example_model, regex_pars = "herd")
#' posterior_interval(example_model, pars = "period2", prob = 0.9)
#' 
posterior_interval <- function(object, prob = 0.5, type = "central",
                    pars = NULL, regex_pars = NULL, ...) {
  if (!is.stanreg(object))
    stop(deparse(substitute(object)), " is not a stanreg object.")
  if (used.optimizing(object))
    STOP_not_optimizing("posterior_interval")
  if (!identical(length(prob), 1L) || prob <= 0 || prob >= 1)
    stop("'prob' should be a single number greater than 0 and less than 1.", 
         call. = FALSE)
  if (!identical(type, "central"))
    stop("Currently the only option for 'type' is 'central'.", call. = FALSE)
  
  mat <- as.matrix.stanreg(object)
  pars <- .collect_pars(object, pars, regex_pars)
  if (!is.null(pars)) {
    mat <- mat[, colnames(mat) %in% pars, drop = FALSE]
  } else {
    pars <- grep("mean_PPD|log-posterior", colnames(mat), invert = TRUE, 
                 value = TRUE)
    mat <- mat[, pars, drop = FALSE]
  }
  alpha <- (1 - prob) / 2
  probs <- c(alpha, 1 - alpha)
  labs <- paste0(100 * probs, "%")
  ci <- t(apply(mat, 2L, quantile, probs = probs))
  return(structure(ci, dimnames = list(colnames(mat), labs)))
}
