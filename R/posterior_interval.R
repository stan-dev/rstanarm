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

#' Posterior uncertainty intervals
#' 
#' For models fit using MCMC (\code{algorithm="sampling"}) or one of the 
#' variational approximations (\code{"meanfield"} or \code{"fullrank"}), the 
#' \code{posterior_interval} function computes Bayesian posterior uncertainty 
#' intervals. These intervals are often referred to as \emph{credible} 
#' intervals, but we use the term \emph{uncertainty} intervals to highlight the 
#' fact that wider intervals correspond to greater uncertainty.
#' 
#' @export
#' @templateVar stanregArg object
#' @template args-stanreg-object
#' @template args-regex-pars
#' @param prob A number between 0 and 1 indicating the desired posterior 
#'   probability mass \eqn{p} to include in the intervals. The default is to 
#'   report 90\% intervals (\code{prob=0.9}) rather than the traditionally used 
#'   95\% (see Details).
#' @param type The type of interval to compute. Currently the only option is 
#'   \code{"central"} (see Details). A central \eqn{100p}\%
#'   interval is defined by the \eqn{\alpha} and \eqn{1 - \alpha} quantiles,
#'   where \eqn{\alpha = (1 - p)/2}.
#' @param pars An optional character vector of parameter names.
#' @param ... Currently ignored.
#' 
#' @return A two-column matrix.
#' 
#' @details
#' \subsection{Interpretation}{
#' Unlike for a frenquentist confidence interval, it is valid to say that, 
#' conditional on the data and model, we believe that with probability \eqn{p} 
#' the value of a parameter is in its \eqn{100p}\% posterior interval. This 
#' intuitive interpretation of Bayesian intervals is often erroneously applied 
#' to frequentist confidence intervals. See Morey et al. (2015) for more details
#' on this issue and the advantages of using Bayesian posterior uncertainty
#' intervals (also known as credible intervals).
#' }
#' \subsection{Default 90\% intervals}{
#' We default to reporting 90\% intervals rather than 95\% intervals for several
#' reasons:
#' \itemize{
#'  \item Computational stability: 90\% intervals are more stable than 95\% 
#'  intervals (for which each end relies on only 2.5\% of the posterior draws).
#'  \item Relation to Type-S errors (Gelman and Carlin, 2014): 95\% of the 
#'  mass in a 90\% central interval is above the lower value and 95\% of the 
#'  mass is below the upper value. For a parameter \eqn{\theta}, it is therefore
#'  easy to see if the posterior probability that \eqn{\theta > 0} (or
#'  \eqn{\theta < 0}) is larger or smaller than 95\%.
#' }
#' Of course, if 95\% intervals are desired they can be computed by specifying
#' \code{prob=0.95}.
#' }
#' \subsection{Types of intervals}{
#' Currently \code{posterior_interval} only computes central intervals because
#' other types of intervals are rarely useful for the models that \pkg{rstanarm}
#' can estimate. Additional possibilities may be provided in future releases as
#' more models become available.
#' }
#' 
#' @seealso \code{\link{confint.stanreg}}, which, for models fit using 
#'   optimization, can be used to compute traditional confidence intervals.
#'   
#' @template reference-gelman-carlin
#' @template reference-morey
#' 
#' @examples 
#' posterior_interval(example_model)
#' posterior_interval(example_model, regex_pars = "herd")
#' posterior_interval(example_model, pars = "period2", prob = 0.5)
#' 
posterior_interval <- function(object, prob = 0.9, type = "central",
                               pars = NULL, regex_pars = NULL, ...) {
  if (!is.stanreg(object))
    stop(deparse(substitute(object)), " is not a stanreg object.", 
         call. = FALSE)
  if (used.optimizing(object))
    STOP_not_optimizing("posterior_interval")
  if (!identical(length(prob), 1L) || prob <= 0 || prob >= 1)
    stop("'prob' should be a single number greater than 0 and less than 1.", 
         call. = FALSE)
  if (!identical(type, "central"))
    stop("Currently the only option for 'type' is 'central'.", 
         call. = FALSE)
  
  mat <- as.matrix.stanreg(object, pars = pars, regex_pars = regex_pars)
  alpha <- (1 - prob) / 2
  probs <- c(alpha, 1 - alpha)
  labs <- paste0(100 * probs, "%")
  ci <- t(apply(mat, 2L, quantile, probs = probs))
  return(structure(ci, dimnames = list(colnames(mat), labs)))
}
