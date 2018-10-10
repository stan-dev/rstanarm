# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2015, 2016, 2017 Trustees of Columbia University
# Copyright (C) 2016, 2017 Sam Brilleman
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
#

#' Graphical checks of the estimated survival function
#' 
#' This function plots the estimated marginal survival function based on draws 
#' from the posterior predictive distribution of the fitted model, 
#' and then overlays a Kaplan-Meier curve based on the observed data.
#' 
#' @importFrom ggplot2 ggplot aes_string geom_step
#' @export
#' @templateVar stanregArg object
#' @templateVar labsArg xlab,ylab
#' @templateVar cigeomArg ci_geom_args
#' @template args-stansurv-stanjm-object
#' @template args-labs
#' @template args-ci-geom-args
#'   
#' @param check The type of plot to show. Currently only "survival" is 
#'   allowed, which compares the estimated marginal survival function 
#'   under the fitted model to the estimated Kaplan-Meier curve based 
#'   on the observed data.
#' @param limits A quoted character string specifying the type of limits to
#'   include in the plot. Can be one of: \code{"ci"} for the Bayesian
#'   posterior uncertainty interval (often known as a credible interval);
#'   or \code{"none"} for no interval limits.
#' @param draws An integer indicating the number of MCMC draws to use to 
#'   to estimate the survival function. The default and maximum number of 
#'   draws is the size of the posterior sample.
#' @param seed An optional \code{\link[=set.seed]{seed}} to use.
#' @param ... Optional arguments passed to 
#'   \code{\link[ggplot2]{geom_line}} and used to control features
#'   of the plotted trajectory.
#' 
#' @return A ggplot object that can be further customized using the
#'   \pkg{ggplot2} package.
#'   
#' @seealso 
#'   \code{\link{posterior_survfit}} for the estimated marginal or
#'   subject-specific survival function based on draws of the model parameters
#'   from the posterior distribution \cr
#'   \code{\link{posterior_predict}} for drawing from the posterior 
#'   predictive distribution for the longitudinal submodel (for 
#'   \code{\link{stan_jm}} models only) \cr
#'   \code{\link{pp_check}} for graphical checks of the longitudinal submodel
#'   (for \code{\link{stan_jm}} models only)
#'    
#' @examples
#' \donttest{
#' if (!exists("example_jm")) example(example_jm)
#' # Compare estimated survival function to Kaplan-Meier curve
#' ps <- ps_check(example_jm)
#' ps + 
#'   ggplot2::scale_color_manual(values = c("red", "black")) + # change colors
#'   ggplot2::scale_size_manual (values = c(0.5, 3)) +         # change line sizes 
#'   ggplot2::scale_fill_manual (values = c(NA, NA))           # remove fill
#' }
#' 
ps_check <- function(object, check = "survival", 
                     limits = c("ci", "none"),
                     draws = NULL, seed = NULL, 
                     xlab = NULL, ylab = NULL,
                     ci_geom_args = NULL, ...) {
  
  if (!requireNamespace("survival"))
    stop("the 'survival' package must be installed to use this function")
  
  if (!any(is.stansurv(object), is.stanjm(object)))
    stop("Object is not a 'stansurv' or 'stanjm' object.")
  
  limits <- match.arg(limits)

  # Obtain standardised survival probabilities for the fitted model
  dat <- posterior_survfit(object, 
                           times       = 0, 
                           extrapolate = TRUE, 
                           standardise = TRUE, 
                           condition   = FALSE, 
                           draws       = draws, 
                           seed        = seed)
  
  # Obtain the response variable for the fitted model
  response <- get_surv(object)
  if (is.null(response)) 
    stop("Bug found: no response variable found in fitted model object.")

  # Obtain the formula for KM curve
  type <- attr(response, "type")
  form <- switch(type,
                 right    = formula(survival::Surv(time, status, type = type) ~ 1),
                 counting = formula(survival::Surv(start, stop,  status, type = type) ~ 1),
                 interval = formula(survival::Surv(time1, time2, status, type = 'interval') ~ 1),
                 interval2= formula(survival::Surv(time1, time2, status, type = 'interval') ~ 1),
                 stop("Bug found: invalid type of survival object."))
  
  # Obtain the KM estimates
  kmfit <- survival::survfit(form, data = data.frame(unclass(response)))
  kmdat <- data.frame(times = kmfit$time, surv = kmfit$surv)
  
  # Plot estimated survival function with KM curve overlaid
  psgraph <- plot.survfit.stanjm(dat, ids = NULL, limits = limits, ...)
  kmgraph <- geom_step(aes_string(x = "times", y = "surv"), kmdat)
  psgraph + kmgraph
}



