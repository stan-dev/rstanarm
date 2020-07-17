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
#' from the posterior predictive distribution of the fitted joint model, and then 
#' overlays the Kaplan-Meier curve based on the observed data.
#' 
#' @export
#' @templateVar stanjmArg object
#' @templateVar labsArg xlab,ylab
#' @templateVar cigeomArg ci_geom_args
#' @template args-stanjm-object
#' @template args-labs
#' @template args-ci-geom-args
#'   
#' @param check The type of plot to show. Currently only "survival" is 
#'   allowed, which compares the estimated marginal survival function under
#'   the joint model to the estimated Kaplan-Meier curve based on the 
#'   observed data.
#' @param limits A quoted character string specifying the type of limits to
#'   include in the plot. Can be one of: \code{"ci"} for the Bayesian
#'   posterior uncertainty interval (often known as a credible interval);
#'   or \code{"none"} for no interval limits.
#' @param draws An integer indicating the number of MCMC draws to use to 
#'   to estimate the survival function. The default and maximum number of 
#'   draws is the size of the posterior sample.
#' @param seed An optional \code{\link[=set.seed]{seed}} to use.
#' @param ... Optional arguments passed to 
#'   \code{\link[ggplot2:geom_path]{geom_line}} and used to control features
#'   of the plotted trajectory.
#' 
#' @return A ggplot object that can be further customized using the
#'   \pkg{ggplot2} package.
#'   
#' @seealso \code{\link{posterior_survfit}} for the estimated marginal or
#'   subject-specific survival function based on draws of the model parameters
#'   from the posterior distribution, 
#'   \code{\link{posterior_predict}} for drawing from the posterior 
#'   predictive distribution for the longitudinal submodel, and 
#'   \code{\link{pp_check}} for graphical checks of the longitudinal submodel.
#'    
#' @examples
#' \donttest{
#' if (!exists("example_jm")) example(example_jm)
#' # Compare estimated survival function to Kaplan-Meier curve
#' ps <- ps_check(example_jm)
#' ps + 
#'  ggplot2::scale_color_manual(values = c("red", "black")) + # change colors
#'  ggplot2::scale_size_manual(values = c(0.5, 3)) + # change line sizes 
#'  ggplot2::scale_fill_manual(values = c(NA, NA)) # remove fill
#' }
#' @importFrom ggplot2 ggplot aes_string geom_step
#' 
ps_check <- function(object, check = "survival", 
                     limits = c("ci", "none"),
                     draws = NULL, seed = NULL, 
                     xlab = NULL, ylab = NULL,
                     ci_geom_args = NULL, ...) {
  if (!requireNamespace("survival"))
    stop("the 'survival' package must be installed to use this function")
  
  validate_stanjm_object(object)
  limits <- match.arg(limits)

  # Predictions for plotting the estimated survival function
  dat <- posterior_survfit(object, standardise = TRUE, 
                           condition = FALSE, 
                           times = 0, extrapolate = TRUE, 
                           draws = draws, seed = seed)
  
  # Estimate KM curve based on response from the event submodel
  form <- reformulate("1", response = formula(object)$Event[[2]])
  coxdat <- object$survmod$mod$y
  if (is.null(coxdat)) 
    stop("Bug found: no response y found in the 'survmod' component of the ", 
         "fitted joint model.")
  resp <- attr(coxdat, "type")
  if (resp == "right") {
    form <- formula(survival::Surv(time, status) ~ 1)
  } else if (resp == "counting") {
    form <- formula(survival::Surv(start, stop, time) ~ 1)
  } else {
    stop("Bug found: only 'right' or 'counting' survival outcomes should ",
         "have been allowed as the response type in the fitted joint model.")
  }
  km <- survival::survfit(form, data = as.data.frame(unclass(coxdat)))
  kmdat <- data.frame(times = km$time, surv = km$surv,
                       lb = km$lower, ub = km$upper)
  
  # Plot estimated survival function with KM curve overlaid
  graph <- plot.survfit.stanjm(dat, ids = NULL, limits = limits, ...)
  kmgraph <- geom_step(data = kmdat, 
                       mapping = aes_string(x = "times", y = "surv"))
  graph + kmgraph
}



