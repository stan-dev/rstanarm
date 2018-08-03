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
#
#' Plot method for stanreg objects
#'
#' The \code{plot} method for \link{stanreg-objects} provides a convenient 
#' interface to the \link[bayesplot]{MCMC} module in the \pkg{\link{bayesplot}} 
#' package for plotting MCMC draws and diagnostics. It is also straightforward 
#' to use the functions from the \pkg{bayesplot} package directly rather than 
#' via the \code{plot} method. Examples of both methods of plotting are given
#' below.
#'
#' @method plot stansurv
#' @export
#' @templateVar stanregArg x
#' @template args-stanreg-object
#' @template args-pars
#' @template args-regex-pars
#' @param plotfun A character string naming the \pkg{bayesplot} 
#'   \link[bayesplot]{MCMC} function to use. The default is to call
#'   \code{\link[bayesplot]{mcmc_intervals}}. \code{plotfun} can be specified
#'   either as the full name of a \pkg{bayesplot} plotting function (e.g.
#'   \code{"mcmc_hist"}) or can be abbreviated to the part of the name following
#'   the \code{"mcmc_"} prefix (e.g. \code{"hist"}). To get the names of all
#'   available MCMC functions see \code{\link[bayesplot]{available_mcmc}}.
#'
#' @param ... Additional arguments to pass to \code{plotfun} for customizing the
#'   plot. These are described on the help pages for the individual plotting 
#'   functions. For example, the arguments accepted for the default
#'   \code{plotfun="intervals"} can be found at
#'   \code{\link[bayesplot]{mcmc_intervals}}.
#'
#' @return Either a ggplot object that can be further customized using the
#'   \pkg{ggplot2} package, or an object created from multiple ggplot objects
#'   (e.g. a gtable object created by \code{\link[gridExtra]{arrangeGrob}}).
#'
#' @seealso
#' \itemize{ 
#'   \item The vignettes in the \pkg{bayesplot} package for many examples.
#'   \item \code{\link[bayesplot]{MCMC-overview}} (\pkg{bayesplot}) for links to
#'   the documentation for all the available plotting functions.
#'   \item \code{\link[bayesplot]{color_scheme_set}} (\pkg{bayesplot}) to change
#'   the color scheme used for plotting.
#'   \item \code{\link{pp_check}} for graphical posterior predictive checks.
#'   \item \code{\link{plot_nonlinear}} for models with nonlinear smooth 
#'   functions fit using \code{\link{stan_gamm4}}.
#' }  
#'
#' @template reference-bayesvis
#' 
#' @examples
#' \donttest{
#' # Use rstanarm example model
#' if (!exists("example_model")) example(example_model)
#' fit <- example_model
#'
#' #####################################
#' ### Intervals and point estimates ###
#' #####################################
#' plot(fit) # same as plot(fit, "intervals"), plot(fit, "mcmc_intervals")
#'
#' p <- plot(fit, pars = "size", regex_pars = "period",
#'           prob = 0.5, prob_outer = 0.9)
#' p + ggplot2::ggtitle("Posterior medians \n with 50% and 90% intervals")
#'
#' # Shaded areas under densities
#' bayesplot::color_scheme_set("brightblue")
#' plot(fit, "areas", regex_pars = "period",
#'      prob = 0.5, prob_outer = 0.9)
#' 
#' # Make the same plot by extracting posterior draws and calling
#' # bayesplot::mcmc_areas directly
#' x <- as.array(fit, regex_pars = "period")
#' bayesplot::mcmc_areas(x, prob = 0.5, prob_outer = 0.9)
#'
#'
#' ##################################
#' ### Histograms & density plots ###
#' ##################################
#' plot_title <- ggplot2::ggtitle("Posterior Distributions")
#' plot(fit, "hist", regex_pars = "period") + plot_title
#' plot(fit, "dens_overlay", pars = "(Intercept)",
#'      regex_pars = "period") + plot_title
#'
#' ####################
#' ### Scatterplots ###
#' ####################
#' bayesplot::color_scheme_set("teal")
#' plot(fit, "scatter", pars = paste0("period", 2:3))
#' plot(fit, "scatter", pars = c("(Intercept)", "size"),
#'      size = 3, alpha = 0.5) +
#'      ggplot2::stat_ellipse(level = 0.9)
#'
#'
#' ####################################################
#' ### Rhat, effective sample size, autocorrelation ###
#' ####################################################
#' bayesplot::color_scheme_set("red")
#' 
#' # rhat
#' plot(fit, "rhat")
#' plot(fit, "rhat_hist")
#' 
#' # ratio of effective sample size to total posterior sample size
#' plot(fit, "neff")
#' plot(fit, "neff_hist")
#' 
#' # autocorrelation by chain
#' plot(fit, "acf", pars = "(Intercept)", regex_pars = "period")
#' plot(fit, "acf_bar", pars = "(Intercept)", regex_pars = "period")
#' 
#' 
#' ##################
#' ### Traceplots ###
#' ##################
#' # NOTE: rstanarm doesn't store the warmup draws (to save space because they
#' # are not so essential for diagnosing the particular models implemented in
#' # rstanarm) so the iterations in the traceplot are post-warmup iterations
#' 
#' bayesplot::color_scheme_set("pink")
#' (trace <- plot(fit, "trace", pars = "(Intercept)"))
#'
#' # change traceplot colors to ggplot defaults or custom values
#' trace + ggplot2::scale_color_discrete()
#' trace + ggplot2::scale_color_manual(values = c("maroon", "skyblue2"))
#'
#' # changing facet layout 
#' plot(fit, "trace", pars = c("(Intercept)", "period2"),
#'      facet_args = list(nrow = 2))
#' # same plot by calling bayesplot::mcmc_trace directly
#' x <- as.array(fit, pars = c("(Intercept)", "period2"))
#' bayesplot::mcmc_trace(x, facet_args = list(nrow = 2))
#'
#'
#' ############
#' ### More ###
#' ############
#'
#' # regex_pars examples
#' plot(fit, regex_pars = "herd:1\\]")
#' plot(fit, regex_pars = "herd:[279]")
#' plot(fit, regex_pars = "herd:[279]|period2")
#' plot(fit, regex_pars = c("herd:[279]", "period2"))
#' }
#'
#' # For graphical posterior predictive checks see
#' # help("pp_check.stanreg")
#'
plot.stansurv <- function(x, plotfun = "basehaz", prob = 0.95, ci = TRUE, 
                          ci_geom_args = NULL, pars = NULL,
                          regex_pars = NULL, ...) {
  
  validate_stansurv_object(x)

  if (!plotfun %in% c("basehaz", "tde"))
    NextMethod("plot")
    
  stanpars <- extract_pars(x)
  has_intercept <- check_for_intercept(x$basehaz)
  
  t_min <- min(x$entrytime)
  t_max <- max(x$exittime)
  times <- seq(t_min, t_max, by = (t_max - t_min) / 1000)
  
  if (plotfun == "basehaz") { 
    
    if (!is.null(pars))
      warning2("'pars' is ignored when plotting the baseline hazard.")
    if (!is.null(regex_pars))
      warning2("'regex_pars' is ignored when plotting the baseline hazard.")
    
    basehaz <- evaluate_basehaz(times, x$basehaz, stanpars$bhcoef, stanpars$alpha)
    basehaz <- median_and_bounds(basehaz, prob, na.rm = TRUE)
    plotdat <- data.frame(times, basehaz)
    
    ylab <- "Baseline hazard rate"
    xlab <- "Time"
    
  } else if (plotfun == "tde") {
    
    smooth_map   <- get_smooth_name(x$s_events, type = "smooth_map")
    smooth_vars  <- get_smooth_name(x$s_events, type = "smooth_vars")
    smooth_coefs <- get_smooth_name(x$s_events, type = "smooth_coefs")
    
    if (is.null(pars))
      pars <- smooth_vars
    if (length(pars) > 1L)
      stop2("Only one variable can be specified in 'pars' .")
    if (!pars %in% smooth_vars)
      stop2("Cannot find variable '", pars, "' amongst the tde terms.")
    
    sel1 <- which(smooth_vars == pars)
    sel2 <- smooth_coefs[smooth_map == sel1]
    
    betas_tf <- stanpars$beta    [, pars, drop = FALSE]
    betas_td <- stanpars$beta_tde[, sel2, drop = FALSE]
    betas    <- cbind(betas_tf, betas_td)

    times__ <- times
    basis   <- eval(parse(text = x$formula$td_basis[sel1]))
    basis   <- add_intercept(basis)
    log_hr  <- linear_predictor(betas, basis)
    plotdat <- median_and_bounds(exp(log_hr), prob, na.rm = TRUE)
    plotdat <- data.frame(times, plotdat)
    
    ylab <- "Hazard ratio"
    xlab <- "Time"
    
  }
  
  geom_defs <- list(color = "black")  # default plot args
  geom_args <- set_geom_args(geom_defs, ...)  
  geom_ylab <- ggplot2::ylab(ylab)
  geom_xlab <- ggplot2::xlab(xlab)
  geom_maps <- list(mapping = aes_string(x = "times", y = "med"))
  geom_base <- ggplot(plotdat) + geom_ylab + geom_xlab + theme_bw()
  geom_plot <- geom_base + do.call("geom_line", c(geom_maps, geom_args))
  if (ci) {
    lim_defs <- list(alpha = 0.3) # default plot args for ci
    lim_args <- c(defaults = list(lim_defs), ci_geom_args)
    lim_args <- do.call("set_geom_args", lim_args)
    lim_maps <- list(mapping = aes_string(x = "times", ymin = "lb", ymax = "ub"))
    lim_plot <- do.call("geom_ribbon", c(lim_maps, lim_args))
  } else {
    lim_plot <- NULL
  }
  return(geom_plot + lim_plot)
  
}