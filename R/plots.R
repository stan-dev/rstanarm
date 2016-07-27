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

#' Plot method for stanreg objects
#' 
#' The \code{plot} method for \link{stanreg-objects} provides a convenient 
#' interface to the \pkg{bayesplot} package's plotting functionality for MCMC 
#' draws. It is also straightforward to use the functions from the 
#' \pkg{bayesplot} package directly rather than via the \code{plot.stanreg} 
#' method. Examples of both methods of plotting are given below.
#' 
#' @method plot stanreg
#' @export
#' @templateVar stanregArg x
#' @template args-stanreg-object
#' @template args-pars
#' @template args-regex-pars
#' @param plotfun A character string naming the \pkg{bayesplot} plotting
#'   function to apply to the stanreg object. See \link[bayesplot]{MCMC-overview} 
#'   for the available plots. Also see the Examples section below. 
#'   
#'   \code{plotfun} can be either the full name of a plotting function (e.g. 
#'   \code{"mcmc_hist"}) or can be abbreviated to the part of the name following
#'   the \code{"mcmc_"} prefix (e.g. \code{"hist"}). The default plot is 
#'   \code{\link[bayesplot]{mcmc_intervals}}, which shows intervals and point 
#'   estimates for all model parameters.
#'   
#' @param ... Additional arguments to pass to \code{plotfun} for customizing the
#'   plot.
#'
#' @return In most cases, a ggplot object (or several) that can be further 
#'   customized using the \pkg{ggplot2} package.
#'   
#' @seealso \code{\link[bayesplot]{MCMC-overview}} (\pkg{bayesplot}) for details
#'   on the individual plotting functions.
#'   
#' @examples
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
#' plot(fit, "areas", regex_pars = "period", 
#'      prob = 0.5, prob_outer = 0.9)
#'      
#' # Make the same plot by extracting posterior draws and calling 
#' # bayesplot::mcmc_areas directly
#' x <- as.array(fit, regex_pars = "period")
#' bayesplot::mcmc_areas(x, prob = 0.5, prob_outer = 0.9)
#' 
#' 
#' ##################
#' ### Traceplots ###
#' ##################
#' # note: rstanarm doesn't store the warmup draws by default 
#' (trace <- plot(fit, "trace", pars = "(Intercept)"))
#' 
#' # change traceplot colors to ggplot defaults or custom values
#' trace + ggplot2::scale_color_discrete()
#' trace + ggplot2::scale_color_manual(values = c("maroon", "skyblue2"))
#' 
#' # changing facet layout
#' plot(fit, "trace", pars = c("(Intercept)", "period2"), 
#'      facet_args = list(nrow = 2))
#' 
#' # same plot by calling bayesplot::mcmc_trace directly
#' x <- as.array(fit, pars = c("(Intercept)", "period2"))
#' bayesplot::mcmc_trace(x, facet_args = list(nrow = 2))
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
#' plot(fit, "scatter", pars = paste0("period", 2:3))
#' plot(fit, "scatter", pars = c("(Intercept)", "size"), 
#'      color = "black", size = 3, alpha = 0.5) +
#'      ggplot2::stat_ellipse(level = 0.9)
#' 
#' ######################################
#' ### Rhat and effective sample size ###
#' ######################################
#' plot(fit, "rhat")
#' plot(fit, "rhat_hist")
#' plot(fit, "neff")
#' plot(fit, "neff_hist")
#' 
#' ############
#' ### More ###
#' ############
#' plot(fit, regex_pars = "herd:1\\]")
#' plot(fit, regex_pars = "herd:[279]")
#' plot(fit, regex_pars = "herd:[279]|period2")
#' plot(fit, regex_pars = c("herd:[279]", "period2"))
#' 
#' # For graphical posterior predictive checks see 
#' # help("pp_check.stanreg")
#' 
#' @importFrom rstan stan_ac stan_diag stan_mcse stan_par quietgg
#' @importFrom ggplot2 ggplot aes_string xlab %+replace% theme
#' 
plot.stanreg <- function(x, plotfun = "intervals", pars = NULL, 
                         regex_pars = NULL, ...) {
  do.call(
    what = set_plotting_fun(plotfun), 
    args = set_plotting_args(x, pars, regex_pars, ..., 
                             plotfun = plotfun)
  )
}

# Check for valid parameters
# @param x stanreg object
# @param pars user specified character vector
check_plotting_pars <- function(x, pars, plotfun = character()) {
  if (used.optimizing(x)) {
    allpars <- c("alpha", "beta", rownames(x$stan_summary))
  } else {
    sim <- x$stanfit@sim
    allpars <- c(sim$pars_oi, sim$fnames_oi)
  }
  m <- which(match(pars, allpars, nomatch = 0) == 0)
  if (length(m) > 0) 
    stop("No parameter ", paste(pars[m], collapse = ', '), 
         call. = FALSE) 
  return(unique(pars))
}

# Prepare argument list to pass to plotting function
# @param x stanreg object
# @param pars, regex_pars user specified pars and regex_pars arguments (can be
#   missing)
# @param ...  additional arguments to pass to the plotting function
# @param plotfun User's 'plotfun' argument
set_plotting_args <- function(x, pars = NULL, regex_pars = NULL, ..., 
                              plotfun = character()) {
  
  plotfun <- mcmc_function_name(plotfun)
  if (grepl("_rhat", plotfun, fixed = TRUE)) {
    rhat <- rhat(x, pars = pars, regex_pars = regex_pars)
    return(list(rhat = rhat, ...))
  }
  if (grepl("_neff", plotfun, fixed = TRUE)) {
    ratio <- neff_ratio(x, pars = pars, regex_pars = regex_pars)
    return(list(ratio = ratio, ...))
  }
  if (!is.null(pars) || !is.null(regex_pars)) {
    pars <- collect_pars(x, pars, regex_pars)
    pars <- allow_special_parnames(x, pars)
  }
  if (needs_chains(plotfun))
    list(x = as.array(x, pars = pars, regex_pars = regex_pars), ...)
  else
    list(x = as.matrix(x, pars = pars, regex_pars = regex_pars), ...)
}

mcmc_function_name <- function(fun) {
  if (fun == "scat") {
    fun <- "scatter"
  } else if (fun == "ess") {
    fun <- "neff"
  }
  
  if (identical(substr(fun, 1, 4), "ppc_"))
    stop(
      "For 'ppc_' functions use the 'pp_check' ", 
      "method instead of 'plot'.", 
      call. = FALSE
    )
  
  if (!identical(substr(fun, 1, 5), "mcmc_"))
    fun <- paste0("mcmc_", fun)
  
  return(fun)
}

# check if a plotting function requires multiple chains
needs_chains <- function(x) {
  nms <- c("trace",
           "trace_highlight",
           "hist_by_chain",
           "dens_overlay",
           "violin",
           "combo")
  x %in% paste0("mcmc_", nms)
}

# Select the correct plotting function
# @param plotfun user specified plotfun argument (can be missing)
set_plotting_fun <- function(plotfun = NULL) {
  if (is.null(plotfun))
    return("mcmc_intervals")
  if (!is.character(plotfun))
    stop("'plotfun' should be a string.", call. = FALSE)
  
  plotfun <- mcmc_function_name(plotfun)
  fun <- try(match.fun(plotfun), silent = TRUE)
  if (!inherits(fun, "try-error"))
    return(fun)
  
  stop(
    "Plotting function ",  plotfun, " not found. ",
    "A valid plotting function is any function from the ",
    "'bayesplot' package beginning with the prefix 'mcmc_'.",
    call. = FALSE
  )
}

# function calling arm::coefplot (only used for models fit using optimization)
stan_plot_opt <- function(x, pars = NULL, varnames = NULL, ...) {
  if (!requireNamespace("arm", quietly = TRUE)) 
    stop("Please install the 'arm' package to use this feature.", 
         call. = FALSE)
  stopifnot(used.optimizing(x))
  coefs <- coef(x)
  sds <- se(x)
  nms <- varnames %ORifNULL% names(x$coefficients)
  if (!is.null(pars)) {
    mark <- NA
    if ("alpha" %in% pars) 
      mark <- c(mark, "(Intercept)")
    if ("beta" %in% pars) 
      mark <- c(mark, setdiff(names(x$coefficients), "(Intercept)"))
    mark <- c(mark, setdiff(pars, c("alpha", "beta")))
    mark <- mark[!is.na(mark)] 
    coefs <- coefs[mark]
    sds <- sds[mark]
    if (is.null(varnames)) 
      nms <- nms[nms %in% mark]
  }
  arm::coefplot.default(coefs = coefs, sds = sds, varnames = nms, ...)
}

#' Pairs method for stanreg objects
#' 
#' @method pairs stanreg
#' @export
#' @param x A stanreg object returned by one of the rstanarm modeling functions.
#' @param ... Arguments to pass to \code{\link[rstan]{pairs.stanfit}}. 
#' 
#' @description See \code{\link[rstan]{pairs.stanfit}} for details.
#' @details See the Details section in \code{\link[rstan]{pairs.stanfit}}.
#' @importFrom graphics pairs
#' @examples
#' if (!exists("example_model")) example(example_model)
#' pairs(example_model, pars = c("(Intercept)", "log-posterior"))
#' 
pairs.stanreg <- function(x, ...) {
  if (!used.sampling(x)) 
    STOP_sampling_only("pairs")
  requireNamespace("rstan")
  requireNamespace("KernSmooth")
  
  if (is.mer(x)) {
    dots <- list(...)
    if (is.null(dots[["pars"]]))
      return(pairs(x$stanfit, pars = names(fixef(x)), ...))
    
    b <- b_names(rownames(x$stan_summary), value = TRUE)
    if (any(dots[["pars"]] %in% b))
      stop("pairs.stanreg does not yet allow group-level parameters in 'pars'.")
  }
  
  pairs(x$stanfit, ...)
}

#' Plots for rstanarm models
#' 
#' Models fit using \code{algorithm='sampling'}, \code{"meanfield"}, or
#' \code{"fullrank"} are compatible with a variety of plotting functions from
#' the \pkg{rstan} package. Each function returns at least one
#' \code{\link[ggplot2]{ggplot}} object that can be customized further using the
#' \pkg{ggplot2} package. The plotting functions described here can be called
#' using the \code{\link[=plot.stanreg]{plot method}} for stanreg objects 
#' without loading the \pkg{rstan} package. For example usage see 
#' \code{\link{plot.stanreg}}.
#' 
#' @name rstanarm-plots
#' 
#' @section Plotting functions:
#' 
#' \describe{
#' \item{Posterior intervals and point
#' estimates}{\code{\link[rstan]{stan_plot}}}
#' \item{Traceplots}{\code{\link[rstan]{stan_trace}}}
#' \item{Histograms}{\code{\link[rstan]{stan_hist}}}
#' \item{Kernel density estimates}{\code{\link[rstan]{stan_dens}}}
#' \item{Scatterplots}{\code{\link[rstan]{stan_scat}}}
#' \item{Diagnostics for Hamiltonian Monte Carlo and the No-U-Turn
#' Sampler}{\code{\link[rstan]{stan_diag}}}
#' \item{Rhat}{\code{\link[rstan]{stan_rhat}}}
#' \item{Ratio of effective sample size to total posterior sample
#' size}{\code{\link[rstan]{stan_ess}}}
#' \item{Ratio of Monte Carlo standard error to posterior standard
#' deviation}{\code{\link[rstan]{stan_mcse}}}
#' \item{Autocorrelation}{\code{\link[rstan]{stan_ac}}}
#' }
#' 
#' @seealso \code{\link{plot.stanreg}} for how to call the \code{plot} method, 
#'   \code{\link{shinystan}} for interactive model exploration,
#'   \code{\link{pp_check}} for graphical posterior predicive checking.
#'
#' @examples
#' # See examples at help("plot.stanreg", package = "rstanarm")
NULL
