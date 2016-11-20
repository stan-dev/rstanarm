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
#' @method plot stanreg
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
#'   The vignettes in the \pkg{bayesplot} package for many examples.
#'   
#'   \code{\link[bayesplot]{MCMC-overview}} (\pkg{bayesplot}) for links to the 
#'   documentation for all the available plotting functions.
#'   
#'   \code{\link[bayesplot]{color_scheme_set}} (\pkg{bayesplot}) to change the
#'   color scheme used for plotting.
#'   
#'   \code{\link{pp_check}} for graphical posterior predictive checks.
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
#' bayesplot::color_scheme_set("brightblue")
#' plot(fit, "areas", regex_pars = "period",
#'      prob = 0.5, prob_outer = 0.9)
#' 
#' \donttest{
#' # Make the same plot by extracting posterior draws and calling
#' # bayesplot::mcmc_areas directly
#' x <- as.array(fit, regex_pars = "period")
#' bayesplot::mcmc_areas(x, prob = 0.5, prob_outer = 0.9)
#' }
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
#' \donttest{
#' # same plot by calling bayesplot::mcmc_trace directly
#' x <- as.array(fit, pars = c("(Intercept)", "period2"))
#' bayesplot::mcmc_trace(x, facet_args = list(nrow = 2))
#' }
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
#'
#'
#' # For graphical posterior predictive checks see
#' # help("pp_check.stanreg")
#'
#' @importFrom ggplot2 ggplot aes_string xlab %+replace% theme
#'
plot.stanreg <- function(x, plotfun = "intervals", pars = NULL,
                         regex_pars = NULL, ...) {
  fun <- set_plotting_fun(plotfun)
  args <- set_plotting_args(x, pars, regex_pars, ..., plotfun = plotfun)
  do.call(fun, args)
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
  if (!used.sampling(x))
    validate_plotfun_for_opt_or_vb(plotfun)

  if (grepl("_nuts", plotfun, fixed = TRUE)) {
    nuts_stuff <- list(x = bayesplot::nuts_params(x), ...)
    if (!grepl("_energy", plotfun))
      nuts_stuff[["lp"]] <- bayesplot::log_posterior(x)
    return(nuts_stuff)
  }
  if (grepl("_rhat", plotfun, fixed = TRUE)) {
    rhat <- bayesplot::rhat(x, pars = pars, regex_pars = regex_pars)
    return(list(rhat = rhat, ...))
  }
  if (grepl("_neff", plotfun, fixed = TRUE)) {
    ratio <- bayesplot::neff_ratio(x, pars = pars, regex_pars = regex_pars)
    return(list(ratio = ratio, ...))
  }
  if (!is.null(pars) || !is.null(regex_pars)) {
    pars <- collect_pars(x, pars, regex_pars)
    pars <- allow_special_parnames(x, pars)
  }
  
  if (!used.sampling(x)) {
    if (!length(pars))
      pars <- NULL
    return(list(x = as.matrix(x, pars = pars), ...))
  }
  
  if (needs_chains(plotfun))
    list(x = as.array(x, pars = pars, regex_pars = regex_pars), ...)
  else
    list(x = as.matrix(x, pars = pars, regex_pars = regex_pars), ...)
}

mcmc_function_name <- function(fun) {
  # to keep backwards compatibility convert old function names
  if (fun == "scat") {
    fun <- "scatter"
  } else if (fun == "ess") {
    fun <- "neff"
  } else if (fun == "ac") {
    fun <- "acf"
  } else if (fun %in% c("diag", "stan_diag")) {
    stop(
      "For NUTS diagnostics, instead of 'stan_diag', ",
      "please specify the name of one of the functions listed at ",
      "help('NUTS', 'bayesplot')",
      call. = FALSE
    )
  }

  if (identical(substr(fun, 1, 4), "ppc_"))
    stop(
      "For 'ppc_' functions use the 'pp_check' ",
      "method instead of 'plot'.",
      call. = FALSE
    )

  if (!identical(substr(fun, 1, 5), "mcmc_"))
    fun <- paste0("mcmc_", fun)
  
  if (!fun %in% bayesplot::available_mcmc())
    stop(
      fun, " is not a valid MCMC function name.",  
      " Use bayesplot::available_mcmc() for a list of available MCMC functions."
    )

  return(fun)
}

# check if a plotting function requires multiple chains
needs_chains <- function(x) {
  nms <- paste0("mcmc_",
    c(
      "trace",
      "trace_highlight",
      "acf",
      "acf_bar",
      "hist_by_chain",
      "dens_overlay",
      "violin",
      "combo"
    )
  )
  mcmc_function_name(x) %in% nms
}

# Select the correct plotting function
# @param plotfun user specified plotfun argument (can be missing)
set_plotting_fun <- function(plotfun = NULL) {
  if (is.null(plotfun))
    return("mcmc_intervals")
  if (!is.character(plotfun))
    stop("'plotfun' should be a string.", call. = FALSE)

  plotfun <- mcmc_function_name(plotfun)
  fun <- try(get(plotfun, pos = asNamespace("bayesplot"), mode = "function"), 
             silent = TRUE)
  if (!inherits(fun, "try-error"))
    return(fun)
  
  stop(
    "Plotting function ",  plotfun, " not found. ",
    "A valid plotting function is any function from the ",
    "'bayesplot' package beginning with the prefix 'mcmc_'.",
    call. = FALSE
  )
}

# check if plotfun is ok to use with vb or optimization
validate_plotfun_for_opt_or_vb <- function(plotfun) {
  plotfun <- mcmc_function_name(plotfun)
  if (needs_chains(plotfun) || 
      grepl("rhat_|neff_|nuts_", plotfun))
    STOP_sampling_only(plotfun)
}

# Check for valid parameters
# @param x stanreg object
# @param pars user specified character vector
# check_plotting_pars <- function(x, pars, plotfun = character()) {
#   if (used.optimizing(x)) {
#     allpars <- c("alpha", "beta", rownames(x$stan_summary))
#   } else {
#     sim <- x$stanfit@sim
#     allpars <- c(sim$pars_oi, sim$fnames_oi)
#   }
#   m <- which(match(pars, allpars, nomatch = 0) == 0)
#   if (length(m) > 0)
#     stop("No parameter ", paste(pars[m], collapse = ', '),
#          call. = FALSE)
#   return(unique(pars))
# }




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
