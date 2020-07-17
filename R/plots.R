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
#' interface to the \link[bayesplot:MCMC-overview]{MCMC} module in the 
#' \pkg{\link{bayesplot}} package for plotting MCMC draws and diagnostics. It is also 
#' straightforward to use the functions from the \pkg{bayesplot} package directly rather than
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
#'   \link[bayesplot:MCMC-overview]{MCMC} function to use. The default is to call
#'   \code{\link[bayesplot:MCMC-intervals]{mcmc_intervals}}. \code{plotfun} can be specified
#'   either as the full name of a \pkg{bayesplot} plotting function (e.g.
#'   \code{"mcmc_hist"}) or can be abbreviated to the part of the name following
#'   the \code{"mcmc_"} prefix (e.g. \code{"hist"}). To get the names of all
#'   available MCMC functions see \code{\link[bayesplot:available_ppc]{available_mcmc}}.
#'
#' @param ... Additional arguments to pass to \code{plotfun} for customizing the
#'   plot. These are described on the help pages for the individual plotting 
#'   functions. For example, the arguments accepted for the default
#'   \code{plotfun="intervals"} can be found at
#'   \code{\link[bayesplot:MCMC-intervals]{mcmc_intervals}}.
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
#'   \item \code{\link[bayesplot:bayesplot-colors]{color_scheme_set}} (\pkg{bayesplot}) to change
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
#' # Ridgelines version of the areas plot
#' bayesplot::mcmc_areas_ridges(x, regex_pars = "period", prob = 0.9)
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
#' @importFrom ggplot2 ggplot aes_string xlab %+replace% theme
#'
plot.stanreg <- function(x, plotfun = "intervals", pars = NULL,
                         regex_pars = NULL, ...) {
  
  if (plotfun %in% c("pairs", "mcmc_pairs"))
    return(pairs.stanreg(x, pars = pars, regex_pars = regex_pars, ...))
  
  fun <- set_plotting_fun(plotfun)
  args <- set_plotting_args(x, pars, regex_pars, ..., plotfun = plotfun)
  do.call(fun, args)
}



# internal for plot.stanreg ----------------------------------------------

# Prepare argument list to pass to plotting function
#
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

  .plotfun_is_type <- function(patt) {
    grepl(pattern = paste0("_", patt), x = plotfun, fixed = TRUE)
  }
  
  if (.plotfun_is_type("nuts")) {
    nuts_stuff <- list(x = bayesplot::nuts_params(x), ...)
    if (!.plotfun_is_type("energy"))
      nuts_stuff[["lp"]] <- bayesplot::log_posterior(x)
    return(nuts_stuff)
  }
  if (.plotfun_is_type("rhat")) {
    rhat <- bayesplot::rhat(x, pars = pars, regex_pars = regex_pars)
    return(list(rhat = rhat, ...))
  }
  if (.plotfun_is_type("neff")) {
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
  
  list(x = as.array(x, pars = pars, regex_pars = regex_pars), ...)
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
  nms <- c(
    "trace",
    "trace_highlight",
    "rank",
    "rank_overlay",
    "acf",
    "acf_bar",
    "hist_by_chain",
    "dens_overlay",
    "violin",
    "combo"
  )
  mcmc_function_name(x) %in% paste0("mcmc_", nms)
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
      grepl("_rhat|_neff|_nuts_", plotfun))
    STOP_sampling_only(plotfun)
}



# pairs method ------------------------------------------------------------
#' Pairs method for stanreg objects
#' 
#' Interface to \pkg{bayesplot}'s \code{\link[bayesplot:MCMC-scatterplots]{mcmc_pairs}} function 
#' for use with \pkg{rstanarm} models. Be careful not to specify too
#' many parameters to include or the plot will be both hard to read and slow to
#' render.
#'
#' @method pairs stanreg
#' @export
#' @importFrom bayesplot pairs_style_np pairs_condition
#' @export pairs_style_np pairs_condition
#' @aliases pairs_style_np pairs_condition
#' 
#' @templateVar stanregArg x
#' @template args-stanreg-object
#' @template args-regex-pars
#' @param pars An optional character vetor of parameter names. All parameters 
#'   are included by default, but for models with more than just a few 
#'   parameters it may be far too many to visualize on a small computer screen 
#'   and also may require substantial computing time.
#' @param condition Same as the \code{condition} argument to 
#'   \code{\link[bayesplot:MCMC-scatterplots]{mcmc_pairs}} except the \emph{default is different}
#'   for \pkg{rstanarm} models. By default, the \code{mcmc_pairs} function in
#'   the \pkg{bayesplot} package plots some of the Markov chains (half, in the
#'   case of an even number of chains) in the panels above the diagonal and the
#'   other half in the panels below the diagonal. However since we know that 
#'   \pkg{rstanarm} models were fit using Stan (which \pkg{bayesplot} doesn't 
#'   assume) we can make the default more useful by splitting the draws 
#'   according to the \code{accept_stat__} diagnostic. The plots below the 
#'   diagonal will contain realizations that are below the median 
#'   \code{accept_stat__} and the plots above the diagonal will contain 
#'   realizations that are above the median \code{accept_stat__}. To change this
#'   behavior see the documentation of the \code{condition} argument at 
#'   \code{\link[bayesplot:MCMC-scatterplots]{mcmc_pairs}}.
#' @param ... Optional arguments passed to 
#'   \code{\link[bayesplot:MCMC-scatterplots]{mcmc_pairs}}. 
#'   The \code{np}, \code{lp}, and \code{max_treedepth} arguments to 
#'   \code{mcmc_pairs} are handled automatically by \pkg{rstanarm} and do not 
#'   need to be specified by the user in \code{...}. The arguments that can be 
#'   specified in \code{...} include \code{transformations}, \code{diag_fun},
#'   \code{off_diag_fun}, \code{diag_args}, \code{off_diag_args},
#'   and \code{np_style}. These arguments are
#'   documented thoroughly on the help page for
#'   \code{\link[bayesplot:MCMC-scatterplots]{mcmc_pairs}}.
#' 
#'   
#' @examples
#' \donttest{
#' if (!exists("example_model")) example(example_model)
#' 
#' bayesplot::color_scheme_set("purple")
#' 
#' # see 'condition' argument above for details on the plots below and 
#' # above the diagonal. default is to split by accept_stat__.
#' pairs(example_model, pars = c("(Intercept)", "log-posterior"))
#' 
#' pairs(
#'   example_model, 
#'   regex_pars = "herd:[2,7,9]", 
#'   diag_fun = "dens",
#'   off_diag_fun = "hex"
#' )
#' }
#' 
#' \donttest{
#' # for demonstration purposes, intentionally fit a model that
#' # will (almost certainly) have some divergences
#' fit <- stan_glm(
#'   mpg ~ ., data = mtcars,
#'   iter = 1000,
#'   # this combo of prior and adapt_delta should lead to some divergences
#'   prior = hs(),
#'   adapt_delta = 0.9,
#'   refresh = 0
#' )
#' 
#' pairs(fit, pars = c("wt", "sigma", "log-posterior"))
#' 
#' pairs(
#'   fit, 
#'   pars = c("wt", "sigma", "log-posterior"), 
#'   transformations = list(sigma = "log"), # show log(sigma) instead of sigma
#'   off_diag_fun = "hex" # use hexagonal heatmaps instead of scatterplots
#' )
#' 
#' 
#' bayesplot::color_scheme_set("brightblue")
#' pairs(
#'   fit, 
#'   pars = c("(Intercept)", "wt", "sigma", "log-posterior"), 
#'   transformations = list(sigma = "log"), 
#'   off_diag_args = list(size = 3/4, alpha = 1/3), # size and transparency of scatterplot points
#'   np_style = pairs_style_np(div_color = "black", div_shape = 2) # color and shape of the divergences
#' )
#' 
#' # Using the condition argument to show divergences above the diagonal 
#' pairs(
#'   fit, 
#'   pars = c("(Intercept)", "wt", "log-posterior"), 
#'   condition = pairs_condition(nuts = "divergent__")
#' )
#' 
#' }
#'
pairs.stanreg <-
  function(x,
           pars = NULL,
           regex_pars = NULL,
           condition = pairs_condition(nuts = "accept_stat__"),
           ...) {
    
    if (!used.sampling(x))
      STOP_sampling_only("pairs")
    
    dots <- list(...)
    ignored_args <- c("np", "lp", "max_treedepth")
    specified <- ignored_args %in% names(dots)
    if (any(specified)) {
      warning(
        "The following arguments were ignored because they are ",
        "specified automatically by rstanarm: ", 
        paste(sQuote(ignored_args[specified]), collapse = ", ")
      )
    }
    
    posterior <- as.array.stanreg(x, pars = pars, regex_pars = regex_pars)
    if (is.null(pars) && is.null(regex_pars)) {
      # include log-posterior by default
      lp_arr <- as.array.stanreg(x, pars = "log-posterior")
      dd <- dim(posterior)
      dn <- dimnames(posterior)
      dd[3] <- dd[3] + 1
      dn$parameters <- c(dn$parameters, "log-posterior")
      tmp <- array(NA, dim = dd, dimnames = dn)
      tmp[,, 1:(dd[3] - 1)] <- posterior
      tmp[,, dd[3]] <- lp_arr
      posterior <- tmp
    }
    posterior <- round(posterior, digits = 12)
    
    bayesplot::mcmc_pairs(
      x = posterior, 
      np = bayesplot::nuts_params(x),  
      lp = bayesplot::log_posterior(x),  
      max_treedepth = .max_treedepth(x),
      condition = condition,
      ...
    )
    
  }


# internal for pairs.stanreg ----------------------------------------------

# @param x stanreg object
.max_treedepth <- function(x) {
  control <- x$stanfit@stan_args[[1]]$control
  if (is.null(control)) {
    max_td <- 10
  } else {
    max_td <- control$max_treedepth
    if (is.null(max_td))
      max_td <- 10
  }
  return(max_td)
}
