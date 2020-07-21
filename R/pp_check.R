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
#' Graphical posterior predictive checks
#' 
#' Interface to the \link[bayesplot:PPC-overview]{PPC} (posterior predictive checking) module
#' in the \pkg{\link{bayesplot}} package, providing various plots comparing the 
#' observed outcome variable \eqn{y} to simulated datasets \eqn{y^{rep}}{yrep} 
#' from the posterior predictive distribution. The \code{pp_check} method for 
#' \link{stanreg-objects} prepares the arguments required for the specified 
#' \pkg{bayesplot} PPC plotting function and then calls that function. It is 
#' also straightforward to use the functions from the \pkg{bayesplot} package 
#' directly rather than via the \code{pp_check} method. Examples of both are
#' given below.
#' 
#' @export
#' @export pp_check
#' @aliases pp_check
#' @method pp_check stanreg
#' @templateVar bdaRef (Ch. 6)
#' @templateVar stanregArg object
#' @template reference-bda
#' @template reference-bayesvis
#' @template args-stanreg-object
#' @param plotfun A character string naming the \pkg{bayesplot} 
#'   \link[bayesplot:PPC-overview]{PPC} function to use. The default is to call
#'   \code{\link[bayesplot:PPC-distributions]{ppc_dens_overlay}}. \code{plotfun} can be specified
#'   either as the full name of a \pkg{bayesplot} plotting function (e.g.
#'   \code{"ppc_hist"}) or can be abbreviated to the part of the name following
#'   the \code{"ppc_"} prefix (e.g. \code{"hist"}). To get the names of all
#'   available PPC functions see \code{\link[bayesplot]{available_ppc}}.
#' @param nreps The number of \eqn{y^{rep}}{yrep} datasets to generate from the 
#'   \link[=posterior_predict]{posterior predictive distribution} and show in
#'   the plots. The default depends on \code{plotfun}. For functions that plot
#'   each \code{yrep} dataset separately (e.g. \code{ppc_hist}), \code{nreps}
#'   defaults to a small value to make the plots readable. For functions that
#'   overlay many \code{yrep} datasets (e.g., \code{ppc_dens_overlay}) a larger
#'   number is used by default, and for other functions (e.g. \code{ppc_stat})
#'   the default is to set \code{nreps} equal to the posterior sample size.
#' @param ... Additonal arguments passed to the \pkg{\link{bayesplot}} function 
#'   called. For many plotting functions \code{...} is optional, however for 
#'   functions that require a \code{group} or \code{x} argument, these arguments
#'   should be specified in \code{...}. If specifying \code{group} and/or 
#'   \code{x}, they can be provided as either strings naming variables (in which
#'   case they are searched for in the model frame) or as vectors containing the
#'   actual values of the variables. See the \strong{Examples} section, below.
#' @param seed An optional \code{\link[=set.seed]{seed}} to pass to 
#'   \code{\link{posterior_predict}}.
#' 
#' @return \code{pp_check} returns a ggplot object that can be further
#'   customized using the \pkg{ggplot2} package.
#' 
#' @note For binomial data, plots of \eqn{y} and \eqn{y^{rep}}{yrep} show the 
#'   proportion of 'successes' rather than the raw count. Also for binomial 
#'   models see \code{\link[bayesplot:PPC-errors]{ppc_error_binned}} for binned residual
#'   plots.
#' 
#' @seealso 
#' \itemize{
#'   \item The vignettes in the \pkg{bayesplot} package for many examples.
#'     Examples of posterior predictive checks can also be found in the
#'     \pkg{rstanarm} vignettes and demos.
#'   \item \code{\link[bayesplot]{PPC-overview}} (\pkg{bayesplot}) for links to 
#'     the documentation for all the available plotting functions.
#'   \item \code{\link{posterior_predict}} for drawing from the posterior 
#'     predictive distribution. 
#'   \item \code{\link[bayesplot:bayesplot-colors]{color_scheme_set}} to change the color scheme 
#'     of the plots.
#' }
#' 
#' @examples 
#' fit <- stan_glmer(
#'   mpg ~ wt + am + (1|cyl), 
#'   data = mtcars, 
#'   iter = 400, # iter and chains small just to keep example quick
#'   chains = 2, 
#'   refresh = 0
#' ) 
#' 
#' # Compare distribution of y to distributions of multiple yrep datasets
#' pp_check(fit)
#' pp_check(fit, plotfun = "boxplot", nreps = 10, notch = FALSE)
#' pp_check(fit, plotfun = "hist", nreps = 3)
#' 
#' \donttest{
#' # Same plot (up to RNG noise) using bayesplot package directly
#' bayesplot::ppc_hist(y = mtcars$mpg, yrep = posterior_predict(fit, draws = 3))
#'
#' # Check histograms of test statistics by level of grouping variable 'cyl'
#' pp_check(fit, plotfun = "stat_grouped", stat = "median", group = "cyl")
#' 
#' # Defining a custom test statistic
#' q25 <- function(y) quantile(y, probs = 0.25) 
#' pp_check(fit, plotfun = "stat_grouped", stat = "q25", group = "cyl")
#' 
#' # Scatterplot of two test statistics
#' pp_check(fit, plotfun = "stat_2d", stat = c("mean", "sd"))
#' 
#' # Scatterplot of y vs. average yrep
#' pp_check(fit, plotfun = "scatter_avg") # y vs. average yrep

#' # Same plot (up to RNG noise) using bayesplot package directly
#' bayesplot::ppc_scatter_avg(y = mtcars$mpg, yrep = posterior_predict(fit))
#' 
#' # Scatterplots of y vs. several individual yrep datasets
#' pp_check(fit, plotfun = "scatter", nreps = 3)
#'
#' # Same plot (up to RNG noise) using bayesplot package directly
#' bayesplot::ppc_scatter(y = mtcars$mpg, yrep = posterior_predict(fit, draws = 3))
#' 
#' # yrep intervals with y points overlaid
#' # by default 1:length(y) used on x-axis but can also specify an x variable
#' pp_check(fit, plotfun = "intervals")
#' pp_check(fit, plotfun = "intervals", x = "wt") + ggplot2::xlab("wt")
#'
#' # Same plot (up to RNG noise) using bayesplot package directly
#' bayesplot::ppc_intervals(y = mtcars$mpg, yrep = posterior_predict(fit), 
#'                          x = mtcars$wt) + ggplot2::xlab("wt")
#' 
#' # predictive errors
#' pp_check(fit, plotfun = "error_hist", nreps = 6)
#' pp_check(fit, plotfun = "error_scatter_avg_vs_x", x = "wt") + 
#'   ggplot2::xlab("wt")
#'   
#' # Example of a PPC for ordinal models (stan_polr)
#' fit2 <- stan_polr(tobgp ~ agegp, data = esoph, method = "probit",
#'                   prior = R2(0.2, "mean"), init_r = 0.1, 
#'                   refresh = 0)
#' pp_check(fit2, plotfun = "bars", nreps = 500, prob = 0.5)
#' pp_check(fit2, plotfun = "bars_grouped", group = esoph$agegp, 
#'          nreps = 500, prob = 0.5)
#' }
pp_check.stanreg <-
  function(object,
           plotfun = "dens_overlay",
           nreps = NULL,
           seed = NULL,
           ...) {
    if (used.optimizing(object))
      STOP_not_optimizing("pp_check")
    
    if (is.stanmvreg(object)) {
      dots <- list(...)
      m <- dots[["m"]]
      if (is.null(m))
        stop("Argument 'm' must be provided for stanmvreg objects.")
    } else m <- NULL
    
    plotfun_name <- .ppc_function_name(plotfun)
    plotfun <- get(plotfun_name, pos = asNamespace("bayesplot"), mode = "function")
    is_binomial_model <- is_binomial_ppc(object, m = m)
    
    y_yrep <-
      .ppc_y_and_yrep(
        object,
        seed = seed,
        nreps = .set_nreps(nreps, fun = plotfun_name),
        binned_resid_plot = isTRUE(plotfun_name == "ppc_error_binned"),
        ...
      )
    args <-
      .ppc_args(
        object,
        y = y_yrep[["y"]],
        yrep = y_yrep[["yrep"]],
        fun = plotfun_name,
        ...
      )
    
    do.call(plotfun, args)
  }




# internal ----------------------------------------------------------------

# check if binomial
is_binomial_ppc <- function(object, ...) {
  if (is_polr(object) && !is_scobit(object)) {
    FALSE
  } else {
    is.binomial(family(object, ...)$family)
  }
}

# prepare y and yrep arguments to bayesplot function
.ppc_y_and_yrep <-
  function(object,
           nreps = NULL,
           seed = NULL,
           binned_resid_plot = FALSE, 
           ...) {
    y <- get_y(object, ...)
    if (binned_resid_plot) {
      yrep <- posterior_epred(object, ...)
      yrep <- yrep[1:nreps, , drop = FALSE]
    } else {
      yrep <- posterior_predict(object, draws = nreps, seed = seed, ...)
    }
    
    if (is_binomial_ppc(object, ...)) { # includes stan_polr's scobit models
      if (NCOL(y) == 2L) {
        trials <- rowSums(y)
        y <- y[, 1L] / trials
        if (!binned_resid_plot)
          yrep <- sweep(yrep, 2L, trials, "/")
      } else if (is.factor(y))
        y <- fac2bin(y)
    } else if (is_polr(object)) { # excluding scobit
      y <- as.integer(y)
      yrep <- polr_yrep_to_numeric(yrep)
    }
    
    nlist(y, yrep)
  }


# prepare 'group' and 'x' variable for certain plots
.ppc_xvar <- .ppc_groupvar <- function(object, var = NULL, ...) {
  if (is.null(var) || !is.character(var))
    return(var)

  mf <- model.frame(object, ...)
  vars <- colnames(mf)
  if (var %in% vars)
    return(mf[, var])
  
  stop("Variable '", var, "' not found in model frame. ")
}

# 
# @param fun user's plotfun argument
.ppc_function_name <- function(fun = character()) {
  if (!length(fun))
    stop("Plotting function not specified.", call. = FALSE)
  
  if (identical(substr(fun, 1, 5), "mcmc_"))
    stop(
      "For 'mcmc_' functions use the 'plot' ",
      "method instead of 'pp_check'.",
      call. = FALSE
    )
  
  if (!identical(substr(fun, 1, 4), "ppc_"))
    fun <- paste0("ppc_", fun)
  
  if (fun == "ppc_loo_pit") {
    warning(
      "'ppc_loo_pit' is deprecated. ", 
      "Use 'ppc_loo_pit_overlay' or 'ppc_loo_pit_qq' instead.", 
      call.=FALSE
    )
    fun <- "ppc_loo_pit_qq"
  }
  if (!fun %in% bayesplot::available_ppc())
    stop(
      fun, " is not a valid PPC function name.",  
      " Use bayesplot::available_ppc() for a list of available PPC functions."
    )
  
  return(fun)
}


# prepare all arguments to pass to bayesplot function
# @param object user's object
# @param y,yrep returned by .ppc_y_and_yrep
# @param fun string returned by .ppc_function_name
# @param ... user's ...
# @return named list
#
.ppc_args <- function(object, y, yrep, fun, ...) {
  funname <- fun
  fun <- match.fun(fun)
  dots <- list(...)
  dots[["y"]] <- as.numeric(y)
  dots[["yrep"]] <- yrep
  argnames <- names(formals(fun))
  
  if (is.stanmvreg(object)) {
    m <- dots[["m"]]
    if (is.null(m))
      stop("Argument 'm' must be provided for stanmvreg objects.")
    dots[["m"]] <- NULL # don't return m as part of bayesplot arguments
  }
  else m <- NULL
  
  if ("group" %in% argnames) {
    groupvar <- dots[["group"]] %ORifNULL% 
        stop("This PPC requires the 'group' argument.", call. = FALSE)
    dots[["group"]] <- .ppc_groupvar(object, groupvar, m = m)
  }
  if ("x" %in% argnames) {
    xvar <- dots[["x"]]  
    if (!is.null(xvar)) {
      dots[["x"]] <- .ppc_xvar(object, xvar, m = m)
    } else {
      if (funname %in% c("ppc_intervals", "ppc_ribbon")) {
        message("'x' not specified in '...'. Using x=1:length(y).")
        dots[["x"]] <- seq_along(y)
      } else {
        stop("This PPC requires the 'x' argument.", call. = FALSE)
      }
    }
  }
  
  if ("psis_object" %in% argnames && is.null(dots[["psis_object"]])) {
    dots[["psis_object"]] <- psis.stanreg(object)
  } else if ("lw" %in% argnames && is.null(dots[["lw"]])) {
    # for LOO predictive checks
    dots[["lw"]] <- weights(psis.stanreg(object))
  }
  
  return(dots)
}

# set default nreps value based on plot
.set_nreps <- function(nreps = NULL, fun = character()) {
  fun <- sub("ppc_", "", fun)
  switch(fun,
    # DISTRIBUTIONS
    "dens_overlay" = nreps %ORifNULL% 50,
    "ecdf_overlay" = nreps %ORifNULL% 50,
    "hist" = nreps %ORifNULL% 8,
    "dens" = nreps %ORifNULL% 8,
    "boxplot" = nreps %ORifNULL% 8,
    "freqpoly" = nreps %ORifNULL% 8,
    "freqpoly_grouped" = nreps %ORifNULL% 3,
    "violin_grouped" = nreps, # NULL ok
    
    # PREDICTIVE ERRORS
    "error_binned" = nreps %ORifNULL% 3,
    "error_hist" = nreps %ORifNULL% 3,
    "error_hist_grouped" = nreps %ORifNULL% 3,
    "error_scatter" = nreps %ORifNULL% 3,
    "error_scatter_avg" = nreps, # NULL ok
    "error_scatter_avg_vs_x" = nreps, # NULL ok
    
    # SCATTERPLOTS
    "scatter" = nreps %ORifNULL% 3, 
    "scatter_avg" = nreps, # NULL ok
    "scatter_avg_grouped" = nreps, # NULL ok
    
    # TEST-STATISTICS
    "stat" = .ignore_nreps(nreps),
    "stat_2d" = .ignore_nreps(nreps),
    "stat_grouped" = .ignore_nreps(nreps),
    "stat_freqpoly" = .ignore_nreps(nreps),
    "stat_freqpoly_grouped" = .ignore_nreps(nreps),
    
    # INTERVALS
    "intervals" = .ignore_nreps(nreps),
    "intervals_grouped" = .ignore_nreps(nreps),
    "ribbon" = .ignore_nreps(nreps),
    "ribbon_grouped" = .ignore_nreps(nreps), 
    
    # DISCRETE ONLY
    "rootogram" = nreps, # NULL ok
    "bars" = nreps, # NULL ok
    "bars_grouped" = nreps, # NULL ok
    
    # LOO PLOTS
    "loo_pit" = .ignore_nreps(nreps),
    "loo_pit_overlay" = .ignore_nreps(nreps),
    "loo_pit_qq" = .ignore_nreps(nreps),
    "loo_intervals" = .ignore_nreps(nreps),
    "loo_ribbon" = .ignore_nreps(nreps),
    
    # otherwise function not found
    stop(
      "Plotting function not supported. ",
      "(If the plotting function is included in the output from ", 
      "bayesplot::available_ppc() then it should be available via pp_check ",
      "and this error is probably a bug.)"
    )
  )
}
.ignore_nreps <- function(nreps) {
  if (!is.null(nreps))
    warning("'nreps' is ignored for this PPC", call. = FALSE)
  return(NULL)
}

# convert a character matrix (returned by posterior_predict for ordinal models) to a 
# numeric matrix
# 
# @param yrep character matrix
polr_yrep_to_numeric <- function(yrep) {
  apply(yrep, 2L, function(x) as.integer(as.factor(x)))
}
