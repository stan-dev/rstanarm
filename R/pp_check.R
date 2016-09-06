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
#' Graphical posterior predictive checks
#' 
#' Interface to the posterior predictive checking functionality in the 
#' \pkg{\link{bayesplot}} package, providing various plots comparing the
#' observed outcome variable \eqn{y} to simulated datasets \eqn{y^{rep}}{yrep}
#' from the posterior predictive distribution.
#' 
#' @export
#' @export pp_check
#' @aliases pp_check
#' @method pp_check stanreg
#' @templateVar bdaRef (Ch. 6)
#' @templateVar stanregArg object
#' @template reference-bda
#' @template args-stanreg-object
#' @param check The type of plot to show. One of \code{"distributions"},
#'   \code{"residuals"}, \code{"scatterplots"}, \code{"test-statistics"},
#'   \code{"vs_x"} (can be abbreviated). See Details for descriptions.
#' @param nreps The number of \eqn{y^{rep}}{yrep} datasets to generate from the 
#'   posterior predictive distribution (\code{\link{posterior_predict}}) and 
#'   show in the plots. The default is \code{nreps=3} for 
#'   \code{check="residuals"} and \code{nreps=8} for 
#'   \code{check="distributions"}. If \code{check="test"} or
#'   \code{check="vs_x"}, \code{nreps} is ignored and the number of simulated
#'   datasets is the number of post-warmup draws from the posterior
#'   distribution. If \code{check="scatter"}, \code{nreps} is not ignored but
#'   defaults to the number of post-warmup draws.
#' @param seed An optional \code{\link[=set.seed]{seed}} to pass to 
#'   \code{\link{posterior_predict}}.
#' @param overlay For \code{check="distributions"} only, should distributions be
#'   plotted as density estimates overlaid in a single plot (\code{TRUE}, the 
#'   default) or as separate histograms (\code{FALSE})?
#' @param test For \code{check="test"} only, a character vector (of length 1 or 
#'   2) naming a single function or a pair of functions. The function(s) should 
#'   take a vector input and return a scalar test statistic. See Details and
#'   Examples.
#' @param group A string naming a grouping variable by which to stratify, or a
#'   factor providing the values of that variable. Specifying \code{group} via a
#'   string is only allowed if the variable is in the \code{\link{model.frame}}.
#'   Not available for all plots.
#' @param x For \code{check="vs_x"} only, either a string naming the variable to
#'   use as the \eqn{x}-variable, or a numeric vector providing the values of 
#'   that variable. Specifying \code{x} via a string is only allowed if
#'   the variable is in the \code{\link{model.frame}}.
#' @param ... Optionally, additonal arguments passed to the 
#'   \pkg{\link{bayesplot}} function called. These are described in the help 
#'   pages for \pkg{bayesplot}, links to which can be found in Details, below.
#' 
#' @return A ggplot object that can be further customized using the 
#'   \pkg{ggplot2} package.
#'   
#' @details Descriptions of the plots corresponding to the different values of 
#' \code{check}:
#' \describe{
#'  \item{\code{distributions}}{
#'    The distributions of \eqn{y} and \code{nreps} simulated
#'    \eqn{y^{rep}}{yrep} datasets.
#'    
#'    \pkg{\link{bayesplot}} reference:
#'    \code{\link[bayesplot]{PPC-distributions}}
#'  } 
#'  \item{\code{residuals}}{
#'    The distributions of residuals computed from \eqn{y} and each of
#'    \code{nreps} simulated datasets. For binomial data, binned residual plots
#'    are generated (similar to \code{binnedplot} in package \pkg{arm}).
#'    
#'    \pkg{\link{bayesplot}} reference:
#'    \code{\link[bayesplot]{PPC-residuals}}
#'  }
#'  \item{\code{scatterplots}}{
#'    If \code{nreps} is \code{NULL} then \eqn{y} is plotted against the average
#'    values of \eqn{y^{rep}}{yrep}, i.e., the points \eqn{(y_n, 
#'    \bar{y}^{rep}_n),\, n = 1, \dots, N}{(y_n, mean(yrep_n)), n = 1,...,N}, 
#'    where each \eqn{y^{rep}_n}{yrep_n} is a vector of length equal to the
#'    number of posterior draws. If \code{nreps} is a (preferably small)
#'    integer, then only \code{nreps} \eqn{y^{rep}}{yrep} datasets are simulated
#'    and they are each plotted separately against \eqn{y}.
#'    
#'    \pkg{\link{bayesplot}} reference:
#'    \code{\link[bayesplot]{PPC-scatterplots}}
#'  }
#'  \item{\code{test-statistics}}{
#'    The distribution of a single test statistic
#'    \eqn{{T(y^{rep})}}{T(yrep)} or a pair of test statistics over the
#'    \code{nreps} simulated datasets. If the \code{test} argument specifies only
#'    one function then the resulting plot is a histogram of
#'    \eqn{{T(y^{rep})}}{T(yrep)} and the value of the test statistic in the 
#'    observed data, \eqn{T(y)}, is shown in the plot as a vertical line. If two 
#'    functions are specified then the plot is a scatterplot and \eqn{T(y)} is 
#'    shown as a large point.
#'    
#'    \pkg{\link{bayesplot}} reference:
#'    \code{\link[bayesplot]{PPC-test-statistics}}
#'  }
#'  \item{\code{vs_x}}{
#'    Medians and central interval estimates of \eqn{y^{rep}}{yrep} by value of 
#'    an "\eqn{x}" variable, with \eqn{y} overlaid.
#'    
#'    \pkg{\link{bayesplot}} reference:
#'    \code{\link[bayesplot]{PPC-vs-x}}
#'  }
#' }
#' 
#' @note For binomial data, plots of \eqn{y} and \eqn{y^{rep}}{yrep} show the
#'   proportion of 'successes' rather than the raw count.
#' 
#' @seealso \code{\link{posterior_predict}} for drawing from the posterior 
#'   predictive distribution. Examples of posterior predictive checks can also 
#'   be found in the \pkg{rstanarm} vignettes and demos.
#'   
#'   \code{\link[bayesplot]{set_color_scheme}} to change the color scheme of the
#'   plots.
#' 
#' @examples 
#' if (!exists("example_model")) example(example_model)
#' 
#' # Compare distribution of y (dark color) to distributions of 
#' # yrep (light color)
#' (pp_dist <- pp_check(example_model, check = "dist"))
#'  
#' # Violin plot of yrep by level of 'herd' grouping variable with 
#' # y points overlaid
#' pp_check(example_model, check = "dist", group = "herd")
#'
#' # Check residuals (default is binned residual plot for binomial 
#' # models, histograms otherwise)
#' pp_check(example_model, check = "resid", nreps = 4)
#'
#' # Check histograms of test statistics
#' pp_check(example_model, check = "test", test = "mean")
#' pp_check(example_model, check = "test", test = "sd")
#' pp_check(example_model, check = "test", test = "mean", group = "herd")
#' 
#' # Scatterplot of two test statistics
#' pp_check(example_model, check = "test", test = c("mean", "sd"))
#' 
#' \donttest{
#' # Scatterplots of y vs. yrep
#' fit <- stan_glm(mpg ~ wt, data = mtcars)
#' pp_check(fit, check = "scatter") # y vs. average yrep
#' pp_check(fit, check = "scatter", nreps = 3) # y vs. a few different yrep datasets 
#' 
#' # yrep "ribbon" vs x (with y points overlaid)
#' pp_check(fit, check = "vs_x", x = mtcars$disp, y_style = "points")
#' pp_check(fit, check = "vs_x", x = "wt") + ggplot2::xlab("wt")
#' 
#' # Defining a function to compute test statistic 
#' roaches$roach100 <- roaches$roach1 / 100
#' fit_pois <- stan_glm(y ~ treatment + roach100 + senior, data = roaches,
#'                      offset = log(exposure2), family = "poisson")
#' fit_nb <- update(fit_pois, family = "neg_binomial_2")
#' 
#' prop0 <- function(y) mean(y == 0) # function to compute proportion of zeros
#' pp_check(fit_pois, check = "test", test = "prop0") # looks bad 
#' pp_check(fit_nb, check = "test", test = "prop0")   # much better
#' }
#' 
pp_check.stanreg <-
  function(object,
           check = "distributions",
           nreps = NULL,
           seed = NULL,
           overlay = TRUE,
           test = "mean",
           group = NULL,
           x = NULL,
           ...) {
    if (used.optimizing(object))
      STOP_not_optimizing("pp_check")
    
    valid_ppcs <- c("distributions", "residuals", "scatterplots", 
                    "test-statistics", "vs_x")
    plotfun <-
      ppc_fun(
        check = match.arg(check, choices = valid_ppcs),
        grouped = !is.null(group),
        nreps = nreps,
        ntests = length(test),
        overlay = isTRUE(overlay),
        binomial_model = is_binomial_ppc(object), 
        has_x = !is.null(x)
      )
    
    y_yrep <-
      ppc_y_and_yrep(
        object,
        seed = seed,
        nreps = set_nreps(nreps, fun = plotfun),
        binned_resid_plot = isTRUE(plotfun == "ppc_resid_binned")
      )
    
    plotargs <-
      ppc_args(
        y = y_yrep[["y"]],
        yrep = y_yrep[["yrep"]],
        group = set_group(object, group),
        x = set_x(object, x),
        fun = plotfun,
        test = test,
        ...
      )
    
    do.call(plotfun, plotargs)
  }

ppc_y_and_yrep <-
  function(object,
           nreps = NULL,
           seed = NULL,
           binned_resid_plot = FALSE) {
    y <- get_y(object)
    if (binned_resid_plot) {
      yrep <- posterior_linpred(object, transform = TRUE)
      yrep <- yrep[1:nreps, , drop = FALSE]
    } else {
      yrep <- posterior_predict(object, draws = nreps, seed = seed)
    }
    
    if (is_binomial_ppc(object)) {
      if (NCOL(y) == 2L) {
        trials <- rowSums(y)
        y <- y[, 1L] / trials
        if (!binned_resid_plot)
          yrep <- sweep(yrep, 2L, trials, "/")
      } else if (is.factor(y))
        y <- fac2bin(y)
    }
    if (is(object, "polr")) {
      y <- as.integer(y)
      yrep <- apply(yrep, 2L, function(x) as.integer(as.factor(x)))
    }
    
    nlist(y, yrep)
  }

ppc_args <-
  function(y, 
           yrep,
           group = NULL,
           fun = character(),
           test = NULL,
           x = NULL,
           ...) {
    args <- nlist(y, yrep, ...)
    if (!is.null(group))
      args$group <- group
    if (!is.null(x))
      args$x <- x
    if (fun == "ppc_resid_binned")
      names(args)[names(args) %in% "yrep"] <- "Ey"
    if (grepl("^ppc_stat", fun))
      args$stat <- test
    
    args
  }

set_x <- set_group <- function(object, group = NULL) {
  if (is.null(group) || !is.character(group))
    return(group)
  
  mf <- model.frame(object)
  vars <- colnames(mf)
  if (group %in% vars)
    return(mf[, group])
  
  stop("Variable '", group, "' not found in model frame. ")
}

is_binomial_ppc <- function(object) {
  if (is(object, "polr") && !is_scobit(object)) {
    FALSE
  } else {
    is.binomial(family(object)$family)
  }
}

ppc_fun <-
  function(check,
           grouped = FALSE,
           nreps = NULL,
           ntests = 1,
           overlay = TRUE,
           binomial_model = FALSE, 
           has_x = FALSE) {
    
    if (check == "distributions") {
      if (grouped) 
        return("ppc_violin_grouped")
      else if (overlay) 
        return("ppc_dens_overlay")
      else 
        return("ppc_hist")
    }
    
    if (check == "residuals") {
      if (grouped)
        warning("'group' is ignored for residuals plots.", call. = FALSE)
      if (binomial_model) 
        return("ppc_resid_binned")
      else 
        return("ppc_resid")  
    }
    
    if (check == "test-statistics") {
      if (ntests > 1) {
        if (grouped)
          warning("'group' is ignored if length(test) > 1.", call. = FALSE)
        return("ppc_stat_2d")
      }
      else if (grouped)
        return("ppc_stat_grouped")
      else
        return("ppc_stat")
    }
    
    if (check == "scatterplots") {
      if (!is.null(nreps)) {
        if (grouped)
          warning("'group' is ignored for scatterplots unless 'nreps' is NULL.", 
                  call. = FALSE)
        return("ppc_scatter" )
      }
      else if (grouped)
        return("ppc_scatter_avg_grouped")
      else
        return("ppc_scatter_avg")
    }
    
    if (check == "vs_x") {
      if (!has_x)
        stop("If 'check' is 'vs_x' then the 'x' argument must be specified.")
      if (grouped) 
        return("ppc_vs_x_grouped")
      else 
        return("ppc_vs_x")
    }
  }

set_nreps <- function(nreps = NULL, fun = character()) {
  fun <- sub("ppc_", "", fun)
  switch(fun,
    # DISTRIBUTIONS
    "dens_overlay" = nreps %ORifNULL% 50,
    "hist" = nreps %ORifNULL% 8,
    "violin_grouped" = nreps, # NULL ok
    
    # RESIDUALS
    "resid" = nreps %ORifNULL% 3,
    "resid_binned" = nreps %ORifNULL% 3,
    
    # SCATTERPLOTS
    "scatter" = nreps %ORifNULL% 3, 
    "scatter_avg" = nreps, # NULL ok
    "scatter_avg_grouped" = nreps, # NULL ok
    
    # TEST-STATISTICS
    "stat" = ignore_nreps(nreps),
    "stat_2d" = ignore_nreps(nreps),
    "stat_grouped" = ignore_nreps(nreps),
    
    # VS X
    "vs_x" = ignore_nreps(nreps),
    "vs_x_grouped" = ignore_nreps(nreps)
  )
}

ignore_nreps <- function(nreps, check = "test") {
  if (!is.null(nreps))
    warning("'nreps' is ignored if check=", check, call. = FALSE)
  return(NULL)
}
