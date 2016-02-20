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
#' For models fit using MCMC or one of the variational approximations, there are
#' a variety of plots that can be generated. For models fit using optimization,
#' the regression coefficients and standard errors are passed to
#' \code{\link[arm]{coefplot}} (\pkg{arm}).
#' 
#' @method plot stanreg
#' @export
#' @templateVar stanregArg x
#' @template args-stanreg-object
#' @template args-pars
#' @template args-regex-pars
#' @param plotfun A character string naming the plotting function to apply to 
#'   the stanreg object. See \code{\link{plots}} for the names and descriptions.
#'   Also see the Examples section below. \code{plotfun} can be either the full 
#'   name of the plotting function (e.g. \code{"stan_hist"}) or can be 
#'   abbreviated to the part of the name following the underscore (e.g. 
#'   \code{"hist"}). The default plot shows intervals and point estimates for 
#'   the coefficients. Note: \code{plotfun} should not be specified for models
#'   fit using \code{algorithm="optimizing"} as there is currently only one
#'   plotting function for these models.
#' @param ... Additional arguments to pass to \code{plotfun} (see
#'   \code{\link{plots}}) or, for models fit using
#'   \code{algorithm="optimizing"}, \code{\link[arm]{coefplot}}.
#'
#' @return In most cases, a ggplot object (or several) that can be further 
#'   customized using the \pkg{ggplot2} package. The exception is for models fit
#'   using \code{"optimizing"} as the estimation algorithm, in which case a plot
#'   is produced but nothing is returned.
#'
#' @seealso \code{\link{plots}} for details on the individual plotting
#'   functions.
#'   
#' @examples
#' # Use rstanarm example model
#' fit <- example_model
#' 
#' # Intervals and point estimates
#' plot(fit) + 
#' ggplot2::ggtitle("Posterior medians \n with 80% and 95% credible intervals")
#' plot(fit, pars = "size", regex_pars = "period", 
#'      ci_level = 0.95, outer_level = 1, show_density = TRUE)
#' 
#' # Traceplot
#' # note: rstanarm doesn't store the warmup draws by default 
#' (trace <- plot(fit, "trace", pars = "(Intercept)"))
#' trace + ggplot2::scale_color_discrete()
#' trace + ggplot2::scale_color_manual(values = c("maroon", "skyblue2"))
#' 
#' # Distributions 
#' plot_title <- ggplot2::ggtitle("Posterior Distributions")
#' plot(fit, "hist", fill = "skyblue", regex_pars = "period") + plot_title
#' plot(fit, "dens", pars = "(Intercept)", regex_pars = "period", 
#'      separate_chains = TRUE, alpha = 1/3) + plot_title
#' 
#' # Scatterplot
#' plot(fit, plotfun = "scat", pars = paste0("period", 2:3))
#' plot(fit, plotfun = "scat", pars = c("(Intercept)", "size"), 
#'      color = "black", size = 5, alpha = 0.2)
#' 
#' # Some diagnostics
#' plot(fit, "rhat")
#' plot(fit, "ess")
#' 
#' # Using regex_pars
#' plot(fit, regex_pars = "period")
#' plot(fit, regex_pars = "herd:1")
#' plot(fit, regex_pars = "herd:1\\]")
#' plot(fit, regex_pars = "herd:[279]")
#' plot(fit, regex_pars = "herd:[279]|period2")
#' plot(fit, regex_pars = c("herd:[279]", "period2"))
#' 
#' # For graphical posterior predictive checks see 
#' # help("pp_check", package = "rstanarm")
#' 
#' @importFrom rstan stan_plot stan_trace stan_scat stan_hist stan_dens stan_ac
#'   stan_diag stan_rhat stan_ess stan_mcse stan_par quietgg
#' 
plot.stanreg <- function(x, plotfun = NULL, pars = NULL, 
                         regex_pars = NULL, ...) {
  args <- set_plotting_args(x, pars, regex_pars, ...)
  fun <- set_plotting_fun(x, plotfun)
  do.call(fun, args)
}


# Check for valid parameters
# @param x stanreg object
# @param pars user specified character vector
check_plotting_pars <- function(x, pars) {
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
set_plotting_args <- function(x, pars = NULL, regex_pars = NULL, ...) {
  args <- list(x, ...)
  pars <- collect_pars(x, pars, regex_pars)
  if (!is.null(pars)) 
    args$pars <- check_plotting_pars(x, pars)
  return(args)
}

# Select the correct plotting function
# @param x stanreg object
# @param plotfun user specified plotfun argument (can be missing)
set_plotting_fun <- function(x, plotfun = NULL) {
  .plotters <- function(x) paste0("stan_", x)
  
  if (used.optimizing(x)) {
    if (!is.null(plotfun)) {
      stop("'plotfun' should not be specified for models fit using ",
           "algorithm='optimizing'.", call. = FALSE)
    } else {
      return("stan_plot_opt")
    }
  } else if (is.null(plotfun)) {
    plotfun <- "stan_plot"
  }
  
  samp_only <- c("ac", "diag", "rhat", "ess", "mcse", "par")
  plotters <- .plotters(c("plot", "trace", "scat", "hist", "dens", samp_only))
  funname <- grep(paste0(plotfun, "$"), plotters, value = TRUE)
  if (used.variational(x) && funname %in% .plotters(samp_only))
    STOP_sampling_only(funname)
  fun <- try(getExportedValue("rstan", funname), silent = TRUE)
  if (inherits(fun, "try-error")) 
    stop("Plotting function not found. See ?rstanarm::plots for valid names.", 
         call. = FALSE)
  
  return(fun)
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
#' pairs(example_model, pars = c("(Intercept)", "log-posterior"))
#' 
pairs.stanreg <- function(x, ...) {
  if (!used.sampling(x)) 
    STOP_sampling_only("pairs")
  requireNamespace("rstan")
  requireNamespace("KernSmooth")
  pairs(x$stanfit, ...)
}

#' Plots for rstanarm models
#' 
#' All \pkg{rstanarm} models fit using \code{algorithm='sampling'},
#' \code{"meanfield"}, or \code{"fullrank"} are compatible with a variety of
#' plotting functions from the \pkg{rstan} package. Each function returns at
#' least one \code{\link[ggplot2]{ggplot}} object that can be customized further
#' using the \pkg{ggplot2} package. The plotting functions described here can be
#' called using the \code{\link[=plot.stanreg]{plot method}} for stanreg objects
#' without loading the \pkg{rstan} package. For example usage see
#' \code{\link{plot.stanreg}}.
#' 
#' @name plots
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
#' @examples
#' # See examples at help("plot.stanreg", package = "rstanarm")
NULL
