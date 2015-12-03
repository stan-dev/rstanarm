# @param x stanreg object
# @param pars user specified character vector
.check_plotting_pars <- function(x, pars) {
  if (used.sampling(x)) {
    sim <- x$stanfit@sim
    allpars <- c(sim$pars_oi, sim$fnames_oi)    
  }
  else if (used.optimizing(x)) {
    allpars <- c("alpha", "beta", rownames(x$stan_summary))
  }
  else stop("Not yet enabled for advi")

  m <- which(match(pars, allpars, nomatch = 0) == 0)
  if (length(m) > 0) stop("no parameter ", paste(pars[m], collapse = ', ')) 
  return(unique(pars))
}

# @param x stanreg object
# @param regex_pars character vector of patterns
.grep_for_pars <- function(x, regex_pars) {
  stopifnot(is.character(regex_pars))
  out <- unlist(lapply(seq_along(regex_pars), function(j) {
    grep(regex_pars[j], rownames(x$stan_summary), value = TRUE) 
  }))
  if (length(out)) return(out) 
  else stop("No matches for regex_pars.")
}


#' Plot method for stanreg objects
#' 
#' For models fit using \code{algorithm="sampling"} there are a variety of
#' \code{\link{plots}} that can be generated. For models fit with
#' \code{algorithm="optimizing"} the regression coefficients and standard errors
#' are passed to \code{\link[arm]{coefplot}} (\pkg{arm}).
#' 
#' @method plot stanreg
#' @export
#' @templateVar stanregArg x
#' @template args-stanreg-object
#' @param plotfun A character string (possibly abbreviated) naming the plotting
#'   function to apply to the stanreg object. See \code{\link{plots}} for the
#'   names and descriptions. The default plot shows intervals and point
#'   estimates for the coefficients. (\strong{Note:} for models fit using
#'   \code{algorithm="optimizing"} the \code{plotfun} argument is ignored as
#'   there is currently only one plotting function for these models.)
#' @param pars An optional character vector of parameter names.
#' @param regex_pars An optional character vector of \link[=grep]{regular
#'   expressions} to use for parameter selection. Can be used in place of
#'   \code{pars} or in addition to \code{pars}. \code{regex_pars}
#'   is ignored for models fit using \code{algorithm="optimizing"}.
#' @param ... Additional arguments to pass to \code{plotfun} (see
#'   \code{\link{plots}}) or, for models fit using
#'   \code{algorithm="optimizing"}, \code{\link[arm]{coefplot}}.
#'
#' @return A ggplot object (or several) that can be further customized using the
#'   \pkg{ggplot2} package. (If \code{x$algorithm="optimizing"} a plot is
#'   produced but nothing is returned.)
#'
#' @seealso \code{\link{plots}} for details on the individual plotting
#'   functions.
#'   
#' @examples 
#' # See examples at help("plots", "rstanarm")
#' 
#' @importFrom rstan stan_plot stan_trace stan_scat stan_hist stan_dens stan_ac
#'   stan_diag stan_rhat stan_ess stan_mcse stan_par quietgg
#' 
plot.stanreg <- function(x, plotfun, pars, regex_pars, ...) {
  args <- list(x, ...)
  if (missing(plotfun)) plotfun <- "plot"
  if (missing(regex_pars)) check_regex <- FALSE
  else {
    check_regex <- !used.optimizing(x)
    if (!check_regex) {
      warning("'regex_pars' ignored for models with algorithm='optimizing'.")
    }
  }
  if (!missing(pars)) {
    pars[pars == "varying"] <- "b"
    if (check_regex) pars <- c(pars, .grep_for_pars(x, regex_pars))
    pars <- .check_plotting_pars(x, pars)
    args$pars <- pars
  }
  else if (check_regex) {
    args$pars <- .check_plotting_pars(x, .grep_for_pars(x, regex_pars))
  }
  
  if (used.optimizing(x)) fun <- "stan_plot_opt"
  else {
    plotters <- paste0("stan_", c("plot", "trace", "scat", "hist", "dens", "ac",
                                  "diag", "rhat", "ess", "mcse", "par"))
    funname <- grep(paste0(plotfun, "$"), plotters, value = TRUE)
    fun <- try(getExportedValue("rstan", funname), silent = TRUE)
    if (inherits(fun, "try-error")) 
      stop("Plotting function not found. See ?rstanarm::plots for valid names.", 
           call. = FALSE)
  }
  do.call(fun, args)
}

# function calling arm::coefplot (only used for models fit using optimization)
stan_plot_opt <- function(x, pars, varnames = NULL, ...) {
  if (!requireNamespace("arm", quietly = TRUE)) 
    stop("Please install the 'arm' package to use this feature.", 
         call. = FALSE)
  stopifnot(used.optimizing(x))
  coefs <- coef(x)
  sds <- se(x)
  nms <- varnames %ORifNULL% names(x$coefficients)
  if (!missing(pars)) {
    mark <- NA
    if ("alpha" %in% pars) mark <- c(mark, "(Intercept)")
    if ("beta" %in% pars) 
      mark <- c(mark, setdiff(names(x$coefficients), "(Intercept)"))
    mark <- c(mark, setdiff(pars, c("alpha", "beta")))
    mark <- mark[!is.na(mark)] 
    coefs <- coefs[mark]
    sds <- sds[mark]
    if (is.null(varnames)) nms <- nms[nms %in% mark]
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
  if (!used.sampling(x)) STOP_sampling_only("pairs")
  requireNamespace("rstan")
  requireNamespace("KernSmooth")
  pairs(x$stanfit, ...)
}

#' Plots
#' 
#' All models fit using \code{algorithm='sampling'} are compatible with a 
#' variety of plotting functions from the \pkg{rstan} package. Each function
#' returns at least one \code{\link[ggplot2]{ggplot}} object that can be
#' customized further using the \pkg{ggplot2} package. The plotting functions
#' described here can also be called using the \code{\link[=plot.stanreg]{plot}}
#' method for stanreg objects without loading the \pkg{rstan} package.
#' 
#' 
#' @name plots
#' 
#' @section Plotting functions:
#' \describe{
#' \item{Posterior intervals and point estimates}{\code{\link[rstan]{stan_plot}}}
#' \item{Traceplots}{\code{\link[rstan]{stan_trace}}}
#' \item{Histograms}{\code{\link[rstan]{stan_hist}}}
#' \item{Kernel density estimates}{\code{\link[rstan]{stan_dens}}}
#' \item{Scatterplots}{\code{\link[rstan]{stan_scat}}}
#' \item{Diagnostics for Hamiltonian Monte Carlo and the No-U-Turn Sampler}{\code{\link[rstan]{stan_diag}}}
#' \item{Rhat}{\code{\link[rstan]{stan_rhat}}}
#' \item{Ratio of effective sample size to total posterior sample size}{\code{\link[rstan]{stan_ess}}}
#' \item{Ratio of Monte Carlo standard error to posterior standard deviation}{\code{\link[rstan]{stan_mcse}}}
#' \item{Autocorrelation}{\code{\link[rstan]{stan_ac}}}
#' }
#' 
#' For graphical posterior predicive checking see \code{\link{ppcheck}}.
#' 
#' @seealso \code{\link{plot.stanreg}}, \code{\link{shinystan}}
#' @examples 
#' # Intervals and point estimates
#' plot(example_model, ci_level = 0.8)
#' plot(example_model, pars = "size", regex_pars = "period", show_density = TRUE)
#' 
#' # Traceplot
#' (trace <- plot(example_model, plotfun = "trace", pars = "(Intercept)"))
#' trace + ggplot2::scale_color_discrete()
#' 
#' # Distributions 
#' plot(example_model, "hist", fill = "skyblue") + ggplot2::ggtitle("Posterior Distributions")
#' plot(example_model, "dens", regex_pars = "period", 
#'      separate_chains = TRUE, alpha = 0.25)
#' 
#' # Scatterplot
#' plot(example_model, plotfun = "scat", pars = paste0("period", 2:3))
#' 
#' # Some diagnostics
#' plot(example_model, "rhat")
#' plot(example_model, "ess")
#' 
#' # Using regex_pars
#' plot(example_model, regex_pars = "herd:1")
#' plot(example_model, regex_pars = "herd:1\\]")
#' plot(example_model, regex_pars = "herd:[279]")
#' plot(example_model, regex_pars = "herd:[279]|period2")
#' plot(example_model, regex_pars = c("herd:[279]", "period2"))
#' 
#' # Posterior predictive checks (see ?ppcheck for examples)
NULL
