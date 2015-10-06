#' Plot method for stanreg objects
#' 
#' For models fit using \code{algorithm="sampling"} there are a variety of
#' \code{\link{plots}} that can be generated. For models fit with
#' \code{algorithm="optimizing"} the regression coefficients and standard errors
#' are passed to \code{\link[arm]{coefplot}} (\pkg{arm}).
#' 
#' @method plot stanreg
#' @export
#' @inheritParams summary.stanreg
#' @param x A stanreg object returned by one of the \pkg{rstanarm} modeling
#'   functions.
#' @param fun A character string naming the plotting function to apply to the 
#'   stanreg object. See \code{\link{plots}} for the names and descriptions. The
#'   default plot shows intervals and point estimates for the coefficients. 
#'   (\strong{Note:} for models fit using \code{algorithm="optimizing"} the
#'   \code{fun} argument is ignored as there is only one plotting function for 
#'   these models.)
#' @param ... Additional arguments to pass to \code{fun} (see
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
#' # See help("plots", "rstanarm")
#' 
plot.stanreg <- function(x, fun = "stan_plot", pars, ...) {
  args <- list(x, ...)
  if (!missing(pars)) {
    pars[pars == "varying"] <- "b"
    args$pars <- pars
  }
  fun <- if (x$algorithm != "optimizing") match.fun(fun) else "stan_plot_opt"
  do.call(fun, args)
}

# function calling arm::coefplot (only used for models fit using optimization)
stan_plot_opt <- function(x, pars, varnames = NULL, ...) {
  if (!requireNamespace("arm", quietly = TRUE)) 
    stop("Please install the 'arm' package to use this feature")
  stopifnot(x$algorithm == "optimizing")
  nms <- varnames %ORifNULL% names(x$coefficients)
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
#' 
#' @examples 
#' \dontrun{
#' (fit <- stan_lm(mpg ~ wt + qsec + am, data = mtcars, prior = R2(0.75), seed = 12345))
#' pairs(fit, pars = c("(Intercept)", "log-posterior"))
#' } 
#' 
pairs.stanreg <- function(x, ...) {
  pairs(x$stanfit, ...)
}

#' Plots
#' 
#' All models fit using \code{algorithm='sampling'} are compatible with a 
#' variety of plotting functions. Each function returns at least one 
#' \code{\link[ggplot2]{ggplot}} object that can be customized further using the
#' \pkg{ggplot2} package. The plotting functions described here can also be
#' called using the \code{\link[=plot.stanreg]{plot}} method for stanreg
#' objects.
#' 
#' 
#' @name plots
#' 
#' @section Plotting functions:
#' \describe{
#' \item{Posterior predictive checks}{\code{\link{ppcheck}}}
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
#' @seealso \code{\link{plot.stanreg}}, \code{\link{shinystan}}
#' 
#' @examples 
#' \dontrun{
#' data("clouds", package = "HSAUR3")
#' f <- rainfall ~ seeding * (sne + cloudcover + prewetness + echomotion) + time
#' fit <- stan_lm(f, data = clouds, prior = R2(location = 0.25), 
#'                cores = 4, seed = 12345)
#'                
#' # stan_plot: posterior intervals and point estimates
#' stan_plot(fit, ci_level = 0.8)
#' stan_plot(fit, pars = c("prewetness", "echomotionstationary"), show_density = TRUE)               
#' 
#' # posterior predictive checks (see ?ppcheck for more details and examples)
#' ppcheck(fit, check = "distributions")
#' ppcheck(fit, check = "distributions", overlay = TRUE)
#' ppcheck(fit, check = "residuals")
#' ppcheck(fit, check = "test", test = sd)
#' 
#' # traceplot
#' (trace <- stan_trace(fit, pars = "(Intercept)"))
#' trace + scale_color_discrete()
#' 
#' # distributions 
#' stan_hist(fit, fill = "skyblue") + ggtitle("Example Plot")
#' stan_dens(fit, pars = c("sne", "cloudcover"), separate_chains = TRUE, alpha = 0.1)
#' 
#' # scatterplot
#' stan_scat(fit, pars = c("sne", "cloudcover"))
#' }
NULL