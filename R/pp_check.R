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
#' An interface to the \pkg{\link{ppcheck}} package for \pkg{rstanarm} models, 
#' providing various plots comparing the observed outcome variable \eqn{y} to
#' simulated datasets \eqn{y^{rep}}{yrep} from the posterior predictive
#' distribution.
#' 
#' @export
#' @method pp_check stanreg
#' @templateVar bdaRef (Ch. 6)
#' @templateVar stanregArg object
#' @template reference-bda
#' @template args-stanreg-object
#' @param check The type of plot (possibly abbreviated) to show. One of 
#'   \code{"distributions"}, \code{"residuals"}, \code{"scatter"}, 
#'   \code{"test"}. See Details for descriptions.
#' @param nreps The number of \eqn{y^{rep}}{yrep} datasets to generate from the
#'   posterior predictive distribution (\code{\link{posterior_predict}}) and
#'   show in the plots. The default is \code{nreps=3} for
#'   \code{check="residuals"} and \code{nreps=8} for
#'   \code{check="distributions"}. If \code{check="test"}, \code{nreps} is
#'   ignored and the number of simulated datasets is the number of post-warmup
#'   draws from the posterior distribution. If \code{check="scatter"},
#'   \code{nreps} is not ignored but defaults to the number of post-warmup
#'   draws.
#' @param seed An optional \code{\link[=set.seed]{seed}} to pass to 
#'   \code{\link{posterior_predict}}.
#' @param overlay For \code{check="distributions"} only, should distributions be
#'   plotted as density estimates overlaid in a single plot (\code{TRUE}, the 
#'   default) or as separate histograms (\code{FALSE})?
#' @param test For \code{check="test"} only, a character vector (of length 1 or 
#'   2) naming a single function or a pair of functions. The function(s) should 
#'   take a vector input and return a scalar test statistic. See Details and
#'   Examples.
#' @param ... Currently unused.
#' 
#' @return A ggplot object that can be further customized using the
#'   \pkg{ggplot2} package.
#'   
#' @details Descriptions of the plots corresponding to the different values of 
#' \code{check}:
#' \describe{
#'  \item{\code{distributions}}{The distributions of \eqn{y} and \code{nreps} 
#'  simulated \eqn{y^{rep}}{yrep} datasets.} 
#'  \item{\code{residuals}}{The distributions of residuals computed from \eqn{y}
#'  and each of \code{nreps} simulated datasets. For binomial data, binned 
#'  residual plots are generated (similar to \code{\link[arm]{binnedplot}}).}
#'  \item{\code{scatter}}{If \code{nreps} is \code{NULL} then \eqn{y} is plotted
#'  against the average values of \eqn{y^{rep}}{yrep}, i.e., the points
#'  \eqn{(y_n, \bar{y}^{rep}_n),\, n = 1, \dots, N}{(y_n, mean(yrep_n)), n =
#'  1,...,N}, where each \eqn{y^{rep}_n}{yrep_n} is a vector of length equal to
#'  the number of posterior draws. If \code{nreps} is a (preferably small)
#'  integer, then only \code{nreps} \eqn{y^{rep}}{yrep} datasets are simulated
#'  and they are each plotted separately against \eqn{y}.}
#'  \item{\code{test}}{The distribution of a single test statistic
#'  \eqn{{T(y^{rep})}}{T(yrep)} or a pair of test statistics over the
#'  \code{nreps} simulated datasets. If the \code{test} argument specifies only
#'  one function then the resulting plot is a histogram of
#'  \eqn{{T(y^{rep})}}{T(yrep)} and the value of the test statistic in the 
#'  observed data, \eqn{T(y)}, is shown in the plot as a vertical line. If two 
#'  functions are specified then the plot is a scatterplot and \eqn{T(y)} is 
#'  shown as a large point.}
#' }
#' 
#' @note For binomial data, plots of \eqn{y} and \eqn{y^{rep}}{yrep} show the
#'   proportion of 'successes' rather than the raw count.
#' 
#' @seealso \code{\link{posterior_predict}} for drawing from the posterior 
#'   predictive distribution. Examples of posterior predictive checks can also
#'   be found in the \pkg{rstanarm} vignettes and demos.
#' 
#' @examples 
#' if (!exists("example_model")) example(example_model)
#' # Compare distribution of y to distributions of yrep
#' (pp_dist <- pp_check(example_model, check = "dist", overlay = TRUE))
#' pp_dist + 
#'  ggplot2::scale_color_manual(values = c("red", "black")) + # change colors
#'  ggplot2::scale_size_manual(values = c(0.5, 3)) + # change line sizes 
#'  ggplot2::scale_fill_manual(values = c(NA, NA)) # remove fill
#'
#' # Check residuals
#' pp_check(example_model, check = "resid", nreps = 6)
#'
#' # Check histograms of test statistics
#' test_mean <- pp_check(example_model, check = "test", test = "mean")
#' test_sd <- pp_check(example_model, check = "test", test = "sd")
#' gridExtra::grid.arrange(test_mean, test_sd, ncol = 2)
#' 
#' # Scatterplot of two test statistics
#' pp_check(example_model, check = "test", test = c("mean", "sd"))
#' 
#' \dontrun{
#' # Scatterplots of y vs. yrep
#' fit <- stan_glm(mpg ~ wt, data = mtcars)
#' pp_check(fit, check = "scatter") # y vs. average yrep
#' pp_check(fit, check = "scatter", nreps = 3) # y vs. a few different yrep datasets 
#' 
#' 
#' # Defining a function to compute test statistic 
#' roaches$roach100 <- roaches$roach1 / 100
#' fit_pois <- stan_glm(y ~ treatment + roach100 + senior, offset = log(exposure2), 
#'                      family = "poisson", data = roaches)
#' fit_nb <- update(fit_pois, family = "neg_binomial_2")
#' 
#' prop0 <- function(y) mean(y == 0) # function to compute proportion of zeros
#' pp_check(fit_pois, check = "test", test = "prop0") # looks bad 
#' pp_check(fit_nb, check = "test", test = "prop0")   # much better
#' }
#' 
#' @importFrom ggplot2 ggplot aes_string xlab %+replace% theme
#' 
pp_check.stanreg <- function(object, check = "distributions", nreps = NULL, 
                     seed = NULL, overlay = TRUE, test = "mean", ...) {
  if (used.optimizing(object)) 
    STOP_not_optimizing("pp_check")
  
  is_binomial <- if (is(object, "polr") && !is_scobit(object))
    FALSE else is.binomial(family(object)$family)
  
  checks <- c("distributions", "residuals", "scatter", "test")
  ppc_fun <- switch(
    match.arg(arg = check, choices = checks),
    'distributions' = if (overlay) 
      "ppc_dens_overlay" else "ppc_hist",
    'residuals' = if (is_binomial) 
      "ppc_resid_binned" else "ppc_resid",
    'test' = if (length(test) > 1) 
      "ppc_stat_2d" else "ppc_stat",
    'scatter' = if (is.null(nreps))  
      "ppc_scatter_average" else "ppc_scatter_multiple"
  )
  sel <- c("ppc_stat", "ppc_stat_2d", "ppc_scatter_average", "ppc_scatter_multiple")
  if (is.null(nreps) && !ppc_fun %in%  sel)
    nreps <- ifelse(ppc_fun == "ppc_dens_overlay", 50, ifelse(ppc_fun == "ppc_hist", 8, 3))
  if (!is.null(nreps) && check == "test") {
    warning("'nreps' is ignored if check='test'", call. = FALSE)
    nreps <- NULL
  }

  y <- get_y(object)
  if (ppc_fun == "ppc_resid_binned") {
    yrep <- posterior_linpred(object, transform = TRUE)
  } else {
    yrep <- posterior_predict(object, draws = nreps, seed = seed)
  }
  
  if (is_binomial) {
    if (NCOL(y) == 2L) {
      trials <- rowSums(y)
      y <- y[, 1L] / trials
      if (!binned_resid)
        yrep <- sweep(yrep, 2L, trials, "/")
    } else if (is.factor(y))
      y <- fac2bin(y)
  }
  if (is(object, "polr")) {
    y <- as.integer(y)
    yrep <- apply(yrep, 2L, function(x) as.integer(as.factor(x)))
  }
  
  ppc_args <- nlist(y, yrep)
  if (ppc_fun == "ppc_resid_binned") {
    ppc_args$yrep <- yrep[1:nreps, , drop = FALSE]
    names(ppc_args)[2] <- "Ey"
  }
  if (ppc_fun %in% c("ppc_stat", "ppc_stat_2d"))
    ppc_args$stat <- test
  
  do.call(ppc_fun, ppc_args)
}
