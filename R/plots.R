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
#' @importFrom ggplot2 ggplot aes_string xlab %+replace% theme
#'
plot.stanreg <- function(x, plotfun = "intervals", pars = NULL,
                         regex_pars = NULL, ...) {
  fun <- set_plotting_fun(plotfun)
  args <- set_plotting_args(x, pars, regex_pars, ..., plotfun = plotfun)
  do.call(fun, args)
}



# internal ----------------------------------------------------------------

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
      grepl("_rhat|_neff|_nuts_", plotfun))
    STOP_sampling_only(plotfun)
}


#' Pairs method for stanreg objects
#' 
#' This is essentially the same as \code{\link[rstan]{pairs.stanfit}} but with a
#' few tweaks for compatibility with \pkg{rstanarm} models. \strong{Be careful
#' not to specify too many parameters to include or the plot will be both hard
#' to read and slow to render.}
#'
#' @method pairs stanreg
#' @export
#' @templateVar stanregArg x
#' @template args-stanreg-object
#' @template args-regex-pars
#' @param pars An optional character vetor of parameter names. All parameters 
#'   are included by default, but for models with more than just a few 
#'   parameters it may be far too many to visualize on a small computer screen 
#'   and also may require substantial computing time.
#' 
#' @inheritParams rstan::pairs.stanfit
#' @param ...,labels,panel,lower.panel,upper.panel,diag.panel Same as in 
#'   \code{\link[graphics]{pairs}} syntactically but see the \strong{Details}
#'   section for different default arguments.
#' @param text.panel,label.pos,cex.labels,font.labels,row1attop,gap Same as in 
#'   \code{\link[graphics]{pairs.default}}.
#' @param log Same as in \code{\link[graphics]{pairs.default}} (but additionally
#'   accepts \code{log=TRUE}), which makes it possible to utilize logarithmic
#'   axes. See the \strong{Details} section.
#' 
#' @details This method differs from the default \code{\link{pairs}} method in 
#'   the following ways. If unspecified, the \code{\link{smoothScatter}} 
#'   function is used for the off-diagonal plots, rather than 
#'   \code{\link{points}}, since the former is more appropriate for visualizing 
#'   thousands of draws from a posterior distribution. Also, if unspecified, 
#'   histograms of the marginal distribution of each quantity are placed on the 
#'   diagonal of the plot, after pooling the chains.
#'   
#'   The draws from the warmup phase are always discarded before plotting.
#'   
#'   By default, the lower (upper) triangle of the plot contains draws with 
#'   below (above) median acceptance probability. Also, if \code{condition} is 
#'   not \code{"divergent__"}, red points will be superimposed onto the smoothed
#'   density plots indicating which (if any) iterations encountered a divergent 
#'   transition. Otherwise, yellow points indicate a transition that hit the 
#'   maximum treedepth rather than terminated its evolution normally.
#'   
#'   You may very well want to specify the \code{log} argument for non-negative 
#'   parameters. Of the various allowed specifications of \code{log}, it is 
#'   probably easiest to specify \code{log = TRUE}, which will utilize 
#'   logarithmic axes for all non-negative quantities.
#'   
#' @examples
#' if (!exists("example_model")) example(example_model)
#' pairs(example_model, pars = c("(Intercept)", "log-posterior"))
#' 
#' \donttest{
#' pairs(example_model, regex_pars = "herd:[279]")
#' }
#'
#' @importFrom graphics hist pairs par points rect smoothScatter text
#' @importFrom rstan get_logposterior
#' 
pairs.stanreg <-
  function(x,
           pars = NULL,
           regex_pars = NULL,
           condition = "accept_stat__",
           ...,
           labels = NULL,
           panel = NULL,
           lower.panel = NULL,
           upper.panel = NULL,
           diag.panel = NULL,
           text.panel = NULL,
           label.pos = 0.5 + 1 / 3,
           cex.labels = NULL,
           font.labels = 1,
           row1attop = TRUE,
           gap = 1,
           log = "") {
    
    if (!used.sampling(x))
      STOP_sampling_only("pairs")
    
    requireNamespace("rstan")
    requireNamespace("KernSmooth")
    
    arr <- as.array.stanreg(x, pars = pars, regex_pars = regex_pars)
    if (is.null(pars) && is.null(regex_pars)) {
      # include log-posterior by default
      lp_arr <- as.array.stanreg(x, pars = "log-posterior")
      dd <- dim(arr)
      dn <- dimnames(arr)
      dd[3] <- dd[3] + 1
      dn$parameters <- c(dn$parameters, "log-posterior")
      tmp <- array(NA, dim = dd, dimnames = dn)
      tmp[,, 1:(dd[3] - 1)] <- arr
      tmp[,, dd[3]] <- lp_arr
      arr <- tmp
    }
    arr <- round(arr, digits = 12)
    x <- x$stanfit
    
    gsp <- rstan::get_sampler_params(x, inc_warmup = FALSE)
    # if ("energy__" %in% colnames(gsp[[1]])) {
    #   dims <- dim(arr)
    #   dims[3] <- dims[3] + 1L
    #   nms <- dimnames(arr)
    #   nms$parameters <- c(nms$parameters, "energy__")
    #   E <- sapply(gsp, FUN = function(y) y[, "energy__"])
    #   arr <- array(c(c(arr), c(E)), dim = dims, dimnames = nms)
    # }
    
    sims <- nrow(arr)
    chains <- ncol(arr)
    varying <- apply(arr, 3, FUN = function(y) length(unique(c(y))) > 1)
    if (any(!varying)) {
      message(
        "The following parameters were dropped because they are constant:\n",
        paste(names(varying)[!varying], collapse = ", ")
      )
      arr <- arr[, , varying, drop = FALSE]
    }
    
    dupes <- duplicated(arr, MARGIN = 3)
    if (any(dupes)) {
      message(
        "The following parameters were dropped because they are duplicative:\n",
        paste(dimnames(arr)[[3]][dupes], collapse = ", ")
      )
      arr <- arr[, , !dupes, drop = FALSE]
    }
    
    tmp <- c(sapply(gsp, FUN = function(y) y[, "divergent__"]))
    divergent__ <- matrix(tmp, nrow = sims * chains, ncol = dim(arr)[3])
    
    max_td <- x@stan_args[[1]]$control
    if (is.null(max_td)) {
      max_td <- 10
    } else {
      max_td <- max_td$max_treedepth
      if (is.null(max_td))
        max_td <- 10
    }
    tmp <- c(sapply(gsp, FUN = function(y) y[, "treedepth__"] > max_td))
    hit <- matrix(tmp, nrow = sims * chains, ncol = dim(arr)[3])
    
    if (is.list(condition)) {
      if (length(condition) != 2)
        stop("If a list, 'condition' must be of length 2.")
      arr <- arr[, c(condition[[1]], condition[[2]]), , drop = FALSE]
      k <- length(condition[[1]])
      mark <- c(rep(TRUE, sims * k), rep(FALSE, sims * length(condition[[2]])))
      
    } else if (is.logical(condition)) {
      
      stopifnot(length(condition) == (sims * chains))
      mark <- !condition
      
    } else if (is.character(condition)) {
      
      condition <- match.arg(condition, several.ok = FALSE, 
                             choices = c("accept_stat__", "stepsize__",
                                         "treedepth__", "n_leapfrog__",
                                         "divergent__", "energy__", "lp__"))
      if (condition == "lp__") {
        mark <- simplify2array(get_logposterior(x, inc_warmup = FALSE))
      } else {
        mark <- sapply(gsp, FUN = function(y) y[, condition])
      }
      if (condition == "divergent__") {
        mark <- as.logical(mark)
      } else {
        mark <- c(mark) >= median(mark)
      }
      if (length(unique(mark)) == 1)
        stop(condition, " is constant so it cannot be used as a condition.")
      
    } else if (!is.null(condition)) {
      
      if (all(condition == as.integer(condition))) {
        arr <- arr[, condition, , drop = FALSE]
        k <- ncol(arr) %/% 2
        mark <- c(rep(FALSE, sims * k), rep(TRUE, sims * (chains - k)))
      } else if (condition > 0 && condition < 1) {
        mark <- rep(1:sims > (condition * sims), times = chains)
      } else {
        stop("If numeric, 'condition' must be an integer (vector) ",
             "or a number between 0 and 1 (exclusive).")
      }
      
    } else {
      
      k <- ncol(arr) %/% 2
      mark <- c(rep(FALSE, sims * k), rep(TRUE, sims * (chains - k)))
    }
    
    x <- apply(arr, MARGIN = "parameters", FUN = function(y) y)
    nc <- ncol(x)
    
    if (isTRUE(log)) {
      xl <- apply(x >= 0, 2, FUN = all)
      if ("log-posterior" %in% names(xl))
        xl["log-posterior"] <- FALSE
    } else if (is.numeric(log)) {
      xl <- log
    } else {
      xl <- grepl("x", log)
    }
    
    if (is.numeric(xl) || any(xl)) {
      x[, xl] <- log(x[, xl])
      colnames(x)[xl] <- paste("log", colnames(x)[xl], sep = "-")
    }
    if (is.null(lower.panel)) {
      if (!is.null(panel)) {
        lower.panel <- panel
      } else {
        lower.panel <- function(x, y, ...) {
          dots <- list(...)
          dots$x <- x[!mark]
          dots$y <- y[!mark]
          if (is.null(mc$nrpoints) &&
              !identical(condition, "divergent__")) {
            dots$nrpoints <- Inf
            dots$col <- ifelse(divergent__[!mark] == 1,
                               "red",
                               ifelse(hit[!mark] == 1, "yellow", NA_character_))
          }
          dots$add <- TRUE
          do.call(smoothScatter, args = dots)
        }
      }
    }
    if (is.null(upper.panel)) {
      if (!is.null(panel)) {
        upper.panel <- panel
      } else {
        upper.panel <- function(x, y, ...) {
          dots <- list(...)
          dots$x <- x[mark]
          dots$y <- y[mark]
          if (is.null(mc$nrpoints) &&
              !identical(condition, "divergent__")) {
            dots$nrpoints <- Inf
            dots$col <- ifelse(divergent__[mark] == 1,
                               "red",
                               ifelse(hit[mark] == 1, "yellow", NA_character_))
          }
          dots$add <- TRUE
          do.call(smoothScatter, args = dots)
        }
      }
    }
    if (is.null(diag.panel))
      diag.panel <- function(x, ...) {
        usr <- par("usr")
        on.exit(par(usr))
        par(usr = c(usr[1:2], 0, 1.5))
        h <- hist(x, plot = FALSE)
        breaks <- h$breaks
        y <- h$counts
        y <- y / max(y)
        nB <- length(breaks)
        rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
      }
    if (is.null(panel))
      panel <- points
    
    if (is.null(text.panel)) {
      textPanel <- function(x = 0.5, y = 0.5, txt, cex, font) {
        text(x, y, txt, cex = cex, font = font)
      }
    } else {
      textPanel <- text.panel
    }
    if (is.null(labels))
      labels <- colnames(x)
    
    mc <- match.call(expand.dots = FALSE)
    mc[1] <- call("pairs")
    mc$x <- x
    mc$labels <- labels
    mc$panel <- panel
    mc$lower.panel <- lower.panel
    mc$upper.panel <- upper.panel
    mc$diag.panel <- diag.panel
    mc$text.panel <- textPanel
    mc$log <- ""
    mc$condition <- NULL
    mc$pars <- NULL
    mc$regex_pars <- NULL
    mc$include <- NULL
    eval(mc)
  }
