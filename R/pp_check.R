# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2015 Trustees of Columbia University
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

#' ppcheck (deprecated)
#' 
#' ppcheck is deprecated and will be removed in the near future. Use
#' \code{\link{pp_check}} instead.
#' 
#' @export
#' @keywords internal
#' @inheritParams pp_check
ppcheck <- function(object,
                    check = "distributions",
                    nreps = NULL, overlay = TRUE, test = "mean", ...) {
  .Deprecated("pp_check")
  mc <- match.call()
  mc[[1L]] <- as.name("pp_check")
  eval(mc, parent.frame())
}

#' Graphical posterior predictive checks
#' 
#' Various plots comparing the observed outcome variable \eqn{y} to simulated 
#' datasets \eqn{yrep} from the \link[=posterior_predict]{posterior predictive
#' distribution}.
#' 
#' @export
#' @templateVar bdaRef (Ch. 6)
#' @templateVar stanregArg object
#' @template reference-bda
#' @template args-stanreg-object
#' @param check The type of plot (possibly abbreviated) to show. One of 
#'   \code{"distributions"}, \code{"residuals"}, \code{"scatter"}, 
#'   \code{"test"}. See Details for descriptions.
#' @param nreps The number of \eqn{yrep} datasets to generate from the posterior
#'   predictive distribution and show in the plots. The default is 
#'   \code{nreps=3} for \code{check="residuals"} and \code{nreps=8} for
#'   \code{check="distributions"}. If \code{check="test"} then \code{nreps} is
#'   ignored and the number of simulated datasets is the number of post-warmup
#'   draws from the posterior distribution. If \code{check="scatter"},
#'   \code{nreps} is not ignored but defaults to the number of post-warmup
#'   draws.
#' @param overlay For \code{check="distributions"} only, should distributions be
#'   plotted as density estimates overlaid in a single plot (\code{TRUE}, the
#'   default) or as separate histograms (\code{FALSE})?
#' @param test For \code{check="test"} only, a character vector (of length 1 or 
#'   2) naming a single function or a pair of functions. The function(s) should 
#'   take a vector input and return a scalar test statistic. See Details.
#' @param ... Optional arguments to geoms to control features of the plots 
#'   (e.g. \code{binwidth} if the plot is a histogram).
#' 
#' @return A ggplot object that can be further customized using the
#'   \pkg{ggplot2} package.
#'   
#' @details Descriptions of the plots corresponding to the different values of 
#' \code{check}:
#' \describe{
#'  \item{\code{distributions}}{The distributions of \eqn{y} and \code{nreps} 
#'  simulated \eqn{yrep} datasets.} 
#'  \item{\code{residuals}}{The distributions of residuals computed from \eqn{y}
#'  and each of \code{nreps} simulated datasets. For binomial data, binned 
#'  residual plots are generated (similar to \code{\link[arm]{binnedplot}}).}
#'  \item{\code{scatter}}{If \code{nreps} is \code{NULL} then \eqn{y} is plotted
#'  against the average values of \eqn{yrep}, i.e., the points \eqn{(y_n, 
#'  \bar{yrep}_n),\, n = 1, \dots, N}{(y_n, mean(yrep_n)), n = 1,...,N},
#'  where each \eqn{yrep_n} is a vector of length equal to the number of
#'  posterior draws. If \code{nreps} is a (preferably small) integer, then only
#'  \code{nreps} \eqn{yrep} datasets are simulated and they are each plotted
#'  separately against \eqn{y}.}
#'  \item{\code{test}}{The distribution of a single test statistic \eqn{T(yrep)}
#'  or a pair of test statistics over the \code{nreps} simulated datasets. If 
#'  the \code{test} argument specifies only one function then the resulting plot
#'  is a histogram of \eqn{T(yrep)} and the value of the test statistic in the
#'  observed data, \eqn{T(y)}, is shown in the plot as a vertical line. If two
#'  functions are specified then the plot is a scatterplot and \eqn{T(y)} is
#'  shown as a large point.}
#' }
#' 
#' @note For binomial data, plots of \eqn{y} and \eqn{yrep} show the proportion
#'   of 'successes' rather than the raw count.
#' 
#' @seealso \code{\link{posterior_predict}} for drawing from the posterior 
#'   predictive distribution. Examples of posterior predictive checks can also
#'   be found in the \pkg{rstanarm} vignettes and demos.
#' 
#' @examples
#' \dontrun{
#' # Scatterplot of y vs. average yrep
#' fit <- stan_glm(mpg ~ wt, data = mtcars)
#' pp_check(fit, check = "scatter")
#' 
#' # Separate scatterplots of y vs. a few different yrep datasets 
#' pp_check(fit, check = "scatter", nreps = 3)
#' }
#' 
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
#' # Define a test function 
#' prop_zero <- function(y) mean(y == 0)
#' pp_check(example_model, check = "test", test = "prop_zero", binwidth = 1/20)
#' 
#' @importFrom ggplot2 ggplot aes_string xlab %+replace% theme
#' 
pp_check <- function(object,
                    check = "distributions",
                    nreps = NULL, overlay = TRUE, test = "mean", ...) {
  if (!is.stanreg(object)) 
    stop(deparse(substitute(object)), " is not a stanreg object")
  if (used.optimizing(object)) 
    STOP_not_optimizing("pp_check")
  
  checks <- c("distributions", "residuals", "scatter", "test") #, "refit")
  fn <- switch(match.arg(arg = check, choices = checks),
               'distributions' = "pp_check_dist",
               'residuals' = "pp_check_resid",
               'test' = "pp_check_stat",
               'scatter' = "pp_check_scatter") #, 'refit' = "pp_check_refit")
  if (is.null(nreps) && !fn %in%  c("pp_check_stat", "pp_check_scatter"))
    nreps <- ifelse(fn == "pp_check_dist", 8, 3)
  if (!is.null(nreps) && fn == "pp_check_stat") {
    warning("'nreps' is ignored if check='test'")
    nreps <- NULL
  }
  
  # if (fn == "pp_check_refit") {
  #   return(pp_check_refit(object, n = nreps, ...))
  # }
  if (fn == "pp_check_resid") {
    if (!is(object, "polr") && is.binomial(object$family$family)) {
      return(pp_check_binned_resid(object, n = nreps, ...))
    }
  }
  
  yrep <- posterior_predict(object, draws = nreps)
  y <- get_y(object)
  if (is(object, "polr")) y <- as.integer(y)
  if (NCOL(y) == 2) {
    trials <- rowSums(y)
    y <- y[, 1] / trials
    yrep <- sweep(yrep, 2, trials, "/")
  }
  if (fn == "pp_check_dist") 
    args <- list(y = y, yrep = yrep, n = nreps, overlay = overlay, ...)
  else if (fn == "pp_check_resid")
    args <- list(y = y, yrep = yrep, n = nreps, ...)
  else if (fn == "pp_check_stat")
    args <- list(y = y, yrep = yrep, test = test, ...)
  else if (fn == "pp_check_scatter")
    args <- list(y = y, yrep = yrep, n = nreps, ...)
  
  return(do.call(fn, args))
}


# pp_check stuff -----------------------------------------------------------
.PP_FILL <- "skyblue"
.PP_DARK <- "skyblue4"
.PP_VLINE_CLR <- "#222222"
.PP_YREP_CLR <- "#487575"
.PP_YREP_FILL <- "#222222"

#' @importFrom ggplot2 ggtitle element_blank element_line element_text
#'   theme_classic
.pp_check_theme <- function(no_y = TRUE) {
  blank <- element_blank()
  thm <- theme_classic() +
    theme(axis.line = element_line(color = "#222222"),
          axis.line.y = if (no_y) blank  else element_line(size = 0.5),
          axis.line.x = element_line(size = 2),
          axis.title = element_text(face = "bold", size = 13),
          strip.background = element_blank(),
          strip.text = element_text(color = "black", face = "bold"),
          legend.position = "none",
          legend.title = element_text(size = 11),
          legend.text = element_text(size = 13),
          plot.title = element_text(size = 18))
  if (no_y) {
    thm <- thm %+replace% theme(axis.text.y = element_blank(),
                                axis.ticks.y = element_blank(),
                                axis.title.y = element_blank())
  }
  thm
}

.set_geom_args <- function(defaults, ...) {
  dots <- list(...)
  if (!length(dots)) return(defaults)
  dot_names <- names(dots)
  def_names <- names(defaults)
  for (j in seq_along(def_names)) {
    if (def_names[j] %in% dot_names)
      defaults[[j]] <- dots[[def_names[j]]]
  }
  
  extras <- setdiff(dot_names, def_names)
  if (length(extras))
  for (j in seq_along(extras)) {
    defaults[[extras[j]]] <- dots[[extras[j]]]
  }
  return(defaults)
}

pp_check_dist <- function(y, yrep, n = 8, overlay = TRUE, ...) {
  fn <- if (overlay) "pp_check_dens" else "pp_check_hist"
  stopifnot(n <= nrow(yrep))
  s <- 1:n
  yrep <- as.data.frame(yrep)
  colnames(yrep) <- paste0("value.", 1:ncol(yrep))
  yrep_melt <- reshape(yrep, direction = "long", v.names = "value", 
                       varying = list(1:ncol(yrep)), ids = paste0('yrep_', s))
  dat <- rbind(yrep_melt, 
               cbind(time = seq_along(y), value = y, id = 'Observed y'))
  rownames(dat) <- NULL
  dat$is_y <- dat$id == "Observed y"
  dat$value <- as.numeric(dat$value)
  graph <- do.call(fn, list(dat = dat, ...))
  return(graph + .pp_check_theme())
}

#' @importFrom ggplot2 geom_histogram facet_wrap facet_grid stat_bin
pp_check_hist <- function(dat, ...) {
  defaults <- list(size = 0.2)
  geom_args <- .set_geom_args(defaults, ...)
  geom_args$mapping <- aes_string(y = "..density..")
  
  ggplot(dat, aes_string(x = 'value', fill = 'is_y', 
                         color = "is_y", size = "is_y")) + 
    do.call("geom_histogram", geom_args) + 
    facet_wrap(~id, scales = "free") + 
    scale_fill_manual(values = c("black", .PP_FILL)) +
    scale_color_manual(values = c(NA, NA)) + 
    scale_size_manual(values = c(NA, NA)) + 
    xlab(NULL)
}

#' @importFrom ggplot2 geom_density scale_alpha_manual scale_size_manual
#'   scale_fill_manual scale_color_manual xlab
pp_check_dens <- function(dat, ...) {
  graph <- ggplot(dat, aes_string(x = 'value', group = 'id',
                                  color = "is_y", fill = "is_y", 
                                  size = 'is_y')) + 
    geom_density(...) + 
    scale_color_manual(values = c("black", .PP_DARK)) +
    scale_fill_manual(values = c(NA, .PP_FILL)) +
    scale_size_manual(values = c(0.2, 1)) +
    xlab("y")
}

#' @importFrom ggplot2 geom_vline annotate
#' @importFrom utils packageVersion
pp_check_stat <- function(y, yrep, test = "mean", ...) {
  vline_color <- .PP_FILL
  fill_color <- "black"
  if (missing(test) || !length(test) || length(test) > 2) 
    stop("'test' should have length 1 or 2.")
  if (!is.character(test)) 
    stop("'test' should be a character vector.")
  
  if (length(test) == 1) {
    defaults <- list(fill = fill_color, na.rm = TRUE)
    geom_args <- .set_geom_args(defaults, ...)
    geom_args$mapping <- aes_string(y = "..density..")
    geom_args$show.legend <- FALSE
    color_scale <-  scale_color_manual(name = "", 
                                       values = c(vline_color, fill_color),
                                       labels = c("T(y)", "T(yrep)"))
    test1 <- match.fun(test)
    T_y <- test1(y)
    T_yrep <- apply(yrep, 1, test1)
    base <- ggplot(data.frame(x = T_yrep), aes_string(x = "x", color = "'A'")) + 
      xlab(paste("Test =", test))
    graph <- base + 
      do.call("geom_histogram", geom_args) + 
      geom_vline(data = data.frame(t = T_y), 
                 aes_string(xintercept = "t", color = "factor(t)"), 
                 size = 2, show.legend = TRUE) +
      color_scale
    thm <- .pp_check_theme() %+replace% theme(legend.position = "right")
    return(graph + thm)
  }
  else { # length(test) == 2
    defaults <- list(shape = 21, color = "black", fill = "black", alpha = 0.75)
    geom_args <- .set_geom_args(defaults, ...)
    
    if (is.character(test[1])) test1 <- match.fun(test[1])
    if (is.character(test[2])) test2 <- match.fun(test[2])
    T_y1 <- test1(y)
    T_y2 <- test2(y)
    T_yrep1 <- apply(yrep, 1, test1)
    T_yrep2 <- apply(yrep, 1, test2)
    base <- ggplot(data.frame(x = T_yrep1, y = T_yrep2), 
                   aes_string(x = "x", y = "y", color = "'A'"))
    graph <- base + 
      do.call("geom_point", geom_args) + 
      annotate("segment", x = c(T_y1, -Inf), xend = c(T_y1, T_y1), 
               y = c(-Inf, T_y2), yend = c(T_y2, T_y2), 
               color = vline_color, linetype = 2) + 
      geom_point(data = data.frame(x = T_y1, y = T_y2), 
                 aes_string(x = "x", y = "y", color = "'B'"), 
                 size = 4) + 
      scale_color_manual(name = "", 
                         values = c('B' = vline_color, 'A' = fill_color),
                         labels = c('B' = "T(y)", 'A' = "T(yrep)")) + 
      labs(x = paste("Test =", test[1]), y = paste("Test =", test[2]))
    
    thm <- .pp_check_theme(no_y = FALSE) %+replace% 
      theme(legend.position = "right")
    return(graph + thm)
  }
}

#' @importFrom ggplot2 labs
pp_check_resid <- function(y, yrep, n = 1, ...) {
  defaults <- list(fill = "black")
  geom_args <- .set_geom_args(defaults, ...)
  geom_args$mapping <- aes_string(y = "..density..")
  stopifnot(n <= nrow(yrep))
  s <- 1:n
  if (n == 1) {
    base <- ggplot(data.frame(x = y - yrep), aes_string(x = "x"))
  } else {
    resids <- as.data.frame(-1 * sweep(yrep, 2, y))
    colnames(resids) <- paste0("r.", 1:ncol(resids))
    resids <- reshape(resids, direction = "long", v.names = "r", 
                      varying = list(1:ncol(resids)), 
                      ids = paste0("resid(yrep_",s,")"))
    base <- ggplot(resids, aes_string(x = "r"))
  }
  graph <- base + do.call("geom_histogram", geom_args)
  if (n == 1)
    graph <- graph + labs(y = NULL, x = paste0("resid(yrep_",s,")"))
  else 
    graph <- graph + labs(y = NULL, x = NULL) + facet_wrap(~id, scales = "free")
  
  return(graph + .pp_check_theme())
}

#' @importFrom ggplot2 geom_hline geom_point geom_path labs facet_wrap
pp_check_binned_resid <- function(object, n = 1, ...) {
  if (!requireNamespace("arm", quietly = TRUE)) 
    stop("This plot requires the 'arm' package (install.packages('arm'))")
  
  binner <- function(rep_id, ey, r, nbins) {
    br <- arm::binned.resids(ey, r, nbins)$binned[, c("xbar", "ybar", "2se")]
    colnames(br) <- c("xbar", "ybar", "se2")
    data.frame(rep = paste0("yrep_", rep_id), br)
  }
  
  dat <- pp_data(object, newdata = NULL)
  stanmat <- as.matrix.stanreg(object)
  sel <- sample(nrow(stanmat), size = n)
  beta <- stanmat[sel, 1:ncol(dat$x), drop = FALSE]
  eta <- linear_predictor(beta, dat$x, dat$offset)
  inverse_link <- linkinv(object)
  Ey <- inverse_link(eta)
  y <- get_y(object)
  if (NCOL(y) == 2) y <- y[, 1] / rowSums(y)
  resids <- sweep(-Ey, MARGIN = 2, STATS = y, "+")
  ny <- length(y)
  stopifnot(ny == ncol(Ey))
  if (ny >= 100) nbins <- floor(sqrt(ny))
  else if (ny > 10 && ny < 100) nbins <- 10
  else nbins <- floor(ny / 2) # if (ny <= 10)
  
  binned <- binner(rep_id = 1, ey = Ey[1, ], r = resids[1, ], nbins)
  if (n > 1) {
    for (i in 2:nrow(resids)) {
      binned <- rbind(binned, binner(rep_id = i, ey = Ey[i, ], 
                                     r = resids[i, ], nbins))
    }
  }
  dots <- list(...)
  line_color <- dots$color %ORifNULL% .PP_FILL
  line_size <- dots$size %ORifNULL% 1
  pt_color <- dots$fill %ORifNULL% .PP_VLINE_CLR
  base <- ggplot(binned, aes_string(x = "xbar"))
  graph <- base + 
    geom_hline(yintercept = 0, linetype = 2) + 
    geom_path(aes_string(y = "se2"), color = line_color, size = line_size) + 
    geom_path(aes_string(y = "-se2"), color = line_color, size = line_size) + 
    geom_point(aes_string(y = "ybar"), shape = 19, color = pt_color) + 
    labs(x = "Expected Values", y = "Average Residual \n (with 2SE bounds)") + 
    ggtitle("Binned Residuals")
  
  if (n > 1) graph <- graph + facet_wrap(~rep, scales = "free")
  
  return(graph + .pp_check_theme(no_y = FALSE))
}


# pp_check_refit <- function(object, n = 1, ...) { # nocov start
#   message("Refitting model using y = yrep...\n")
#   yrep <- as.vector(posterior_predict(object, draws = 1))
#   mf <- model.frame(object)
#   if (is(object, "polr")) fam <- "polr"
#   else fam <- object$family$family
#   
#   if (!is.binomial(fam)) {
#     if (is(object, "polr"))
#       yrep <- factor(yrep, labels = levels(get_y(object)), ordered = TRUE)
#     mf[[1]] <- yrep
#     refit <- update(object, data = mf)
#   }
#   else {
#     y <- get_y(object)
#     if (NCOL(y) == 2) {
#       new_f <- update.formula(formula(object), cbind(yrep_1s, yrep_0s) ~ .)
#       mf2 <- data.frame(yrep_1s = yrep, yrep_0s = rowSums(y) - yrep, mf[, -1])
#       refit <- update(object, formula = new_f, 
#                       data = get_all_vars(new_f, data = mf2))
#     } 
#     else {
#       if (NCOL(y) == 1 && !all(y %in% c(0, 1)))
#         yrep <- yrep / object$weights
#       mf[[1]] <- yrep
#       refit <- update(object, data = mf)
#     }
#   }
#   
#   pp1 <- posterior_predict(object, draws = n)
#   pp2 <- posterior_predict(refit, draws = n)
#   if (is.binomial(fam)) {
#     if (NCOL(y) == 2) {
#       trials <- rowSums(y)
#       pp1 <- sweep(pp1, 2, trials, "/")
#       pp2 <- sweep(pp2, 2, trials, "/")
#     }
#   }
#   varying <- list(1:ncol(pp1))
#   pp1 <- reshape(as.data.frame(pp1), direction = "long", v.names = "value", 
#                  varying = varying)[, c("value", "id")]
#   pp2 <- reshape(as.data.frame(pp2), direction = "long", v.names = "value", 
#                  varying = varying)[, c("value", "id")]
#   dat <- cbind(rbind(pp1, pp2), 
#                model = rep(c("Model", "Checking model"), each = nrow(pp1)))
#   
#   defaults <- list(size = 0.2)
#   geom_args <- .set_geom_args(defaults, ...)
#   geom_args$mapping <- aes_string(y = "..density..")
#   clr_vals <- c("black", .PP_FILL)
#   base <- ggplot(dat, aes_string(x = 'value', fill = "model", color = "model"))
#   graph <- base +
#     do.call("geom_histogram", geom_args) + 
#     scale_fill_manual("", values = clr_vals) +
#     scale_color_manual("", values = clr_vals) + 
#     facet_grid(model ~ id, scales = "fixed") + 
#     xlab("yrep")
#   thm <- .pp_check_theme() %+replace%
#     theme(strip.text = element_blank(), legend.position = "right")
#   
#   return(graph + thm)
# } # nocov end

#' @importFrom ggplot2 geom_abline
pp_check_scatter <- function(y, yrep, n = NULL, ...){
  # If n is null the avg yrep is compared to y, otherwise n single yrep datasets
  # are compared to y in separate facets.
  
  defaults <- list(shape = 21, fill = .PP_FILL, color = "black", 
                   size = 2.5, alpha = 1)
  geom_args <- .set_geom_args(defaults, ...)

  if (is.null(n)) {
    avg_yrep <- colMeans(yrep)
    dat <- data.frame(x = y, y = avg_yrep, z = abs(y - avg_yrep))
    graph <- ggplot(dat, aes_string("x", "y")) + 
      geom_abline(intercept = 0, slope = 1, linetype = 2) +
      do.call("geom_point", geom_args) +
      labs(x = "y", y = "Average yrep")
  }
  else {
    stopifnot(n <= nrow(yrep))
    s <- 1:n
    yrep <- as.data.frame(yrep)
    colnames(yrep) <- paste0("value.", 1:ncol(yrep))
    yrep_melt <- reshape(yrep, direction = "long", v.names = "value", 
                         varying = list(1:ncol(yrep)), ids = paste0('yrep_', s))
    yrep_melt$y <- rep(y, each = n)
    graph <- ggplot(yrep_melt, aes_string(x = "y", y = "value")) + 
      geom_abline(intercept = 0, slope = 1, linetype = 2) +
      do.call("geom_point", geom_args) +
      labs(x = "y", y = "yrep") + 
      facet_wrap(~id, scales = "free")
  }
  return(graph + .pp_check_theme(no_y = FALSE))
}
