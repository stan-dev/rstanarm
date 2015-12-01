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
#' @param check The type of plot (possibly abbreviated) to show. See Details.
#' @param nreps The number of \eqn{yrep} datasets to generate from the posterior
#'   predictive distribution and show in the plots. The default is
#'   \code{nreps=3} for \code{check='residuals'} and \code{check='refit'}, and
#'   \code{nreps=8} for \code{check='distributions'}. If \code{check='test'}
#'   then \code{nreps} is ignored and the number of simulated datasets is the
#'   number of post-warmup draws from the posterior distribution.
#' @param overlay For \code{check='distributions'} only, should distributions be
#'   plotted as density estimates overlaid in a single plot (\code{TRUE}, the
#'   default) or as separate histograms (\code{FALSE})?
#' @param test For \code{check='test'} only, a character vector (of length 1 or 
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
#'  \item{\code{test}}{The distribution of a single test statistic \eqn{T(yrep)}
#'  or a pair of test statistics over the \code{nreps} simulated datasets. If 
#'  the \code{test} argument specifies only one function then the resulting plot
#'  is a histogram of \eqn{T(yrep)} and the value of the test statistic in the
#'  observed data, \eqn{T(y)}, is shown in the plot as a vertical line. If two
#'  functions are specified then the plot is a scatterplot and \eqn{T(y)} is
#'  shown as a large point.}
#'  \item{\code{residuals}}{The distributions of residuals computed from
#'  \eqn{y} and each of \code{nreps} simulated datasets. For binomial data,
#'  binned residual plots are generated.}
#'  \item{\code{refit}}{First a \emph{checking model} is fit using the same 
#'  predictors but replacing the outcome \eqn{y} with a single realization of
#'  \eqn{yrep} from the posterior predictive distribution. Then \code{nreps}
#'  datasets are simulated from both the original model and the checking model
#'  and plotted.}
#' }
#' 
#' @note For binomial data with a number of trials greater than 1
#'   (i.e. not Bernoulli), plots of \eqn{y} and \eqn{yrep} show the
#'   proportion of 'successes' rather than the raw count.
#' 
#' @seealso \code{\link{posterior_predict}} for drawing from the posterior 
#'   predictive distribution. Examples of posterior predictive checks can also
#'   be found in the \pkg{rstanarm} vignettes and demos.
#' 
#' @examples
#' # Compare distribution of y to distributions of yrep
#' (pp_dist <- ppcheck(example_model, check = "distributions", overlay = TRUE))
#' pp_dist + 
#'  ggplot2::scale_color_manual(values = c("red", "black")) + # change colors
#'  ggplot2::scale_size_manual(values = c(0.5, 3)) + # change line sizes 
#'  ggplot2::scale_fill_manual(values = c(NA, NA)) # remove fill
#'
#' # Check residuals
#' ppcheck(example_model, check = "residuals", nreps = 3, fill = "blue") + 
#'   ggplot2::ggtitle("Residuals (y - yrep)")
#'
#' # Check histograms of test statistics
#' test_mean <- ppcheck(example_model, check = "test", test = "mean")
#' test_sd <- ppcheck(example_model, check = "test", test = "sd")
#' gridExtra::grid.arrange(test_mean, test_sd, ncol = 2)
#' # Scatterplot of two test statistics
#' ppcheck(example_model, check = "test", test = c("mean", "sd"))
#' 
#' # Refit using yrep and compare posterior predictive distributions of 
#' # original model and checking model
#' \dontrun{
#' ppcheck(example_model, check = "refit", nreps = 3)
#' }
#' 
#' @importFrom ggplot2 ggplot aes_string xlab %+replace% theme
#' 
ppcheck <- function(object,
                    check = c("distributions", "test", "residuals", "refit"),
                    nreps = NULL, overlay = TRUE, test = "mean", ...) {
  if (!is.stanreg(object))
    stop(deparse(substitute(object)), " is not a stanreg object")
  if (!used.sampling(object)) 
    STOP_sampling_only("ppcheck")
  fn <- switch(match.arg(check), 
               'distributions' = "ppcheck_dist",
               'residuals' = "ppcheck_resid",
               'test' = "ppcheck_stat",
               'refit' = "ppcheck_refit")
  if (is.null(nreps) && fn != "ppcheck_stat")
    nreps <- ifelse(fn == "ppcheck_dist", 8, 3) 
  if (fn == "ppcheck_refit") {
    graph <- ppcheck_refit(object, n = nreps, ...) + 
      .ppcheck_theme() +
      theme(strip.text = element_blank(), legend.position = "right")
    return(graph)
  }
  if (fn == "ppcheck_resid") {
    if (!is(object, "polr") && is.binomial(object$family$family)) {
    graph <- ppcheck_binned_resid(object, n = nreps, ...) + 
      ggtitle("Binned Residuals") + 
      .ppcheck_theme(no_y = FALSE)
    return(graph)
    }
  }
  thm <- .ppcheck_theme()
  yrep <- posterior_predict(object, draws = nreps)
  y <- get_y(object)
  if (is(object, "polr")) 
    y <- as.integer(y)
  if (NCOL(y) == 2) {
    trials <- rowSums(y)
    y <- y[, 1] / trials
    yrep <- sweep(yrep, 2, trials, "/")
  }
  if (fn == "ppcheck_dist") 
    args <- list(y = y, yrep = yrep, n = nreps, overlay = overlay, ...)
  else if (fn == "ppcheck_resid")
    args <- list(y = y, yrep = yrep, n = nreps, ...)
  else if (fn == "ppcheck_stat")
    args <- list(y = y, yrep = yrep, test = test, ...)
  
  graph <- do.call(fn, args)
  if (fn == "ppcheck_stat") {
    if (length(test) == 2) thm <- .ppcheck_theme(no_y = FALSE)
    return(graph + thm %+replace% theme(legend.position = "right"))
  }
  else return(graph + thm)
}


# ppcheck stuff -----------------------------------------------------------
.PP_FILL <- "skyblue"
.PP_DARK <- "skyblue4"
.PP_VLINE_CLR <- "#222222"
.PP_YREP_CLR <- "#487575"
.PP_YREP_FILL <- "#222222"

#' @importFrom ggplot2 ggtitle element_blank element_line element_text
#'   theme_classic
.ppcheck_theme <- function(no_y = TRUE) {
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

ppcheck_dist <- function(y, yrep, n = 8, overlay = TRUE, ...) {
  fn <- if (overlay) "ppcheck_dens" else "ppcheck_hist"
  stopifnot(n <= nrow(yrep))
  s <- 1:n
  yrep <- as.data.frame(yrep)
  colnames(yrep) <- paste0("value.", 1:ncol(yrep))
  yrep_melt <- reshape(yrep, direction = "long", v.names = "value", 
                       varying = list(1:ncol(yrep)), ids = paste0('rep_', s))
  dat <- rbind(yrep_melt, cbind(time = seq_along(y), value = y, id = 'Observed'))
  rownames(dat) <- NULL
  dat$is_y <- dat$id == "Observed"
  dat$value <- as.numeric(dat$value)
  do.call(fn, list(dat = dat, ...))
}

#' @importFrom ggplot2 geom_histogram facet_wrap facet_grid stat_bin
ppcheck_hist <- function(dat, ...) {
  ggplot(dat, aes_string(x = 'value', fill = 'is_y', color = "is_y", size = "is_y")) + 
    geom_histogram(aes_string(y="..density.."), size = .2, ...) + 
    facet_wrap(~id, scales = "free") + 
    scale_fill_manual(values = c("black", .PP_FILL)) +
    scale_color_manual(values = c(NA, NA)) + 
    scale_size_manual(values = c(NA, NA)) + 
    xlab(NULL)
}

#' @importFrom ggplot2 geom_density scale_alpha_manual scale_size_manual
#'   scale_fill_manual scale_color_manual xlab
ppcheck_dens <- function(dat, ...) {
  # dat$id <- factor(dat$id, levels = unique(dat$id))
  graph <- ggplot(dat, aes_string(x = 'value', group = 'id',
                         color = "is_y", fill = "is_y", size = 'is_y')) + 
    geom_density(...) + 
    scale_color_manual(values = c("black", .PP_DARK)) +
    scale_fill_manual(values = c(NA, .PP_FILL)) +
    scale_size_manual(values = c(0.2, 1)) +
    xlab("y")
}

#' @importFrom ggplot2 geom_vline annotate
#' @importFrom utils packageVersion
ppcheck_stat <- function(y, yrep, test = "mean", ...) {
  vline_color <- .PP_FILL
  fill_color <- "black"
  if (missing(test) || !length(test) || length(test) > 2) 
    stop("'test' should have length 1 or 2.")
  if (!is.character(test)) 
    stop("'test' should be a character vector.")
  if (length(test) == 1) {
    test1 <- match.fun(test)
    T_y <- test1(y)
    T_yrep <- apply(yrep, 1, test1)
    graph <- ggplot(data.frame(x = T_yrep), aes_string(x = "x", color = "'A'"))
    graph <- graph + xlab(paste("Test =", test))
    if (packageVersion("ggplot2") < "1.1.0") {
      graph + 
        geom_histogram(aes_string(y = "..count../sum(..count..)"), 
                       fill = fill_color, show_guide = FALSE, na.rm = TRUE, ...)  +
        geom_vline(data = data.frame(t = T_y), 
                   aes_string(xintercept = "t", color = "factor(t)"), 
                   size = 2, show_guide = TRUE) +
        scale_color_manual(name = "", 
                           values = c(vline_color, fill_color),
                           labels = c("T(y)", "T(yrep)"))
    } else {
      graph + 
        geom_histogram(aes_string(y = "..count../sum(..count..)"), 
                       fill = fill_color, show.legend = FALSE, na.rm = TRUE, ...)  +
        geom_vline(data = data.frame(t = T_y), 
                   aes_string(xintercept = "t", color = "factor(t)"), 
                   size = 2, show.legend = TRUE) +
        scale_color_manual(name = "", 
                           values = c(vline_color, fill_color),
                           labels = c("T(y)", "T(yrep)"))
    }
  }
  else { # length(test) == 2
    if (is.character(test[1])) test1 <- match.fun(test[1])
    if (is.character(test[2])) test2 <- match.fun(test[2])
    T_y1 <- test1(y)
    T_y2 <- test2(y)
    T_yrep1 <- apply(yrep, 1, test1)
    T_yrep2 <- apply(yrep, 1, test2)
    graph <- ggplot(data.frame(x = T_yrep1, y = T_yrep2), 
                    aes_string(x = "x", y = "y", color = "'A'"))
    graph + 
      geom_point(...) + 
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
  }
}

#' @importFrom ggplot2 labs
ppcheck_resid <- function(y, yrep, n = 1, ...) {
  stopifnot(n <= nrow(yrep))
  s <- 1:n
  if (n == 1) {
    base <- ggplot(data.frame(x = y - yrep), aes_string(x = "x"))
  } else {
    resids <- as.data.frame(-1 * sweep(yrep, 2, y))
    colnames(resids) <- paste0("r.", 1:ncol(resids))
    resids <- reshape((resids), direction = "long", v.names = "r", 
                      varying = list(1:ncol(resids)), ids = paste0('rep_', s))
    base <- ggplot(resids, aes_string(x = "r"))
  }
  graph <- base + geom_histogram(aes_string(y="..count../sum(..count..)"), ...)
  
  if (n == 1) 
    graph + labs(y = NULL, x = paste0("resids(yrep_",s,")"))
  else 
    graph + labs(y = NULL, x = NULL) + facet_wrap(~id, scales = "free")
}

#' @importFrom ggplot2 geom_hline geom_point geom_path labs facet_wrap
ppcheck_binned_resid <- function(object, n = 1, ...) {
  if (!requireNamespace("arm", quietly = TRUE)) 
    stop("This plot requires the 'arm' package (install.packages('arm'))")
  dat <- .pp_data(object, newdata = NULL)
  stanmat <- if (used.sampling(object))
    as.matrix(object$stanfit) else stop("MLE not implemented yet")
  sel <- sample(nrow(stanmat), size = n)
  beta <- stanmat[sel, 1:ncol(dat$x), drop=FALSE]
  eta <- linear_predictor(beta, dat$x, dat$offset)
  Ey <- family(object)$linkinv(eta)
  y <- get_y(object)
  if (NCOL(y) == 2) y <- y[, 1] / rowSums(y)
  resids <- sweep(-Ey, MARGIN = 2, STATS = y, "+")
  ny <- length(y)
  stopifnot(ny == ncol(Ey))
  if (ny >= 100) nbins <- floor(sqrt(ny))
  else if (ny > 10 && ny < 100) nbins <- 10
  else nbins <- floor(ny / 2) # if (ny <= 10)
  
  binner <- function(rep_id, ey, r, nbins) {
    br <- arm::binned.resids(ey, r, nbins)$binned[, c("xbar", "ybar", "2se")]
    colnames(br) <- c("xbar", "ybar", "se2")
    data.frame(rep = paste0("rep_", sel[rep_id]), br)
  }
  binned <- binner(rep_id = 1, ey = Ey[1, ], r = resids[1, ], nbins)
  if (n > 1) {
    for (i in 2:nrow(resids)) {
      binned <- rbind(binned, 
                      binner(rep_id = i, ey = Ey[i, ], r = resids[i, ], nbins))
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
    labs(x = "Expected Values", y = "Average Residual \n (with 2SE bounds)")
  
  if (n == 1) graph
  else graph + facet_wrap(~rep, scales = "free")
}

ppcheck_refit <- function(object, n = 1, ...) {
  message("Refitting model using y = yrep...\n")
  yrep <- as.vector(posterior_predict(object, draws = 1))
  mf <- model.frame(object)
  if (is(object, "polr")) fam <- "polr"
  else fam <- object$family$family
  
  if (!is.binomial(fam)) {
    if (is(object, "polr"))
      yrep <- factor(yrep, labels = levels(get_y(object)), ordered = TRUE)
    mf[[1]] <- yrep
    refit <- update(object, data = mf)
  }
  else {
    y <- get_y(object)
    if (NCOL(y) == 2) {
      new_f <- update.formula(formula(object), cbind(yrep_1s, yrep_0s) ~ .)
      mf2 <- data.frame(yrep_1s = yrep, yrep_0s = rowSums(y) - yrep, mf[, -1])
      refit <- update(object, formula = new_f, data = get_all_vars(new_f, data = mf2))
    } 
    else {
      if (NCOL(y) == 1 && !all(y %in% c(0, 1)))
        yrep <- yrep / object$weights
      mf[[1]] <- yrep
      refit <- update(object, data = mf)
    }
  }
  
  pp1 <- posterior_predict(object, draws = n)
  pp2 <- posterior_predict(refit, draws = n)
  if (is.binomial(fam)) {
    if (NCOL(y) == 2) {
      trials <- rowSums(y)
      pp1 <- sweep(pp1, 2, trials, "/")
      pp2 <- sweep(pp2, 2, trials, "/")
    }
  }

  varying <- list(1:ncol(pp1))
  pp1 <- reshape(as.data.frame(pp1), direction = "long", v.names = "value", 
                 varying = varying)[, c("value", "id")]
  pp2 <- reshape(as.data.frame(pp2), direction = "long", v.names = "value", 
                 varying = varying)[, c("value", "id")]
  dat <- cbind(rbind(pp1, pp2), model = rep(c("Model", "Checking model"), each = nrow(pp1)))
  clr_vals <- c("black", .PP_FILL)
  ggplot(dat, aes_string(x = 'value', fill = "model", color = "model")) + 
    geom_histogram(aes_string(y = "..density.."), size = .2, ...) +
    scale_fill_manual("", values = clr_vals) +
    scale_color_manual("", values = clr_vals) + 
    facet_grid(model ~ id, scales = "fixed") + 
    xlab("yrep")
}

