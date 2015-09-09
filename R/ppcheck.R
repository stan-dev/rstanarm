#' Graphical posterior predictive checks
#' 
#' Various plots comparing the observed outcome variable \code{y} to simulated
#' datasets \code{yrep} from the posterior predictive distribution.
#' 
#' @export
#' @inheritParams stanreg-methods
#' @param check The type of plot (possibly abbreviated) to show.
#' @param nreps The number of datasets to generate from the posterior predictive
#'   distribution and show in the plots. The default is \code{1} for 
#'   \code{check='residuals'} and \code{8} for \code{check='distributions'}. If 
#'   \code{check='test statistics'} then \code{nreps} is ignored and the number
#'   of simulated datasets is the number of post-warmup draws from the posterior
#'   distribution.
#' @param overlay For \code{check="distributions"} only, should distributions be
#'   plotted separately (\code{FALSE}) or overlaid in a single plot
#'   (\code{TRUE})? For other values of \code{check} this is ignored.
#' @param test For \code{check="test statistics"} only, \code{test} should be a 
#'   function that computes the desired test statistic. It can be the name of a 
#'   function as a character string (e.g., \code{test = 'mean'}) or a function 
#'   object (e.g., \code{test = sd}, \code{test = function(x) mean(x == 0)}, 
#'   etc.). The resulting plot will show the distribution of \code{test(yrep)} 
#'   over the \code{yrep} datasets. The value of \code{test(y)} be shown in the
#'   plot as a vertical line.
#' @param ... Optional arguments to geoms to control features of the plots.
#' 
#' @return A ggplot object that can be further customized using the
#'   \pkg{ggplot2} package.
#' 
#' @seealso \code{\link{posterior_predict}} for drawing from the posterior 
#'   predictive distribution. Examples of posterior predictive checking can also
#'   be found in the \pkg{rstanarm} vignettes and demos.
#' 
#' @examples
#' fit <- stan_glm(mpg ~ wt + cyl, data = mtcars, chains = 1)
#' 
#' # Compare distribution of y (mpg) to simulated datasets from the model
#' ppcheck(fit, check = "dist")
#' ppcheck(fit, check = "dist", nreps = 15)
#' ppcheck(fit, check = "dist", overlay = TRUE)
#' 
#' # Check histograms of test statistics
#' library(gridExtra)
#' test_mean <- ppcheck(fit, check = "test", test = 'mean')
#' test_sd <- ppcheck(fit, check = "test", test = 'sd')
#' grid.arrange(test_mean, test_sd, ncol = 2)
#' 
#' q25 <- function(x) quantile(x, 0.25)
#' ppcheck(fit, check = "test", test = q25)
#' @importFrom ggplot2 xlab %+replace% theme
ppcheck <- function(object,
                    check = c("distributions", "residuals", "test statistics"),
                    nreps = NULL, overlay = FALSE, test = 'mean',
                    ...) {
  if (!is(object, "stanreg"))
    stop(deparse(substitute(object)), " is not a stanreg object")
  fn <- switch(match.arg(check), 
               'distributions' = "ppcheck_dist",
               'residuals' = "ppcheck_resid",
               'test statistics' = "ppcheck_stat")
  yrep <- posterior_predict(object)
  y <- if (!is.null(object$y)) 
    object$y else model.response(model.frame(object))
  
  stopifnot(NCOL(y) == 1L, is.matrix(yrep), length(y) == ncol(yrep))
  if (is.null(nreps) && fn != "ppcheck_stat") {
    nreps <- ifelse(fn == "ppcheck_resid", 1, 8)
  }
  thm <- .ppcheck_theme()
  args <- list(y = y, yrep = yrep, n = nreps, overlay = overlay, test = test, 
               ...)
  graph <- do.call(fn, args)
  if (fn == "ppcheck_stat") {
    test_lab <- as.character(match.call()[["test"]])
    graph + 
      xlab(paste("Test = ", test_lab)) + 
      thm %+replace% theme(legend.position = "right")
  }
  else graph + thm
}


# ppcheck stuff -----------------------------------------------------------
#' @importFrom ggplot2 element_blank element_line element_text theme_classic
.ppcheck_theme <- function() {
  theme_classic() +
    theme(axis.line = element_line(color = "#222222"),
          axis.line.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          axis.line.x = element_line(size = 3),
          axis.title = element_text(face = "bold", size = 13),
          strip.background = element_blank(),
          strip.text = element_text(color = "black", face = "bold"),
          legend.position = "none",
          legend.title = element_text(size = 11),
          legend.text = element_text(size = 13),
          plot.title = element_text(size = 18)) 
}

.PP_FILL <- "skyblue"
.PP_VLINE_CLR <- "#222222"
.PP_YREP_CLR <- "#487575"
.PP_YREP_FILL <- "#222222"

ppcheck_dist <- function(y, yrep, n = 8, overlay = FALSE, ...) {
  fn <- if (overlay) "ppcheck_dens" else "ppcheck_hist"
  stopifnot(n <= nrow(yrep))
  s <- sample.int(nrow(yrep), n)
  yrep <- as.data.frame(yrep[s, ])
  colnames(yrep) <- paste0("value.", 1:ncol(yrep))
  yrep_melt <- reshape(yrep, direction = "long", v.names = "value", 
                       varying = list(1:ncol(yrep)), ids = paste0('yrep.', s))
  dat <- rbind(yrep_melt, cbind(time = seq_along(y), value = y, id = 'y'))
  rownames(dat) <- NULL
  dat$is_y <- dat$id == "y"
  dat$value <- as.numeric(dat$value)
  do.call(fn, list(dat=dat, n=n, ...))
}

#' @importFrom ggplot2 aes_string facet_wrap ggplot scale_fill_manual stat_bin
ppcheck_hist <- function(dat, ...) {
  ggplot(dat, aes_string(x = 'value', fill = 'is_y')) + 
    stat_bin(aes_string(y="..density.."), size = .2, ...) + 
    facet_wrap(~id, scales = "free") + 
    scale_fill_manual(values = c("black", .PP_FILL)) +
    xlab(NULL)
}

#' @importFrom ggplot2 geom_density scale_size_manual
ppcheck_dens <- function(dat, ...) {
  dat$id <- factor(dat$id, levels = unique(dat$id))
  ggplot(dat, aes_string(x = 'value', group = 'id', 
                         fill = "is_y", size = 'is_y')) + 
    geom_density(alpha = 0.4, ...) + 
    scale_fill_manual(values = c(NA, .PP_FILL)) +
    scale_size_manual(values = c(.5, 1.5)) +
    xlab("y")
}

#' @importFrom ggplot2 geom_vline scale_color_manual
ppcheck_stat <- function(y, yrep, test, ...) {
  if (is.character(test))
    test <- match.fun(test)
  T_y <- test(y)
  T_yrep <- apply(yrep, 1L, test)
  dots <- list(...)
  vline_color <- if ("color" %in% names(dots)) dots$color else .PP_FILL
  fill_color <- if ("fill" %in% names(dots)) dots$fill else "black"
  graph <- ggplot(data.frame(x = T_yrep), aes_string(x = "x", color = "'A'")) +
    stat_bin(aes_string(y = "..count../sum(..count..)"), 
             fill = fill_color, show_guide = FALSE, ...) 
  graph + 
    geom_vline(data = data.frame(t = T_y), 
               aes_string(xintercept = "t", color = "factor(t)"), 
               size = 2, show_guide = TRUE) +
    scale_color_manual(name = "", 
                       values = c(vline_color, fill_color),
                       labels = c("T(y)", "T(yrep)"))
}

#' @importFrom ggplot2 labs
ppcheck_resid <- function(y, yrep, n = 1, ...) {
  stopifnot(n <= nrow(yrep))
  s <- sample.int(nrow(yrep), n)
  yrep <- yrep[s, ]
  if (n == 1) {
    base <- ggplot(data.frame(x = y - yrep), aes_string(x = "x"))
  } else {
    resids <- as.data.frame(-1 * sweep(yrep, 2, y))
    colnames(resids) <- paste0("r.", 1:ncol(resids))
    resids <- reshape((resids), direction = "long", v.names = "r", 
                      varying = list(1:ncol(resids)), ids = paste0('yrep_', s))
    base <- ggplot(resids, aes_string(x = "r"))
  }
  graph <- base + stat_bin(aes_string(y="..count../sum(..count..)"), ...)
  
  if (n == 1) 
    graph + labs(y = NULL, x = paste0("resids(yrep_",s,")"))
  else 
    graph + labs(y = NULL, x = NULL) + facet_wrap(~id, scales = "free")
}
