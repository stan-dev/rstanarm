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
#'   (\code{TRUE}, the default)? For other values of \code{check} this is
#'   ignored.
#' @param test For \code{check="test statistics"} only, \code{test} should be a 
#'   function that computes the desired test statistic. It can be the name of a 
#'   function as a character string (e.g., \code{test = 'mean'}) or a function 
#'   object (e.g., \code{test = sd}, \code{test = function(x) mean(x == 0)}, 
#'   etc.). The resulting plot will show the distribution of \code{test(yrep)} 
#'   over the \code{yrep} datasets. The value of \code{test(y)} be shown in the
#'   plot as a vertical line.
#' @param ... Optional arguments to geoms to control features of the plots. For 
#'   the \code{'distributions'} plots, don't use \code{...} to specify the 
#'   aesthetics (\code{color}, \code{fill}, \code{size}, and \code{alpha}) as
#'   they are mapped to a binary variable so that \code{y} and \code{yrep} can
#'   be styled separately. To change these aesthetics add a
#'   \code{\link[ggplot2]{discrete_scale}} (e.g.
#'   \code{\link[ggplot2]{scale_fill_manual}}) with two values to the ggplot
#'   object returned by \code{ppcheck}. See Examples.
#' 
#' @return A ggplot object that can be further customized using the
#'   \pkg{ggplot2} package.
#' 
#' @seealso \code{\link{posterior_predict}} for drawing from the posterior 
#'   predictive distribution. Examples of posterior predictive checking can also
#'   be found in the \pkg{rstanarm} vignettes and demos.
#' 
#' @examples
#' \dontrun{
#' options(mc.cores = parallel::detectCores())
#' fit <- stan_glm(mpg ~ wt + cyl, data = mtcars)
#' 
#' # Compare distribution of y (mpg) to simulated datasets 
#' # from the model (posterior predictive distribution)
#' (pp_dist <- ppcheck(fit, check = "dist", nreps = 20))
#' 
#' pp_dist + 
#'  scale_color_manual(values = c("red", "black")) + # change line colors
#'  scale_size_manual(values = c(0.1, 2)) + # change line sizes 
#'  scale_fill_manual(values = c(NA, NA)) # remove fill
#'  
#' # Show the distributions as separate histograms instead of overlaid  
#' (pp_dist_sep <- ppcheck(fit, check = "dist", overlay = FALSE)) 
#' pp_dist_sep + scale_fill_manual(values = c("blue", "red")) # change fill colors
#' 
#' # Check histograms of test statistics
#' library(gridExtra)
#' test_mean <- ppcheck(fit, check = "test", test = 'mean')
#' test_sd <- ppcheck(fit, check = "test", test = 'sd')
#' grid.arrange(test_mean, test_sd, ncol = 2)
#' 
#' q25 <- function(x) quantile(x, 0.25)
#' ppcheck(fit, check = "test", test = q25)
#' 
#' # Residuals
#' ppcheck(fit, check = "resid", nreps = 3) + ggtitle("Residuals y - yrep")
#' 
#' # For logistic regressions binned residual plots are generated instead.
#' # See the Examples for the stan_glm function 
#' # help("stan_glm", package = "rstanarm")
#' }
#' @importFrom ggplot2 xlab %+replace% theme
#' 
ppcheck <- function(object,
                    check = c("distributions", "residuals", "test statistics"),
                    nreps = NULL, overlay = TRUE, test = 'mean',
                    ...) {
  if (!is(object, "stanreg"))
    stop(deparse(substitute(object)), " is not a stanreg object")
  
  fn <- switch(match.arg(check), 
               'distributions' = "ppcheck_dist",
               'residuals' = "ppcheck_resid",
               'test statistics' = "ppcheck_stat")
  if (is.null(nreps) && fn != "ppcheck_stat") {
    nreps <- ifelse(fn == "ppcheck_resid", 1, 8)
  }
  if (fn == "ppcheck_resid" && family(object) == "binomial") {
    graph <- pp_check_binned_resid(object, n = nreps, ...) + 
      .ppcheck_theme(no_y = FALSE) + 
      ggtitle("Binned Residual Plot")
    return(graph)
  } 
  thm <- .ppcheck_theme()
  yrep <- posterior_predict(object)
  y <- get_y(object)
  if (NCOL(y) == 2L) {
    trials <- rowSums(y)
    y <- y[, 1L] / trials
    yrep <- sweep(yrep, 2, trials, "/")
  }
  args <- list(y = y, yrep = yrep, n = nreps, overlay = overlay, test = test, 
               ...)
  graph <- do.call(fn, args)
  if (fn == "ppcheck_stat") {
    test_lab <- as.character(match.call()[["test"]])
    if (!length(test_lab)) test_lab <- formals(ppcheck)$test
    graph + 
      xlab(paste("Test = ", test_lab)) + 
      thm %+replace% theme(legend.position = "right")
  }
  else graph + thm
}


# ppcheck stuff -----------------------------------------------------------
.PP_FILL <- "skyblue"
.PP_DARK <- "skyblue4"
.PP_VLINE_CLR <- "#222222"
.PP_YREP_CLR <- "#487575"
.PP_YREP_FILL <- "#222222"

#' @importFrom ggplot2 ggtitle element_blank element_line element_text theme_classic
.ppcheck_theme <- function(no_y = TRUE) {
  blank <- element_blank()
  thm <- theme_classic() +
    theme(axis.line = element_line(color = "#222222"),
          axis.line.y = if (no_y) blank  else element_line(size = 0.5),
          axis.line.x = element_line(size = 3),
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

ppcheck_dist <- function(y, yrep, n = 8, overlay = FALSE, ...) {
  fn <- if (overlay) "ppcheck_dens" else "ppcheck_hist"
  stopifnot(n <= nrow(yrep))
  s <- sample.int(nrow(yrep), n)
  yrep <- as.data.frame(yrep[s, ])
  colnames(yrep) <- paste0("value.", 1:ncol(yrep))
  yrep_melt <- reshape(yrep, direction = "long", v.names = "value", 
                       varying = list(1:ncol(yrep)), ids = paste0('rep_', s))
  dat <- rbind(yrep_melt, cbind(time = seq_along(y), value = y, id = 'Observed'))
  rownames(dat) <- NULL
  dat$is_y <- dat$id == "Observed"
  dat$value <- as.numeric(dat$value)
  do.call(fn, list(dat=dat, n=n, ...))
}

#' @importFrom ggplot2 aes_string facet_wrap ggplot stat_bin
ppcheck_hist <- function(dat, ...) {
  ggplot(dat, aes_string(x = 'value', fill = 'is_y', color = "is_y", size = "is_y")) + 
    stat_bin(aes_string(y="..density.."), size = .2, ...) + 
    facet_wrap(~id, scales = "free") + 
    scale_fill_manual(values = c("black", .PP_FILL)) +
    scale_color_manual(values = c(NA, NA)) + 
    scale_size_manual(values = c(NA, NA)) + 
    xlab(NULL)
}

#' @importFrom ggplot2 geom_density xlab scale_size_manual scale_fill_manual scale_color_manual
ppcheck_dens <- function(dat, ...) {
  # dat$id <- factor(dat$id, levels = unique(dat$id))
  dots <- list(...)
  ggplot(dat, aes_string(x = 'value', group = 'id',
                         color = "is_y", fill = "is_y", size = 'is_y',
                         alpha = "is_y")) + 
    geom_density(...) + 
    scale_alpha_manual(values = c(0, dots$alpha %ORifNULL% 1)) +
    scale_color_manual(values = c("black", .PP_DARK)) +
    scale_fill_manual(values = c("black", .PP_FILL)) +
    scale_size_manual(values = c(0.25, 1)) +
    xlab("y")
}

#' @importFrom ggplot2 geom_vline
ppcheck_stat <- function(y, yrep, test = "mean", ...) {
  if (is.character(test))
    test <- match.fun(test)
  T_y <- test(y)
  T_yrep <- apply(yrep, 1L, test)
  dots <- list(...)
  vline_color <- dots$color %ORifNULL% .PP_FILL
  fill_color <- dots$fill %ORifNULL% "black"
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
                      varying = list(1:ncol(resids)), ids = paste0('rep_', s))
    base <- ggplot(resids, aes_string(x = "r"))
  }
  graph <- base + stat_bin(aes_string(y="..count../sum(..count..)"), ...)
  
  if (n == 1) 
    graph + labs(y = NULL, x = paste0("resids(yrep_",s,")"))
  else 
    graph + labs(y = NULL, x = NULL) + facet_wrap(~id, scales = "free")
}

#' @importFrom ggplot2 geom_hline geom_point geom_path labs facet_wrap
pp_check_binned_resid <- function(object, n = 1, ...) {
  if (!requireNamespace("arm", quietly = TRUE)) 
    stop("This plot requires the 'arm' package (install.packages('arm'))")
  dat <- .pp_data(object, newdata = NULL)
  stanmat <- if (object$algorithm == "sampling") 
    as.matrix(object$stanfit) else stop("MLE not implemented yet")
  sel <- sample(nrow(stanmat), size = n)
  beta <- stanmat[sel, 1:ncol(dat$x), drop=FALSE]
  eta <- linear_predictor(beta, dat$x, dat$offset)
  Ey <- family(object)$linkinv(eta)
  y <- get_y(object)
  if (NCOL(y) == 2L) y <- y[, 1L] / rowSums(y)
  resids <- sweep(-Ey, MARGIN = 2L, STATS = y, "+")
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
    geom_path(aes_string(y = "se2"), color = line_color, size = line_size, ...) + 
    geom_path(aes_string(y = "-se2"), color = line_color, size = line_size, ...) + 
    geom_point(aes_string(y = "ybar"), shape = 19, color = pt_color, ...) + 
    labs(x = "Expected Values", y = "Average Residual \n (with 2SE bounds)")
  
  if (n == 1) graph
  else graph + facet_wrap(~rep, scales = "free")
}
