#' Juxtapose prior and posterior
#' 
#' Plot medians and central intervals comparing parameter draws from the prior 
#' and posterior distributions.
#' 
#' @export
#' @templateVar stanregArg object
#' @template args-stanreg-object
#' @inheritParams summary.stanreg
#' @param group_by_parameter Should estimates be grouped together by parameter
#'   (\code{TRUE}) or by posterior and prior (\code{FALSE}, the default)?
#' @param color_by How should the estimates be colored? Use \code{"parameter"} 
#'   to color by parameter name, \code{"vs"} to color the prior one color and 
#'   the posterior another, and \code{"none"} to use no color.
#' @param prob A number \eqn{p \in (0,1)}{p (0 < p < 1)} indicating the desired 
#'   posterior probability mass to include in the (central posterior) interval 
#'   estimates displayed in the plot.
#' @param facet_args A named list of arguments passed to
#'   \code{\link[ggplot2]{facet_wrap}} (other than the \code{facets} argument),
#'   e.g., \code{nrow} or \code{ncol} to change the layout, \code{scales} to 
#'   allow axis scales to vary across facets, etc. See Examples.
#' @param ... Arguments (other than \code{color}) passed to 
#'   \code{\link[ggplot2]{geom_pointrange}} to control the appearance of the 
#'   plotted intervals.
#'   
#' @return A ggplot object that can be further customized using the 
#'   \pkg{ggplot2} package. Except when \code{color_by="none"}, a variable is 
#'   also mapped to the color aesthetic and it is therefore also possible to 
#'   change the default colors by adding one of the various discrete color 
#'   scales available in \code{ggplot2}
#'   (\code{\link[ggplot2]{scale_color_manual}}, 
#'   \code{\link[ggplot2]{scale_color_brewer}}, etc.). See Examples.
#'   
#' @examples 
#' posterior_vs_prior(example_model)
#' posterior_vs_prior(example_model, pars = NULL)
#' posterior_vs_prior(example_model, group_by_parameter = TRUE, 
#'                    facet_args = list(scales = "free"))
#' posterior_vs_prior(example_model, pars = "varying", color_by = "none")
#' 
#' # assign to object and customize with functions from ggplot2
#' (gg <- posterior_vs_prior(example_model, prob = 0.5))
#' gg + ggplot2::coord_flip()
#' gg + ggplot2::scale_color_manual(values = c("orange", "purple", "blue", "green"))
#' gg + ggplot2::scale_color_grey()
#' gg + 
#'  ggplot2::scale_color_brewer() + 
#'  ggplot2::theme(panel.background = element_rect(fill = "gray30"))
#' 
#' @importFrom ggplot2 geom_pointrange facet_wrap aes_string labs
#' 
posterior_vs_prior <- function(object,
                               pars = "beta",
                               color_by = c("parameter", "vs", "none"),
                               group_by_parameter = FALSE,
                               prob = 0.9,
                               facet_args = list(),
                               ...) {
  if (!is.stanreg(object))
    stop(deparse(substitute(object)), " is not a stanreg object.")
  stopifnot(isTRUE(prob > 0 && prob < 1))
  group_by <- if (group_by_parameter) "parameter" else "model"
  color_by <- switch(match.arg(color_by), 
                     parameter = "parameter", 
                     vs = "model", 
                     none = NA)
  
  objects <- list(Posterior = object, Prior = update(object, prior_PD = TRUE))
  plot_data <- stack_estimates(objects, pars = pars, prob = prob)
  if (!length(facet_args)) {
    facet_args <- list(facets = group_by)
  } else {
    facet_args$facets <- group_by
  }
  aes_args <- list(y = "estimate", ymin = "lb", ymax = "ub", 
                   x = if (group_by == "parameter") "model" else "parameter")
  if (!is.na(color_by))
    aes_args$color <- color_by
  
  graph <- ggplot(plot_data, mapping = do.call("aes_string", aes_args)) + 
    geom_pointrange(...) + 
    do.call("facet_wrap", facet_args) +
    labs(x = NULL, y = NULL) +
    pp_check_theme(no_y = FALSE) +
    theme(panel.grid.major.y = element_line(size = 0.1, color = "gray"))
  
  if (group_by == "parameter")
    return(graph)
  
  # clean up x-axis labels a bit if tick labels are parameter names
  # (user can override this after plot is created if need be, 
  # but this makes the default a bit nicer if many parameters)
  graph + 
    theme(axis.text.x = element_text(angle = -30, hjust = 0)) +
    scale_x_discrete(labels = abbreviate(plot_data$parameter, 12, 
                                         method = "both.sides", 
                                         dot = TRUE))
}


stack_estimates <- function(models = list(), pars = NULL, prob = NULL) {
  mnames <- names(models)
  if (is.null(mnames)) {
    mnames <- paste0("model_", seq_along(models))
  } else {
    has_name <- nzchar(mnames)
    if (!all(has_name))
      stop("Either all or none of the elements in 'models' should be named.")
  }
  
  alpha <- (1 - prob) / 2
  probs <- sort(c(0.5, alpha, 1 - alpha))
  labs <- c(paste0(100 * probs, "%"))
  ests <- lapply(models, function(x) {
    s <- summary(x, pars = pars, probs = probs)
    if (is.null(pars))
      s <- s[!rownames(s) %in% c("log-posterior", "mean_PPD"), ]
    s[, labs, drop = FALSE]
  })
  est_column <- function(list_of_matrices, col) {
    x <- sapply(list_of_matrices, function(x) x[, col])
    if (is.list(x)) unlist(x) else as.vector(x)
  }
  data.frame(
    model = rep(mnames, times = sapply(ests, nrow)),
    parameter = unlist(lapply(ests, rownames)),
    estimate = est_column(ests, labs[2]),
    lb = est_column(ests, labs[1]),
    ub = est_column(ests, labs[3]))
}
