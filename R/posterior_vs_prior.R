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
#'   the posterior another, and \code{"none"} to use no color. Except when 
#'   \code{color_by="none"}, a variable is mapped to the color 
#'   \code{\link[ggplot2]{aes}}thetic and it is therefore also possible to
#'   change the default colors by adding one of the various discrete color
#'   scales available in \code{ggplot2} 
#'   (\code{\link[ggplot2]{scale_color_manual}}, 
#'   \code{\link[ggplot2]{scale_color_brewer}}, etc.). See Examples.
#' @param prob A number \eqn{p \in (0,1)}{p (0 < p < 1)} indicating the desired 
#'   posterior probability mass to include in the (central posterior) interval 
#'   estimates displayed in the plot. The default is \eqn{0.9}.
#' @param facet_args A named list of arguments passed to
#'   \code{\link[ggplot2]{facet_wrap}} (other than the \code{facets} argument),
#'   e.g., \code{nrow} or \code{ncol} to change the layout, \code{scales} to 
#'   allow axis scales to vary across facets, etc. See Examples.
#' @param ... Arguments (other than \code{color}) passed to 
#'   \code{\link[ggplot2]{geom_pointrange}} to control the appearance of the 
#'   plotted intervals.
#'   
#' @return A ggplot object that can be further customized using the 
#'   \pkg{ggplot2} package.
#'   
#' @examples 
#' # display non-varying (i.e. not group-level) coefficients
#' posterior_vs_prior(example_model, pars = "beta")
#' 
#' # show group-level (varying) parameters and group by parameter
#' posterior_vs_prior(example_model, pars = "varying",
#'                    group_by_parameter = TRUE, color_by = "vs")
#' \dontrun{
#' # group by parameter and allow axis scales to vary across facets
#' posterior_vs_prior(example_model, regex_pars = "period",
#'                    group_by_parameter = TRUE, color_by = "none",
#'                    facet_args = list(scales = "free"))
#' 
#' # assign to object and customize with functions from ggplot2
#' (gg <- posterior_vs_prior(example_model, pars = c("beta", "varying"), prob = 0.8))
#' gg + 
#'  ggplot2::scale_color_brewer() + 
#'  ggplot2::theme(panel.background = ggplot2::element_rect(fill = "gray30"))
#' 
#'  ggplot2::geom_hline(yintercept = 0, size = 0.3, linetype = 3) + 
#'  ggplot2::coord_flip() + 
#'  ggplot2::ggtitle("Comparing the prior and posterior")
#' gg + 
#'  ggplot2::theme(panel.background = element_rect(fill = "gray30"), 
#'                 axis.text.x = element_blank())
#'                 
#'                 
#' # compare very wide and very narrow priors using roaches example
#' # (see help(roaches, "rstanarm") for info on the dataset)
#' roaches$roach100 <- roaches$roach1 / 100
#' wide_prior <- normal(0, 10)
#' narrow_prior <- normal(0, 0.1)
#' fit_pois_wide_prior <- stan_glm(y ~ treatment + roach100 + senior, 
#'                                 offset = log(exposure2), 
#'                                 family = "poisson", data = roaches, 
#'                                 prior = wide_prior)
#' posterior_vs_prior(fit_pois_wide_prior, pars = "beta", prob = 0.5, 
#'                    group_by_parameter = TRUE, color_by = "vs", 
#'                    facet_args = list(scales = "free"))
#'                    
#' fit_pois_narrow_prior <- update(fit_pois_wide_prior, prior = narrow_prior)
#' posterior_vs_prior(fit_pois_narrow_prior, pars = "beta", prob = 0.5, 
#'                    group_by_parameter = TRUE, color_by = "vs", 
#'                    facet_args = list(scales = "free"))
#' }
#' 
#' @importFrom ggplot2 geom_pointrange facet_wrap aes_string labs
#'   scale_x_discrete element_line element_text
#' 
posterior_vs_prior <- function(object, 
                               pars = NULL, regex_pars = NULL, prob = 0.9, 
                               color_by = c("parameter", "vs", "none"),
                               group_by_parameter = FALSE,
                               facet_args = list(),
                               ...) {
  if (!is.stanreg(object))
    stop(deparse(substitute(object)), " is not a stanreg object.")
  if (!used.sampling(object))
    STOP_sampling_only("posterior_vs_prior")
  stopifnot(isTRUE(prob > 0 && prob < 1))
  
  # stuff needed for ggplot
  color_by <- switch(match.arg(color_by), 
                     parameter = "parameter", 
                     vs = "model", 
                     none = NA)
  if (group_by_parameter) {
    group_by <- "parameter" 
    xvar <- "model"
  } else {
    group_by <- "model"
    xvar <- "parameter"
  }
  aes_args <- list(x = xvar, y = "estimate", ymin = "lb", ymax = "ub")
  if (!is.na(color_by))
    aes_args$color <- color_by
  if (!length(facet_args)) {
    facet_args <- list(facets = group_by)
  } else {
    facet_args$facets <- group_by
  }

  # draw from prior distribution and prepare plot data
  message("\nDrawing from prior...")
  Prior <- suppressWarnings(update(object, prior_PD = TRUE, 
                                   refresh = -1, chains = 2))
  objects <- nlist(Prior, Posterior = object)
  plot_data <- stack_estimates(objects, prob = prob,
                               pars = pars, regex_pars = regex_pars)
  

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


stack_estimates <- function(models = list(), pars = NULL, regex_pars = NULL, 
                            prob = NULL) {
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
    s <- summary(x, pars = pars, regex_pars = regex_pars, probs = probs)
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
