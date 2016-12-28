#' Predictive performance
#' 
#' @name predictive-performance
#' 
#' @templateVar stanregArg object
#' @template args-stanreg-object
#' @param psis Should PSIS-weighted estimates also be plotted? Defaults to
#'   \code{TRUE}.
#' @param ... Currently ignored.
#' 
#' @examples 
#' head(wells) # included with rstanarm
#' wells200 <- wells[1:200, ] # use only part of the data to emphasize the effect
#' wells200$dist100 <- wells200$dist / 100 # rescale 
#' 
#' fit <- stan_glm(
#'  switch ~ dist100 + arsenic + assoc + educ, 
#'  data = wells200, 
#'  family = binomial("logit"),
#'  prior = student_t(7, 0, 10), 
#'  prior_intercept = student_t(7, 0, 10)
#' )
#' roc_plot(fit)
#' roc_plot(fit, psis = TRUE)
#' 
NULL

#' @rdname predictive-performance
#' @export 
roc_plot <- function(object, psis = FALSE, ...) {
  validate_bernoulli_model(object)
  df <- roc_data(object, psis = psis)
  roc_curve(df)
}

#' @rdname predictive-performance
#' @export 
#' @param balanced Should balanced or unbalanced classification accuracy be
#'   reported? The default is \code{FALSE} (unbalanced).
posterior_classification <- function(object, balanced = FALSE, psis = FALSE, ...) {
  validate_bernoulli_model(object)
  pp <- postpred_probs(object, psis = psis)
  y <- validate_glm_outcome_support(get_y(object), family = binomial())
  if (!balanced) {
    p_gt_half <- pp > 0.5
    m <- mean(xor(p_gt_half, y))
  } else {
    y0 <- y == 0
    y1 <- y == 1
    p_gt_half_0 <- pp[y0] > 0.5
    p_gt_half_1 <- pp[y1] > 0.5
    m0 <- mean(xor(p_gt_half_0, y[y0]))
    m1 <- mean(xor(p_gt_half_1, y[y1]))
    m <- (m0 + m1) / 2
  }
  structure(m, 
            class = c("class_acc", class(m)), 
            balanced = isTRUE(balanced), 
            psis = isTRUE(psis))
}

#' @export
print.class_acc <- function(x, digits, ...) {
  atts <- attributes(x)
  txt <- paste0(if (atts$psis) "PSIS-LOO" else "Posterior", 
               if (atts$balanced) " balanced", 
               " classification accuracy:\n")
  x2 <- x
  if (!missing(digits)) {
    x2 <- format(round(x2, digits = digits), nsmall = digits)
  }
  cat(txt, x2)
  invisible(x)
}



# internal ----------------------------------------------------------------
validate_bernoulli_model <- function(x) {
  validate_stanreg_object(x)
  if (!is_binomial_model(x) || NCOL(get_y(x)) == 2)
    stop("Only available for Bernoulli models.", call. = FALSE)
}

# @param x stanreg object
# @param psis see above
postpred_probs <- function(x, psis = FALSE) {
  pps <- posterior_linpred(x, transform = TRUE)
  if (!psis) 
    return(colMeans(pps))
  
  message("Computing PSIS-weighted estimates...")
  psis_object <- loo::psislw(llfun = ll_fun(x), llargs = ll_args(x))
  colSums(pps * exp(psis_object$lw_smooth))
}


# @param x stanreg object
roc_data <- function(x, psis = TRUE) {
  y <- validate_glm_outcome_support(get_y(x), family = binomial())
  
  # posterior predictive probabilities
  pp <- postpred_probs(x)
  y_pp <- y[order(pp, decreasing=TRUE)]
  
  # true and false positive rates
  out <- data.frame(tpr = cumsum(y_pp) / sum(y_pp),
                    fpr = cumsum(!y_pp) / sum(!y_pp))
  
  # PSIS-LOO predictive probabilities
  if (psis) {
    ploo <- postpred_probs(x, psis = TRUE)
    y_loo <- y[order(ploo, decreasing=TRUE)]
    out$tpr_loo <- cumsum(y_loo) / sum(y_loo)
    out$fpr_loo <- cumsum(!y_loo) / sum(!y_loo)
  }
  
  return(out)
}

#' @importFrom ggplot2 geom_abline geom_line scale_x_reverse scale_color_manual
roc_curve <- function(roc_df) {
  plot_data <- data.frame(
    x = 100 * (1 - roc_df$fpr),
    y = 100 * roc_df$tpr
  )
  
  compare_to_psis <- "tpr_loo" %in% colnames(roc_df)
  if (compare_to_psis) {
    plot_data$x_loo <- 100 * (1 - roc_df$fpr_loo)
    plot_data$y_loo <- 100 * roc_df$tpr_loo
  }
  
  clrs <- bayesplot::color_scheme_get()
  graph <- ggplot(plot_data) +
    geom_abline(
      intercept = 100,
      slope = 1,
      size = 0.5,
      linetype = 2,
      color = "gray30"
    ) +
    geom_line(aes_(x = ~ x, y = ~ y, color = "Posterior ROC"))
  
  if (compare_to_psis) {
    graph <- graph + 
      geom_line(aes_(x = ~ x_loo, y = ~ y_loo, color = "PSIS-LOO ROC"))
  }
  
  graph + 
    scale_color_manual(
      name = "", 
      values = c("Posterior ROC" = clrs[[5]], "PSIS-LOO ROC" = clrs[[2]])
    ) + 
    labs(x = "Specificity (%)", y = "Sensitivity (%)") +
    scale_x_reverse() +
    bayesplot::theme_default()
}


area_approx <- function(x, y) {
  0.5 * sum(diff(x) * (head(y, -1) + tail(y, -1)))
}
