#' Predictive performance
#' 
#' In-sample predictive performance is an overoptimistic estimate of 
#' out-of-sample predictive performance. To better estimate the predictive
#' performance for new (not yet observed) data we can use use
#' \code{\link[=loo.stanreg]{leave-one-out cross-validation}}.
#' 
#' @name predictive-performance
#' 
#' @templateVar stanregArg object
#' @template args-stanreg-object
#' @param psis Should PSIS-weighted estimates of out-of-sample performance also 
#'   be plotted? Recommended but defaults to \code{FALSE} because it may be time
#'   consuming for very large datasets.
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
#' 
#' # too optimistic 
#' class_accuracy(fit) 
#' class_accuracy(fit, balanced = TRUE) 
#' 
#' # approximate out-of-sample performance
#' class_accuracy(fit, psis = TRUE) 
#' class_accuracy(fit, psis = TRUE, balanced = TRUE) 
#' 
#' # posterior ROC plot
#' roc_plot(fit)
#' 
#' # compare posterior ROC to PSIS-LOO ROC
#' roc_plot(fit, psis = TRUE, print_auc = TRUE)
#' 
NULL

#' @rdname predictive-performance
#' @export 
#' @param balanced Should balanced (\code{TRUE}) or unbalanced (\code{FALSE}) 
#'   classification error be reported? The default is unbalanced (\code{FALSE}).
#'   
#' @return \code{class_error} returns a scalar between 0 and 1, which has a
#'   custom print method. \code{class_accuracy} is a simple wrapper around
#'   \code{class_error} that returns \code{1 - class_error(object)}.
#' 
class_error <-
  function(object,
           balanced = FALSE,
           psis = FALSE,
           ...) {
    validate_bernoulli_model(object)
    y <- validate_glm_outcome_support(get_y(object), family = binomial())
    pp <- postpred_probs(object, psis = psis)
    
    if (!balanced) {
      y_pred <- as.integer(pp >= 0.5)
      err <- mean(xor(y_pred, y))
    } else {
      y1 <- y == 1
      y1_pred <- as.integer(pp[y1] >= 0.5)
      y0_pred <- as.integer(pp[!y1] >= 0.5)
      e1 <- mean(xor(y1_pred, y[y1]))
      e0 <- mean(xor(y0_pred, y[!y1]))
      err <- (e0 + e1) / 2
    }
    structure(
      err,
      class = c("classification_performance", class(err)),
      type = "error",
      balanced = isTRUE(balanced),
      psis = isTRUE(psis)
    )
  }

#' @rdname predictive-performance
#' @export
class_accuracy <-
  function(object,
           balanced = FALSE,
           psis = FALSE,
           ...) {
    acc <- 1 - class_error(object, balanced, psis, ...)
    structure(acc, type = "accuracy")
  }

#' @export
print.classification_performance <- function(x, digits = 3, ...) {
  atts <- attributes(x)
  txt <- paste0(if (atts$psis) "PSIS-LOO" else "Posterior", 
               if (atts$balanced) " balanced", 
               " classification ", atts$type, ":\n")

  x2 <- format(round(x, digits = digits), nsmall = digits)
  cat(txt, x2)
  invisible(x)
}


#' @rdname predictive-performance
#' @export 
#' @param print_auc For \code{roc_plot}, should the area under the curve (AUC)
#'   be printed in the legend? Defaults to \code{FALSE}.
#' @param size For \code{roc_plot}, passed to \code{\link[ggplot2]{geom_line}}.
#' 
#' @return \code{roc_plot} returns a ggplot object.
#' 
roc_plot <- function(object, psis = FALSE, print_auc = FALSE, size = 1, ...) {
  validate_bernoulli_model(object)
  df <- roc_data(object, psis = psis)
  roc_curve(df, print_auc = print_auc, size = size)
}


# internal ----------------------------------------------------------------

# @param x stanreg object
validate_bernoulli_model <- function(x) {
  validate_stanreg_object(x)
  if (!is_binomial_model(x) || NCOL(get_y(x)) == 2)
    stop("Only available for Bernoulli models.", call. = FALSE)
}

# @param x stanreg object
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
  
  structure(out, psis = isTRUE(psis))
}

# @param roc_df Object returned by roc_data
# @param print_auc,size See above.
#' @importFrom ggplot2 geom_abline geom_line scale_x_reverse scale_color_manual geom_blank
roc_curve <- function(roc_df, print_auc = FALSE, size = 1) {
  compare_to_psis <- isTRUE(attr(roc_df, "psis"))
  graph <- ggplot(roc_df) +
    geom_abline(
      intercept = 0,
      slope = 1,
      size = 0.5,
      linetype = 2,
      color = "gray30"
    ) +
    ggplot2::geom_line(aes_(x = ~ fpr, y = ~ tpr, color = "Posterior ROC"), 
              size = size)
  
  scale_labels <- "Posterior ROC"
  if (print_auc)
    scale_labels <- 
      paste0(scale_labels, "\n  AUC: ", 
             round(auc_approx(roc_df$fpr, roc_df$tpr), 2), "\n")
    
  if (compare_to_psis) {
    scale_labels <- c(scale_labels, "PSIS-LOO ROC")
    graph <- graph + 
      geom_line(aes_(x = ~ fpr_loo, y = ~ tpr_loo, color = "PSIS-LOO ROC"), 
                size = 0.75 * size)
    if (print_auc)
      scale_labels[2] <- 
        paste0(scale_labels[2], "\n  AUC: ", 
               round(auc_approx(roc_df$fpr_loo, roc_df$tpr_loo), 2), "\n")
  }
  
  scheme <- bayesplot::color_scheme_get()
  graph + 
    scale_color_manual(
      name = "", 
      labels = scale_labels,
      values = c("Posterior ROC" = scheme[[5]], "PSIS-LOO ROC" = scheme[[2]])
    ) + 
    labs(
      x = "False positive rate \n(1-specificity)", 
      y = "True positive rate \n(sensitivity)"
    ) +
    bayesplot::theme_default()
}

auc_approx <- function(x, y) {
  0.5 * sum(diff(x) * (head(y, -1) + tail(y, -1)))
}
