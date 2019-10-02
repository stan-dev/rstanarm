#' Compute weighted expectations using LOO
#' 
#' These functions are wrappers around the \code{\link[loo]{E_loo}} function
#' (\pkg{loo} package) that provide compatibility for \pkg{rstanarm} models.
#' 
#' @export
#' @aliases loo_predict loo_linpred loo_predictive_interval
#' 
#' @template reference-loo
#' @template reference-bayesvis
#' @templateVar stanregArg object
#' @template args-stanreg-object
#' @param psis_object An object returned by \code{\link[loo]{psis}}. If missing 
#' then \code{psis} will be run internally, which may be time consuming
#' for models fit to very large datasets.
#' @param ... Currently unused.
#' @inheritParams loo::E_loo
#'   
#' @return A list with elements \code{value} and \code{pareto_k}. 
#'   
#'   For \code{loo_predict} and \code{loo_linpred} the value component is a 
#'   vector with one element per observation. 
#'   
#'   For \code{loo_predictive_interval} the \code{value} component is a matrix
#'   with one row per observation and two columns (like
#'   \code{\link{predictive_interval}}). \code{loo_predictive_interval(..., prob
#'   = p)} is equivalent to \code{loo_predict(..., type = "quantile", probs =
#'   c(a, 1-a))} with \code{a = (1 - p)/2}, except it transposes the result and
#'   adds informative column names.
#'   
#'   See \code{\link[loo]{E_loo}} and \code{\link[loo]{pareto-k-diagnostic}} for
#'   details on the \code{pareto_k} diagnostic.
#'   
#' @examples
#' \dontrun{
#' if (!exists("example_model")) example(example_model)
#' 
#' # optionally, log-weights can be pre-computed and reused
#' psis_result <- loo::psis(log_ratios = -log_lik(example_model))
#' 
#' loo_probs <- loo_linpred(example_model, type = "mean", transform = TRUE, psis_object = psis_result)
#' str(loo_probs)
#' 
#' loo_pred_var <- loo_predict(example_model, type = "var", psis_object = psis_result)
#' str(loo_pred_var)
#' 
#' loo_pred_ints <- loo_predictive_interval(example_model, prob = 0.8, psis_object = psis_result)
#' str(loo_pred_ints)
#' }
#' 
loo_predict.stanreg <-
  function(object, 
           type = c("mean", "var", "quantile"), 
           probs = 0.5,
           ...,
           psis_object = NULL) {
    
    if ("lw" %in% names(list(...))) {
      stop(
        "Due to changes in the 'loo' package, the 'lw' argument ", 
        "is no longer supported. Use the 'psis_object' argument instead."
      )
    }
    
    type <- match.arg(type)
    log_ratios <- -log_lik(object)
    if (is.null(psis_object)) {
      message("Running PSIS to compute weights...")
      r_eff <- loo::relative_eff(exp(-log_ratios), chain_id = chain_id_for_loo(object))
      psis_object <- loo::psis(log_ratios, r_eff = r_eff)
    }
    
    preds <- posterior_predict(object)
    if (is_polr(object) && !is_scobit(object)) {
      preds <- polr_yrep_to_numeric(preds)
    }
    
    loo::E_loo(
      x = preds,
      psis_object = psis_object,
      type = type,
      probs = probs,
      log_ratios = log_ratios
    )
  }

#' @rdname loo_predict.stanreg
#' @export
#' @param transform Passed to \code{\link{posterior_linpred}}.
#'    
loo_linpred.stanreg <-
  function(object,
           type = c("mean", "var", "quantile"),
           probs = 0.5,
           transform = FALSE,
           ..., 
           psis_object = NULL) {
    
    if ("lw" %in% names(list(...))) {
      stop(
        "Due to changes in the 'loo' package, the 'lw' argument ", 
        "is no longer supported. Use the 'psis_object' argument instead."
      )
    }
    
    type <- match.arg(type)
    log_ratios <- -log_lik(object)
    if (is.null(psis_object)) {
      message("Running PSIS to compute weights...")
      r_eff <- loo::relative_eff(exp(-log_ratios), chain_id = chain_id_for_loo(object))
      psis_object <- loo::psis(log_ratios, r_eff = r_eff)
    }
    
    type <- match.arg(type)
    linpreds <- posterior_linpred(object, transform = transform)
    
    loo::E_loo(
      x = linpreds,
      psis_object = psis_object,
      type = type,
      probs = probs,
      log_ratios = log_ratios
    )
  }


#' @rdname loo_predict.stanreg
#' @export
#' @param prob For \code{loo_predictive_interval}, a scalar in \eqn{(0,1)}
#'   indicating the desired probability mass to include in the intervals. The
#'   default is \code{prob=0.9} (\eqn{90}\% intervals).
loo_predictive_interval.stanreg <-
  function(object,
           prob = 0.9,
           ...,
           psis_object = NULL) {
    stopifnot(length(prob) == 1)
    alpha <- (1 - prob) / 2
    probs <- c(alpha, 1 - alpha)
    labs <- paste0(100 * probs, "%")
    E_loo_result <-
      loo_predict.stanreg(object,
                          type = "quantile",
                          probs = probs,
                          psis_object = psis_object,
                          ...)
    intervals <- E_loo_result$value
    rownames(intervals) <- labs
    intervals <- t(intervals)
    list(value = intervals, pareto_k = E_loo_result$pareto_k)
  }



# internal ----------------------------------------------------------------

psis.stanreg <- function(object, ...) {
  message("Running PSIS to compute weights...")
  ll <- log_lik(object)
  r_eff <- loo::relative_eff(exp(ll), chain_id = chain_id_for_loo(object))
  loo::psis(-ll, r_eff = r_eff, ...)
}

