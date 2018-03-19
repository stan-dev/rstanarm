#' Compute weighted expectations using LOO
#' 
#' These functions are wrappers around the \code{\link[loo]{E_loo}} function 
#' (\pkg{loo} package).
#' 
#' @export
#' @aliases loo_predict loo_linpred loo_predictive_interval
#' 
#' @template reference-loo
#' @templateVar stanregArg object
#' @template args-stanreg-object
#' @param lw Deprecated. Use \code{psis_object} instead.
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
#' head(loo_predictive_interval(example_model, prob = 0.8, cores = 2))
#' 
#' # optionally, log-weights can be pre-computed and reused
#' psis_result <- loo::psis(log_ratios = -log_lik(example_model))
#' 
#' loo_linpred(example_model, type = "mean", transform = TRUE, psis_object = psis_result)
#' loo_predict(example_model, type = "var", psis_object = psis_result)
#' loo_predictive_interval(example_model, prob = 0.8, psis_object = psis_result)
#' }
#' 
loo_predict.stanreg <-
  function(object, 
           type = c("mean", "var", "quantile"), 
           probs = 0.5,
           ...,
           psis_object = NULL,
           lw = NULL) {
    
    if (!is.null(lw)) {
      warning(
        "The 'lw' argument is now deprectated. ", 
        "Use the 'psis_object' argument instead."
      )
      
      if (!is.null(psis_object)) {
        warning("Using 'psis_object' instead of 'lw'.")
      } else {
        return(
          loo_predict_old(
            object = object,
            type = type,
            probs = probs,
            ...,
            lw = lw
          )
        )
      }
    } 
    
    log_ratios <- -log_lik(object)
    
    if (is.null(lw) && is.null(psis_object)) {
      message("Running PSIS to compute weights...")
      r_eff <- loo::relative_eff(exp(-log_ratios), chain_id = chain_id_for_loo(object))
      psis_object <- loo::psis(log_ratios, r_eff = r_eff)
    }
    
    
    type <- match.arg(type)
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
           psis_object = NULL,
           lw = NULL) {
    
    if (!is.null(lw)) {
      warning(
        "The 'lw' argument is now deprectated. ", 
        "Use the 'psis_object' argument instead."
      )
      
      if (!is.null(psis_object)) {
        warning("Using 'psis_object' instead of 'lw'.")
      } else {
        return(
          loo_linpred_old(
            object = object,
            type = type,
            probs = probs,
            transform = transform,
            ...,
            lw = lw
          )
        )
      }
    } 
    
    log_ratios <- -log_lik(object)
    
    if (is.null(lw) && is.null(psis_object)) {
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
           psis_object = NULL,
           lw = NULL) {
    stopifnot(length(prob) == 1)
    alpha <- (1 - prob) / 2
    probs <- c(alpha, 1 - alpha)
    labs <- paste0(100 * probs, "%")
    E_loo_result <-
      loo_predict.stanreg(object,
                          type = "quantile",
                          probs = probs,
                          psis_object = psis_object,
                          lw = lw,
                          ...)
    intervals <- E_loo_result$value
    rownames(intervals) <- labs
    intervals <- t(intervals)
    list(value = intervals, pareto_k = E_loo_result$pareto_k)
  }

# internal ----------------------------------------------------------------

# for backwards compatibility if lw argument and not psis_object is provided
loo_predict_old <- function(object,
                            type = c("mean", "var", "quantile"),
                            probs = 0.5,
                            ...,
                            lw) {
  type <- match.arg(type)
  lwts <- loo_weights(object, lw, log = TRUE, ...)
  preds <- posterior_predict(object)
  if (is_polr(object) && !is_scobit(object))
    preds <- polr_yrep_to_numeric(preds)
  
  loo::E_loo(
    x = preds,
    lw = lwts,
    type = type,
    probs = probs
  )
}
loo_linpred_old <- function(object,
                            type = c("mean", "var", "quantile"),
                            probs = 0.5,
                            transform = FALSE,
                            ...,
                            lw) {
  type <- match.arg(type)
  lwts <- loo_weights(object, lw, log = TRUE, ...)
  linpreds <- posterior_linpred(object, transform = transform)
  
  loo::E_loo(
    x = linpreds,
    lw = lwts,
    type = type,
    probs = probs
  )
}

# @param object,lw,... Same as above.
# @param log If FALSE (default) the weights are exponentiated before returning
# @return A matrix.
loo_weights <- function(object, lw, log = FALSE, ...) {
  if (!missing(lw)) {
    stopifnot(is.matrix(lw))
  } else {
    message("Running PSIS to compute weights...")
    psis <- loo::psislw(llfun = ll_fun(object), llargs = ll_args(object), ...)
    lw <- psis[["lw_smooth"]]
  }
  if (log) 
    return(lw) 
  
  exp(lw)
}
