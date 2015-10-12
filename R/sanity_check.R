#' Perform a sanity check
#' 
#' Refit a model using the same predictors but taking the outcome vector 
#' \code{y} to be a draw from the \link[=posterior_predict]{posterior predictive
#' distribution}. It can be useful to compare the estimates from the sanity
#' check model to the estimates obtained from fitting the model to the real
#' data. See Examples.
#' 
#' @export
#' @inheritParams stanreg-methods
#' @param ... Arguments to or from other methods.
#' @return The stanreg object resulting from fitting the model to the simulated
#'   data.
#'   
#' @seealso \code{\link{posterior_predict}}
#' @examples 
#' \dontrun{
#' (fit <- stan_glm(mpg ~ wt, data = mtcars, prior = normal(0,1), 
#'                  prior_intercept = normal(0,1)))
#' (fit_check <- sanity_check(fit))
#'
#' library(gridExtra)
#' grid.arrange(plot(fit), plot(fit_check))
#' 
#' estimates <- cbind(fit = coef(fit), check = coef(fit_check), 
#'                    se_fit = se(fit), se_check = se(fit_check))
#' round(estimates)                    
#' diffs <- summary(fit) - summary(fit_check)
#' diffs[, colnames(diffs) %in% c("mean", "sd", "50%")]
#' }  
#' 
sanity_check <- function(object, ...) {
  if (object$algorithm != "sampling") {
    message("Only available for models fit using MCMC (algorithm = 'sampling').")
    return(invisible(NULL))
  }
  yrep <- as.vector(posterior_predict(object, draws = 1))
  mf <- model.frame(object)
  if (is.binomial(object$family$family)) {
    y <- get_y(object)
    if (NCOL(y) == 1L && !all(y %in% c(0, 1)))
      yrep <- yrep / object$weights
    else 
      yrep <- cbind(yrep, rowSums(y) - yrep)
  }
  mf[[1L]] <- yrep
  update(object, data = mf, ...)
}
