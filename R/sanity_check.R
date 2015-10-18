#' Perform a sanity check
#' 
#' Refit a model using the same predictors but replacing the outcome \code{y} 
#' with data simulated from the \link[=posterior_predict]{posterior predictive 
#' distribution}. It can be useful to compare the resulting estimates to the 
#' estimates obtained from fitting the model to the observed data. See Examples.
#' 
#' @export
#' @inheritParams stanreg-methods
#' @param ... Arguments to or from other methods.
#' @return The \code{\link[=stanreg-objects]{stanreg}} object resulting from
#'   fitting the model to the simulated data.
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
#' diffs <- round(summary(fit) - summary(fit_check), 1)
#' diffs[, c("mean", "sd", "50%")]
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
    if (NCOL(y) == 2L) 
      yrep <- cbind(yrep, rowSums(y) - yrep)
  }
  mf[[1L]] <- yrep
  update(object, data = mf, ...)
}
