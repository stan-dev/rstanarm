#' Using the ShinyStan GUI with stanreg objects
#' 
#' The \code{\link[shinystan]{launch_shinystan}} function will accept a 
#' \link[=stanreg-objects]{'stanreg' object} as input. Currently, almost any
#' model fit using one of \pkg{rstanarm}'s model-fitting functions can be
#' used with ShinyStan. The only exception is that ShinyStan does not currently
#' support \pkg{rstanarm} models fit using \code{algorithm='optimizing'}.
#' 
#' See the \pkg{\link[=shinystan-package]{shinystan}} package documentation for
#' more information.
#' 
#' @name shinystan
#' @aliases launch_shinystan
#' @importFrom shinystan launch_shinystan
#'   
#' @examples
#' \dontrun{
#' # Fit a model
#' fit <- stan_glm(mpg ~ wt, data = mtcars, seed = 12345)   
#' 
#' # Launch the ShinyStan app (saving resulting shinystan object as fit_sso)
#' fit_sso <- launch_shinystan(fit) 
#' 
#' # Launch the ShinyStan app (without saving shinystan object)
#' launch_shinystan(fit) 
#' }
#' 
#' 
NULL