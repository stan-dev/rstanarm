#' Using the ShinyStan GUI to explore \pkg{rstanarm} models
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
#' @seealso \code{\link[shinystan]{launch_shinystan}}
#' @importFrom shinystan launch_shinystan
#'   
#' @examples
#' \dontrun{
#' 
#' # Fit a model
#' fit <- stan_glm(mpg ~ wt, data = mtcars, seed = 12345)   
#' 
#' # Launch the ShinyStan app
#' launch_shinystan(fit) 
#' 
#' # Launch and also save shinystan object as fit_sso
#' fit_sso <- launch_shinystan(fit) 
#' }
#' 
#' 
NULL