#' Standard errors
#' 
#' @name se
#' @keywords internal
#' @export
#' @param object The object to use. 
#' @param parm A character vector of parameter names.
#' @return estimated standard errors for the parameters named in \code{parm}.
#' 
se <- function(object, parm) UseMethod("se")


#' Pointwise log-likelihood
#' 
#' @name log_lik
#' @keywords internal
#' @export
#' @param object The object to use. 
#' @return pointwise log-likelihood matrix
#'
log_lik <- function(object) UseMethod("log_lik")
