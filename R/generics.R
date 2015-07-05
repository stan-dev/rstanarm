#' Standard errors
#' 
#' @name se
#' @keywords internal
#' @export
#' @param object The object to use. 
#' @param parm A character vector of parameter names.
#' 
se <- function(object, parm) UseMethod("se")


#' Pointwise log-likelihood
#' 
#' @name log_lik
#' @keywords internal
#' @export
#' @param object The object to use. 
#'
log_lik <- function(object) UseMethod("log_lik")
