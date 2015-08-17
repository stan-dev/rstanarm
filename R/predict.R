#' Predict method for stanreg objects
#' 
#' @export
#' 
#' @inheritParams stanreg-methods
#' @param ... Ignored.
#' @param newdata Optionally, a data frame in which to look for variables with 
#'   which to predict. If omitted, the model matrix is used.
#' @param type The type of prediction. The default \code{'link'} is on the scale
#'   of the linear predictors; the alternative \code{'response'} is on the scale
#'   of the response variable.
#' @param se.fit A logical scalar indicating if standard errors should be 
#'   returned. The default is \code{FALSE}.
#'   
#' @return A vector if \code{se.fit} is \code{FALSE} and a list if \code{se.fit}
#'   is \code{TRUE}.
#'
#' @seealso \code{\link{posterior_predict}}, \code{\link[stats]{predict.glm}}
#' 

predict.stanreg <- function(object, ..., newdata = NULL, 
                            type = c("link", "response"), se.fit = FALSE) {
  
  type <- match.arg(type)
  if (!se.fit && is.null(newdata)) {
    if (type == "link") return(object$linear.predictors)
    else return(object$fitted.values)
  }
  mcmc <- object$algorithm == "sampling"
  if (mcmc && type == "response")
    stop("type='response' should not be used for models estimated by MCMC.",
         "Use posterior_predict() to draw from the posterior predictive distribution.",
         call. = FALSE)
  dat <- .pp_data(object, newdata)
  stanmat <- if (mcmc) as.matrix(object$stanfit) else stop("MLE not implemented yet")
  beta <- stanmat[, 1:ncol(dat$x)]
  eta <- linear_predictor(beta, dat$x, dat$offset)
  fit <- colMeans(eta)
  se.fit <- apply(eta, 2L, sd)
  if (type == "response") { 
    stop("MLE not implemented yet")
  }
  nlist(fit, se.fit)
}

.pp_data <- function(object, newdata = NULL) {
  if (is.null(newdata)) {
    x <- model.matrix(object, data = object$data) 
    offset <- if (is.null(object$offset)) rep(0, nrow(x)) else object$offset
    return(nlist(x, offset))
  }
  tt <- terms(object)
  Terms <- delete.response(tt)
  m <- model.frame(Terms, newdata, xlev = object$xlevels)
  if (!is.null(cl <- attr(Terms, "dataClasses"))) 
    .checkMFClasses(cl, m)
  x <- model.matrix(Terms, m, contrasts.arg = object$contrasts)
  offset <- rep(0, nrow(x))
  if (!is.null(off.num <- attr(tt, "offset"))) 
    for (i in off.num) {
      offset <- offset + eval(attr(tt, "variables")[[i + 1]], newdata)
    }
  if (!is.null(object$call$offset)) 
    offset <- offset + eval(object$call$offset, newdata)
  nlist(x, offset)
}
