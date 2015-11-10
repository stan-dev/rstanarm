#' Predict method for stanreg objects
#' 
#' @export
#' @templateVar stanregArg object
#' @template args-stanreg-object
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
  if (used.sampling(object) && type == "response")
    stop("type='response' should not be used for models estimated by MCMC.",
         "Use posterior_predict() to draw from the posterior predictive distribution.",
         call. = FALSE)
  
  dat <- pp_data(object, newdata)
  stanmat <- if (used.sampling(object))
    as.matrix.stanreg(object) else stop("MLE not implemented yet")
  beta <- stanmat[, 1:ncol(dat$x)]
  eta <- linear_predictor(beta, dat$x, dat$offset)
  if (is(object, "lmerMod")) {
    b <- stanmat[, .bnames(colnames(stanmat)), drop = FALSE]
    eta <- eta + linear_predictor(b, dat$z)
  }
  fit <- colMeans(eta)
  if (type == "response") { 
    stop("MLE not implemented yet")
  }
  if (!se.fit) return(fit)
  else {
    se.fit <- apply(eta, 2L, sd)
    nlist(fit, se.fit) 
  }
}
