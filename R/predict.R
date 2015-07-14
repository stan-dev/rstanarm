#' Predict method for stan_lm and stan_glm models
#' 
#' @export
#' @param object a model fit with \code{\link{stan_lm}} or 
#'   \code{\link{stan_glm}}.
#' @param ... ignored.
#' @param newdata optionally, a data frame in which to look for variables with 
#'   which to predict. If omitted, the model matrix is used.
#' @param type the type of prediction. The default \code{'link'} is on the scale
#'   of the linear predictors; the alternative \code{'response'} is on the scale
#'   of the response variable.
#' @param se.fit logical switch indicating if standard errors should be
#'   returned. Ignored if \code{output} argument is set to \code{'sims'}.
#' @param output the type of output returned. If \code{'point'} then the output 
#'   is similar to \code{\link[stats]{predict.glm}}. If \code{'sims'} then the 
#'   output is a matrix containing either simulations from the posterior
#'   predictive distribution ( if \code{type='response'}) or the
#'   posterior distribution of the linear predictor (if \code{type='link'}).
#'   
#' @return If \code{output='sims'}, a matrix. If \code{output='point'}, a 
#'   vector is returned if \code{se.fit} is \code{FALSE} and a list is returned 
#'   if \code{se.fit} is \code{TRUE}.
#'

predict.stanreg <- function(object, ..., newdata = NULL, type = c("link", "response"),
                            se.fit = FALSE, output = c("point", "sims")) {
  
  type <- match.arg(type)
  output <- match.arg(output)
  family <- object$family
  famname <- family$family
  ppfun <- paste0(".pp_", famname)

  if ((!se.fit && output == "point") && is.null(newdata)) {
    if (type == "link") return(object$linear.predictors)
    else return(object$fitted.values)
  }
  dat <- .pp_data(object, newdata)
  stanmat <- as.matrix(object$stanfit)
  beta <- stanmat[, 1:ncol(dat$x)]
  eta <- linear_predictor(beta, dat$x, dat$offset)
  
  if (type == "link") {
    if (output == "sims") return(eta)
    else { # output == "point"
      if (!se.fit) return(colMeans(eta))
      else return(list(fit = colMeans(eta), se.fit = apply(eta, 2, sd)))
    }
  }
  
  ppargs <- list(mu = family$linkinv(eta))
  if (famname == "gaussian")
    ppargs$sigma <- stanmat[, "sigma"]
  if (famname == "binomial") {
    y <- if (!is.null(object$y)) 
      object$y else model.response(model.frame(object))
    ppargs$trials <- if (NCOL(y) == 2L) rowSums(y) else rep(1, NROW(y))
  }
  yrep <- do.call(ppfun, ppargs) # yrep is niter x ndata matrix
  
  if (output == "sims") return(yrep)
  else {
    if (famname == "binomial" && !all(ppargs$trials == 1))
      yrep <- yrep / ppargs$trials
    if (!se.fit) return(colMeans(yrep))
    else return(list(fit = colMeans(yrep), se.fit = apply(yrep, 2, sd))) 
  }
}

.pp_gaussian <- function(mu, sigma) {
  t(sapply(1:nrow(mu), function(s) {
    rnorm(ncol(mu), mu[s,], sigma[s])
  }))
}
.pp_poisson <- function(mu) {
  t(sapply(1:nrow(mu), function(s) {
    rpois(ncol(mu), mu[s,])
  }))
}
.pp_binomial <- function(mu, trials) {
  t(sapply(1:nrow(mu), function(s) {
    rbinom(ncol(mu), size = trials, prob = mu[s,])
  }))
}
.pp_data <- function(object, newdata = NULL) {
  if (is.null(newdata)) {
    x <- model.matrix(object) 
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
