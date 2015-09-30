stan_control <- list(adapt_delta = 0.95, max_treedepth = 15L)

is.binomial <- function(x) x == "binomial"
is.gaussian <- function(x) x == "gaussian"
is.gamma <- function(x) x == "Gamma"
is.ig <- function(x) x == "inverse.gaussian"
is.nb <- function(x) x == "neg_binomial_2"
is.poisson <- function(x) x == "poisson"

`%ORifNULL%` <- function(a, b) {
  if (is.null(a)) b else a
}

`%ORifINF%` <- function(a, b) {
  if (a == Inf) b else a
}

maybe_broadcast <- function(x, n) {
  # if x has no length replicate 0 n times, 
  # if x has length 1 replicate x n times
  # else return x itself
  if (!length(x)) rep(0, times = n)
  else if (length(x) == 1L) rep(x, times = n)
  else x
}

nlist <- function(...) {
  # named lists
  m <- match.call()
  out <- list(...)
  no_names <- is.null(names(out))
  has_name <- if (no_names)
    FALSE else nzchar(names(out))
  if (all(has_name))
    return(out)
  nms <- as.character(m)[-1]
  if (no_names) {
    names(out) <- nms
  } else {
    names(out)[!has_name] <- nms[!has_name]
  }
  out
}

validate_parameter_value <- function(x) {
  # check for positive scale or df parameter
  nm <- deparse(substitute(x))
  if (!is.null(x)) {
    if (!is.numeric(x)) stop(nm, " should be NULL or numeric")
    if (any(x <= 0)) stop(nm, " should be positive")
  }
  invisible(TRUE)
}

set_prior_scale <- function(scale, default, link) {
  stopifnot(is.numeric(default), is.character(link))
  if (is.null(scale))
    scale <- default
  if (link == "probit")
    scale * dnorm(0) / dlogis(0)
  else 
    scale
}

# create linear predictor vector from x and point estimates for beta, or linear
# predictor matrix from x and full posterior sample of beta
linear_predictor <- function(beta, x, offset = NULL) {
  UseMethod("linear_predictor")
}
linear_predictor.default <- function(beta, x, offset = NULL) {
  eta <- as.vector(if (NCOL(x) == 1L) x * beta else x %*% beta)
  if (!length(offset)) eta
  else eta + offset
}
linear_predictor.matrix <- function(beta, x, offset = NULL) {
  if (NCOL(beta) == 1L) 
    beta <- as.matrix(beta)
  eta <- beta %*% t(x)
  if (!length(offset)) eta
  else sweep(eta, MARGIN = 2L, offset, `+`)
}


#' Extracy X, Y or Z
#' 
#' @keywords internal
#' @export
#' @param object object
get_y <- function(object) UseMethod("get_y")
#' @rdname get_y
#' @export
get_x <- function(object) UseMethod("get_x")
#' @rdname get_y
#' @export
get_z <- function(object) UseMethod("get_z")

#' @export
get_y.default <- function(object) {
  object$y %ORifNULL% model.response(model.frame(object))
}
#' @export
get_x.default <- function(object) {
  object$x %ORifNULL% model.matrix(object)
}
#' @export
get_x.lmerMod <- function(object) {
  object$glmod$X %ORifNULL% stop("X not found")
}
#' @export
get_z.lmerMod <- function(object) {
  Zt <- object$glmod$reTrms$Zt %ORifNULL% stop("Z not found")
  t(as.matrix(Zt))
}
