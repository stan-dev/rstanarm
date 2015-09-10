`%ORifNULL%` <- function(a, b) {
  if (is.null(a)) b else a
}

get_y <- function(object) {
  object$y %ORifNULL% model.response(model.frame(object))
}
get_x <- function(object) {
  object$x %ORifNULL% model.matrix(object)
}

na_replace <- function(x, replacement) {
  # if x is NA return replacement, else return x itself
  if (is.na(x)) 
    replacement 
  else 
    x
}

maybe_broadcast <- function(x, n) {
  # if x has length <= 1 replicate it n times, else return x itself
  if (length(x) == 0)
    rep(0, times = n)
  else if (length(x) == 1L) 
    rep(x, times = n)
  else 
    x
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
  if (!is.null(x) & any(x <= 0)) {
    nm <- deparse(substitute(x))
    stop(paste(nm, "should be positive"))
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
  if (length(offset) == 0) eta
  else eta + offset
}
linear_predictor.matrix <- function(beta, x, offset = NULL) {
  if (NCOL(beta) == 1L) 
    beta <- as.matrix(beta)
  eta <- beta %*% t(x)
  if (length(offset) == 0) eta
  else sweep(eta, MARGIN = 2L, offset, `+`)
}