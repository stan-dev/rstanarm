# Default 'control' argument for stan() if none specified by user
# value of adapt_delta depends on prior
default_stan_control <- function(prior, adapt_delta = NULL, max_treedepth = 15L) {
  if (is.null(adapt_delta)) {
    adapt_delta <- switch(prior$dist, 
                          "R2" = 0.99,
                          "hs" = 0.99,
                          "hs_plus" = 0.99,
                          "t" = if (any(prior$df <= 2)) 0.99 else 0.95,
                          0.95) # default
  }
  nlist(adapt_delta, max_treedepth)
}

# Prepares a list of arguments to use with rstan::sampling via do.call()
# @param object The stanfit object to use for sampling
# @param user_dots The contents of ... from the user's call to stan_*
# @param user_adapt_delta The value for adapt_delta specified by the user
# @param prior Prior distribution list
# @param ... Other arguments to sampling not coming from user_dots (e.g. data,
#   pars, init, etc.)
set_sampling_args <- function(object, prior, user_dots = list(), 
                              user_adapt_delta = NULL, ...) {
  args <- list(object = object, ...)
  unms <- names(user_dots)
  for (j in seq_along(user_dots)) args[[unms[j]]] <- user_dots[[j]]
  if ("control" %in% unms && !is.null(user_adapt_delta)) {
    args$control$adapt_delta <- user_adapt_delta
  }
  else args$control <- 
    default_stan_control(prior = prior, adapt_delta = user_adapt_delta)
  return(args)
}

# Test if an object is a stanreg object
is.stanreg <- function(x) inherits(x, "stanreg")

# Test for a given family
# @param x Character vector (probably x = family(fit)$family)
is.binomial <- function(x) x == "binomial"
is.gaussian <- function(x) x == "gaussian"
is.gamma <- function(x) x == "Gamma"
is.ig <- function(x) x == "inverse.gaussian"
is.nb <- function(x) x == "neg_binomial_2"
is.poisson <- function(x) x == "poisson"

# Grep for "b" parameters (ranef)
# @param x Character vector (often rownames(fit$stan_summary))
.bnames <- function(x, ...) grep("^b\\[", x, ...)

# Use 50% or Median depending on algorithm
# @param algorithm String naming the algorithm (probably fit$algorithm)
.select_median <- function(algorithm) {
  switch(algorithm, 
         sampling = "50%",
         optimizing = "Median",
         stop("Incorrect algorithm name"))
}

# If a is NULL (and Inf, respectively) return b, otherwise just return a
# @param a,b Objects
`%ORifNULL%` <- function(a, b) {
  if (is.null(a)) b else a
}
`%ORifINF%` <- function(a, b) {
  if (a == Inf) b else a
}

# If x has no length replicate 0 n times, x has length 1 replicate x n times, 
# otherwise return x itself
maybe_broadcast <- function(x, n) {
  if (!length(x)) rep(0, times = n)
  else if (length(x) == 1L) rep(x, times = n)
  else x
}

# Create a named list using specified names or, if names are omitted, using the
# names of the objects in ...
nlist <- function(...) {
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

# Check for positive scale or df parameter (NULL ok)
validate_parameter_value <- function(x) {
  nm <- deparse(substitute(x))
  if (!is.null(x)) {
    if (!is.numeric(x)) stop(nm, " should be NULL or numeric")
    if (any(x <= 0)) stop(nm, " should be positive")
  }
  invisible(TRUE)
}

# Check and set scale parameters for priors
set_prior_scale <- function(scale, default, link) {
  stopifnot(is.numeric(default), is.character(link))
  if (is.null(scale))
    scale <- default
  if (link == "probit")
    scale * dnorm(0) / dlogis(0)
  else 
    scale
}

# Create linear predictor vector from x and point estimates for beta, or linear 
# predictor matrix from x and full posterior sample of beta
# @param beta A vector or matrix or parameter estimates
# @param x Predictor matrix
# @param offset Optional offset vector
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


#' Extract X, Y or Z
#' 
#' @keywords internal
#' @export
#' @param object A stanreg object.
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
