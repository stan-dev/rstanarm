# Default 'control' argument for stan() if none specified by user,
# value of adapt_delta depends on prior
default_stan_control <- function(prior, adapt_delta = NULL, 
                                 max_treedepth = 15L) {
  if (is.null(prior)) adapt_delta <- 0.95
  else if (is.null(adapt_delta)) {
    adapt_delta <- switch(prior$dist, 
                          "R2" = 0.99,
                          "hs" = 0.99,
                          "hs_plus" = 0.99,
                          "t" = if (any(prior$df <= 2)) 0.99 else 0.95,
                          0.95) # default
  }
  return(nlist(adapt_delta, max_treedepth))
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
  for (j in seq_along(user_dots)) {
    args[[unms[j]]] <- user_dots[[j]]
  }
  defaults <- default_stan_control(prior = prior, 
                                   adapt_delta = user_adapt_delta)
  
  if (!"control" %in% unms) { # no user-specified 'control' argument
    args$control <- defaults
  }
  else { # user specifies a 'control' argument
    if (!is.null(user_adapt_delta)) { 
      # if user specified adapt_delta argument to stan_* then 
      # set control$adapt_delta to user-specified value
      args$control$adapt_delta <- user_adapt_delta
    } else {
      # use default adapt_delta for the user's chosen prior
      args$control$adapt_delta <- defaults$adapt_delta
    }
    if (is.null(args$control$max_treedepth)) {
      # if user's 'control' has no max_treedepth set it to rstanarm default
      args$control$max_treedepth <- defaults$max_treedepth
    }
  }
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

# Test for a given estimation method
# @param x a stanreg object
used.optimizing <- function(x) {
  stopifnot(is.stanreg(x))
  x$algorithm == "optimizing"
}
used.sampling <- function(x) {
  stopifnot(is.stanreg(x))
  x$algorithm == "sampling"
}
# used.variational <- function(x) {
#   stopifnot(is.stanreg(x))
#   x$algorithm %in% c("meanfield", "fullrank")
# }

# Consistent error message to use when something is only available if
# algorithm="sampling"
# @param what An optional message to prepend to the default message.
STOP_sampling_only <- function(what) {
  msg <- "only available for models fit using MCMC (algorithm='sampling')."
  if (!missing(what)) msg <- paste(what, msg)
  stop(msg, call. = FALSE)
}

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
  has_name <- if (no_names) FALSE else nzchar(names(out))
  if (all(has_name)) return(out)
  nms <- as.character(m)[-1]
  if (no_names) names(out) <- nms
  else names(out)[!has_name] <- nms[!has_name]
  return(out)
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
  if (is.null(scale)) scale <- default
  if (link == "probit")
    return(scale * dnorm(0) / dlogis(0))
  else 
    return(scale)
}

# Make prior.info list
# @param user_call the user's call, i.e. match.call(expand.dots = TRUE)
# @param function_formals formal arguments of stan_* function, i.e. formals()
get_prior_info <- function(user_call, function_formals) {
  user <- grep("prior", names(user_call), value = TRUE)
  default <- setdiff(grep("prior", names(function_formals), value = TRUE), user)
  U <- length(user)
  D <- length(default)
  priors <- list()
  for (j in 1:(U+D)) {
    if (j <= U) priors[[user[j]]] <- eval(user_call[[user[j]]])
    else priors[[default[j-U]]] <- eval(function_formals[[default[j-U]]])
  }
  return(priors)
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
#' @templateVar stanregArg object
#' @template args-stanreg-object
#' @return For \code{get_x} and \code{get_z}, a matrix. For \code{get_y}, either
#'   a vector or a matrix, depending on how the response variable was specified.
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

linkinv <- function(x, ...) UseMethod("linkinv")
linkinv.stanreg <- function(x, ...) {
  if (is(x, "polr")) 
    return(polr_linkinv(x))
  else 
    return(family(x)$linkinv)
}
linkinv.family <- function(x, ...) {
  return(x$linkinv)
}
linkinv.character <- function(x, ...) {
  stopifnot(length(x) == 1)
  return(polr_linkinv(x))
}

# Make inverse link function for stan_polr models
# @param x A stanreg object or character scalar giving the "method"
polr_linkinv <- function(x) {
  if (is.stanreg(x) && is(x, "polr")) method <- x$method
  else if (is.character(x) && length(x) == 1) method <- x
  else stop("'x' should be a stanreg object created by stan_polr or a single string.")
  
  if (method == "logistic") make.link("logit")$linkinv
  else if (method == "loglog") pgumbel
  else make.link(method)$linkinv
}
