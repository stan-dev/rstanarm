# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2015 Trustees of Columbia University
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

# Set arguments for sampling 
#
# Prepare a list of arguments to use with \code{rstan::sampling} via
# \code{do.call}.
#
# @param object The stanfit object to use for sampling.
# @param user_dots The contents of \code{...} from the user's call to
#   the \code{stan_*} modeling function.
# @param user_adapt_delta The value for \code{adapt_delta} specified by the
#   user.
# @param prior Prior distribution list (can be NULL).
# @param ... Other arguments to \code{\link[rstan]{sampling}} not coming from
#   \code{user_dots} (e.g. \code{data}, \code{pars}, \code{init}, etc.)
# @return A list of arguments to use for the \code{args} argument for 
#   \code{do.call(sampling, args)}.
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

# Default control arguments for sampling
# 
# Called by set_sampling_args to set the default 'control' argument for
# \code{rstan::sampling} if none specified by user. This allows the value of
# \code{adapt_delta} to depend on the prior.
# 
# @param prior Prior distribution list (can be NULL).
# @param adapt_delta User's \code{adapt_delta} argument.
# @param max_treedepth Default for \code{max_treedepth}.
# @return A list with \code{adapt_delta} and \code{max_treedepth}.
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

# Test if an object is a stanreg object
#
# @param x The object to test. 
is.stanreg <- function(x) inherits(x, "stanreg")

# Test for a given family
#
# @param x A character vector (probably x = family(fit)$family)
is.binomial <- function(x) x == "binomial"
is.gaussian <- function(x) x == "gaussian"
is.gamma <- function(x) x == "Gamma"
is.ig <- function(x) x == "inverse.gaussian"
is.nb <- function(x) x == "neg_binomial_2"
is.poisson <- function(x) x == "poisson"

# Test for a given estimation method
#
# @param x A stanreg object.
used.optimizing <- function(x) {
  stopifnot(is.stanreg(x))
  x$algorithm == "optimizing"
}
used.sampling <- function(x) {
  stopifnot(is.stanreg(x))
  x$algorithm == "sampling"
}
used.variational <- function(x) {
  stopifnot(is.stanreg(x))
  x$algorithm %in% c("meanfield", "fullrank")
}

# Consistent error message to use when something is only available for 
# models fit using MCMC
#
# @param what An optional message to prepend to the default message.
STOP_sampling_only <- function(what) {
  msg <- "only available for models fit using MCMC (algorithm='sampling')."
  if (!missing(what)) msg <- paste(what, msg)
  stop(msg, call. = FALSE)
}

# Consistent error message to use when something is only available for models
# fit using MCMC or VB but not optimization
# 
# @param what An optional message to prepend to the default message.
STOP_not_optimizing <- function(what) {
  msg <- "not available for models fit using algorithm='optimizing'."
  if (!missing(what)) msg <- paste(what, msg)
  stop(msg, call. = FALSE)
}

# Check weights argument
# 
# @param w The \code{weights} argument specified by user or the result of
#   calling \code{model.weights} on a model frame.
# @return If no error is thrown then \code{w} is returned.
validate_weights <- function(w) {
  if (missing(w) || is.null(w)) return(double(0))
  else {
    if (!is.numeric(w)) stop("'weights' must be a numeric vector.")
    if (any(w < 0)) stop("Negative weights are not allowed.")
  }
  return(w)
}

# Check offset argument
#
# @param o The \code{offset} argument specified by user or the result of calling
#   \code{model.offset} on a model frame.
# @param y The result of calling \code{model.response} on a model frame.
# @return If no error is thrown then \code{o} is returned.
validate_offset <- function(o, y) {
  if (is.null(o)) return(double(0))
  if (length(o) != NROW(y))
    stop(gettextf("Number of offsets is %d but should be %d (number of observations)",
                  length(o), NROW(y)), domain = NA)
  return(o)
}


# Check family argument
#
# @param f The \code{family} argument specified by user (or the default).
# @return If no error is thrown, then either \code{f} itself is returned (if
#   already a family) or the family object created from \code{f} is returned (if
#   \code{f} is a string or function).
validate_family <- function(f) {
  if (is.character(f)) f <- get(f, mode = "function", envir = parent.frame(2))
  if (is.function(f)) f <- f()
  if (!is(f, "family")) stop("'family' must be a family.")
  return(f)
}

# Check if any variables in a model frame are constants
# @param mf A model frame or model matrix
# @return If no constant variables are found mf is returned, otherwise an error
#   is thrown.
check_constant_vars <- function(mf) {
  lu <- function(x) length(unique(x))
  nocheck <- c("(weights)", "(offset)", "(Intercept)")
  sel <- !colnames(mf) %in% nocheck
  is_constant <- apply(mf[, sel, drop = FALSE], 2, lu) == 1
  if (any(is_constant)) 
    stop("Constant variable(s) found: ", 
         paste(names(is_constant)[is_constant], collapse = ", "), 
         call. = FALSE)
  return(mf)
}


# Grep for "b" parameters (ranef)
#
# @param x Character vector (often rownames(fit$stan_summary))
# @param ... Passed to grep
.bnames <- function(x, ...) {
  grep("^b\\[", x, ...)
}

# Get the correct column name to use for selecting the median
#
# @param algorithm String naming the estimation algorithm (probably
#   \code{fit$algorithm}).
# @return Either \code{"50%"} or \code{"Median"} depending on \code{algorithm}.
.select_median <- function(algorithm) {
  switch(algorithm, 
         sampling = "50%",
         meanfield = "50%",
         fullrank = "50%",
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

# Maybe broadcast 
#
# @param x A vector or scalar.
# @param n Number of replications to possibly make. 
# @return If \code{x} has no length the \code{0} replicated \code{n} times is
#   returned. If \code{x} has length 1, the \code{x} replicated \code{n} times
#   is returned. Otherwise \code{x} itself is returned.
maybe_broadcast <- function(x, n) {
  if (!length(x)) rep(0, times = n)
  else if (length(x) == 1L) rep(x, times = n)
  else x
}

# Create a named list using specified names or, if names are omitted, using the
# names of the objects in the list
#
# @param ... Objects to include in the list. 
# @return A named list.
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
#
# @param x The value to check.
# @return Either an error is thrown or \code{TRUE} is returned invisibly.
validate_parameter_value <- function(x) {
  nm <- deparse(substitute(x))
  if (!is.null(x)) {
    if (!is.numeric(x)) stop(nm, " should be NULL or numeric")
    if (any(x <= 0)) stop(nm, " should be positive")
  }
  return(invisible(TRUE))
}

# Check and set scale parameters for priors
#
# @param scale Value of scale parameter (can be NULL).
# @param default Default value to use if \code{scale} is NULL.
# @param link String naming the link function.
# @return If a probit link is being used, \code{scale} (or \code{default} if
#   \code{scale} is NULL) is scaled by \code{dnorm(0) / dlogis(0)}. Otherwise
#   either \code{scale} or \code{default} is returned.
set_prior_scale <- function(scale, default, link) {
  stopifnot(is.numeric(default), is.character(link))
  if (is.null(scale)) scale <- default
  if (link == "probit")
    return(scale * dnorm(0) / dlogis(0))
  else 
    return(scale)
}

# Make prior.info list
# @param user_call The user's call, i.e. match.call(expand.dots = TRUE).
# @param function_formals Formal arguments of the stan_* function, i.e.
#   formals().
# @return A list containing information about the prior distributions and
#   options used.
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


# Methods for creating linear predictor
#
# Make linear predictor vector from x and point estimates for beta, or linear 
# predictor matrix from x and full posterior sample of beta.
#
# @param beta A vector or matrix or parameter estimates.
# @param x Predictor matrix.
# @param offset Optional offset vector.
# @return A vector or matrix.
linear_predictor <- function(beta, x, offset = NULL) {
  UseMethod("linear_predictor")
}
linear_predictor.default <- function(beta, x, offset = NULL) {
  eta <- as.vector(if (NCOL(x) == 1L) x * beta else x %*% beta)
  if (!length(offset)) return(eta)
  else return(eta + offset)
}
linear_predictor.matrix <- function(beta, x, offset = NULL) {
  if (NCOL(beta) == 1L) 
    beta <- as.matrix(beta)
  eta <- beta %*% t(x)
  if (!length(offset)) return(eta)
  else return(sweep(eta, MARGIN = 2L, offset, `+`))
}


#' Extract X, Y or Z from a stanreg object
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


# Get inverse link function
#
# @param x A stanreg object, family object, or string. 
# @param ... Ignored. 
# @return The inverse link function associated with x.
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
#
# @param x A stanreg object or character scalar giving the "method".
# @return The inverse link function associated with x.
polr_linkinv <- function(x) {
  if (is.stanreg(x) && is(x, "polr")) method <- x$method
  else if (is.character(x) && length(x) == 1) method <- x
  else stop("'x' should be a stanreg object created by stan_polr ", 
            "or a single string.")
  
  if (method == "logistic") make.link("logit")$linkinv
  else if (method == "loglog") pgumbel
  else make.link(method)$linkinv
}
