# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2015, 2016 Trustees of Columbia University
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
  
  if (!"control" %in% unms) { 
    # no user-specified 'control' argument
    args$control <- defaults
  } else { 
    # user specifies a 'control' argument
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
  args$save_warmup <- FALSE
  
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
  if (is.null(prior)) {
    adapt_delta <- 0.95
  } else if (is.null(adapt_delta)) {
    adapt_delta <- switch(prior$dist, 
                          "R2" = 0.99,
                          "hs" = 0.99,
                          "hs_plus" = 0.99,
                          "t" = if (any(prior$df <= 2)) 0.99 else 0.95,
                          0.95) # default
  }
  nlist(adapt_delta, max_treedepth)
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

# Test if stanreg object used stan_(g)lmer
#
# @param x A stanreg object.
is.mer <- function(x) {
  stopifnot(is.stanreg(x))
  check1 <- is(x, "lmerMod")
  check2 <- !is.null(x$glmod)
  if (check1 && !check2) {
    stop("Bug found. 'x' has class 'lmerMod' but no 'glmod' component.")
  } else if (!check1 && check2) {
    stop("Bug found. 'x' has 'glmod' component but not class 'lmerMod'.")
  }
  isTRUE(check1 && check2)
}

# Consistent error message to use when something is only available for 
# models fit using MCMC
#
# @param what An optional message to prepend to the default message.
STOP_sampling_only <- function(what) {
  msg <- "only available for models fit using MCMC (algorithm='sampling')."
  if (!missing(what)) 
    msg <- paste(what, msg)
  stop(msg, call. = FALSE)
}

# Consistent error message to use when something is only available for models
# fit using MCMC or VB but not optimization
# 
# @param what An optional message to prepend to the default message.
STOP_not_optimizing <- function(what) {
  msg <- "not available for models fit using algorithm='optimizing'."
  if (!missing(what)) 
    msg <- paste(what, msg)
  stop(msg, call. = FALSE)
}

# Message to issue when mean-field selected but 'QR=FALSE'. 
msg_meanfieldQR <- function() {
  message("Setting 'QR' to TRUE can often be helpful when ", 
          "using the 'meanfield' algorithm.",
          "\nSee the documentation for the 'QR' argument.")
}

# Issue warning if high rhat values
# 
# @param rhats Vector of rhat values.
# @param threshold Threshold value. If any rhat values are above threshold a 
#   warning is issued.
check_rhats <- function(rhats, threshold = 1.1) {
  if (any(rhats > threshold, na.rm = TRUE)) 
    warning("Markov chains did not converge! Do not analyze results!", 
            call. = FALSE, noBreaks. = TRUE)
}

# If y is a 1D array keep any names but convert to vector (used in stan_glm)
#
# @param y Result of calling model.response
array1D_check <- function(y) {
  if (length(dim(y)) == 1L) {
    nms <- rownames(y)
    dim(y) <- NULL
    if (!is.null(nms)) 
      names(y) <- nms
  }
  return(y)
}


# Check for a binomial model with Y given as proportion of successes and weights 
# given as total number of trials
# 
binom_y_prop <- function(y, family, weights) {
  if (!is.binomial(family$family)) 
    return(FALSE)

  yprop <- NCOL(y) == 1L && 
    is.numeric(y) && 
    any(y > 0 & y < 1) && 
    !any(y < 0 | y > 1)
  if (!yprop)
    return(FALSE)
  
  wtrials <- !identical(weights, double(0)) && 
    all(weights > 0) && 
    all(abs(weights - round(weights)) < .Machine$double.eps^0.5)
  isTRUE(wtrials)
}

# Convert 2-level factor to 0/1
fac2bin <- function(y) {
  if (!is.factor(y)) 
    stop("Bug found: non-factor as input to fac2bin.", 
         call. = FALSE)
  if (!identical(nlevels(y), 2L)) 
    stop("Bug found: factor with nlevels != 2 as input to fac2bin.", 
         call. = FALSE)
  as.integer(y != levels(y)[1L])
}

# Check weights argument
# 
# @param w The \code{weights} argument specified by user or the result of
#   calling \code{model.weights} on a model frame.
# @return If no error is thrown then \code{w} is returned.
validate_weights <- function(w) {
  if (missing(w) || is.null(w)) {
    w <- double(0)
  } else {
    if (!is.numeric(w)) 
      stop("'weights' must be a numeric vector.", 
           call. = FALSE)
    if (any(w < 0)) 
      stop("Negative weights are not allowed.", 
           call. = FALSE)
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
  if (is.null(o)) {
    o <- double(0)
  } else {
    if (length(o) != NROW(y))
      stop(gettextf("Number of offsets is %d but should be %d (number of observations)",
                    length(o), NROW(y)), domain = NA, call. = FALSE)
  }
  return(o)
}


# Check family argument
#
# @param f The \code{family} argument specified by user (or the default).
# @return If no error is thrown, then either \code{f} itself is returned (if
#   already a family) or the family object created from \code{f} is returned (if
#   \code{f} is a string or function).
validate_family <- function(f) {
  if (is.character(f)) 
    f <- get(f, mode = "function", envir = parent.frame(2))
  if (is.function(f)) 
    f <- f()
  if (!is(f, "family")) 
    stop("'family' must be a family.", call. = FALSE)
  
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
b_names <- function(x, ...) {
  grep("^b\\[", x, ...)
}

# Get the correct column name to use for selecting the median
#
# @param algorithm String naming the estimation algorithm (probably
#   \code{fit$algorithm}).
# @return Either \code{"50%"} or \code{"Median"} depending on \code{algorithm}.
select_median <- function(algorithm) {
  switch(algorithm, 
         sampling = "50%",
         meanfield = "50%",
         fullrank = "50%",
         optimizing = "Median",
         stop("Bug found (incorrect algorithm name passed to select_median)", 
              call. = FALSE))
}

# Regex parameter selection
#
# @param x stanreg object
# @param regex_pars Character vector of patterns
grep_for_pars <- function(x, regex_pars) {
  if (used.optimizing(x)) {
    warning("'regex_pars' ignored for models fit using algorithm='optimizing'.",
            call. = FALSE)
    return(NULL)
  }
  stopifnot(is.character(regex_pars))
  out <- unlist(lapply(seq_along(regex_pars), function(j) {
    grep(regex_pars[j], rownames(x$stan_summary), value = TRUE) 
  }))
  if (!length(out))
    stop("No matches for 'regex_pars'.", call. = FALSE)
  
  return(out)
}

# Combine pars and regex_pars
#
# @param x stanreg object
# @param pars Character vector of parameter names
# @param regex_pars Character vector of patterns
collect_pars <- function(x, pars = NULL, regex_pars = NULL) {
  if (is.null(pars) && is.null(regex_pars)) 
    return(NULL)
  if (!is.null(pars)) 
    pars[pars == "varying"] <- "b"
  if (!is.null(regex_pars)) 
    pars <- c(pars, grep_for_pars(x, regex_pars))
  unique(pars)
}

# Get the posterior sample size
#
# @param x A stanreg object
# @return NULL if used.optimizing(x), otherwise the posterior sample size
posterior_sample_size <- function(x) {
  stopifnot(is.stanreg(x))
  if (used.optimizing(x)) 
    return(NULL)
  pss <- x$stanfit@sim$n_save
  if (used.variational(x))
    return(pss)
  sum(pss - x$stanfit@sim$warmup2)
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
  if (!length(x)) {
    rep(0, times = n)
  } else if (length(x) == 1L) {
    rep(x, times = n)
  } else {
    x
  }
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
  if (all(has_name)) 
    return(out)
  nms <- as.character(m)[-1L]
  if (no_names) {
    names(out) <- nms
  } else {
    names(out)[!has_name] <- nms[!has_name]
  } 
  
  return(out)
}

# Check for positive scale or df parameter (NULL ok)
#
# @param x The value to check.
# @return Either an error is thrown or \code{TRUE} is returned invisibly.
validate_parameter_value <- function(x) {
  nm <- deparse(substitute(x))
  if (!is.null(x)) {
    if (!is.numeric(x)) 
      stop(nm, " should be NULL or numeric", call. = FALSE)
    if (any(x <= 0)) 
      stop(nm, " should be positive", call. = FALSE)
  }
  invisible(TRUE)
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
  if (is.null(scale)) 
    scale <- default
  if (link == "probit")
    scale <- scale * dnorm(0) / dlogis(0)
  
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
  default <- setdiff(grep("prior", names(function_formals), value = TRUE), 
                     user)
  U <- length(user)
  D <- length(default)
  priors <- list()
  for (j in 1:(U + D)) {
    if (j <= U) {
      priors[[user[j]]] <- eval(user_call[[user[j]]])
    } else {
      priors[[default[j-U]]] <- eval(function_formals[[default[j-U]]])
    } 
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
  if (length(offset))
    eta <- eta + offset
  
  return(eta)
}
linear_predictor.matrix <- function(beta, x, offset = NULL) {
  if (NCOL(beta) == 1L) 
    beta <- as.matrix(beta)
  eta <- beta %*% t(x)
  if (length(offset)) 
    eta <- sweep(eta, 2L, offset, `+`)

  return(eta)
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
  object[["y"]] %ORifNULL% model.response(model.frame(object))
}
#' @export
get_x.default <- function(object) {
  object[["x"]] %ORifNULL% model.matrix(object)
}
#' @export
get_x.lmerMod <- function(object) {
  object$glmod$X %ORifNULL% stop("X not found")
}
#' @export
get_z.lmerMod <- function(object) {
  Zt <- object$glmod$reTrms$Zt %ORifNULL% stop("Z not found")
  t(Zt)
}


# Get inverse link function
#
# @param x A stanreg object, family object, or string. 
# @param ... Ignored. 
# @return The inverse link function associated with x.
linkinv <- function(x, ...) UseMethod("linkinv")
linkinv.stanreg <- function(x, ...) {
  if (is(x, "polr")) polr_linkinv(x) else family(x)$linkinv
}
linkinv.family <- function(x, ...) {
  x$linkinv
}
linkinv.character <- function(x, ...) {
  stopifnot(length(x) == 1)
  polr_linkinv(x)
}

# Make inverse link function for stan_polr models, neglecting any
# exponent in the scobit case
#
# @param x A stanreg object or character scalar giving the "method".
# @return The inverse link function associated with x.
polr_linkinv <- function(x) {
  if (is.stanreg(x) && is(x, "polr")) {
    method <- x$method
  } else if (is.character(x) && length(x) == 1L) {
    method <- x
  } else {
    stop("'x' should be a stanreg object created by stan_polr ", 
         "or a single string.")
  }
  if (method == "logistic") 
    method <- "logit"
  
  if (method == "loglog") {
    pgumbel
  } else {
    make.link(method)$linkinv
  } 
}
