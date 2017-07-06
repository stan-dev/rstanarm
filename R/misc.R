# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2015, 2016, 2017 Trustees of Columbia University
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
  if (!length(prior)) {
    if (is.null(adapt_delta)) adapt_delta <- 0.95
  } else if (is.null(adapt_delta)) {
    adapt_delta <- switch(prior$dist, 
                          "R2" = 0.99,
                          "hs" = 0.99,
                          "hs_plus" = 0.99,
                          "lasso" = 0.99,
                          "product_normal" = 0.99,
                          0.95) # default
  }
  nlist(adapt_delta, max_treedepth)
}

# Test if an object is a stanreg object
#
# @param x The object to test. 
is.stanreg <- function(x) inherits(x, "stanreg")

# Throw error if object isn't a stanreg object
# 
# @param x The object to test.
validate_stanreg_object <- function(x, call. = FALSE) {
  if (!is.stanreg(x))
    stop("Object is not a stanreg object.", call. = call.) 
}

# Test for a given family
#
# @param x A character vector (probably x = family(fit)$family)
is.binomial <- function(x) x == "binomial"
is.gaussian <- function(x) x == "gaussian"
is.gamma <- function(x) x == "Gamma"
is.ig <- function(x) x == "inverse.gaussian"
is.nb <- function(x) x == "neg_binomial_2"
is.poisson <- function(x) x == "poisson"
is.beta <- function(x) x == "beta" || x == "Beta regression"

# test if a stanreg object has class polr 
is_polr <- function(object) {
  inherits(object, "polr")
}

# test if a stanreg object is a scobit model
is_scobit <- function(object) {
  validate_stanreg_object(object)
  if (!is(object, "polr")) return(FALSE)
  return("alpha" %in% rownames(object$stan_summary))
}

# Test for a given estimation method
#
# @param x A stanreg object.
used.optimizing <- function(x) {
  x$algorithm == "optimizing"
}
used.sampling <- function(x) {
  x$algorithm == "sampling"
}
used.variational <- function(x) {
  x$algorithm %in% c("meanfield", "fullrank")
}

# Test if stanreg object used stan_(g)lmer
#
# @param x A stanreg object.
is.mer <- function(x) {
  stopifnot(is.stanreg(x))
  check1 <- inherits(x, "lmerMod")
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

# Message to issue when fitting model with ADVI but 'QR=FALSE'. 
recommend_QR_for_vb <- function() {
  message(
    "Setting 'QR' to TRUE can often be helpful when using ", 
    "one of the variational inference algorithms. ", 
    "See the documentation for the 'QR' argument."
  )
}

# Issue warning if high rhat values
# 
# @param rhats Vector of rhat values.
# @param threshold Threshold value. If any rhat values are above threshold a 
#   warning is issued.
check_rhats <- function(rhats, threshold = 1.1, check_lp = FALSE) {
  if (!check_lp)
    rhats <- rhats[!names(rhats) %in% c("lp__", "log-posterior")]
  
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


# Check for glmer syntax in formulas for non-glmer models
#
# @param f The model \code{formula}.
# @return Nothing is returned but an error might be thrown
validate_glm_formula <- function(f) {
  if (any(grepl("\\|", f)))
    stop("Using '|' in model formula not allowed. ",
         "Maybe you meant to use 'stan_(g)lmer'?", call. = FALSE)
}


# Check if any variables in a model frame are constants
# @param mf A model frame or model matrix
# @return If no constant variables are found mf is returned, otherwise an error
#   is thrown.
check_constant_vars <- function(mf) {
  # don't check if columns are constant for binomial
  mf1 <- if (NCOL(mf[, 1]) == 2) mf[, -1, drop=FALSE] else mf
  
  lu <- function(x) length(unique(x))
  nocheck <- c("(weights)", "(offset)", "(Intercept)")
  sel <- !colnames(mf1) %in% nocheck
  is_constant <- apply(mf1[, sel, drop=FALSE], 2, lu) == 1
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

# Return names of the last dimension in a matrix/array (e.g. colnames if matrix)
#
# @param x A matrix or array
last_dimnames <- function(x) {
  ndim <- length(dim(x))
  dimnames(x)[[ndim]]
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
  validate_stanreg_object(x)
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
  validate_stanreg_object(x)
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


# Check and set scale parameters for priors
#
# @param scale Value of scale parameter (can be NULL).
# @param default Default value to use if \code{scale} is NULL.
# @param link String naming the link function or NULL.
# @return If a probit link is being used, \code{scale} (or \code{default} if
#   \code{scale} is NULL) is scaled by \code{dnorm(0) / dlogis(0)}. Otherwise
#   either \code{scale} or \code{default} is returned.
set_prior_scale <- function(scale, default, link) {
  stopifnot(is.numeric(default), is.character(link) || is.null(link))
  if (is.null(scale)) 
    scale <- default
  if (isTRUE(link == "probit"))
    scale <- scale * dnorm(0) / dlogis(0)
  
  return(scale)
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
#' @param ... Other arguments passed to methods. For a \code{stanmvreg} object
#'   this can be an integer \code{m} specifying the submodel.
#' @return For \code{get_x} and \code{get_z}, a matrix. For \code{get_y}, either
#'   a vector or a matrix, depending on how the response variable was specified.
get_y <- function(object, ...) UseMethod("get_y")
#' @rdname get_y
#' @export
get_x <- function(object, ...) UseMethod("get_x")
#' @rdname get_y
#' @export
get_z <- function(object, ...) UseMethod("get_z")

#' @export
get_y.default <- function(object, ...) {
  object[["y"]] %ORifNULL% model.response(model.frame(object))
}
#' @export
get_x.default <- function(object, ...) {
  object[["x"]] %ORifNULL% model.matrix(object)
}
#' @export
get_x.gamm4 <- function(object, ...) {
  as.matrix(object[["x"]])
}
#' @export
get_x.lmerMod <- function(object, ...) {
  object$glmod$X %ORifNULL% stop("X not found")
}
#' @export
get_z.lmerMod <- function(object, ...) {
  Zt <- object$glmod$reTrms$Zt %ORifNULL% stop("Z not found")
  t(Zt)
}
#' @export
get_y.stanmvreg <- function(object, m = NULL, ...) {
  ret <- lapply(object$glmod, lme4::getME, "y") %ORifNULL% stop("y not found")
  stub <- get_stub(object)
  if (!is.null(m)) ret[[m]] else list_nms(ret, stub = stub)
}
#' @export
get_x.stanmvreg <- function(object, m = NULL, ...) {
  ret <- lapply(object$glmod, lme4::getME, "X") %ORifNULL% stop("X not found")
  stub <- get_stub(object)
  if (!is.null(m)) ret[[m]] else list_nms(ret, stub = stub)
}
#' @export
get_z.stanmvreg <- function(object, m = NULL, ...) {
  ret <- lapply(object$glmod, lme4::getME, "Z") %ORifNULL% stop("Z not found")
  stub <- get_stub(object)
  if (!is.null(m)) ret[[m]] else list_nms(ret, stub = stub)
}

# Get inverse link function
#
# @param x A stanreg object, family object, or string. 
# @param ... Other arguments passed to methods. For a \code{stanmvreg} object
#   this can be an integer \code{m} specifying the submodel.
# @return The inverse link function associated with x.
linkinv <- function(x, ...) UseMethod("linkinv")
linkinv.stanreg <- function(x, ...) {
  if (is(x, "polr")) polr_linkinv(x) else family(x)$linkinv
}
linkinv.stanmvreg <- function(x, m = NULL, ...) {
  ret <- lapply(family(x), `[[`, "linkinv")
  stub <- get_stub(x)
  if (!is.null(m)) ret[[m]] else list_nms(ret, stub = stub)
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
  if (is.null(method) || method == "logistic") 
    method <- "logit"
  
  if (method == "loglog")
    return(pgumbel)
  
  make.link(method)$linkinv
}

# Wrapper for rstan::summary
# @param stanfit A stanfit object created using rstan::sampling or rstan::vb
# @return A matrix of summary stats
make_stan_summary <- function(stanfit) {
  levs <- c(0.5, 0.8, 0.95)
  qq <- (1 - levs) / 2
  probs <- sort(c(0.5, qq, 1 - qq))
  rstan::summary(stanfit, probs = probs, digits = 10)$summary  
}

check_reTrms <- function(reTrms) {
  stopifnot(is.list(reTrms))
  nms <- names(reTrms$cnms)
  dupes <- duplicated(nms)
  for (i in which(dupes)) {
    original <- reTrms$cnms[[nms[i]]]
    dupe <- reTrms$cnms[[i]]
    overlap <- dupe %in% original
    if (any(overlap))
      stop("rstanarm does not permit formulas with duplicate group-specific terms.\n", 
           "In this case ", nms[i], " is used as a grouping factor multiple times and\n",
           dupe[overlap], " is included multiple times.\n", 
           "Consider using || or -1 in your formulas to prevent this from happening.")
  }
  return(invisible(NULL))
}

#' @importFrom lme4 glmerControl
make_glmerControl <- function(...) {
  glmerControl(check.nlev.gtreq.5 = "ignore",
               check.nlev.gtr.1 = "stop",
               check.nobs.vs.rankZ = "ignore",
               check.nobs.vs.nlev = "ignore",
               check.nobs.vs.nRE = "ignore", ...)  
}

# Check if a fitted model (stanreg object) has weights
# 
# @param x stanreg object
# @return Logical. Only TRUE if x$weights has positive length and the elements
#   of x$weights are not all the same.
#
model_has_weights <- function(x) {
  wts <- x[["weights"]]
  if (!length(wts)) {
    FALSE
  } else if (all(wts == wts[1])) {
    FALSE
  } else {
    TRUE
  }
}


# Validate newdata argument for posterior_predict, log_lik, etc.
#
# Doesn't check if the correct variables are included (that's done in pp_data),
# just that newdata is either NULL or a data frame with no missing values.
# @param x User's 'newdata' argument
# @return Either NULL or a data frame
#
validate_newdata <- function(x) {
  if (is.null(x))
    return(NULL)
  if (!is.data.frame(x))
    stop("If 'newdata' is specified it must be a data frame.", call. = FALSE)
  if (any(is.na(x)))
    stop("NAs are not allowed in 'newdata'.", call. = FALSE)

  as.data.frame(x)
}

# Check that a stanfit object (or list returned by rstan::optimizing) is valid
#
check_stanfit <- function(x) {
  if (is.list(x)) {
    if (!all(c("par", "value") %in% names(x)))
      stop("Invalid object produced please report bug")
  }
  else {
    stopifnot(is(x, "stanfit"))
    if (x@mode != 0)
      stop("Invalid stanfit object produced please report bug")
  }
  return(TRUE)
}

# Validate data argument
#
# Make sure that, if specified, data is a data frame.
# 
# @param data User's data argument
# @param if_missing Object to return if data is missing/null
# @return If no error is thrown, data itself is returned if not missing/null, 
#   otherwise if_missing is returned.
# 
validate_data <- function(data, if_missing = NULL) {
  if (missing(data) || is.null(data)) {
    warn_data_arg_missing()
    return(if_missing)
  }
  if (!is.data.frame(data))
    stop("'data' must be a data frame.", call. = FALSE)
  
  return(data)
}

# Throw a warning if 'data' argument to modeling function is missing
warn_data_arg_missing <- function() {
  warning(
    "Omitting the 'data' argument is not recommended ",
    "and may not be allowed in future versions of rstanarm. ", 
    "Some post-estimation functions (in particular 'update', 'loo', 'kfold') ", 
    "are not guaranteed to work properly unless 'data' is specified as a data frame.",
    call. = FALSE
  )
}


#---------------------- for stan_jm only -----------------------------

# Test if family object corresponds to a linear mixed model
#
# @param x A family object
is.lmer <- function(x) {
  if (!is(x, "family"))
    stop("x should be a family object.", call. = FALSE)
  isTRUE((x$family == "gaussian") && (x$link == "identity"))
}

# Split a 2D array into nsplits subarrays, returned as a list
#
# @param x A 2D array or matrix
# @param nsplits An integer, the number of subarrays or submatrices
# @param bycol A logical, if TRUE then the subarrays are generated by
#   splitting the columns of x
# @return A list of nsplits arrays or matrices
array2list <- function(x, nsplits, bycol = TRUE) {
  len <- if (bycol) ncol(x) else nrow(x)
  len_k <- len %/% nsplits 
  if (!len == (len_k * nsplits))
    stop("Dividing x by nsplits does not result in an integer.")
  lapply(1:nsplits, function(k) {
    if (bycol) x[, (k-1) * len_k + 1:len_k, drop = FALSE] else
      x[(k-1) * len_k + 1:len_k, , drop = FALSE]})
}

# Convert a standardised quadrature node to an unstandardised value based on 
# the specified integral limits
#
# @param x An unstandardised quadrature node
# @param a The lower limit(s) of the integral, possibly a vector
# @param b The upper limit(s) of the integral, possibly a vector
unstandardise_quadpoints <- function(x, a, b) {
  if (!identical(length(x), 1L) || !is.numeric(x))
    stop("'x' should be a single numeric value.", call. = FALSE)
  if (!all(is.numeric(a), is.numeric(b)))
    stop("'a' and 'b' should be numeric.", call. = FALSE)
  if (!length(a) %in% c(1L, length(b)))
    stop("'a' and 'b' should be vectors of length 1, or, be the same length.", call. = FALSE)
  if (any((b - a) < 0))
    stop("The upper limits for the integral ('b' values) should be greater than ",
         "the corresponding lower limits for the integral ('a' values).", call. = FALSE)
  ((b - a) / 2) * x + ((b + a) / 2)
}

# Convert a standardised quadrature weight to an unstandardised value based on 
# the specified integral limits
#
# @param x An unstandardised quadrature weight
# @param a The lower limit(s) of the integral, possibly a vector
# @param b The upper limit(s) of the integral, possibly a vector
unstandardise_quadweights <- function(x, a, b) {
  if (!identical(length(x), 1L) || !is.numeric(x))
    stop("'x' should be a single numeric value.", call. = FALSE)
  if (!all(is.numeric(a), is.numeric(b)))
    stop("'a' and 'b' should be numeric.", call. = FALSE)
  if (!length(a) %in% c(1L, length(b)))
    stop("'a' and 'b' should be vectors of length 1, or, be the same length.", call. = FALSE)
  if (any((b - a) < 0))
    stop("The upper limits for the integral ('b' values) should be greater than ",
         "the corresponding lower limits for the integral ('a' values).", call. = FALSE)
  ((b - a) / 2) * x
}

# Test if object is stanmvreg class
#
# @param x An object to be tested.
is.stanmvreg <- function(x) {
  is(x, "stanmvreg")
}

# Test if object is a joint longitudinal and survival model
#
# @param x An object to be tested.
is.jm <- function(x) {
  x$stan_function == "stan_jm"
}

# Test if object contains a multivariate GLM
#
# @param x An object to be tested.
is.mvmer <- function(x) {
  x$stan_function %in% c("stan_mvmer", "stan_jm")
}

# Test if object contains a survival model
#
# @param x An object to be tested.
is.surv <- function(x) {
  x$stan_function %in% c("stan_jm")
}

# Throw error if object isn't a stanmvreg object
# 
# @param x The object to test.
validate_stanmvreg_object <- function(x, call. = FALSE) {
  if (!is.stanmvreg(x))
    stop("Object is not a stanmvreg object.", call. = call.) 
}

# Throw error if parameter isn't a positive scalar
#
# @param x The object to test.
validate_positive_scalar <- function(x, not_greater_than = NULL) {
  nm <- deparse(substitute(x))
  if (is.null(x))
    stop(nm, " cannot be NULL", call. = FALSE)
  if (!is.numeric(x))
    stop(nm, " should be numeric", call. = FALSE)
  if (any(x <= 0)) 
    stop(nm, " should be postive", call. = FALSE)
  if (!is.null(not_greater_than)) {
    if (!is.numeric(not_greater_than) || (not_greater_than <= 0))
      stop("'not_greater_than' should be numeric and postive")
    if (!all(x <= not_greater_than))
      stop(nm, " should less than or equal to ", not_greater_than, call. = FALSE)
  }
}

# Return a list with the median and prob% CrI bounds for each column of a 
# matrix or 2D array
#
# @param x A matrix or 2D array
# @param prob Value between 0 and 1 indicating the desired width of the CrI
median_and_bounds <- function(x, prob) {
  if (!any(is.matrix(x), is.array(x)))
    stop("x should be a matrix or 2D array.")
  med <- apply(x, 2, median)
  lb  <- apply(x, 2, quantile, (1 - prob)/2)
  ub  <- apply(x, 2, quantile, (1 + prob)/2)
  nlist(med, lb, ub)
}

# Return the stub for variable names from one submodel of a stan_jm model
#
# @param m An integer specifying the number of the longitudinal submodel or
#   a character string specifying the submodel (e.g. "Long1", "Event", etc)
# @param stub A character string to prefix to m, if m is supplied as an integer
get_m_stub <- function(m, stub = "Long") {
  if (is.null(m)) {
    return(NULL)
  } else if (is.numeric(m)) {
    return(paste0(stub, m, "|"))
  } else if (is.character(m)) {
    return(paste0(m, "|"))
  }
}

# Return the appropriate stub for variable names
#
# @param object A stanmvreg object
get_stub <- function(object) {
  if (is.jm(object)) "Long" else "y"  
} 

# Separates a names object into separate parts based on the longitudinal, 
# event, or association parameters.
# 
# @param x Character vector (often rownames(fit$stan_summary))
# @param M An integer specifying the number of longitudinal submodels.
# @param stub The character string used at the start of the names of variables
#   in the longitudinal/GLM submodels
# @param ... Arguments passed to grep
# @return A list with x separated out into those names corresponding
#   to parameters from the M longitudinal submodels, the event submodel
#   or association parameters.
collect_nms <- function(x, M, stub = "Long", ...) {
  ppd <- grep(paste0("^", stub, ".{1}\\|mean_PPD"), x, ...)      
  y <- lapply(1:M, function(m) grep(mod2rx(m, stub = stub), x, ...))
  y_extra <- lapply(1:M, function(m) 
    c(grep(paste0("^", stub, m, "\\|sigma"), x, ...),
      grep(paste0("^", stub, m, "\\|shape"), x, ...),
      grep(paste0("^", stub, m, "\\|lambda"), x, ...),
      grep(paste0("^", stub, m, "\\|reciprocal_dispersion"), x, ...)))             
  y <- lapply(1:M, function(m) setdiff(y[[m]], c(y_extra[[m]], ppd[m])))
  e <- grep(mod2rx("^Event"), x, ...)     
  e_extra <- c(grep("^Event\\|weibull-shape|^Event\\|basehaz-coef", x, ...))         
  e <- setdiff(e, e_extra)
  a <- grep(mod2rx("^Assoc"), x, ...)
  b <- b_names(x, ...)
  y_b <- lapply(1:M, function(m) b_names_M(x, m, stub = stub, ...))
  alpha <- grep("^.{5}\\|\\(Intercept\\)", x, ...)      
  alpha <- c(alpha, grep("^", stub, ".{1}\\|\\(Intercept\\)", x, ...))      
  beta <- setdiff(c(unlist(y), e, a), alpha)  
  nlist(y, y_extra, y_b, e, e_extra, a, b, alpha, beta, ppd) 
}

# Grep for "b" parameters (ranef), can optionally be specified
# for a specific longitudinal submodel
#
# @param x Character vector (often rownames(fit$stan_summary))
# @param submodel Optional integer specifying which long submodel
# @param ... Passed to grep
b_names_M <- function(x, submodel = NULL, stub = "Long", ...) {
  if (is.null(submodel)) {
    grep("^b\\[", x, ...)
  } else {
    grep(paste0("^b\\[", stub, submodel, "\\|"), x, ...)
  }
}

# Grep for regression coefs (fixef), can optionally be specified
# for a specific submodel
#
# @param x Character vector (often rownames(fit$stan_summary))
# @param submodel Character vector specifying which submodels
#   to obtain the coef names for. Can be "Long", "Event", "Assoc", or 
#   an integer specifying a specific longitudinal submodel. Specifying 
#   NULL selects all submodels.
# @param ... Passed to grep
beta_names <- function(x, submodel = NULL, ...) {
  if (is.null(submodel)) {
    rxlist <- c(mod2rx("^Long"), mod2rx("^Event"), mod2rx("^Assoc"))
  } else {
    rxlist <- c()
    if ("Long" %in% submodel) rxlist <- c(rxlist, mod2rx("^Long"))
    if ("Event" %in% submodel) rxlist <- c(rxlist, mod2rx("^Event"))
    if ("Assoc" %in% submodel) rxlist <- c(rxlist, mod2rx("^Assoc"))
    miss <- setdiff(submodel, c("Long", "Event", "Assoc"))
    if (length(miss)) rxlist <- c(rxlist, sapply(miss, mod2rx))
  }
  unlist(lapply(rxlist, function(y) grep(y, x, ...)))
}

# Converts "Long", "Event" or "Assoc" to the regular expression
# used at the start of variable names for the fitted joint model
#
# @param x The submodel for which the regular expression should be
#   obtained. Can be "Long", "Event", "Assoc", or an integer specifying
#   a specific longitudinal submodel.
mod2rx <- function(x, stub = "Long") {
  if (x == "^Long") {
    c("^Long[1-9]\\|")
  } else if (x == "^Event") {
    c("^Event\\|")
  } else if (x == "^Assoc") {
    c("^Assoc\\|")
  } else if (x == "Long") {
    c("Long[1-9]\\|")
  } else if (x == "Event") {
    c("Event\\|")
  } else if (x == "Assoc") {
    c("Assoc\\|")
  } else if (x == "^y") {
    c("^y[1-9]\\|")
  } else if (x == "y") {
    c("y[1-9]\\|")
  } else {
    paste0("^", stub, x, "\\|")
  }   
}

# Return the number of longitudinal submodels
#
# @param object A stanmvreg object
get_M <- function(object) {
  validate_stanmvreg_object(object)
  return(object$n_markers)
}

# Supplies names for the output list returned by most stanmvreg methods
#
# @param object The list object to which the names are to be applied
# @param M The number of longitudinal/GLM submodels. If NULL then the number of
#   longitudinal/GLM submodels is assumed to be equal to the length of object.
# @param stub The character string to use at the start of the names for
#   list items related to the longitudinal/GLM submodels
list_nms <- function(object, M = NULL, stub = "Long") {
  if (!is.list(object)) 
    stop("'object' argument should be a list")
  if (is.null(M)) M <- length(object)
  nms <- paste0(stub, 1:M)
  if (length(object) > M) nms <- c(nms, "Event")
  names(object) <- nms
  object
}

# Removes the submodel identifying text (e.g. "Long1|", "Event|", etc 
# from variable names
#
# @param x Character vector (often rownames(fit$stan_summary)) from which
#   the stub should be removed
rm_stub <- function(x) {
  x <- gsub(mod2rx("y"), "", x)
  x <- gsub(mod2rx("Long"), "", x)
  x <- gsub(mod2rx("Event"), "", x)
}

# Removes a specified character string from the names of an
# object (for example, a matched call)
#
# @param x The matched call
# @param string The character string to be removed
strip_nms <- function(x, string) {
  names(x) <- gsub(string, "", names(x))
  x
}

# Check argument contains one of the allowed options
check_submodelopt2 <- function(x) {
  if (!x %in% c("long", "event"))
    stop("submodel option must be 'long' or 'event'") 
}
check_submodelopt3 <- function(x) {
  if (!x %in% c("long", "event", "both"))
    stop("submodel option must be 'long', 'event' or 'both'") 
}

# Error message when not specifying an argument required for stanmvreg objects
#
# @param arg The argument
STOP_arg_required_for_stanmvreg <- function(arg) {
  nm <- deparse(substitute(arg))
  msg <- paste0("Argument '", nm, "' required for stanmvreg objects.")
  stop(msg, call. = FALSE)
}

# Error message when a function is not yet implemented for stanmvreg objects
#
# @param what A character string naming the function not yet implemented
STOP_if_stanmvreg <- function(what) {
  msg <- "not yet implemented for stanmvreg objects."
  if (!missing(what)) 
    msg <- paste(what, msg)
  stop(msg, call. = FALSE)
}

# Error message when a function is not yet implemented for stan_mvmer models
#
# @param what An optional message to prepend to the default message.
STOP_stan_mvmer <- function(what) {
  msg <- "is not yet implemented for models fit using stan_mvmer."
  if (!missing(what)) 
    msg <- paste(what, msg)
  stop(msg, call. = FALSE)
}

# Consistent error message to use when something that is only available for 
# models fit using stan_jm
#
# @param what An optional message to prepend to the default message.
STOP_jm_only <- function(what) {
  msg <- "can only be used with stan_jm models."
  if (!missing(what)) 
    msg <- paste(what, msg)
  stop(msg, call. = FALSE)
}

# Consistent error message when binomial models with greater than
# one trial are not allowed
#
STOP_binomial <- function() {
  stop("Binomial models with number of trials greater than one ",
       "are not allowed (i.e. only bernoulli models are allowed).", 
       call. = FALSE)
}

# Error message when a required variable is missing from the data frame
#
# @param var The name of the variable that could not be found
STOP_no_var <- function(var) {
  stop("Variable '", var, "' cannot be found in the data frame.", call. = FALSE)
}

# Check if individuals in ids argument were also used in model estimation
#
# @param object A stanmvreg object
# @param ids A vector of ids appearing in the pp data
# @param m Integer specifying which submodel to get the estimation IDs from
# @return A logical. TRUE indicates their are new ids in the prediction data,
#   while FALSE indicates all ids in the prediction data were used in fitting
#   the model. This return is used to determine whether to draw new b pars.
check_pp_ids <- function(object, ids, m = 1) {
  ids2 <- unique(model.frame(object, m = m)[[object$id_var]])
  if (any(ids %in% ids2))
    warning("Some of the IDs in the 'newdata' correspond to individuals in the ",
            "estimation dataset. Please be sure you want to obtain subject-",
            "specific predictions using the estimated random effects for those ",
            "individuals. If you instead meant to marginalise over the distribution ",
            "of the random effects (for posterior_predict or posterior_traj), or ",
            "to draw new random effects conditional on outcome data provided in ",
            "the 'newdata' arguments (for posterior_survfit), then please make ",
            "sure the ID values do not correspond to individuals in the ",
            "estimation dataset.", immediate. = TRUE)
  if (!all(ids %in% ids2)) TRUE else FALSE
}

# Validate newdataLong and newdataEvent arguments
#
# @param object A stanmvreg object
# @param newdataLong A data frame, or a list of data frames
# @param newdataEvent A data frame
# @param duplicate_ok A logical. If FALSE then only one row per individual is
#   allowed in the newdataEvent data frame
# @return A list of validated data frames
validate_newdatas <- function(object, newdataLong = NULL, newdataEvent = NULL,
                              duplicate_ok = FALSE) {
  validate_stanmvreg_object(object)
  id_var <- object$id_var
  newdatas <- list()
  if (!is.null(newdataLong)) {
    if (!is(newdataLong, "list"))
      newdataLong <- rep(list(newdataLong), get_M(object))
    newdatas <- c(newdatas, newdataLong)
  }
  if (!is.null(newdataEvent)) {
    if (!duplicate_ok && any(duplicated(newdataEvent[[id_var]])))
      stop("'newdataEvent' should only contain one row per individual, since ",
           "time varying covariates are not allowed in the prediction data.")
    newdatas <- c(newdatas, list(Event = newdataEvent))
  }
  if (length(newdatas)) {
    idvar_check <- sapply(newdatas, function(x) id_var %in% colnames(x)) 
    if (!all(idvar_check)) 
      STOP_no_var(id_var)
    ids <- lapply(newdatas, function(x) unique(x[[id_var]]))
    sorted_ids <- lapply(ids, sort)
    if (!length(unique(sorted_ids)) == 1L) 
      stop("The same subject ids should appear in each new data frame.")
    if (!length(unique(ids)) == 1L) 
      stop("The subject ids should be ordered the same in each new data frame.")  
    newdatas <- lapply(newdatas, validate_newdata)
    return(newdatas)
  } else return(NULL)
}

# Return data frames only including the specified subset of individuals
#
# @param object A stanmvreg object
# @param data A data frame, or a list of data frames
# @param ids A vector of ids indicating which individuals to keep
# @return A data frame, or a list of data frames, depending on the input
subset_ids <- function(object, data, ids) {
  validate_stanmvreg_object(object)
  id_var <- object$id_var
  is_list <- is(data, "list")
  if (!is_list) data <- list(data)
  is_df <- sapply(data, is.data.frame)
  if (!all(is_df)) stop("'data' should be a data frame, or list of data frames.")
  data <- lapply(data, function(x) {
    if (!id_var %in% colnames(x)) STOP_no_var(id_var)
    sel <- which(!ids %in% x[[id_var]])
    if (length(sel)) 
      stop("The following 'ids' do not appear in the data: ", 
           paste(ids[[sel]], collapse = ", "))
    x[x[[id_var]] %in% ids, , drop = FALSE]
  })
  if (is_list) return(data) else return(data[[1]])
}

# Return an array or list with the time sequence used for posterior predictions
#
# @param increments An integer with the number of increments (time points) at
#   which to predict the outcome for each individual
# @param t0,t1 Numeric vectors giving the start and end times across which to
#   generate prediction times
# @param simplify Logical specifying whether to return each increment as a 
#   column of an array (TRUE) or as an element of a list (FALSE) 
get_time_seq <- function(increments, t0, t1, simplify = TRUE) {
  val <- sapply(0:(increments - 1), function(x, t0, t1) {
    t0 + (t1 - t0) * (x / (increments - 1))
  }, t0 = t0, t1 = t1, simplify = simplify)
  if (simplify && is.vector(val)) {
    # need to transform if there is only one individual
    val <- t(val)
    rownames(val) <- if (!is.null(names(t0))) names(t0) else 
      if (!is.null(names(t1))) names(t1) else NULL
  }
  return(val)
}

# Extract parameters from stanmat and return as a list
# 
# @param object A stanmvreg object
# @param stanmat A matrix of posterior draws, may be provided if the desired 
#   stanmat is only a subset of the draws from as.matrix(object$stanfit)
# @return A named list
extract_pars <- function(object, stanmat = NULL, means = FALSE) {
  validate_stanmvreg_object(object)
  M <- get_M(object)
  if (is.null(stanmat)) 
    stanmat <- as.matrix(object$stanfit)
  if (means) 
    stanmat <- t(colMeans(stanmat)) # return posterior means
  nms   <- collect_nms(colnames(stanmat), M, stub = get_stub(object))
  beta  <- lapply(1:M, function(m) stanmat[, nms$y[[m]], drop = FALSE])
  ebeta <- stanmat[, nms$e, drop = FALSE]
  abeta <- stanmat[, nms$a, drop = FALSE]
  bhcoef <- stanmat[, nms$e_extra, drop = FALSE]
  b     <- lapply(1:M, function(m) stanmat[, nms$y_b[[m]], drop = FALSE])
  nlist(beta, ebeta, abeta, bhcoef, b, stanmat)
}

# Return a data frame for each submodel that:
# (1) only includes variables used in the model formula
# (2) only includes rows contained in the actual glmod/coxmod model frames
# (3) ensures the ID variable is included in every model frame
#
# @param object A stanmvreg object
# @param m Integer specifying which submodel
get_model_data <- function(object, m = NULL) {
  validate_stanmvreg_object(object)
  forms <- formula(object)
  RHS <- paste(deparse(forms$Event[[3L]]), "+", object$id_var)
  forms$Event <- reformulate(RHS, response = forms$Event[[2L]])
  datas <- c(object$dataLong, list(object$dataEvent))
  row_nms <- lapply(model.frame(object), rownames)
  mfs <- mapply(function(w, x, y) {
    get_all_vars(w, x)[y, , drop = FALSE]
  }, w = forms, x = datas, y = row_nms, SIMPLIFY = FALSE)
  if (is.null(m)) return(mfs) else return(mfs[[m]])
}
