# Part of the rstanjm package
# Copyright (C) 2015, 2016 Trustees of Columbia University
# Copyright (C) 2016 Sam Brilleman
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

# Test if family object corresponds to a linear mixed model
#
# @param x A family object
is.lmer <- function(x) {
  if (!is(x, "family"))
    stop("x should be a family object.", call. = FALSE)
  isTRUE((x$family == "gaussian") && (x$link == "identity"))
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
  if (!identical(length(a), length(b)))
    stop("'a' and 'b' should be vectors of the same length.", call. = FALSE)
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
  if (!identical(length(a), length(b)))
    stop("'a' and 'b' should be vectors of the same length.", call. = FALSE)
  if (any((b - a) < 0))
    stop("The upper limits for the integral ('b' values) should be greater than ",
         "the corresponding lower limits for the integral ('a' values).", call. = FALSE)
  ((b - a) / 2) * x
}

# Test if object is stanjm class
#
# @param x An object to be tested.
is.stanjm <- function(x) {
  is(x, "stanjm")
}

# Throw error if object isn't a stanjm object
# 
# @param x The object to test.
validate_stanjm_object <- function(x, call. = FALSE) {
  if (!is.stanjm(x))
    stop("Object is not a stanjm object.", call. = call.) 
}

# Throw error if parameter isn't a positive scalar
#
# @param x The object to test.
validate_positive_scalar <- function(x) {
  nm <- deparse(substitute(x))
  if (is.null(x))
    stop(nm, " cannot be NULL", call. = FALSE)
  if (!is.numeric(x))
    stop(nm, " should be numeric", call. = FALSE)
  if (any(x <= 0)) 
    stop(nm, " should be postive", call. = FALSE)
  invisible(TRUE)
}

# Separates a names object into separate parts based on the longitudinal, 
# event, or association parameters.
# 
# @param x Character vector (often rownames(fit$stan_summary))
# @param M An integer specifying the number of longitudinal submodels.
# @param ... Arguments passed to grep
# @return A list with x separated out into those names corresponding
#   to parameters from the M longitudinal submodels, the event submodel
#   or association parameters.
collect_nms <- function(x, M, ...) {
  y <- lapply(1:M, function(m) grep(mod2rx(m), x, ...))
  y_extra <- lapply(1:M, function(m) 
    c(grep(paste0("^Long", m, "\\|sigma"), x, ...),
      grep(paste0("^Long", m, "\\|shape"), x, ...),
      grep(paste0("^Long", m, "\\|lambda"), x, ...),
      grep(paste0("^Long", m, "\\|reciprocal_dispersion"), x, ...)))             
  y <- lapply(1:M, function(m) setdiff(y[[m]], y_extra[[m]]))
  e <- grep(mod2rx("^Event"), x, ...)     
  e_extra <- c(grep("^Event\\|weibull-shape|^Event\\|basehaz-coef", x, ...))         
  e <- setdiff(e, e_extra)
  a <- grep(mod2rx("^Assoc"), x, ...)
  b <- b_names(x, ...)
  y_b <- lapply(1:M, function(m) b_names_M(x, m, ...))
  alpha <- grep("^.{5}\\|\\(Intercept\\)", x, ...)      
  beta <- setdiff(c(unlist(y), e, a), alpha)  
  nlist(y, y_extra, y_b, e, e_extra, a, b, alpha, beta) 
}

# Grep for "b" parameters (ranef), can optionally be specified
# for a specific longitudinal submodel
#
# @param x Character vector (often rownames(fit$stan_summary))
# @param submodel Optional integer specifying which long submodel
# @param ... Passed to grep
b_names_M <- function(x, submodel = NULL, ...) {
  if (is.null(submodel)) {
    grep("^b\\[", x, ...)
  } else {
    grep(paste0("^b\\[Long", submodel, "\\|"), x, ...)
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
mod2rx <- function(x) {
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
  } else {
    paste0("^Long", x, "\\|")
  }   
}

# Supplies names for the output list returned by most stanjm methods
#
# @param x The list object to which the names are to be applied
# @param M The number of longitudinal submodels
list_nms <- function(object, M) {
  if (!is.list(object)) stop("'object' argument should be a list")
  nms <- paste0("Long", 1:M)
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

# Get inverse link function
#
# @param x A stanjm object or family object. 
# @param ... Ignored. 
# @return The inverse link function associated with x.
linkinv.stanjm <- function(x, ...) {
  lapply(family(x), function(y) y$linkinv)
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