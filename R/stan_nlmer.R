# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2016 Trustees of Columbia University
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

#' Bayesian nonlinear models with group-specific terms via Stan
#' 
#' \if{html}{\figure{stanlogo.png}{options: width="25px" alt="http://mc-stan.org/about/logo/"}}
#' Bayesian inference for NLMMs with group-specific coefficients that have 
#' unknown covariance matrices with flexible priors.
#'
#' @export
#' @templateVar fun stan_nlmer
#' @templateVar pkg lme4
#' @templateVar pkgfun nlmer
#' @template return-stanreg-object
#' @template see-also
#' @template args-dots
#' @template args-prior_aux
#' @template args-priors
#' @template args-prior_PD
#' @template args-algorithm
#' @template args-adapt_delta
#' @template args-sparse
#' @template args-QR
#' 
#' @param formula,data Same as for \code{\link[lme4]{nlmer}}. \emph{We strongly
#'   advise against omitting the \code{data} argument}. Unless \code{data} is
#'   specified (and is a data frame) many post-estimation functions (including
#'   \code{update}, \code{loo}, \code{kfold}) are not guaranteed to work
#'   properly.
#' @param subset,weights,offset Same as \code{\link[stats]{glm}}.
#' @param na.action,contrasts Same as \code{\link[stats]{glm}}, but rarely 
#'   specified.
#' @param prior_covariance Cannot be \code{NULL}; see \code{\link{decov}} for
#'   more information about the default arguments.
#'
#' @details The \code{stan_nlmer} function is similar in syntax to 
#'   \code{\link[lme4]{nlmer}} but rather than performing (approximate) maximum 
#'   marginal likelihood estimation, Bayesian estimation is by default performed
#'   via MCMC. The Bayesian model adds independent priors on the "coefficients"
#'   --- which are really intercepts --- in the same way as 
#'   \code{\link{stan_nlmer}} and priors on the terms of a decomposition of the 
#'   covariance matrices of the group-specific parameters. See
#'   \code{\link{priors}} for more information about the priors.
#'   
#'   The supported transformation functions are limited to the named 
#'   "self-starting" functions in the \pkg{stats} library:
#'   \code{\link[stats]{SSasymp}}, \code{\link[stats]{SSasympOff}},
#'   \code{\link[stats]{SSasympOrig}}, \code{\link[stats]{SSbiexp}},
#'   \code{\link[stats]{SSfol}}, \code{\link[stats]{SSfpl}},
#'   \code{\link[stats]{SSgompertz}}, \code{\link[stats]{SSlogis}},
#'   \code{\link[stats]{SSmicmen}}, and \code{\link[stats]{SSweibull}}.
#'   
#'   
#' @seealso The vignette for \code{stan_glmer}, which also discusses 
#'   \code{stan_nlmer} models. \url{http://mc-stan.org/rstanarm/articles/}
#'   
#' @examples
#' \donttest{
#' data("Orange", package = "datasets")
#' Orange$circumference <- Orange$circumference / 100
#' Orange$age <- Orange$age / 100
#' fit <- stan_nlmer(
#'   circumference ~ SSlogis(age, Asym, xmid, scal) ~ Asym|Tree, 
#'   data = Orange, 
#'   # for speed only
#'   chains = 1, 
#'   iter = 1000
#'  ) 
#' print(fit)
#' posterior_interval(fit)
#' plot(fit, regex_pars = "b\\[")
#' }
#' @importFrom lme4 nlformula
#' @importFrom stats getInitial
stan_nlmer <-
  function(formula,
           data = NULL,
           subset,
           weights,
           na.action,
           offset,
           contrasts = NULL,
           ...,
           prior = normal(autoscale=TRUE),
           prior_aux = exponential(autoscale=TRUE),
           prior_covariance = decov(),
           prior_PD = FALSE,
           algorithm = c("sampling", "meanfield", "fullrank"),
           adapt_delta = NULL,
           QR = FALSE,
           sparse = FALSE) {

  if (!has_outcome_variable(formula[[2]])) {
    stop("LHS of formula must be specified.")
  }
  f <- as.character(formula[-3])
  SSfunctions <- grep("^SS[[:lower:]]+", ls("package:stats"), value = TRUE) 
  SSfun <- sapply(SSfunctions, function(ss) 
    grepl(paste0(ss, "("), x = f[2], fixed = TRUE))
  if (!any(SSfun)) {
    stop("'stan_nlmer' requires a named self-starting nonlinear function.")
  }
  SSfun <- which(SSfun)
  SSfun_char <- names(SSfun)
  
  mc <- match.call(expand.dots = FALSE)
  mc$prior <- mc$prior_aux <- mc$prior_covariance <- mc$prior_PD <-
    mc$algorithm <- mc$adapt_delta <- mc$QR <- mc$sparse <- NULL
  mc$start <-
    unlist(getInitial(
      object = as.formula(f[-1]),
      data = data,
      control = list(maxiter = 0, warnOnly = TRUE)
    ))
  
  nlf <- nlformula(mc)
  X <- nlf$X
  y <- nlf$respMod$y
  weights <- nlf$respMod$weights
  offset  <- nlf$respMod$offset
  
  nlf$reTrms$SSfun <- SSfun
  nlf$reTrms$decov <- prior_covariance
  
  nlf_inputs <- parse_nlf_inputs(nlf$respMod)
  if (SSfun_char == "SSfol") {
    nlf$reTrms$Dose <- nlf$frame[[nlf_inputs[2]]]
    nlf$reTrms$input <- nlf$frame[[nlf_inputs[3]]]
  } else {
    nlf$reTrms$input <- nlf$frame[[nlf_inputs[2]]]
  }
  
  
  algorithm <- match.arg(algorithm)
  stanfit <- stan_glm.fit(x = X, y = y, family = gaussian(link = "identity"),
                          weights = weights, offset = offset,
                          prior = prior, prior_intercept = NULL,
                          prior_aux = prior_aux, prior_PD = prior_PD, 
                          algorithm = algorithm, adapt_delta = adapt_delta,
                          group = nlf$reTrms, QR = QR, sparse = sparse, ...)
  if (algorithm != "optimizing" && !is(stanfit, "stanfit")) {
    return(stanfit)
  }
  
  if (SSfun_char == "SSfpl") { # SSfun = 6
    stanfit@sim$samples <- lapply(stanfit@sim$samples, FUN = function(x) {
      x[[4L]] <- exp(x[[4L]])
      return(x)
    })
  } else if (SSfun_char == "SSlogis") { # SSfun = 8
    stanfit@sim$samples <- lapply(stanfit@sim$samples, FUN = function(x) {
      x[[3L]] <- exp(x[[3L]])
      return(x)
    })
  }
  
  Z <- pad_reTrms(Ztlist = nlf$reTrms$Ztlist, cnms = nlf$reTrms$cnms, 
                  flist = nlf$reTrms$flist)$Z
  colnames(Z) <- b_names(names(stanfit), value = TRUE)

  fit <- nlist(stanfit, 
               family = make_nlf_family(SSfun_char, nlf), 
               formula, offset, weights, 
               x = cbind(X, Z), y = y, data, call = match.call(), terms = NULL, 
               model = NULL, na.action = na.omit, contrasts, algorithm, 
               glmod = nlf, stan_function = "stan_nlmer")
  out <- stanreg(fit)
  class(out) <- c(class(out), "nlmerMod", "lmerMod")
  return(out)
}


# internal ----------------------------------------------------------------

# @param respMod The respMod slot of the object returned by nlformula
# @return A character vector, the first element of which is the name of the SS
#   function and the rest of the elements are the names of the arguments to the
#   SS function
parse_nlf_inputs <- function(respMod) {
  inputs <- as.character(respMod$nlmod[2])
  inputs <- sub("(", ",", inputs, fixed = TRUE)
  inputs <- sub(")", "", inputs, fixed = TRUE)
  scan(
    text = inputs,
    what = character(),
    sep = ",",
    strip.white = TRUE,
    quiet = TRUE
  )
}

# Make family object 
# 
# @param SSfun_char SS function name as a string
# @param nlf Object returned by nlformula
# @return A family object
make_nlf_family <- function(SSfun_char, nlf) {
  g <- gaussian(link = "identity")
  g$link <- paste("inv", SSfun_char, sep = "_")
  g$linkinv <- function(eta, arg1, arg2 = NULL, FUN = SSfun_char) {
    if (is.matrix(eta)) {
      len <- length(arg1)
      nargs <- ncol(eta) / len
      SSargs <- lapply(1:nargs, FUN = function(i) {
        start <- 1 + (i - 1) * len
        end <- i * len
        t(eta[, start:end, drop = FALSE])
      })
      if (is.null(arg2)) SSargs <- c(list(arg1), SSargs)
      else SSargs <- c(list(arg1, arg2), SSargs)
    } else {
      SSargs <- as.data.frame(matrix(eta, nrow = length(arg1))) 
      if (is.null(arg2)) SSargs <- cbind(arg1, SSargs)
      else SSargs <- cbind(arg1, arg2, SSargs)
    }
    names(SSargs) <- names(formals(FUN))
    if (FUN == "SSbiexp") 
      SSargs$A1 <- SSargs$A1 + exp(SSargs$A2)
    
    do.call(FUN, args = SSargs)
  }
  
  nlf_inputs <- parse_nlf_inputs(nlf$respMod)
  if (SSfun_char == "SSfol") {
    formals(g$linkinv)$arg1 <- nlf$frame[[nlf_inputs[2]]]
    formals(g$linkinv)$arg2 <- nlf$frame[[nlf_inputs[3]]]
  } else {
    formals(g$linkinv)$arg1 <- nlf$frame[[nlf_inputs[2]]]
  }
  
  g$linkfun  <- function(mu) stop("'linkfun' should not have been called")
  g$variance <- function(mu) stop("'variance' should not have been called")
  g$mu.eta   <- function(mu) stop("'mu.eta' should not have been called")
  return(g)
}
