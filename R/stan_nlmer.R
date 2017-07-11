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
#' 
#' @param formula,data Same as for \code{\link[lme4]{nlmer}}.
#' @param subset,weights,offset Same as \code{\link[stats]{glm}}.
#' @param na.action,contrasts Same as \code{\link[stats]{glm}}, but rarely 
#'   specified.
#' @param prior_covariance Cannot be \code{NULL}; see \code{\link{decov}} for
#'   more information about the default arguments.
#'
#' @details The \code{stan_nlmer} function is similar in syntax to 
#'   \code{\link[lme4]{nlmer}} but rather than performing (approximate) maximum 
#'   marginal likelihood estimation, Bayesian estimation is by default
#'   performed via MCMC. The Bayesian model adds independent priors on the 
#'   "coefficients" --- which are really intercepts --- in the same way as 
#'   \code{\link{stan_nlmer}} and priors on the terms of a decomposition of the 
#'   covariance matrices of the group-specific parameters. See \code{\link{priors}} 
#'   for more information about the priors.
#'   
#'   The supported transformation functions are limited to the named 
#'   "self-starting" functions in the stats library: \code{\link[stats]{SSasymp}},
#'   \code{\link[stats]{SSasympOff}}, \code{\link[stats]{SSasympOrig}},
#'   \code{\link[stats]{SSbiexp}}, \code{\link[stats]{SSfol}}, 
#'   \code{\link[stats]{SSfpl}}, \code{\link[stats]{SSgompertz}},
#'   \code{\link[stats]{SSlogis}}, \code{\link[stats]{SSmicmen}}, and
#'   \code{\link[stats]{SSweibull}}.
#' @examples
#' \donttest{
#' data("Orange", package = "datasets")
#' Orange$circumference <- Orange$circumference / 100
#' Orange$age <- Orange$age / 100
#' stan_nlmer(circumference ~ SSlogis(age, Asym, xmid, scal) ~ Asym|Tree, 
#'            data = Orange)
#' }
#' @importFrom lme4 nlformula
#' @importFrom stats getInitial
stan_nlmer <- function (formula, data = NULL, subset, weights, na.action, offset, 
                        contrasts = NULL, ..., 
                        prior = normal(), prior_aux = cauchy(0, 5),
                        prior_covariance = decov(), prior_PD = FALSE, 
                        algorithm = c("sampling", "meanfield", "fullrank"), 
                        adapt_delta = NULL, sparse = FALSE) {
  f <- as.character(formula[-3])
  SSfunctions <- grep("^SS[[:lower:]]+", ls("package:stats"), value = TRUE)
  SSfun <- sapply(SSfunctions, function(ss) 
    grepl(paste0(ss, "("), x = f[2], fixed = TRUE))
  if (any(SSfun)) SSfun <- which(SSfun)
  else stop("'stan_nlmer' requires a named self-starting nonlinear function")
  mc <- match.call(expand.dots = FALSE)
  mc$start <- unlist(getInitial(as.formula(f), data, 
                                control = list(maxiter = 0, warnOnly = TRUE)))
  nlf <- nlformula(mc)
  y <- nlf$respMod$y
  X <- nlf$X
  nlf$reTrms$SSfun <- SSfun
  inputs <- as.character(nlf$respMod$nlmod[2])
  inputs <- sub("(", ",", inputs, fixed = TRUE)
  inputs <- sub(")", "", inputs, fixed = TRUE)
  inputs <- scan(text = inputs, what = character(), sep = ",", 
                 strip.white = TRUE, quiet = TRUE)
  if (SSfun == 5) {
    nlf$reTrms$Dose <- nlf$frame[[inputs[2]]]
    nlf$reTrms$input <- nlf$frame[[inputs[3]]]
  }
  else nlf$reTrms$input <- nlf$frame[[inputs[2]]]
  
  nlf$reTrms$decov <- prior_covariance
  algorithm <- match.arg(algorithm)
  g <- gaussian(link = "identity")
  weights <- nlf$respMod$weights
  offset  <- nlf$respMod$offset
  stanfit <- stan_glm.fit(x = X, y = y, family = g,
                          weights = weights, offset = offset,
                          prior = prior, prior_intercept = NULL,
                          prior_aux = prior_aux, prior_PD = prior_PD, 
                          algorithm = algorithm, adapt_delta = adapt_delta,
                          group = nlf$reTrms, QR = FALSE, sparse = sparse, ...)
  if (SSfun == 6L) {
    stanfit@sim$samples <- lapply(stanfit@sim$samples, FUN = function(x) {
      x[[4L]] <- exp(x[[4L]])
      return(x)
    })
  }
  if (SSfun == 8L) {
    stanfit@sim$samples <- lapply(stanfit@sim$samples, FUN = function(x) {
      x[[3L]] <- exp(x[[3L]])
      return(x)
    })
  }
  Z <- pad_reTrms(Ztlist = nlf$reTrms$Ztlist, cnms = nlf$reTrms$cnms, 
                  flist = nlf$reTrms$flist)$Z
  colnames(Z) <- b_names(names(stanfit), value = TRUE)
  g$link <- paste("inv", SSfunctions[SSfun], sep = "_")
  g$linkinv <- function(eta, arg1, arg2 = NULL, FUN = SSfunctions[SSfun]) {
    if (is.matrix(eta)) {
      len <- length(arg1)
      nargs <- ncol(eta) / len
      SSargs <- lapply(1:nargs, FUN = function(i) {
        start <- 1 + (i - 1) * len
        end <- i * len
        t(eta[,start:end, drop = FALSE])
      })
      if (is.null(arg2)) SSargs <- c(list(arg1), SSargs)
      else SSargs <- c(list(arg1, arg2), SSargs)
    }
    else {
      SSargs <- as.data.frame(matrix(eta, nrow = length(arg1))) 
      if (is.null(arg2)) SSargs <- cbind(arg1, SSargs)
      else SSargs <- cbind(arg1, arg2, SSargs)
    }
    names(SSargs) <- names(formals(FUN))
    if (FUN == "SSbiexp") SSargs$A1 <- SSargs$A1 + exp(SSargs$A2)
    do.call(FUN, args = SSargs)
  }
  if (SSfun == 5) {
    formals(g$linkinv)$arg1 <- nlf$frame[[inputs[2]]]
    formals(g$linkinv)$arg2 <- nlf$frame[[inputs[3]]]
  }
  else formals(g$linkinv)$arg1 <- nlf$frame[[inputs[2]]]
  g$linkfun  <- function(mu) stop("'linkfun' should not have been called")
  g$variance <- function(mu) stop("'variance' should not have been called")
  g$mu.eta   <- function(mu) stop("'mu.eta' should not have been called")
  fit <- nlist(stanfit, family = g, formula, offset, weights, 
               x = if (getRversion() < "3.2.0") cBind(X, Z) else cbind2(X, Z), 
               y = y, data, call = match.call(), terms = NULL, model = NULL, 
               na.action = na.omit, contrasts, algorithm, glmod = nlf, 
               stan_function = "stan_nlmer")
  out <- stanreg(fit)
  class(out) <- c(class(out), "nlmerMod", "lmerMod")
  return(out)
}
