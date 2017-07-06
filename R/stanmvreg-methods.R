# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2015, 2016, 2017 Trustees of Columbia University
# Copyright (C) 2016, 2017 Sam Brilleman
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

#' Methods for stanmvreg objects
#' 
#' S3 methods for \link[=stanmvreg-objects]{stanmvreg} objects. There are also 
#' several methods (listed in \strong{See Also}, below) with their own
#' individual help pages. 
#' The main difference between these methods and the 
#' \link[=stanreg-methods]{stanreg} methods is that the methods described here
#' generally include an additional argument \code{m} which allows the user to
#' specify which submodel they wish to return the result for. If the argument
#' \code{m} is set to \code{NULL} then the result will generally be a named list
#' with each element of the list containing the result for one of the submodels.
#' 
#' @name stanmvreg-methods
#' 
#' @templateVar stanmvregArg object,x
#' @templateVar mArg m
#' @template args-stanmvreg-object
#' @template args-m
#' @template args-remove-stub
#' @param ... Ignored, except by the \code{update} method. See
#'   \code{\link{update}}.
#' 
#' @details Most of these methods are similar to the methods defined for objects
#'   of class 'lm', 'glm', 'glmer', etc. However there are a few exceptions:
#'   
#' \describe{
#' \item{\code{coef}}{
#' Medians are used for point estimates. See the \emph{Point estimates} section
#' in \code{\link{print.stanmvreg}} for more details. \code{coef} returns a list 
#' equal to the length of the number of submodels. The first
#' elements of the list are the coefficients from each of the fitted longitudinal
#' submodels and are the same layout as those returned by \code{coef} method of the
#' \pkg{lme4} package, that is, the sum of the random and fixed effects coefficients 
#' for each explanatory variable for each level of each grouping factor. The final
#' element of the returned list is a vector of fixed effect coefficients from the
#' event submodel. 
#' }
#' \item{\code{se}}{
#' The \code{se} function returns standard errors based on 
#' \code{\link{mad}}. See the \emph{Uncertainty estimates} section in
#' \code{\link{print.stanmvreg}} for more details.
#' }
#' \item{\code{confint}}{
#' Not supplied, since the \code{\link{posterior_interval}} function should 
#' be used instead to compute Bayesian uncertainty intervals.
#' }
#' \item{\code{residuals}}{
#' Residuals are \emph{always} of type \code{"response"} (not \code{"deviance"}
#' residuals or any other type).
#' }
#' \item{\code{log_lik}}{
#' Returns the \eqn{S} by \eqn{N} pointwise log-likelihood matrix,
#' where \eqn{S} is the size of the posterior sample and \eqn{N} is the number
#' of individuals in the fitted model. The likelihood for a single individual 
#' is therefore the sum of the likelihood contributions from their observed
#' longitudinal measurements and their event time data.
#' Note: we use \code{log_lik} rather than defining a
#' \code{\link[stats]{logLik}} method because (in addition to the conceptual
#' difference) the documentation for \code{logLik} states that the return value
#' will be a single number, whereas we return a matrix.
#' }
#' }
#' 
#' @seealso
#' Other S3 methods for stanmvreg objects, which have separate documentation, 
#' including \code{\link{print.stanmvreg}}, and \code{\link{summary.stanmvreg}}.
#' 
#' Also \code{\link{posterior_interval}} for an alternative to \code{confint}, 
#' and \code{posterior_predict}, \code{posterior_traj} and 
#' \code{posterior_survfit} for predictions based on the fitted joint model.
#' 
NULL


#' @rdname stanmvreg-methods
#' @export
#'    
coef.stanmvreg <- function(object, m = NULL, ...) {
  M <- get_M(object)
  if (length(list(...))) 
    warning("Arguments named \"", paste(names(list(...)), collapse = ", "), 
            "\" ignored.", call. = FALSE)
  fef <- lapply(fixef(object), function(x) data.frame(rbind(x), check.names = FALSE))
  ref <- ranef(object)
  refnames <- lapply(ref, function(x) unlist(lapply(x, colnames)))
  missnames <- lapply(1:M, function(m) setdiff(refnames[[m]], names(fef[[m]])))
  nmiss <- sapply(missnames, length)
  if (any(nmiss > 0)) for (x in 1:M) {
    if (nmiss[x] > 0) {
      fillvars <- setNames(data.frame(rbind(rep(0, nmiss[x]))), missnames[[x]])
      fef[[x]] <- cbind(fillvars, fef[[x]])
    }
  }
  val <- lapply(1:M, function(m) 
    lapply(ref[[m]], function(x) fef[[m]][rep.int(1L, nrow(x)), , drop = FALSE]))
  for (x in 1:M) {  # loop over number of markers
    for (i in seq(a = val[[x]])) {  # loop over number of grouping factors
      refi <- ref[[x]][[i]]
      row.names(val[[x]][[i]]) <- row.names(refi)
      nmsi <- colnames(refi)
      if (!all(nmsi %in% names(fef[[x]]))) 
        stop("Unable to align random and fixed effects.", call. = FALSE)
      for (nm in nmsi) 
        val[[x]][[i]][[nm]] <- val[[x]][[i]][[nm]] + refi[, nm]
    }
  }
  val <- lapply(val, function(x) structure(x, class = "coef.mer"))
  if (is.jm(object))
    val <- c(val, list(fixef(object)$Event))        
  if (is.null(m)) list_nms(val, M, stub = get_stub(object)) else val[[m]]       
}

#' @rdname stanmvreg-methods
#' @export
#' 
fitted.stanmvreg <- function(object, m = NULL, ...)  {
  M <- get_M(object)
  stub <- get_stub(object)
  if (is.null(m)) 
    list_nms(object$fitted.values, M, stub = stub) else object$fitted.values[[m]]
}

#' @rdname stanmvreg-methods
#' @export 
residuals.stanmvreg <- function(object, m = NULL, ...) {
  M <- get_M(object)
  stub <- get_stub(object)
  if (is.null(m)) 
    list_nms(object$residuals, M, stub = stub) else object$residuals[[m]]
}

#' @rdname stanmvreg-methods
#' @export
se.stanmvreg <- function(object, m = NULL, ...) {
  M <- get_M(object)
  stub <- get_stub(object)
  if (is.null(m)) list_nms(object$ses, M, stub = stub) else object$ses[[m]]
}

#' @rdname stanmvreg-methods
#' @export
#' @param fixed.only A logical specifying whether to only retain the fixed effect
#'   part of the longitudinal submodel formulas
#' @param random.only A logical specifying whether to only retain the random effect
#'   part of the longitudinal submodel formulas  
formula.stanmvreg <- function (x, fixed.only = FALSE, random.only = FALSE, m = NULL, ...) {
  if (missing(fixed.only) && random.only) 
    fixed.only <- FALSE
  if (fixed.only && random.only) 
    stop("'fixed.only' and 'random.only' can't both be TRUE.", call. = FALSE)
  
  M <- get_M(x)
  fr <- lapply(x$glmod, model.frame) 
  form <- lapply(x$formula, as.formula, ...) # includes Event submodel formula
  glmod_form <- lapply(seq(M), function(i) attr(fr[[i]], "formula")) 
  if (!is.null(glmod_form)) form[1:M] <- glmod_form[1:M]
  if (any(fixed.only, random.only)) {
    if (fixed.only) {
      for (i in 1:M)
        form[[i]][[length(form[[i]])]] <- lme4::nobars(form[[i]][[length(form[[i]])]])
    }
    if (random.only)
      for (i in 1:M)
        form[[i]] <- justRE(form[[i]], response = TRUE)
  }
  if (is.null(m)) return(list_nms(form, M, stub = get_stub(x))) else return(form[[m]])
}

#' terms method for stanmvreg objects
#' @export
#' @keywords internal
#' @templateVar mArg m
#' @template args-m
#' @param x,fixed.only,random.only,... See lme4:::terms.merMod.
#' 
terms.stanmvreg <- function(x, fixed.only = TRUE, random.only = FALSE, m = NULL, ...) {
  if (!is.stanmvreg(x))
    return(NextMethod("terms"))
  if (missing(fixed.only) && random.only) 
    fixed.only <- FALSE
  if (fixed.only && random.only) 
    stop("'fixed.only' and 'random.only' can't both be TRUE.", call. = FALSE)
  Terms <- list()
  if (is.mvmer(x)) {
    M <- get_M(x)
    fr <- lapply(x$glmod, model.frame)
    Terms[1:M] <- lapply(fr, function(i) attr(i, "terms"))
    if (fixed.only) {
      Terms <- lapply(seq(M), function(i) {
        Terms <- terms.formula(formula(x, fixed.only = TRUE, m = i))
        attr(Terms, "predvars") <- attr(terms(fr[[i]]), "predvars.fixed")
        Terms
      })     
    } 
    if (random.only) {
      Terms <- lapply(seq(M), function(i) {
        Terms <- terms.formula(lme4::subbars(formula.stanmvreg(x, random.only = TRUE, m = i)))
        attr(Terms, "predvars") <- attr(terms(fr[[i]]), "predvars.random")
        Terms
      })      
    }
    list_nms(Terms, M, stub = get_stub(x))
  }
  if (is.surv(x)) {
    Terms$Event <- terms(x$coxmod)
  }
  if (is.null(m)) Terms else Terms[[m]]
}

#' @rdname stanmvreg-methods
#' @export
#' @method update stanmvreg
#' @param formula. An updated formula for the model. For a multivariate model  
#'   \code{formula.} should be a list of formulas, as described for the 
#'   \code{formula} argument in \code{\link{stan_mvmer}}.
#' @param formulaLong.,formulaEvent. An updated formula for the longitudinal
#'   or event submodel, when \code{object} was estimated using 
#'   \code{\link{stan_jm}}. For a multivariate joint model \code{formulaLong.} 
#'   should be a list of formulas, as described for the \code{formulaLong}
#'   argument in \code{\link{stan_jm}}.
#' @param evaluate See \code{\link[stats]{update}}.
#'
update.stanmvreg <- function(object, formula., formulaLong., 
                             formulaEvent., ..., evaluate = TRUE) {
  call <- getCall(object)
  M <- get_M(object)
  if (is.null(call)) 
    stop("'object' does not contain a 'call' component.", call. = FALSE)
  extras <- match.call(expand.dots = FALSE)$...
  fm <- formula(object)
  if (!missing(formula.)) {
    if (is.jm(object))
      stop("'formula.' should not be specified for joint models estimated ",
           "using stan_jm. Specify 'formulaLong.' and 'formulaEvent' instead.")
    if (M > 1) {
      if (!is.list(formula.))
        stop("To update the formula for a multivariate model ",
             "'formula.' should be a list of formula objects. Use ",
             "'~ .' if you do not wish to alter the formula for one or ",
             "more of the submodels.", call. = FALSE)
      if (length(formula.) != M)
        stop(paste0("The list provided in 'formula.' appears to be the ",
                    "incorrect length; should be length ", M), call. = FALSE)     
    } else {
      if (!is.list(formula.)) formula. <- list(formula.)
    }
    fm_mvmer <- lapply(1:M, function(m) 
      update.formula(fm[[m]], formula.[[m]]))
    names(fm_mvmer) <- NULL
    fm_mvmer <- as.call(c(quote(list), fm_mvmer))
    call$formula <- fm_mvmer
  }  
  if (!missing(formulaLong.)) {
    if (!is.jm(object))
      stop("'formulaLong.' should only be specified for joint models estimated ",
           "using stan_jm. Specify 'formula.' instead.")
    if (M > 1) {
      if (!is.list(formulaLong.))
        stop("To update the formula for a multivariate joint model ",
             "'formulaLong.' should be a list of formula objects. Use ",
             "'~ .' if you do not wish to alter the formula for one or ",
             "more of the longitudinal submodels.", call. = FALSE)
      if (length(formulaLong.) != M)
        stop(paste0("The list provided in 'formulaLong.' appears to be the ",
             "incorrect length; should be length ", M), call. = FALSE)     
    } else {
      if (!is.list(formulaLong.)) formulaLong. <- list(formulaLong.)
    }
    fm_long <- lapply(1:M, function(m) 
      update.formula(fm[[m]], formulaLong.[[m]]))
    names(fm_long) <- NULL
    fm_long <- as.call(c(quote(list), fm_long))
    call$formulaLong <- fm_long
  }
  if (!missing(formulaEvent.)) {
    if (!is.jm(object))
      stop("'formulaEvent.' should only be specified for joint models estimated ",
           "using stan_jm.")
    call$formulaEvent <- update.formula(fm[[length(fm)]], formulaEvent.)  
  }
  if (length(extras)) {
    existing <- !is.na(match(names(extras), names(call)))
    for (a in names(extras)[existing]) 
      call[[a]] <- extras[[a]]
    if (any(!existing)) {
      call <- c(as.list(call), extras[!existing])
      call <- as.call(call)
    }
  }
  
  if (!evaluate) 
    return(call)
  
  # do this like lme4 update.merMod instead of update.default
  ff <- environment(formula(object))
  pf <- parent.frame()
  sf <- sys.frames()[[1L]]
  tryCatch(eval(call, envir = ff),
           error = function(e) {
             tryCatch(eval(call, envir = sf),
                      error = function(e) {
                        eval(call, pf)
                      })
           })
}

#' @rdname stanmvreg-methods
#' @export
#' @export fixef
#' @importFrom lme4 fixef
#' 
fixef.stanmvreg <- function(object, m = NULL, remove_stub = TRUE, ...) {
  M <- get_M(object)
  coefs <- object$coefficients
  coefs <- lapply(coefs, function(x) x[b_names(names(x), invert = TRUE)])
  if (remove_stub) {
    for (i in 1:length(coefs)) names(coefs[[i]]) <- rm_stub(names(coefs[[i]]))
  }
  if (is.null(m)) list_nms(coefs, M, stub = get_stub(object)) else coefs[[m]]
}

#' @rdname stanmvreg-methods
#' @export
#' @export ngrps
#' @importFrom lme4 ngrps
#' 
ngrps.stanmvreg <- function(object, ...) {
  object$n_grps  
}

#' @rdname stanmvreg-methods
#' @export
#' @export ranef
#' @importFrom lme4 ranef
#'
ranef.stanmvreg <- function(object, m = NULL, ...) {
  M <- get_M(object)
  stub <- get_stub(object)
  all_names <- if (used.optimizing(object))
    rownames(object$stan_summary) else object$stanfit@sim$fnames_oi
  ans_list <- lapply(1:M, function(x) { 
    sel <- b_names_M(all_names, x, stub = stub)
    ans <- object$stan_summary[sel, select_median(object$algorithm)]
    # avoid returning the extra levels that were included
    ans <- ans[!grepl("_NEW_", names(ans), fixed = TRUE)]
    fl <- .flist(object, m = x) 
    levs <- lapply(fl, levels)
    asgn <- attr(fl, "assign")
    cnms <- .cnms(object, m = x) 
    nc <- vapply(cnms, length, 1L)
    nb <- nc * vapply(levs, length, 1L)[asgn]
    nbseq <- rep.int(seq_along(nb), nb)
    ml <- split(ans, nbseq)
    for (i in seq_along(ml)) {
      ml[[i]] <- matrix(ml[[i]], ncol = nc[i], byrow = TRUE, 
                        dimnames = list(NULL, cnms[[i]]))
    }
    ans <- lapply(seq_along(fl), function(i) {
      data.frame(do.call(cbind, ml[asgn == i]), row.names = levs[[i]], 
                 check.names = FALSE)
    })
    names(ans) <- names(fl)
    class(ans) <- c("ranef.mer")
    ans
  })
  if (is.null(m)) list_nms(ans_list, M, stub = get_stub(object)) else ans_list[[m]]
}

#' @rdname stanmvreg-methods
#' @export
#' @export sigma
#' @rawNamespace if(getRversion()>='3.3.0') importFrom(stats, sigma) else
#'   importFrom(lme4,sigma)
#'
sigma.stanmvreg <- function(object, m = NULL, ...) {
  stub <- if (is.null(m)) "Long[1-9]" else if 
    (is.numeric(m)) paste0("Long", m) else if (is.character(m)) m else
      stop("Could not reconcile 'm' argument.")
  nms <- grep("^", stub, "\\|sigma", rownames(object$stan_summary), value = TRUE)
  if (!length(nms)) 
    return(1)
  sigma <- object$stan_summary[nms, select_median(object$algorithm)]
  if (length(sigma) > 1L) {
    new_nms <- gsub("\\|sigma", "", nms)
    names(sigma) <- new_nms
  }
  return(sigma)
}


# Exported but doc kept internal ----------------------------------------------

#' family method for stanmvreg objects
#'
#' @keywords internal
#' @export
#' @templateVar mArg m
#' @template args-m
#' @param object,... See \code{\link[stats]{family}}.
family.stanmvreg <- function(object, m = NULL, ...) {
  M <- get_M(object)
  stub <- get_stub(object)
  if (!is.null(m)) object$family[[m]] else 
    list_nms(object$family, M , stub = stub)
}

#' model.frame method for stanmvreg objects
#' 
#' @keywords internal
#' @export
#' @templateVar mArg m
#' @template args-m
#' @param formula,... See \code{\link[stats]{model.frame}}.
#' @param fixed.only See \code{\link[lme4]{model.frame.merMod}}.
#' 
model.frame.stanmvreg <- function(formula, fixed.only = FALSE, m = NULL, ...) {
  if (is.stanmvreg(formula)) {
    M <- get_M(formula)
    fr <- lapply(formula$glmod, model.frame)
    if (fixed.only) {
      fr <- lapply(seq(M), function(i) {
        ff <- formula(formula, fixed.only = TRUE, m = i)
        vars <- rownames(attr(terms.formula(ff), "factors"))
        fr[[i]][vars]
      })
    }
    fr$Event <- formula$coxmod_stuff$model_frame
    if (is.null(m)) 
      return(list_nms(fr, M, stub = get_stub(formula))) else return(fr[[m]])
  } 
  NextMethod("model.frame")
}


# internal ----------------------------------------------------------------

.stanmvreg_check <- function(object) {
  if (!is.stanmvreg(object))
    stop("This method is for stanmvreg objects only.", call. = FALSE)
}
.cnms.stanmvreg <- function(object, m = NULL, remove_stub = FALSE, ...) {
  .stanmvreg_check(object)
  cnms <- if (is.null(m)) object$cnms else object$glmod[[m]]@cnms
  if (remove_stub) lapply(cnms, rm_stub) else cnms
}
.flist.stanmvreg <- function(object, m = NULL, ...) {
  .stanmvreg_check(object)
  if (is.null(m)) {
    stop("'m = NULL' cannot currently be handled by .flist.stanmvreg method.")
  } else as.list(object$glmod[[m]]@flist)
}
.p <- function(object) {
  .stanmvreg_check(object)
  sapply(object$cnms, length)
}

