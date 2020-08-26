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

#' Methods for stanreg objects
#' 
#' The methods documented on this page are actually some of the least important 
#' methods defined for \link[=stanreg-objects]{stanreg} objects. The most 
#' important methods are documented separately, each with its own page. Links to
#' those pages are provided in the \strong{See Also} section, below.
#' 
#' @name stanreg-methods
#' @aliases VarCorr fixef ranef ngrps sigma nsamples
#' 
#' @templateVar stanregArg object,x
#' @template args-stanreg-object
#' @param ... Ignored, except by the \code{update} method. See
#'   \code{\link{update}}.
#' 
#' @details The methods documented on this page are similar to the methods 
#'   defined for objects of class 'lm', 'glm', 'glmer', etc. However there are a
#'   few key differences:
#'   
#' \describe{
#' \item{\code{residuals}}{
#' Residuals are \emph{always} of type \code{"response"} (not \code{"deviance"}
#' residuals or any other type). However, in the case of \code{\link{stan_polr}}
#' with more than two response categories, the residuals are the difference 
#' between the latent utility and its linear predictor.
#' }
#' \item{\code{coef}}{
#' Medians are used for point estimates. See the \emph{Point estimates} section
#' in \code{\link{print.stanreg}} for more details.
#' }
#' \item{\code{se}}{
#' The \code{se} function returns standard errors based on 
#' \code{\link{mad}}. See the \emph{Uncertainty estimates} section in
#' \code{\link{print.stanreg}} for more details.
#' }
#' \item{\code{confint}}{
#' For models fit using optimization, confidence intervals are returned via a 
#' call to \code{\link[stats:confint]{confint.default}}. If \code{algorithm} is 
#' \code{"sampling"}, \code{"meanfield"}, or \code{"fullrank"}, the
#' \code{confint} will throw an error because the
#' \code{\link{posterior_interval}} function should be used to compute Bayesian 
#' uncertainty intervals.
#' }
#' \item{\code{nsamples}}{
#' The number of draws from the posterior distribution obtained
#' }
#' }
#' 
#' @seealso 
#' \itemize{
#'  \item The \code{\link[=print.stanreg]{print}},
#'    \code{\link[=summary.stanreg]{summary}}, and \code{\link{prior_summary}} 
#'    methods for stanreg objects for information on the fitted model.
#'  \item \code{\link{launch_shinystan}} to use the ShinyStan GUI to explore a
#'    fitted \pkg{rstanarm} model.
#'  \item The \code{\link[=plot.stanreg]{plot}} method to plot estimates and
#'    diagnostics.
#'  \item The \code{\link{pp_check}} method for graphical posterior predictive
#'    checking.
#'  \item The \code{\link{posterior_predict}} and \code{\link{predictive_error}}
#'    methods for predictions and predictive errors.
#'  \item The \code{\link{posterior_interval}} and \code{\link{predictive_interval}}
#'    methods for uncertainty intervals for model parameters and predictions.
#'  \item The \code{\link[=loo.stanreg]{loo}}, \code{\link{kfold}}, and
#'  \code{\link{log_lik}} methods for leave-one-out or K-fold cross-validation, 
#'    model comparison, and computing the log-likelihood of (possibly new) data.
#'  \item The \code{\link[=as.matrix.stanreg]{as.matrix}}, \code{as.data.frame}, 
#'    and \code{as.array} methods to access posterior draws.
#' }
#' 
NULL

#' @rdname stanreg-methods
#' @export
coef.stanreg <- function(object, ...) {
  if (is.mer(object)) 
    return(coef_mer(object, ...))
  
  object$coefficients
}

#' @rdname stanreg-methods
#' @export
#' @param parm For \code{confint}, an optional character vector of parameter
#'   names.
#' @param level For \code{confint}, a scalar between \eqn{0} and \eqn{1}
#'   indicating the confidence level to use.
#'
confint.stanreg <- function(object, parm, level = 0.95, ...) {
  if (!used.optimizing(object)) {
    stop("For models fit using MCMC or a variational approximation please use ", 
         "posterior_interval() to obtain Bayesian interval estimates.", 
         call. = FALSE)
  }
  confint.default(object, parm, level, ...)
}

#' @rdname stanreg-methods
#' @export
fitted.stanreg <- function(object, ...)  {
  object$fitted.values
}

#' @rdname stanreg-methods
#' @export 
nobs.stanreg <- function(object, ...) {
  nrow(model.frame(object))
}

#' @rdname stanreg-methods
#' @export 
residuals.stanreg <- function(object, ...) {
  object$residuals
}

#' Extract standard errors
#' 
#' Generic function for extracting standard errors from fitted models.
#' 
#' @export
#' @keywords internal
#' @param object A fitted model object.
#' @param ... Arguments to methods.
#' @return Standard errors of model parameters.
#' @seealso \code{\link{se.stanreg}}
#' 
se <- function(object, ...) UseMethod("se")

#' @rdname stanreg-methods
#' @export
se.stanreg <- function(object, ...) {
  object$ses
}

#' @rdname stanreg-methods
#' @export
#' @method update stanreg
#' @param formula.,evaluate See \code{\link[stats]{update}}.
#'
update.stanreg <- function(object, formula., ..., evaluate = TRUE) {
  call <- getCall(object)
  if (is.null(call)) 
    stop("'object' does not contain a 'call' component.", call. = FALSE)
  extras <- match.call(expand.dots = FALSE)$...
  if (!missing(formula.)) 
    call$formula <- update.formula(formula(object), formula.)
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

#' @rdname stanreg-methods
#' @export 
#' @param correlation For \code{vcov}, if \code{FALSE} (the default) the
#'   covariance matrix is returned. If \code{TRUE}, the correlation matrix is
#'   returned instead.
#'
vcov.stanreg <- function(object, correlation = FALSE, ...) {
  out <- object$covmat
  if (!correlation) return(out)
  cov2cor(out)
}


#' @rdname stanreg-methods
#' @export
#' @export fixef
#' @importFrom lme4 fixef
#' 
fixef.stanreg <- function(object, ...) {
  coefs <- object$coefficients
  coefs[b_names(names(coefs), invert = TRUE)]
}

#' @rdname stanreg-methods
#' @export
#' @export ngrps
#' @importFrom lme4 ngrps
#' 
ngrps.stanreg <- function(object, ...) {
  vapply(.flist(object), nlevels, 1)  
}

#' @rdname stanreg-methods
#' @export
#' @export nsamples
#' @importFrom rstantools nsamples
nsamples.stanreg <- function(object, ...) {
  posterior_sample_size(object)
}

#' @rdname stanreg-methods
#' @export
#' @export ranef
#' @importFrom lme4 ranef
#' 
ranef.stanreg <- function(object, ...) {
  .glmer_check(object)
  point_estimates <- object$stan_summary[, select_median(object$algorithm)]
  out <- ranef_template(object)
  group_vars <- names(out)
  
  for (j in seq_along(out)) {
    tmp <- out[[j]]
    pars <- colnames(tmp) 
    levs <- rownames(tmp)
    levs <- gsub(" ", "_", levs) 
    for (p in seq_along(pars)) {
      stan_pars <- paste0("b[", pars[p], " ", group_vars[j],  ":", levs, "]")
      tmp[[pars[p]]] <- unname(point_estimates[stan_pars])
    }
    out[[j]] <- tmp
  }
  out
}

# Call lme4 to get the right structure for ranef objects
#' @importFrom lme4 lmerControl glmerControl nlmerControl lmer glmer nlmer
ranef_template <- function(object) {
  stan_fun <- object$stan_function %ORifNULL% "stan_glmer"
  
  if (stan_fun != "stan_gamm4") {
    new_formula <- formula(object)
  } else {
    # remove the part of the formula with s() terms just so we can call lme4
    # to get the ranef template without error
    new_formula_rhs <- as.character(object$call$random)[2]
    new_formula_lhs <- as.character(formula(object))[2]
    new_formula <- as.formula(paste(new_formula_lhs, "~", new_formula_rhs))
  }
  
  if (stan_fun != "stan_nlmer" && 
      (is.gaussian(object$family$family) || is.beta(object$family$family))) {
    stan_fun <- "stan_lmer"
  }
  lme4_fun <- switch(
    stan_fun,
    "stan_lmer" = "lmer",
    "stan_nlmer" = "nlmer",
    "glmer" # for stan_glmer, stan_glmer.nb, stan_gamm4 (unless gaussian or beta)
  )
  cntrl_args <- list(optimizer = "Nelder_Mead", optCtrl = list(maxfun = 1))
  if (lme4_fun != "nlmer") { # nlmerControl doesn't allow these
    cntrl_args$check.conv.grad <- "ignore"
    cntrl_args$check.conv.singular <- "ignore"
    cntrl_args$check.conv.hess <- "ignore"
    cntrl_args$check.nlev.gtreq.5 <- "ignore"
    cntrl_args$check.nobs.vs.rankZ <- "ignore"
    cntrl_args$check.nobs.vs.nlev <- "ignore"
    cntrl_args$check.nobs.vs.nRE <- "ignore"
    if (lme4_fun == "glmer") {
      cntrl_args$check.response.not.const <- "ignore"
    }
  }
  
  cntrl <- do.call(paste0(lme4_fun, "Control"), cntrl_args)
  
  fit_args <- list(
    formula = new_formula,
    data = object$data,
    control = cntrl
  )
  
  if (lme4_fun == "nlmer") { # create starting values to avoid error
    fit_args$start <- unlist(getInitial(
      object = as.formula(as.character(formula(object))[2]),
      data = object$data,
      control = list(maxiter = 0, warnOnly = TRUE)
    ))
  }
  
  family <- family(object)
  fam <- family$family
  if (!(fam %in% c("gaussian", "beta"))) {
    if (fam == "neg_binomial_2") {
      family <- stats::poisson()
    } else if (fam == "beta_binomial") {
      family <- stats::binomial()
    } else if (fam == "binomial" && family$link == "clogit") {
      family <- stats::binomial()
    }
    fit_args$family <- family
  }
  
  lme4_fit <- suppressWarnings(do.call(lme4_fun, args = fit_args))
  ranef(lme4_fit)
}


#' @rdname stanreg-methods
#' @export
#' @export sigma
#' @rawNamespace if(getRversion()>='3.3.0') importFrom(stats, sigma) else
#'   importFrom(lme4,sigma)
#'
sigma.stanreg <- function(object, ...) {
  if (!("sigma" %in% rownames(object$stan_summary))) 
    return(1)
  
  object$stan_summary["sigma", select_median(object$algorithm)]
}

#' @rdname stanreg-methods
#' @param sigma Ignored (included for compatibility with
#'   \code{\link[nlme]{VarCorr}}).
#' @export
#' @export VarCorr
#' @importFrom nlme VarCorr
#' @importFrom stats cov2cor
VarCorr.stanreg <- function(x, sigma = 1, ...) {
  dots <- list(...) # used to pass stanmat with a single draw for posterior_survfit
  mat <- if ("stanmat" %in% names(dots)) as.matrix(dots$stanmat) else as.matrix(x)
  cnms <- .cnms(x)
  useSc <- "sigma" %in% colnames(mat)
  if (useSc) sc <- mat[,"sigma"] else sc <- 1
  Sigma <- colMeans(mat[,grepl("^Sigma\\[", colnames(mat)), drop = FALSE])
  nc <- vapply(cnms, FUN = length, FUN.VALUE = 1L)
  nms <- names(cnms)
  ncseq <- seq_along(nc)
  if (length(Sigma) == sum(nc * nc)) { # stanfit contains all Sigma entries
    spt <- split(Sigma, rep.int(ncseq, nc * nc))
    ans <- lapply(ncseq, function(i) {
      Sigma <- matrix(0, nc[i], nc[i])
      Sigma[,] <- spt[[i]]
      rownames(Sigma) <- colnames(Sigma) <- cnms[[i]]
      stddev <- sqrt(diag(Sigma))
      corr <- cov2cor(Sigma)
      structure(Sigma, stddev = stddev, correlation = corr)
    })       
  } else { # stanfit contains lower tri Sigma entries
    spt <- split(Sigma, rep.int(ncseq, (nc * (nc + 1)) / 2))
    ans <- lapply(ncseq, function(i) {
      Sigma <- matrix(0, nc[i], nc[i])
      Sigma[lower.tri(Sigma, diag = TRUE)] <- spt[[i]]
      Sigma <- Sigma + t(Sigma)
      diag(Sigma) <- diag(Sigma) / 2
      rownames(Sigma) <- colnames(Sigma) <- cnms[[i]]
      stddev <- sqrt(diag(Sigma))
      corr <- cov2cor(Sigma)
      structure(Sigma, stddev = stddev, correlation = corr)
    })    
  }
  names(ans) <- nms
  structure(ans, sc = mean(sc), useSc = useSc, class = "VarCorr.merMod")
}

# Exported but doc kept internal ----------------------------------------------

#' family method for stanreg objects
#'
#' @keywords internal
#' @export
#' @param object,... See \code{\link[stats]{family}}.
family.stanreg <- function(object, ...) object$family

#' model.frame method for stanreg objects
#' 
#' @keywords internal
#' @export
#' @param formula,... See \code{\link[stats]{model.frame}}.
#' @param fixed.only See \code{\link[lme4:merMod-class]{model.frame.merMod}}.
#' 
model.frame.stanreg <- function(formula, fixed.only = FALSE, ...) {
  if (is.mer(formula)) {
    fr <- formula$glmod$fr
    if (fixed.only) {
      ff <- formula(formula, fixed.only = TRUE)
      vars <- rownames(attr(terms.formula(ff), "factors"))
      fr <- fr[vars]
    }
    return(fr)
  }
  
  NextMethod("model.frame")
}

#' model.matrix method for stanreg objects
#' 
#' @keywords internal
#' @export
#' @param object,... See \code{\link[stats]{model.matrix}}.
#' 
model.matrix.stanreg <- function(object, ...) {
  if (inherits(object, "gamm4")) return(object$jam$X)
  if (is.mer(object)) return(object$glmod$X)
    
  NextMethod("model.matrix")
}

#' formula method for stanreg objects
#' 
#' @keywords internal
#' @export
#' @param x A stanreg object.
#' @param ... Can contain \code{fixed.only} and \code{random.only} arguments 
#'   that both default to \code{FALSE}.
#' 
formula.stanreg <- function(x, ..., m = NULL) {
  if (is.mer(x) && !isTRUE(x$stan_function == "stan_gamm4")) return(formula_mer(x, ...))
  x$formula
}

#' terms method for stanreg objects
#' @export
#' @keywords internal
#' @param x,fixed.only,random.only,... See lme4:::terms.merMod.
#' 
terms.stanreg <- function(x, ..., fixed.only = TRUE, random.only = FALSE) {
  if (!is.mer(x))
    return(NextMethod("terms"))
  
  fr <- x$glmod$fr
  if (missing(fixed.only) && random.only) 
    fixed.only <- FALSE
  if (fixed.only && random.only) 
    stop("'fixed.only' and 'random.only' can't both be TRUE.", call. = FALSE)
  
  Terms <- attr(fr, "terms")
  if (fixed.only) {
    Terms <- terms.formula(formula(x, fixed.only = TRUE))
    attr(Terms, "predvars") <- attr(terms(fr), "predvars.fixed")
  } 
  if (random.only) {
    Terms <- terms.formula(lme4::subbars(formula.stanreg(x, random.only = TRUE)))
    attr(Terms, "predvars") <- attr(terms(fr), "predvars.random")
  }
  
  return(Terms)
}



# internal ----------------------------------------------------------------
.glmer_check <- function(object) {
  if (!is.mer(object))
    stop("This method is for stan_glmer and stan_lmer models only.", 
         call. = FALSE)
}
.cnms <- function(object, ...) UseMethod(".cnms")
.cnms.stanreg <- function(object, ...) {
  .glmer_check(object)
  object$glmod$reTrms$cnms
}
.flist <- function(object, ...) UseMethod(".flist")
.flist.stanreg <- function(object, ...) {
  .glmer_check(object)
  as.list(object$glmod$reTrms$flist)
}

coef_mer <- function(object, ...) {
  if (length(list(...))) 
    warning("Arguments named \"", paste(names(list(...)), collapse = ", "), 
            "\" ignored.", call. = FALSE)
  fef <- data.frame(rbind(fixef(object)), check.names = FALSE)
  ref <- ranef(object)
  refnames <- unlist(lapply(ref, colnames))
  missnames <- setdiff(refnames, names(fef))
  nmiss <- length(missnames)
  if (nmiss > 0) {
    fillvars <- setNames(data.frame(rbind(rep(0, nmiss))), missnames)
    fef <- cbind(fillvars, fef)
  }
  val <- lapply(ref, function(x) fef[rep.int(1L, nrow(x)), , drop = FALSE])
  for (i in seq(a = val)) {
    refi <- ref[[i]]
    row.names(val[[i]]) <- row.names(refi)
    nmsi <- colnames(refi)
    if (!all(nmsi %in% names(fef))) 
      stop("Unable to align random and fixed effects.", call. = FALSE)
    for (nm in nmsi) 
      val[[i]][[nm]] <- val[[i]][[nm]] + refi[, nm]
  }
  structure(val, class = "coef.mer")
}

justRE <- function(f, response = FALSE) {
  response <- if (response && length(f) == 3) f[[2]] else NULL
  reformulate(paste0("(", vapply(lme4::findbars(f), 
                                 function(x) paste(deparse(x, 500L), 
                                                   collapse = " "), 
                                 ""), ")"), 
              response = response)
}
formula_mer <- function (x, fixed.only = FALSE, random.only = FALSE, ...) {
  if (missing(fixed.only) && random.only) 
    fixed.only <- FALSE
  if (fixed.only && random.only) 
    stop("'fixed.only' and 'random.only' can't both be TRUE.", call. = FALSE)
  
  fr <- x$glmod$fr
  if (is.null(form <- attr(fr, "formula"))) {
    if (!grepl("lmer$", deparse(getCall(x)[[1L]]))) 
      stop("Can't find formula stored in model frame or call.", call. = FALSE)
    form <- as.formula(formula(getCall(x), ...))
  }
  if (fixed.only) {
    form <- attr(fr, "formula")
    form[[length(form)]] <- lme4::nobars(form[[length(form)]])
  }
  if (random.only)
    form <- justRE(form, response = TRUE)
  
  return(form)
}
