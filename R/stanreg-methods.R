#' Methods for stanreg objects
#' 
#' S3 methods for \link[=stanreg-objects]{stanreg} objects.
#' 
#' @name stanreg-methods
#' @aliases VarCorr fixef ranef ngrps
#' 
#' @templateVar stanregArg object,x
#' @template args-stanreg-object
#' @param ... Ignored.
#' @param parm For \code{confint}, a character vector of parameter names.
#' @param level For \code{confint}, a number between 0 and 1 indicating the
#'   confidence level to use.
#' @param correlation For \code{vcov}, if \code{FALSE} (the default) the
#'   covariance matrix is return. If \code{TRUE}, the correlation matrix is
#'   returned instead.
#' 
#' @details Most of these methods are similar to the methods defined for objects
#'   of class 'lm', 'glm', 'glmer', etc. However there are a few exceptions:
#'   
#' \describe{
#' \item{\code{confint}}{
#' If \code{algorithm="sampling"}, \code{confint} returns 
#' Bayesian \emph{credible} intervals based on posterior quantiles. For
#' \code{"meanfield"} and \code{"fullrank"}, approximate credible intervals are
#' computed from the variational approximation to the posterior. For 
#' \code{"optimizing"}, confidence intervals are returned via a call to
#' \code{\link[stats]{confint.default}}.
#' }
#' 
#' \item{\code{log_lik}}{
#' For \code{algorithm="sampling"} only, the \code{log_lik} function returns the
#' \eqn{S} by \eqn{N} pointwise log-likelihood matrix, where \eqn{S} is the size
#' of the posterior sample and \eqn{N} is the number of data points.
#' }
#' 
#' \item{\code{residuals}}{
#' Residuals are \emph{always} of type \code{'response'} 
#' (not \code{'deviance'} residuals or any other type). 
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
#' }
#' 
#' @seealso Some other S3 methods for stanreg objects have separate 
#'   documentation, including \code{\link{as.matrix.stanreg}},
#'   \code{\link{plot.stanreg}}, \code{\link{print.stanreg}},
#'   \code{\link{summary.stanreg}}.
#' 
NULL

#' @rdname stanreg-methods
#' @export
coef.stanreg <- function(object, ...) {
  if (is(object, "lmerMod")) .mermod_coef(object, ...)
  else object$coefficients
}

#' @rdname stanreg-methods
#' @export
confint.stanreg <- function(object, parm, level = 0.95, ...) {
  if (used.optimizing(object))
    return(confint.default(object, parm, level, ...))
  if (missing(parm)) mat <- as.matrix.stanreg(object)
  else {
    parm[parm == "varying"] <- "b"
    mat <- as.matrix.stanreg(object, pars = parm)
  }
  alpha <- (1 - level) / 2
  t(apply(mat, 2, FUN = quantile, probs = c(alpha, 1 - alpha)))
}

#' @rdname stanreg-methods
#' @export
fitted.stanreg <- function(object, ...)  {
  object$fitted.values
}

#' Extract pointwise log-likelihood matrix
#' 
#' @export
#' @param object object
#' @keywords internal
log_lik <- function(object, ...) UseMethod("log_lik")

#' @rdname stanreg-methods
#' @export
log_lik.stanreg <- function(object, ...) {
  if (!used.sampling(object)) 
    STOP_sampling_only("Pointwise log-likelihood matrix")
  fun <- .llfun(object$family)
  args <- .llargs(object)
  sapply(seq_len(args$N), function(i) {
    as.vector(fun(i = i, data = args$data[i,, drop=FALSE], draws = args$draws)) 
  })
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
#' 
se <- function(object, ...) UseMethod("se")

#' @rdname stanreg-methods
#' @export
se.stanreg <- function(object, ...) {
  object$ses
}

#' @rdname stanreg-methods
#' @export 
vcov.stanreg <- function(object, correlation = FALSE, ...) {
  if (!is(object, "lmerMod")) out <- object$covmat
  else {
    sel <- seq_along(fixef(object))
    out <- object$covmat[sel, sel, drop=FALSE]
  }
  if (correlation) out <- cov2cor(out)
  return(out)
}


.glmer_check <- function(object) {
  if (is.null(object$glmod)) {
    stop("This method is for stan_glmer and stan_lmer models only.")
  }
}
.cnms <- function(object) {
  .glmer_check(object)
  object$glmod$reTrms$cnms
}
.flist <- function(object) {
  .glmer_check(object)
  as.list(object$glmod$reTrms$flist)
}

.mermod_coef <- function(object, ...) {
  if (length(list(...))) 
    warning("arguments named \"", paste(names(list(...)), 
                                        collapse = ", "), "\" ignored")
  fef <- data.frame(rbind(fixef(object)), check.names = FALSE)
  ref <- ranef(object)
  refnames <- unlist(lapply(ref, colnames))
  nmiss <- length(missnames <- setdiff(refnames, names(fef)))
  if (nmiss > 0) {
    fillvars <- setNames(data.frame(rbind(rep(0, nmiss))), 
                         missnames)
    fef <- cbind(fillvars, fef)
  }
  val <- lapply(ref, function(x) fef[rep.int(1L, nrow(x)), , drop = FALSE])
  for (i in seq(a = val)) {
    refi <- ref[[i]]
    row.names(val[[i]]) <- row.names(refi)
    nmsi <- colnames(refi)
    if (!all(nmsi %in% names(fef))) 
      stop("unable to align random and fixed effects")
    for (nm in nmsi) val[[i]][[nm]] <- val[[i]][[nm]] + refi[, nm]
  }
  class(val) <- "coef.mer"
  val
}

#' @rdname stanreg-methods
#' @export
#' @export fixef
#' @importFrom lme4 fixef
#' 
fixef.stanreg <- function(object, ...) {
  coefs <- object$coefficients
  coefs[.bnames(names(coefs), invert = TRUE)]
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
#' @export ranef
#' @importFrom lme4 ranef
#' 
ranef.stanreg <- function(object, ...) {
  if (used.optimizing(object)) 
    sel <- .bnames(rownames(object$stan_summary))
  else sel <- .bnames(object$stanfit@sim$fnames_oi)
  ans <- object$stan_summary[sel, .select_median(object$algorithm)]
  # avoid returning the extra levels that were included
  ans <- ans[!grepl("_NEW_", names(ans), fixed = TRUE)]
  fl <- .flist(object)
  levs <- lapply(fl, levels)
  asgn <- attr(fl, "assign")
  cnms <- .cnms(object)
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
  class(ans) <- "ranef.mer"
  ans
}

#' Extract residual standard deviation
#' 
#' @export
#' @keywords internal
#' @param object Fitted model object.
#' @param ... Arguments to methods.
sigma <- function(object, ...) UseMethod("sigma")

#' @rdname stanreg-methods
#' @export
#' @method sigma stanreg
sigma.stanreg <- function(object, ...) {
  if (!("sigma" %in% rownames(object$stan_summary))) return(1)
  else object$stan_summary["sigma", .select_median(object$algorithm)]
}

#' @rdname stanreg-methods
#' @param sigma Ignored scalar.
#' @param rdig Ignored integer.
#' @export
#' @export VarCorr
#' @importFrom lme4 VarCorr mkVarCorr
VarCorr.stanreg <- function(x, sigma = 1, rdig = 3) {
  cnms <- .cnms(x)
  means <- get_posterior_mean(x$stanfit)
  means <- means[,ncol(means)]
  theta <- means[grepl("^theta_L", names(means))]
  sc <- sigma.stanreg(x)
  out <- lme4::mkVarCorr(sc = sc, cnms = cnms, 
                         nc = vapply(cnms, FUN = length, FUN.VALUE = 1L),
                         theta = theta / sc, nms = names(cnms))
  structure(out, useSc = sc != 1, class = "VarCorr.merMod")
}



# Exported but kept internal ----------------------------------------------

#' formula method for stanreg objects
#' 
#' @keywords internal
#' @export
#' @param x A stanreg object.
#' @param ... Ignored. 
#' 
formula.stanreg <- function(x, ...) x$formula

#' model.frame method for stanreg objects
#' 
#' @keywords internal
#' @export
#' @param formula,... See \code{\link[stats]{model.frame}}.
model.frame.stanreg <- function(formula, ...) {
  if (is(formula, "lmerMod")) return(formula$glmod$fr)
  else NextMethod("model.frame")
}

#' model.matrix method for stanreg objects
#' 
#' @keywords internal
#' @export
#' @param object,... See \code{\link[stats]{model.matrix}}.
model.matrix.stanreg <- function(object, ...) {
  if (is(object, "lmerMod")) object$glmod$X
  else NextMethod("model.matrix")
}
