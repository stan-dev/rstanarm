#' Methods
#' 
#' Methods for \link[=stanreg-objects]{stanreg} objects.
#' 
#' @name stanreg-methods
#' @aliases VarCorr fixef ranef ngrps
#' 
#' @templateVar stanregArg object,x
#' @template args-stanreg-object
#' @param ... Ignored.
#' @param parm A character vector of parameter names.
#' @param level The confidence level to use.
#' 
#' @details Most of these methods are similar to the methods defined for objects
#'   of class 'lm', 'glm', 'glmer', etc. However there are a few exceptions:
#'   
#' \itemize{
#' \item \code{confint} If \code{algorithm='sampling'}, \code{confint} returns
#' Bayesian \emph{credible} intervals based on posterior quantiles. If 
#' \code{algorithm='optimizing'}, confidence intervals are returned via a
#' call to \code{\link[stats]{confint.default}}.
#' \item \code{log_lik} If \code{algorithm='sampling'}, \code{log_lik} returns 
#' the \eqn{S} by \eqn{N} pointwise log-likelihood matrix, where \eqn{S} is the 
#' size of the posterior sample and \eqn{N} is the number of data points. There 
#' is no pointwise log-likelihood matrix for models fit using
#' \code{algorithm='optimizing'}.
#' \item \code{residuals} Residuals are \emph{always} of type \code{'response'} 
#' and never \code{'deviance'} residuals or any other type.
#' \item \code{coef} If \code{algorithm='sampling'}, posterior medians are used 
#' as point estimates. If \code{algorithm='optimizing'}, the coefficients are
#' also medians, but they are computed from 1000 draws from the asymptotic
#' sampling distribution of the parameters.
#' \item \code{se} The \code{se} function returns standard errors, which are 
#' proportional to the median absolute deviation (\code{\link[stats]{mad}}) from
#' the posterior median (if \code{algorithm='sampling'}) or the median of the
#' asymptotic sampling distribution (if \code{algorithm='optimizing'}).
#' }
#' 
#' @seealso \code{\link{summary.stanreg}}, \code{\link{as.matrix.stanreg}}
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

#' Pointwise log-likelihood matrix
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

#' Standard errors
#' @export
#' @keywords internal
#' @param object object
se <- function(object, ...) UseMethod("se")

#' @rdname stanreg-methods
#' @export
se.stanreg <- function(object, ...) {
  object$ses
}

#' @rdname stanreg-methods
#' @export 
vcov.stanreg <- function(object, ...) {
  if (is(object, "lmerMod")) {
    sel <- seq_along(fixef(object))
    object$covmat[sel, sel, drop=FALSE]
  }
  else object$covmat
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

#' Residual standard deviation
#' @export
#' @keywords internal
#' @param object object
sigma <- function(object, ...) UseMethod("sigma")

#' @rdname stanreg-methods
#' @export
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
  theta <- means[grepl("^theta_L\\[[[:digit:]]+\\]", names(means))]
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
  if (is(formula, "lmerMod")) {
    fit <- formula
    model.frame(lme4::subbars(fit$formula), fit$data)
  }
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
