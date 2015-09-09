#' Methods
#' 
#' Methods for \link[=stanreg-objects]{'stanreg' objects}.
#' 
#' @name stanreg-methods
#' 
#' @param object A fitted model object returned by one of the \pkg{rstanarm} 
#'   modeling functions. This will be a list with class 'stanreg' as well as at 
#'   least one of 'lm', 'glm', 'polr', 'lmerMod', or 'aov'.
#' @param ... Ignored.
#' @param parm A character vector of parameter names.
#' @param level The confidence level to use.
#' @details The \code{se} method returns standard errors and the \code{log_lik} 
#'   method returns the pointwise log-likelihood matrix. Unlike 
#'   \code{\link[stats]{residuals.glm}}, residuals are of type \code{'response'}
#'   not \code{'deviance'}.
#' @note For the \code{sigma}, \code{fixef}, \code{ranef}, and \code{VarCorr}
#'   methods, \code{object} must be a model fit using \code{stan_lmer} or
#'   \code{stan_glmer}.
#'
#' @seealso \code{\link{stanreg-objects}}
#' 
NULL

#' @rdname stanreg-methods
#' @export 
residuals.stanreg <- function(object, ...) {
  object$residuals
}

#' @rdname stanreg-methods
#' @export 
vcov.stanreg <- function(object, ...) {
  object$covmat
}

#' @rdname stanreg-methods
#' @export
confint.stanreg <- function (object, parm = NULL, level = 0.95, ...) {
  mat <- as.matrix(object$stanfit)
  if (!is.null(parm)) mat <- mat[,parm,drop=FALSE]
  alpha <- (1 - level) / 2
  t(apply(mat, 2, FUN = quantile, probs = c(alpha, 1 - alpha)))
}


#' @rdname stanreg-methods
#' @export
fitted.stanreg <- function(object, ...)  {
  object$fitted.values
}

#' @rdname stanreg-methods
#' @export
se <- function(object) UseMethod("se")

#' @rdname stanreg-methods
#' @export
se.stanreg <- function(object) {
  object$ses
}

# Compute pointwise log-likelihood matrix
log_lik <- function(object, ...) UseMethod("log_lik")

#' @rdname stanreg-methods
#' @export
log_lik.stanreg <- function(object, ...) {
  if (object$algorithm != "sampling")
    stop("Only available for MCMC.", call. = FALSE)
  fun <- .llfun(object)
  args <- .llargs(object)
  sapply(seq_len(args$N), function(i) {
    as.vector(fun(i = i, data = args$data, draws = args$draws)) 
  })
}

#' @rdname stanreg-methods
#' @export
#' 
coef.stanreg <- function(object, ...) {
  if (is(object, "lmerMod")) .mermod_coef(object, ...)
  else object$coefficients
}



.glmer_check <- function(object) {
  if (!is(object, "lmerMod")) {
    message("This method is for stan_glmer and stan_lmer models only.")
    invisible(FALSE)
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
  val <- lapply(ref, function(x) fef[rep.int(1L, nrow(x)), 
                                     , drop = FALSE])
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
#' @importFrom lme4 sigma
#' 
sigma.stanreg <- function(object, ...) {
  if (!("sigma" %in% rownames(object$stan_summary))) {
    warning("sigma not found. This method is only for Gaussian models.", 
            call. = FALSE)
    invisible(NULL)
  }
  else object$stan_summary["sigma", "mean"]
}

#' @rdname stanreg-methods
#' @export
#' @importFrom lme4 VarCorr
#' 
VarCorr.stanreg <- function(object, ...) {
  cnms <- .cnms(object)
  nms <- names(cnms)
  stan_nms <- grep("^var\\[", object$stanfit@sim$fnames_oi, value = TRUE)
  out <- lapply(seq_along(nms), function(j) {
    patt <- paste0("\\|", nms[j], "\\]")
    sel <- grep(patt, stan_nms, value = TRUE)
    vc <- cov(as.matrix(object$stanfit, pars = sel))
    colnames(vc) <- rownames(vc) <- cnms[[j]]
    structure(vc, stddev = sqrt(diag(vc)), correlation = cov2cor(vc))
  })
  names(out) <- nms
  # return object printable using lmer's method for VarCorr objects
  gaus <- family(object)$family == "gaussian"
  out <- structure(out, sc = if (gaus) sigma(object) else NULL, 
                   useSc = gaus, class = "VarCorr.merMod")
  out
}


#' @rdname stanreg-methods
#' @export
#' @importFrom lme4 fixef
#' 
fixef.stanreg <- function(object, ...) {
  coefs <- object$coefficients
  coefs[grep("^b\\[", names(coefs), invert = TRUE)]
}

#' @rdname stanreg-methods
#' @export
#' @importFrom lme4 ranef
#' 
ranef.stanreg <- function(object, ...) {
  if (object$algorithm == "optimizing")
    sel <- grep("^b\\[", rownames(object$stan_summary))
  else sel <- grep("^b\\[", object$stanfit@sim$fnames_oi)
  ans <- object$stan_summary[sel, 1]
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
#   stopifnot(is(whichel, "character"))
#   whchL <- names(ans) %in% whichel
#   ans <- ans[whchL]
  class(ans) <- "ranef.mer"
  ans
}

#' @rdname stanreg-methods
#' @export
#' @importFrom lme4 ngrps
#' 
ngrps.stanreg <- function(object, ...) {
  vapply(.flist(object), nlevels, 1)  
}

