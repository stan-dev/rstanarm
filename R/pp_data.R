# This file is part of rstanarm.
# Copyright 1995-2007 R Core Development Team
# Copyright 2015 Douglas Bates, Martin Maechler, Ben Bolker, Steve Walker
# Copyright 2015 Stan Development Team
#
# rstanarm is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# rstanarm is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with rstanarm.  If not, see <http://www.gnu.org/licenses/>.

pp_data <- function(object, newdata = NULL, ...) {
  if (is(object, "lmerMod")) .pp_data_mer(object, newdata, ...)
  else .pp_data(object, newdata, ...)
}

.pp_data <- function(object, newdata = NULL, ...) {
  if (is.null(newdata)) {
    x <- get_x(object)
    offset <- object$offset %ORifNULL% rep(0, nrow(x))
    return(nlist(x, offset))
  }
  tt <- terms(object)
  Terms <- delete.response(tt)
  m <- model.frame(Terms, newdata, xlev = object$xlevels)
  if (!is.null(cl <- attr(Terms, "dataClasses"))) 
    .checkMFClasses(cl, m)
  x <- model.matrix(Terms, m, contrasts.arg = object$contrasts)
  offset <- rep(0, nrow(x))
  if (!is.null(off.num <- attr(tt, "offset"))) 
    for (i in off.num) {
      offset <- offset + eval(attr(tt, "variables")[[i + 1]], newdata)
    }
  if (!is.null(object$call$offset)) 
    offset <- offset + eval(object$call$offset, newdata)
  nlist(x, offset)
}

.pp_data_mer <- function(object, newdata, ...) {
  x <- .pp_data_mer_x(object, newdata, ...)
  z <- .pp_data_mer_z(object, newdata, ...)
  if (!is.null(z)) {
    K <- ncol(x)
    x <- cbind(x, z)
    if (!is.null(attr(z, "NEW_cols")))
      attr(x, "NEW_cols") <- attr(z, "NEW_cols") + K
  }
  return(nlist(x, offset = object$offset))
}

.pp_data_mer_x <- function(object, newdata, ...) {
  x <- get_x(object)
  if (is.null(newdata)) return(x)
  form <- attr(object$glmod$fr, "formula")
  L <- length(form)
  form[[L]] <- lme4::nobars(form[[L]])
  RHS <- formula(substitute(~R, list(R = form[[L]])))
  Terms <- terms(object)
  mf <- model.frame(object)
  ff <- formula(form)
  vars <- rownames(attr(terms.formula(ff), "factors"))
  mf <- mf[vars]
  isFac <- vapply(mf, is.factor, FUN.VALUE = TRUE)
  isFac[attr(Terms, "response")] <- FALSE
  orig_levs <- if (length(isFac) == 0) 
    NULL else lapply(mf[isFac], levels)
  mfnew <- model.frame(delete.response(Terms), newdata, xlev = orig_levs)
  x <- model.matrix(RHS, data = mfnew, contrasts.arg = attr(x, "contrasts"))
  return(x)
}

# based on lme4:::mkNewReTrms
.pp_data_mer_z <- function(object, newdata, re.form = NULL, na.action = na.pass,
                           allow.new.levels = FALSE) {
  NAcheck <- !is.null(re.form) && !is(re.form, "formula") && is.na(re.form)
  fmla0check <- is(re.form, "formula") && length(re.form) == 2 && identical(re.form[[2]], 0)
  if (NAcheck || fmla0check) return(NULL)
  if (is.null(newdata) && is.null(re.form)) return(get_z(object))
  else if (is.null(newdata)) {
    rfd <- mfnew <- model.frame(object)
  } else {
    mfnew <- model.frame(delete.response(terms(object, fixed.only=TRUE)),
                         newdata, na.action=na.action)
    newdata.NA <- newdata
    if (!is.null(fixed.na.action <- attr(mfnew,"na.action"))) {
      newdata.NA <- newdata.NA[-fixed.na.action,]
    }
    tt <- delete.response(terms(object, random.only=TRUE))
    rfd <- model.frame(tt,newdata.NA,na.action=na.pass)
    if (!is.null(fixed.na.action))
      attr(rfd,"na.action") <- fixed.na.action
  }
  if (is.null(re.form)) 
    re.form <- justRE(formula(object))
  if (!inherits(re.form, "formula"))
    stop("'re.form' must be NULL, NA, or a formula.")
  if (length(fit.na.action <- attr(mfnew,"na.action")) > 0) {
    newdata <- newdata[-fit.na.action,]
  }
  ## note: mkReTrms automatically *drops* unused levels
  ReTrms <- lme4::mkReTrms(lme4::findbars(re.form[[2]]), rfd)
  if (!allow.new.levels && any(vapply(ReTrms$flist, anyNA, NA)))
    stop("NAs are not allowed in prediction data",
         " for grouping variables unless 'allow.new.levels' is TRUE.")
  ns.re <- names(re <- ranef(object))
  nRnms <- names(Rcnms <- ReTrms$cnms)
  if (!all(nRnms %in% ns.re))
    stop("Grouping factors specified in re.form that were not present in original model.")
  new_levels <- lapply(ReTrms$flist, function(x) levels(factor(x)))
  ## fill in/delete levels as appropriate
  re_x <- Map(function(r,n) levelfun(r,n,allow.new.levels=allow.new.levels),
              re[names(new_levels)], new_levels)
  ## pick out random effects values that correspond to
  ##  random effects incorporated in re.form ...
  ## NB: Need integer indexing, as nRnms can be duplicated: (age|Subj) + (sex|Subj) :
  re_new <- lapply(seq_along(nRnms), function(i) {
    rname <- nRnms[i]
    if (!all(Rcnms[[i]] %in% names(re[[rname]])))
      stop("Terms specified in re.form that were not present in original model")
    re_x[[rname]][,Rcnms[[i]]]
  })
  names(re_new) <- nRnms
  NEW_cols <- which(is.na(unlist(lapply(re_new, t))))
  if (!length(NEW_cols)) NEW_cols <- NULL
  else {
    NEW_vars <- lapply(nRnms, function(x) {
      grep(paste0("^", x, "[1-9]"), names(NEW_cols))
      })
    names(NEW_vars) <- nRnms
    for (j in seq_along(NEW_vars)) {
      if (length(NEW_vars[[j]])) names(NEW_cols)[NEW_vars[[j]]] <- nRnms[j]
    } 
  }
  z <- structure(t(as.matrix(ReTrms$Zt)), 
                 na.action = attr(mfnew, "na.action"), 
                 NEW_cols = NEW_cols)
  return(z)
}

#copied from lme4:::levelfun except use NAs instead of 0s in matrix
levelfun <- function(x, nl.n, allow.new.levels = FALSE) {
  if (!all(nl.n %in% rownames(x))) {
    if (!allow.new.levels) stop("new levels detected in newdata")
    newx <- as.data.frame(matrix(NA, # use NA instead of zero
                                 nrow = length(nl.n), ncol = ncol(x), 
                                 dimnames = list(nl.n, names(x))))
    newx[rownames(x), ] <- x
    x <- newx
  }
  if (!all(r.inn <- rownames(x) %in% nl.n)) {
    x <- x[r.inn, , drop = FALSE]
  }
  return(x)
}

