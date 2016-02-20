# Part of the rstanarm package for estimating model parameters
# Copyright 1995-2007 R Core Development Team
# Copyright 2015 Douglas Bates, Martin Maechler, Ben Bolker, Steve Walker
# Copyright (C) 2015, 2016 Trustees of Columbia University
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

pp_data <- function(object, newdata = NULL, re.form = NULL, ...) {
  if (is.mer(object)) .pp_data_mer(object, newdata, re.form, ...)
  else .pp_data(object, newdata, ...)
}

# for models not fit using stan_(g)lmer or stan_gamm4
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

# for models fit using stan_(g)lmer or stan_gamm4
.pp_data_mer <- function(object, newdata, re.form, ...) {
  x <- .pp_data_mer_x(object, newdata, ...)
  z <- .pp_data_mer_z(object, newdata, re.form, ...)
  return(nlist(x, offset = object$offset, Zt = z$Zt, Z_names = z$Z_names))
}

# the functions below are heavily based on a combination of 
# lme4:::predict.merMod and lme4:::mkNewReTrms, although they do also have 
# substantial modifications
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

.pp_data_mer_z <- function(object, newdata, re.form = NULL,
                           allow.new.levels = TRUE, na.action = na.pass) {
  NAcheck <- !is.null(re.form) && !is(re.form, "formula") && is.na(re.form)
  fmla0check <- (is(re.form, "formula") && 
                   length(re.form) == 2 && 
                   identical(re.form[[2]], 0))
  if (NAcheck || fmla0check) return(list())
  if (is.null(newdata) && is.null(re.form)) {
    Z <- get_z(object)
    return(list(Zt = t(Z)))
  }
  else if (is.null(newdata)) {
    rfd <- mfnew <- model.frame(object)
  } else {
    if ("gam" %in% names(object))
      stop("'posterior_predict' with non-NULL 're.form' not yet supported ", 
           "for models estimated via 'stan_gamm4'")
    mfnew <- model.frame(delete.response(terms(object, fixed.only = TRUE)),
                         newdata, na.action = na.action)
    newdata.NA <- newdata
    if (!is.null(fixed.na.action <- attr(mfnew,"na.action"))) {
      newdata.NA <- newdata.NA[-fixed.na.action,]
    }
    tt <- delete.response(terms(object, random.only = TRUE))
    rfd <- model.frame(tt, newdata.NA, na.action = na.pass)
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
  ReTrms <- lme4::mkReTrms(lme4::findbars(re.form[[2]]), rfd)
  if (!allow.new.levels && any(vapply(ReTrms$flist, anyNA, NA)))
    stop("NAs are not allowed in prediction data",
         " for grouping variables unless 'allow.new.levels' is TRUE.")
  ns.re <- names(re <- ranef(object))
  nRnms <- names(Rcnms <- ReTrms$cnms)
  if (!all(nRnms %in% ns.re))
    stop("Grouping factors specified in re.form that were not present in original model.")
  new_levels <- lapply(ReTrms$flist, function(x) levels(factor(x)))
  Zt <- ReTrms$Zt
  p <- sapply(ReTrms$cnms, FUN = length)
  l <- sapply(attr(ReTrms$flist, "assign"), function(i) 
    nlevels(ReTrms$flist[[i]]))
  t <- length(p)
  group_nms <- names(ReTrms$cnms)
  Z_names <- character()
  for (i in seq_along(ReTrms$cnms)) {
    # if you change this, change it in stan_glm.fit() as well
    nm <- group_nms[i]
    nms_i <- paste(ReTrms$cnms[[i]], group_nms[i])
    if (length(nms_i) == 1) {
      Z_names <- c(Z_names, paste0(nms_i, ":", levels(ReTrms$flist[[nm]])))
    } else {
      Z_names <- c(Z_names, c(t(sapply(nms_i, paste0, ":", new_levels[[nm]]))))
    }
  }
  z <- nlist(Zt = ReTrms$Zt, Z_names)
  return(z)
}
