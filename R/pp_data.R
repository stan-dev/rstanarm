# Part of the rstanarm package for estimating model parameters
# Copyright 2015 Douglas Bates, Martin Maechler, Ben Bolker, Steve Walker
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

pp_data <-
  function(object,
           newdata = NULL,
           re.form = NULL,
           offset = NULL,
           m = NULL,
           ...) {
    validate_stanreg_object(object)
    if (is.mer(object)) {
      out <- .pp_data_mer(object, newdata = newdata,
                          re.form = re.form, m = m, ...)
      if (!is.null(offset)) out$offset <- offset
      return(out)
    }
    else
      .pp_data(object, newdata = newdata, offset = offset, ...)
  }

# for models without lme4 structure
.pp_data <- function(object, newdata = NULL, offset = NULL, ...) {
  if (is(object, "gamm4")) {
    if (is.null(newdata))   x <- predict(object$jam, type = "lpmatrix")
    else x <- predict(object$jam, newdata = newdata, type = "lpmatrix")
    if (is.null(offset)) 
      offset <- object$offset %ORifNULL% rep(0, nrow(x))
    return(nlist(x, offset))
  }
  if (is.null(newdata)) {
    x <- get_x(object)
    if (is.null(offset)) 
      offset <- object$offset %ORifNULL% rep(0, nrow(x))
    return(nlist(x, offset))
  }

  offset <- .pp_data_offset(object, newdata, offset)
  Terms <- delete.response(terms(object))
  m <- model.frame(Terms, newdata, xlev = object$xlevels)
  if (!is.null(cl <- attr(Terms, "dataClasses"))) 
    .checkMFClasses(cl, m)
  x <- model.matrix(Terms, m, contrasts.arg = object$contrasts)
  if (is(object, "polr") && !is_scobit(object)) 
    x <- x[,colnames(x) != "(Intercept)", drop = FALSE]
  
  nlist(x, offset)
}


# for models fit using stan_(g)lmer or stan_gamm4
.pp_data_mer <- function(object, newdata, re.form, m = NULL, ...) {
  if (is(object, "gamm4")) {
    if (is.null(newdata))   x <- predict(object$jam, type = "lpmatrix")
    else x <- predict(object$jam, newdata = newdata, type = "lpmatrix")
    if (is.null(re.form)) {
      re.form <- as.formula(object$call$random)
      if (length(re.form) == 0) re.form <- NA
      z <- .pp_data_mer_z(object, newdata, re.form, ...)
    }
    else z <- .pp_data_mer_z(object, newdata, re.form, ...)
  } else {
    x <- .pp_data_mer_x(object, newdata, m = m, ...)
    z <- .pp_data_mer_z(object, newdata, re.form, m = m, ...)
  }
  offset <- model.offset(model.frame(object, m = m))
  if (!missing(newdata) && (!is.null(offset) || !is.null(object$call$offset))) {
    offset <- try(eval(object$call$offset, newdata), silent = TRUE)
    if (!is.numeric(offset)) offset <- NULL
  }
  return(nlist(x, offset = offset, Zt = z$Zt, Z_names = z$Z_names))
}


# the functions below are heavily based on a combination of 
# lme4:::predict.merMod and lme4:::mkNewReTrms, although they do also have 
# substantial modifications
.pp_data_mer_x <- function(object, newdata, m = NULL, ...) {
  x <- get_x(object, m = m)
  if (is.null(newdata)) return(x)
  form <- if (is.null(m)) attr(object$glmod$fr, "formula") else 
    formula(object, m = m)
  L <- length(form)
  form[[L]] <- lme4::nobars(form[[L]])
  RHS <- formula(substitute(~R, list(R = form[[L]])))
  Terms <- terms(object, m = m)
  mf <- model.frame(object, m = m)
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
                           allow.new.levels = TRUE, na.action = na.pass, 
                           m = NULL, ...) {
  NAcheck <- !is.null(re.form) && !is(re.form, "formula") && is.na(re.form)
  fmla0check <- (is(re.form, "formula") && 
                   length(re.form) == 2 && 
                   identical(re.form[[2]], 0))
  if (NAcheck || fmla0check) return(list())
  if (is.null(newdata) && is.null(re.form)) {
    Z <- get_z(object, m = m)
    return(list(Zt = t(Z)))
  }
  else if (is.null(newdata)) {
    rfd <- mfnew <- model.frame(object, m = m)
  } 
  else if (inherits(object, "gamm4")) {
    if (is.null(newdata))   x <- predict(object$jam, type = "lpmatrix")
    else x <- predict(object$jam, newdata = newdata, type = "lpmatrix")
    NAs <- apply(is.na(x), 1, any)
    rfd <- mfnew <- newdata[!NAs,]
    attr(rfd,"na.action") <- "na.omit"
  } else {
    terms_fixed <- delete.response(terms(object, fixed.only = TRUE, m = m))
    mfnew <- model.frame(terms_fixed, newdata, na.action = na.action)
    newdata.NA <- newdata
    if (!is.null(fixed.na.action <- attr(mfnew,"na.action"))) {
      newdata.NA <- newdata.NA[-fixed.na.action,]
    }
    tt <- delete.response(terms(object, random.only = TRUE, m = m))
    rfd <- model.frame(tt, newdata.NA, na.action = na.pass)
    if (!is.null(fixed.na.action))
      attr(rfd,"na.action") <- fixed.na.action
  }
  if (is.null(re.form)) 
    re.form <- justRE(formula(object, m = m))
  if (!inherits(re.form, "formula"))
    stop("'re.form' must be NULL, NA, or a formula.")
  if (length(fit.na.action <- attr(mfnew,"na.action")) > 0) {
    newdata <- newdata[-fit.na.action,]
  }
  ReTrms <- lme4::mkReTrms(lme4::findbars(re.form[[2]]), rfd)
  if (!allow.new.levels && any(vapply(ReTrms$flist, anyNA, NA)))
    stop("NAs are not allowed in prediction data",
         " for grouping variables unless 'allow.new.levels' is TRUE.")
  ns.re <- names(re <- ranef(object, m = m))
  nRnms <- names(Rcnms <- ReTrms$cnms)
  if (!all(nRnms %in% ns.re))
    stop("Grouping factors specified in re.form that were not present in original model.")
  new_levels <- lapply(ReTrms$flist, function(x) levels(factor(x)))
  Zt <- ReTrms$Zt
  Z_names <- make_b_nms(ReTrms, m = m)
  z <- nlist(Zt = ReTrms$Zt, Z_names)
  return(z)
}



# handle offsets ----------------------------------------------------------
null_or_zero <- function(x) {
  isTRUE(is.null(x) || all(x == 0))
}

.pp_data_offset <- function(object, newdata = NULL, offset = NULL) {
  if (is.null(newdata)) {
    # get offset from model object (should be null if no offset)
    if (is.null(offset)) 
      offset <- object$offset %ORifNULL% model.offset(model.frame(object))
  } else {
    if (!is.null(offset))
      stopifnot(length(offset) == nrow(newdata))
    else {
      # if newdata specified but not offset then confirm that model wasn't fit
      # with an offset (warning, not error)
      if (!is.null(object$call$offset) || 
          !null_or_zero(object$offset) || 
          !null_or_zero(model.offset(model.frame(object)))) {
        warning(
          "'offset' argument is NULL but it looks like you estimated ", 
          "the model using an offset term.", 
          call. = FALSE
        )
      }
      offset <- rep(0, nrow(newdata))
    }
  }
  return(offset)
}
