# This file is part of rstanarm.
# Copyright 1995-2007 R Core Development Team
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

# .pp_data_mer <- function(object, newdata = NULL, ...) {
#   if (is.null(newdata)) {
#     X <- get_x(object)
#     Z <- get_z(object)
#   } else {
#     if (!is.data.frame(newdata))
#       stop("newdata should be a data.frame")
#     if (any(is.na(newdata))) 
#       stop("NAs not allowed in newdata")
#     
#     # check levels of grouping variables in newdata
#     levs <- lapply(.flist(object), levels)
#     grps <- names(levs)
#     has_new_levels <- sapply(seq_along(levs), function(j) {
#       new_levs <- unique(newdata[, grps[j]])
#       out <- c(!all(new_levs %in% levs[[j]]))
#       names(out) <- grps[j]
#       out
#     })
#     if (any(has_new_levels)) { # FIXME (allow new levels)
#       stop("New levels found in grouping variable(s) ", 
#            paste(grps[has_new_levels], collapse = ", "))
#     }
#     
#     browser()
#     # get X and Z matrices
#     fr <- object$glmod$fr # original model frame
#     m <- model.frame(attr(fr, "terms"), newdata)
#     control <- glmerControl(check.nlev.gtreq.5 = "ignore",
#                             check.nlev.gtr.1 = "stop",
#                             check.nobs.vs.rankZ = "ignore",
#                             check.nobs.vs.nlev = "ignore",
#                             check.nobs.vs.nRE = "ignore")
#     # frX <- fr[, -1, drop = FALSE] # drop response
#     # newdata <- newdata[, colnames(frX)] # select columns in correct order
#     # keep <- 1:nrow(newdata)
#     # fr2 <- rbind(data.frame(newdata, y = 0), data.frame(frX, y = 0))
#     # newf <- update.formula(formula(object), y ~ .)
#     # glF <- glFormula(newf, data = fr2)
#     glF <- glFormula(formula(object), data = m, control = control)
#     X <- glF$X[keep,, drop=FALSE]
#     Z <- t(as.matrix(glF$reTrms$Zt))[keep,, drop=FALSE]
#   }
#   # combine X and Z matrices for posterior_predict
#   x <- cbind(X, Z)
#   nlist(x, offset = object$offset)
# }

.pp_data_mer <- function(object, newdata, ...) {
  X <- .pp_data_mer_x(object, newdata, ...)
  Z <- .pp_data_mer_z(object, newdata, ...)
  if (is.null(Z)) return(list(x = X, offset = object$offset))
  else {
    x <- cbind(X, Z)
    if (!is.null(attr(Z, "NEW_cols"))) {
      attr(x, "NEW_cols") <- attr(Z, "NEW_cols") + ncol(X)
    }
    return(list(x = x, offset = object$offset))
  }
}
.pp_data_mer_x <- function(object, newdata, ...) {
  X <- get_x(object)
  if (is.null(newdata)) return(X)
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
  X <- model.matrix(RHS, data = mfnew, contrasts.arg = attr(X, "contrasts"))
  return(X)
}

.pp_data_mer_z <- function(object, newdata, re.form=NULL, na.action=na.pass,
                           allow.new.levels=FALSE) {
  NAcheck <- !is.null(re.form) && !is(re.form, "formula") && is.na(re.form)
  fmla0check <- is(re.form, "formula") && length(re.form) == 2 && identical(re.form[[2]], 0)
  if (NAcheck || fmla0check) return(NULL)
  
  ## construct (fixed) model frame in order to find out whether there are
  ## missing data/what to do about them
  ## need rfd to inherit appropriate na.action; need grouping
  ## variables as well as any covariates that are included
  ## in RE terms
  ## FIXME: mfnew is new data frame, rfd is processed new data
  ##        why do we need both/what is each doing/how do they differ?
  ##        rfd is *only* used in mkReTrms
  ##        mfnew is *only* used for its na.action attribute (!) [fixed only]
  ##        using model.frame would mess up matrix-valued predictors (GH #201)
  if (is.null(newdata)) {
    rfd <- mfnew <- model.frame(object)
  } else {
    mfnew <- model.frame(delete.response(terms(object,fixed.only=TRUE)),
                         newdata, na.action=na.action)
    ## make sure we pass na.action with new data
    ## it would be nice to do something more principled like
    ## rfd <- model.frame(~.,newdata,na.action=na.action)
    ## but this adds complexities (stored terms, formula, etc.)
    ## that mess things up later on ...
    ## rfd <- na.action(get_all_vars(delete.response(terms(object,fixed.only=FALSE)), newdata))
    old <- FALSE
    if (old) {
      rfd <- na.action(newdata)
      if (is.null(attr(rfd,"na.action")))
        attr(rfd,"na.action") <- na.action
    } else {
      newdata.NA <- newdata
      if (!is.null(fixed.na.action <- attr(mfnew,"na.action"))) {
        newdata.NA <- newdata.NA[-fixed.na.action,]
      }
      tt <- delete.response(terms(object, random.only=TRUE))
      ## need to let NAs in RE components go through -- they're handled downstream
      rfd <- model.frame(tt,newdata.NA,na.action=na.pass)
      if (!is.null(fixed.na.action))
        attr(rfd,"na.action") <- fixed.na.action
    }
    ##
    ## ## need terms to preserve info about spline/orthog polynomial bases
    ## attr(rfd,"terms") <- terms(object)
    ## ## ... but variables list messes things up; can we fix it?
    ## vlist <- lapply(all.vars(terms(object)), as.name)
    ## attr(attr(rfd,"terms"),"variables") <-  as.call(c(quote(list), vlist))
    ##
    ## take out variables that appear *only* in fixed effects
    ## all.v <- all.vars(delete.response(terms(object,fixed.only=FALSE)))
    ## ran.v <- vapply(findbars(formula(object)),all.vars,"")
    ## fix.v <- all.vars(delete.response(terms(object,fixed.only=TRUE)))
    ## rfd <- model.frame(delete.response(terms(object,fixed.only=FALSE)),
    ## newdata,na.action=na.action)
  }
  if (is.null(re.form)) 
    re.form <- justRE(formula(object))
  if (!inherits(re.form, "formula"))
    stop("'re.form' must be NULL, NA, or a formula.")
  
  ## DROP values with NAs in fixed effects
  if (length(fit.na.action <- attr(mfnew,"na.action")) > 0) {
    newdata <- newdata[-fit.na.action,]
  }
  ## note: mkReTrms automatically *drops* unused levels
  ReTrms <- lme4::mkReTrms(lme4::findbars(re.form[[2]]), rfd)
  ## update Lambdat (ugh, better way to do this?)
  ReTrms <- within(ReTrms, Lambdat@x <- unname(object$glmod$reTrms$theta[Lind]))
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
      stop("random effects specified in re.form that were not present in original model")
    re_x[[rname]][,Rcnms[[i]]]
  })
  names(re_new) <- nRnms
  re_new <- unlist(lapply(re_new, t))  ## must TRANSPOSE RE matrices before unlisting
  NEW_cols <- which(is.na(re_new), useNames = TRUE)
  if (!length(NEW_cols)) NEW_cols <- NULL
  else {
    NEW_vars <- lapply(nRnms, function(x) grep(paste0("^", x, "[1-9]"), names(NEW_cols)))
    names(NEW_vars) <- nRnms
    for (j in seq_along(NEW_vars)) {
      if (length(NEW_vars[[j]])) names(NEW_cols)[NEW_vars[[j]]] <- nRnms[j]
    } 
  }
  # nRnms[NEW_]
  ## FIXME? use vapply(re_new, t, FUN_VALUE=????)
  Zt <- ReTrms$Zt
  Z <- t(as.matrix(Zt))
  attr(Z, "na.action") <- attr(mfnew, "na.action")
  if (!is.null(NEW_cols)) attr(Z, "NEW_cols") <- NEW_cols
  return(Z)
}


# (!is.null(re.form) && !is(re.form, "formula") && is.na(re.form)) || 
#   (is(re.form, "formula") && length(re.form) == 2 && identical(re.form[[2]], 
#                                                                0))

#copied from lme4:::levelfun
levelfun <- function(x, nl.n, allow.new.levels = FALSE) {
  if (!all(nl.n %in% rownames(x))) {
    if (!allow.new.levels) 
      stop("new levels detected in newdata")
    newx <- as.data.frame(matrix(NA, nrow = length(nl.n), 
                                 ncol = ncol(x), dimnames = list(nl.n, names(x))))
    newx[rownames(x), ] <- x
    x <- newx
  }
  if (!all(r.inn <- rownames(x) %in% nl.n)) {
    x <- x[r.inn, , drop = FALSE]
  }
  return(x)
}

