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

.pp_data_mer <- function(object, newdata = NULL, ...) {
  if (is.null(newdata)) {
    x <- get_x(object)
    z <- get_z(object)
  } else {
    if (!is.data.frame(newdata))
      stop("newdata should be a data.frame")
    if (any(is.na(newdata))) 
      stop("NAs not allowed in newdata")
    fr <- object$glmod$fr # original model frame
#     notfound <- setdiff(colnames(newdata), colnames(fr))
#     if (length(notfound)) {
#       notfound <- paste(notfound, collapse = ", ")
#       stop("Variable(s) ", notfound, " in newdata but not original formula.")
#     }
    # check levels of grouping variables in newdata
    levs <- lapply(.flist(object), levels)
    grps <- names(levs)
    has_new_levels <- sapply(seq_along(levs), function(j) {
      new_levs <- unique(newdata[, grps[j]])
      !all(new_levs %in% levs[[j]])
    })
    if (any(has_new_levels)) {
      stop("New levels found in grouping variable(s) ", 
           paste(grps[has_new_levels], collapse = ", "))
    }
    
    # get X and Z matrices
    frX <- fr[, -1, drop = FALSE]
    newdata <- newdata[, colnames(frX)] # make sure columns in same order
    keep <- 1:nrow(newdata)
    fr2 <- rbind(data.frame(newdata, y = 0), data.frame(frX, y = 0))
    newf <- update.formula(formula(object), y ~ .)
    glF <- glFormula(newf, data = fr2)
    x <- glF$X[keep,, drop=FALSE]
    z <- t(as.matrix(glF$reTrms$Zt))[keep,, drop=FALSE]
  }
  nlist(x, z, offset = object$offset)
}

