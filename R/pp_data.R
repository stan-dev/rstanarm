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

.pp_data <- function(object, newdata = NULL) {
  if (is.null(newdata)) {
    x <- model.matrix(object, data = object$data) 
    offset <- if (is.null(object$offset)) rep(0, nrow(x)) else object$offset
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
