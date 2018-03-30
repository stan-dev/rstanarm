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

#' Create a list of fitted model objects
#' 
#' @export
#' @param ... One or more fitted model objects.
#' @return A list of class \code{"stanreg_list"}, \code{"stanmvreg_list"}, or
#'   \code{"stanjm_list"}, containing the fitted model objects and metadata
#'   stored as attributes.
#'   
#' @seealso \code{\link{loo_model_weights}} for example usage of
#'   \code{stanreg_list}.
#' 
stanreg_list <- function(...) {
  mods <- list(...)
  mod_nms <- match.call(expand.dots = FALSE)$...
  mods <- setNames(mods, sapply(mod_nms, deparse))
  stanreg_list_generator(mods, model_class = "stanreg")
}


#' @rdname stanreg_list
#' @export
stanmvreg_list <- function(...) {
  mods <- list(...)
  mod_nms <- match.call(expand.dots = FALSE)$...
  mods <- setNames(mods, sapply(mod_nms, deparse))
  stanreg_list_generator(mods, model_class = "stanmvreg")
}

#' @rdname stanreg_list
#' @export
stanjm_list <- function(...) {
  mods <- list(...)
  mod_nms <- match.call(expand.dots = FALSE)$...
  mods <- setNames(mods, sapply(mod_nms, deparse))
  stanreg_list_generator(mods, model_class = "stanjm")
}


# internal ----------------------------------------------------------------

stanreg_list_generator <- function(mods, model_class = "stanreg") {
  if (!length(mods)) {
    stop("At least one model must be provided.", call. = FALSE)
  }
  is_model_class <- get(paste0("is.", model_class), mode = "function")
  if (!all(sapply(mods, is_model_class))) {
    stop("All objects in '...' must be ", model_class, " objects.", call. = FALSE)
  }

  structure(
    mods, 
    class = paste0(model_class, "_list"),
    names = names(mods), 
    families = sapply(mods, function(x) {
      fam <- family(x)
      if (!is.character(fam)) fam <- fam$family
      return(fam)
    })
  )
}

#' @export
names.stanreg_list <- function(x) {
  attr(x, "names")
}

#' @export
names.stanmvreg_list <- function(x) {
  attr(x, "names")
}

#' @export
names.stanjm_list <- function(x) {
  attr(x, "names")
}
