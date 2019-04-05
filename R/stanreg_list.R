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

#' Create lists of fitted model objects, combine them, or append new models to
#' existing lists of models.
#' 
#' @export
#' @param ... Objects to combine into a \code{"stanreg_list"},
#'   \code{"stanmvreg_list"}, or \code{"stanjm_list"}. Can be fitted model
#'   objects, existing \code{"stan*_list"} objects to combine, or one existing
#'   \code{"stan*_list"} object followed by fitted model objects to append to
#'   the list.
#' @param model_names Optionally, a character vector of model names. If not
#'   specified then the names are inferred from the name of the objects passed
#'   in via \code{...}. These model names are used, for example, when printing
#'   the results of the \code{loo_compare.stanreg_list} and
#'   \code{loo_model_weights.stanreg_list} methods.
#' @return A list of class \code{"stanreg_list"}, \code{"stanmvreg_list"}, or
#'   \code{"stanjm_list"}, containing the fitted model objects and some metadata
#'   stored as attributes.
#'   
#' @seealso \code{\link{loo_model_weights}} for usage of \code{stanreg_list}.
#' 
stanreg_list <- function(..., model_names = NULL) {
  mods <- list(...)
  names(mods) <- stanreg_list_names(
    model_names, 
    n_models = length(mods), 
    call_dots = match.call(expand.dots = FALSE)$...
  )
  .stanreg_list(mods, model_class = "stanreg")
}

#' @rdname stanreg_list
#' @export
stanmvreg_list <- function(..., model_names = NULL) {
  mods <- list(...)
  names(mods) <- stanreg_list_names(
    model_names, 
    n_models = length(mods), 
    call_dots = match.call(expand.dots = FALSE)$...
  )
  .stanreg_list(mods, model_class = "stanmvreg")
}

#' @rdname stanreg_list
#' @export
stanjm_list <- function(..., model_names = NULL) {
  mods <- list(...)
  names(mods) <- stanreg_list_names(
    model_names, 
    n_models = length(mods), 
    call_dots = match.call(expand.dots = FALSE)$...
  )
  .stanreg_list(mods, model_class = "stanjm")
}


#' @export
names.stanreg_list <- function(x) {
  attr(x, "names")
}


#' @rdname stanreg_list
#' @export
#' @param x The object to print.
print.stanreg_list <- function(x, ...) {
  cl <- class(x)
  if (length(cl) > 1) {
    cl <- cl[1]
  }
  cat(cl, " with ", length(x), " models: \n\n")
  df <- data.frame(
    name = attr(x, "names"), 
    family = unname(attr(x, "families")), 
    formula = sapply(x, function(y) formula_string(formula(y))), 
    row.names = seq_along(x)
  )
  print(df, right = FALSE, ...)
  invisible(x)
}

# internal ----------------------------------------------------------------

#' Create, combine, or append new models to a stanreg_list, stanmvreg_list, or
#' stanjm_list object.
#' 
#' @noRd
#' @param mods List of objects to combine. Can be fitted model objects (stanreg,
#'   stanmvreg, stanjm) or stan*_list objects.
#' @param model_class The type of objects to allow. 
#' @return A stanreg_list, stanmvreg_list, or stanjm_list with one component per
#'   model and attributes containing various metadata about the models.
#'   
.stanreg_list <- function(mods, model_class = c("stanreg", "stanmvreg", "stanjm")) {
  stopifnot(length(mods) >= 1, is.list(mods))
  model_class <- match.arg(model_class)
  is_stanreg_list <- sapply(mods, is.stanreg_list)

  if (!any(is_stanreg_list)) {
    .stopifnot_valid_objects(mods, valid_for = "create", model_class = model_class)
    out <- stanreg_list_create(mods, model_class = model_class)
  } else if (all(is_stanreg_list)) {
    .stopifnot_valid_objects(mods, valid_for = "combine", model_class = model_class)
    out <- stanreg_list_combine(mods, model_class = model_class)
  } else {
    .stopifnot_valid_objects(mods, valid_for = "append", model_class = model_class)
    out <- stanreg_list_append(base_list = mods[[1]], mods = mods[-1], 
                               model_class = model_class)
  }
  
  # set model_name attributes of loo/waic/kfold objects to stanreg_list names
  out <- rename_loos.stanreg_list(out)
  
  return(out)
}


#' Create a stanreg_list from list of fitted model objects
#'
#' @noRd 
#' @param mods List of fitted model objects. 
#' @param model_class What type of list is it? ('stanreg', 'stanmvreg', 'stanjm')
#' @return A stanreg_list object
stanreg_list_create <- function(mods, model_class) {
  list_class <- unique(c(paste0(model_class, "_list"), "stanreg_list"))
  structure(mods, 
    class = list_class,
    names = names(mods), 
    families = stanreg_list_families(mods)
  )
}

#' Combine existing stanreg_list objects
#'
#' @noRd 
#' @param lists List of stanreg_list objects.
#' @param model_class What type of list is it? ('stanreg', 'stanmvreg', 'stanjm')
#' @return A stanreg_list object
#' 
stanreg_list_combine <- function(lists, model_class) {
  N_models_per_list <- sapply(lists, length)
  N_models <- sum(N_models_per_list)

  classes <- lapply(lists, class)
  classes <- sapply(classes, function(x) x[1])
  if (!all(classes == classes[1])) {
    stop("Can't combine ", classes[1], " with ", 
         paste(unique(classes[-1]), collapse = ", "))
  }
  
  new_names <- unlist(lapply(lists, attr, "names", exact = TRUE), use.names = FALSE)
  new_families <- unlist(lapply(lists, attr, "families", exact = TRUE), use.names = FALSE)
  
  new_list <- vector(mode = "list", length = N_models)
  pos <- 1
  for (j in seq_along(lists)) {
    for (m in seq_len(N_models_per_list[j])) {
      new_list[[pos]] <- lists[[j]][[m]]
      pos <- pos + 1
    }
  }

  structure(
    new_list, 
    class = unique(c(paste0(model_class, "_list"), "stanreg_list")),
    names = new_names,
    families = new_families
  )
}

#' Append new models to an existing stanreg_list object
#'
#' @noRd 
#' @param base_list The existing stanreg_list to append the new models to.
#' @param mods List of fitted model objects to append to the existing list.
#' @param model_class What type of list is it? ('stanreg', 'stanmvreg', 'stanjm')
#' @return A stanreg_list object
#' 
stanreg_list_append <- function(base_list, mods, model_class) {
  new_list <- stanreg_list_create(mods, model_class = model_class)
  stanreg_list_combine(list(base_list, new_list), model_class = model_class)
}


is.stanreg_list <- function(x) inherits(x, "stanreg_list")
is.stanmvreg_list <- function(x) is.stanreg_list(x) && inherits(x, "stanmvreg_list")
is.stanjm_list <- function(x) is.stanreg_list(x) && inherits(x, "stanjm_list")

.stopifnot_valid_objects <-
  function(mods,
           valid_for = c("create", "combine", "append"),
           model_class = c("stanreg", "stanmvreg", "stanjm")) {
    valid_for <- match.arg(valid_for)
    model_class <- match.arg(model_class)
    list_class <- paste0(model_class, "_list")
    error_msg <- paste0(
      "For ", list_class,"() objects in '...' must: ", 
      "\n(1) all be ", model_class, " objects, or",
      "\n(2) all be ", list_class, " objects, or", 
      "\n(3) be one ", list_class, " object followed by all ", 
      model_class, " objects"
    )
    
    is_model_class <- sapply(mods, FUN = match.fun(paste0("is.", model_class)))
    is_list_class <- sapply(mods, FUN = match.fun(paste0("is.", list_class)))
    
    throw_error <-
      (valid_for == "create" &&
         !all(is_model_class)) ||
      (valid_for == "combine" &&
         !all(is_list_class)) ||
      (valid_for == "append" &&
         !(is_list_class[1] && all(is_model_class[-1])))
    
    if (throw_error) {
      stop(error_msg, call.  = FALSE)
    }
  }



#' Determine names of the models in a stanreg_list 
#' @noRd 
#' @param user_model_names Either NULL or user-specified model_names argument
#' @param n_models The number of models in the stanreg_list
#' @param call_dots The result of match.call(expand.dots = FALSE)$...
#' @return Either the user-specified model names or names inferred from the
#'   names of the fitted model objects passed to '...'.
#'   
stanreg_list_names <- function(user_model_names, n_models, call_dots) {
  if (!is.null(user_model_names)) {
    stopifnot(is.character(user_model_names))
    if (length(user_model_names) != n_models) {
      stop("Length of 'model_names' must be the same as the number of models.")
    }
    nms <- user_model_names
  } else {
    nms <- sapply(call_dots, FUN = deparse) 
  }
  return(nms)
}

#' Determine the families of the models in a stanreg_list 
#' @noRd 
#' @param mods List of fitted model objects
#' @return Character vector of family names
#' 
stanreg_list_families <- function(mods) {
  fams <- sapply(mods, FUN = function(x) {
    fam <- family(x)
    if (!is.character(fam)) fam <- fam$family
    return(fam)
  })
  unname(fams)
}



# loo/waic/kfold objects created by rstanarm have a model_name attribute. 
# when a stanreg_list is created those attributes should be changed to match 
# the names of the models used for the stanreg_list in case user has specified 
# the model_names argument
rename_loos <- function(x,...) UseMethod("rename_loos")

# Change model_name attributes of a loo/waic/kfold object stored in a stanreg object,
rename_loos.stanreg <- function(x, new_model_name,...) {
  for (criterion in c("loo", "waic", "kfold")) {
     if (!is.null(x[[criterion]])) {
       attr(x[[criterion]], "model_name") <- new_model_name
     }
  }
  return(x)
}

# Change model_name attributes of loo/waic/kfold objects to correspond to 
# model names used for stanreg_list
rename_loos.stanreg_list <- function(x, ...) {
  for (j in seq_along(x)) {
    x[[j]] <- rename_loos.stanreg(x[[j]], new_model_name = names(x)[j])
  }
  return(x)
}

