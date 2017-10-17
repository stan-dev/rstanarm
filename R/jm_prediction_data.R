# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2015, 2016, 2017 Trustees of Columbia University
# Copyright (C) 2016, 2017 Sam Brilleman
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

# Return the design matrices required for evaluating the linear predictor or
# log-likelihood in post-estimation functions for a \code{stan_jm} model
#
# @param object A stanmvreg object
# @param newdataLong A data frame or list of data frames with the new 
#   covariate data for the longitudinal submodel
# @param newdataEvent A data frame with the new covariate data for the
#   event submodel
# @param ids An optional vector of subject IDs specifying which individuals
#   should be included in the returned design matrices.
# @param etimes An optional vector of times at which the event submodel
#   design matrices should be evaluated (also used to determine the 
#   quadrature times). If NULL then times are taken to be the eventimes in
#   the fitted object (if newdataEvent is NULL) or in newdataEvent.
# @param long_parts,event_parts A logical specifying whether to return the
#   design matrices for the longitudinal and/or event submodels.
# @return A named list (with components M, Npat, ndL, ndE, yX, tZt, 
#   yZnames, eXq, assoc_parts) 
jm_data <- function(object, newdataLong = NULL, newdataEvent = NULL, 
                    ids = NULL, etimes = NULL, long_parts = TRUE, 
                    event_parts = TRUE) {
  M <- get_M(object)
  id_var   <- object$id_var
  time_var <- object$time_var
  
  if (!is.null(newdataLong) || !is.null(newdataEvent))
    newdatas <- validate_newdatas(object, newdataLong, newdataEvent)
  
  # prediction data for longitudinal submodels
  ndL <- if (is.null(newdataLong)) 
    get_model_data(object)[1:M] else newdatas[1:M]
  
  # prediction data for event submodel
  ndE <- if (is.null(newdataEvent)) 
    get_model_data(object)[["Event"]] else newdatas[["Event"]]   
  
  # possibly subset
  if (!is.null(ids)) {
    ndL <- subset_ids(object, ndL, ids)
    ndE <- subset_ids(object, ndE, ids)
  }
  id_list <- unique(ndE[[id_var]]) # unique subject id list
  
  # evaluate the last known survival time and status
  if (!is.null(newdataEvent) && is.null(etimes)) {
    # prediction data for the event submodel was provided but  
    # no event times were explicitly specified by the user, so
    # they must be evaluated using the data frame
    surv <- eval(formula(object, m = "Event")[[2L]], ndE)
    etimes  <- unclass(surv)[,"time"]
    estatus <- unclass(surv)[,"status"]    
  } else if (is.null(etimes)) {
    # if no prediction data was provided then event times are 
    # taken from the fitted model
    etimes  <- object$eventtime[as.character(id_list)]
    estatus <- object$status[as.character(id_list)]
  } else { 
    # otherwise, event times ('etimes') are only directly specified for dynamic   
    # predictions via posterior_survfit in which case the 'etimes' correspond 
    # to the last known survival time and therefore we assume everyone has survived
    # up to that point (ie, set estatus = 0 for all individuals), this is true 
    # even if there is an event indicated in the data supplied by the user.
    estatus <- rep(0, length(etimes))
  }
  res <- nlist(M, Npat = length(id_list), ndL, ndE)
  
  if (long_parts && event_parts) 
    lapply(ndL, function(x) {
      if (!time_var %in% colnames(x)) STOP_no_var(time_var)
      mt <- tapply(x[[time_var]], factor(x[[id_var]]), max)
      if (any(mt > etimes))
        stop("There appears to be observation times in the longitudinal data that ",
             "are later than the event time specified in the 'etimes' argument.")      
    }) 
  
  # response and design matrices for longitudinal submodels
  if (long_parts) {
    y <- lapply(1:M, function(m) eval(formula(object, m = m)[[2L]], ndL[[m]]))
    ydat <- lapply(1:M, function(m) pp_data(object, ndL[[m]], m = m))
    yX <- fetch(ydat, "x")
    yZt <- fetch(ydat, "Zt")
    yZ_names <- fetch(ydat, "Z_names")
    flist <- lapply(ndL, function(x) factor(x[[id_var]]))
    res <- c(res, nlist(y, yX, yZt, yZ_names, flist))
  }
  
  # design matrices for event submodel and association structure
  if (event_parts) {
    qnodes <- object$qnodes
    qq <- get_quadpoints(qnodes)
    qtimes <- uapply(qq$points,  unstandardise_qpts, 0, etimes)
    qwts   <- uapply(qq$weights, unstandardise_qwts, 0, etimes)
    starttime <- deparse(formula(object, m = "Event")[[2L]][[2L]])
    edat <- prepare_data_table(ndE, id_var, time_var = starttime)
    times <- c(etimes, qtimes) # times used to design event submodel matrices
    edat <- rolling_merge(edat, ids = rep(id_list, qnodes + 1), times = times)
    eXq  <- .pp_data_mer_x(object, newdata = edat, m = "Event")
    assoc_parts <- lapply(1:M, function(m) {
      ymf <- ndL[[m]]
      grp_stuff <- object$grp_stuff[[m]]
      if (grp_stuff$has_grp) {
        grp_stuff <- get_extra_grp_info( # update grp_info with new data
          grp_stuff, flist = ymf, id_var = id_var, 
          qnodes = qnodes, grp_assoc = object$grp_assoc)
        grp_var <- grp_stuff$grp_var
        ymf[[grp_var]] <- factor(ymf[[grp_var]])
        ymf <- data.table::data.table(ymf, key = c(id_var, grp_var, time_var))
      } else {
        ymf <- data.table::data.table(ymf, key = c(id_var, time_var))
      }
      make_assoc_parts(
        ymf, assoc = object$assoc[,m], id_var = id_var, time_var = time_var, 
        ids = id_list, times = times, grp_stuff = grp_stuff,
        use_function = pp_data, object = object, m = m)
    })
    assoc_attr <- nlist(.Data = assoc_parts, qnodes, qtimes, qwts, etimes, estatus)
    assoc_parts <- do.call("structure", assoc_attr)
    res <- c(res, nlist(eXq, assoc_parts))
  }
  
  return(res)
}


# Return a data frame for each submodel that:
# (1) only includes variables used in the model formula
# (2) only includes rows contained in the glmod/coxmod model frames
# (3) ensures that additional variables that are required
#     such as the ID variable or variables used in the 
#     interaction-type association structures, are included.
#
# It is necessary to drop unneeded variables though so that 
# errors are not encountered if the original data contained 
# NA values for variables unrelated to the model formula.
# We generate a data frame here for in-sample predictions 
# rather than using a model frame, since some quantities will
# need to be recalculated at quadrature points etc, for example
# in posterior_survfit.
#
# @param object A stanmvreg object.
# @param m Integer specifying which submodel to get the
#   prediction data frame for.
# @return A data frame or list of data frames with all the
#   (unevaluated) variables required for predictions.
get_model_data <- function(object, m = NULL) {
  validate_stanmvreg_object(object)
  M <- get_M(object)
  terms <- terms(object, fixed.only = FALSE)
  
  # identify variables to add to the terms objects
  if (is.jm(object)) {
    extra_vars <- lapply(1:M, function(m) {
      # for each submodel loop over the four possible assoc  
      # interaction formulas and collect any variables used
      forms_m <- object$assoc["which_formulas",][[m]]
      uapply(forms_m, function(x) {
        if (length(x)) {
          rownames(attr(terms.formula(x), "factors")) 
        } else NULL
      })
    })
    # also ensure that id_var is in the event data
    extra_vars$Event <- object$id_var
    
    if (!identical(length(terms), length(extra_vars)))
      stop2("Bug found: terms and extra_vars should be same length.")
    
    # add the extra variables to the terms formula for each submodel
    terms <- xapply(terms, extra_vars, FUN = function(x, y) {
      lhs <- x[[2L]]
      rhs <- deparse(x[[3L]], 500L)
      if (!is.null(y))
        rhs <- c(rhs, y)
      reformulate(rhs, response = lhs)
    })
    
    datas <- c(object$dataLong, list(object$dataEvent))
  } else {
    datas <- object$data
  }
  
  # identify rows that were in the model frame
  row_nms <- lapply(model.frame(object), rownames)
  
  # drop rows and variables not required for predictions
  mfs <- xapply(w = terms, x = datas, y = row_nms,
                FUN = function(w, x, y) 
                  get_all_vars(w, x)[y, , drop = FALSE])
  
  mfs <- list_nms(mfs, M, stub = get_stub(object))
  if (is.null(m)) mfs else mfs[[m]]
}

