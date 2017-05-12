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

# Return design matrices required for evaluating the linear predictor or
# log-likelihood of the longitudinal/event submodels in a stan_jm model
#
# @param object A stanjm object
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
  ndL <- if (is.null(newdataLong)) 
    get_model_data(object)[1:M] else newdatas[1:M]
  ndE <- if (is.null(newdataEvent)) 
    get_model_data(object)[["Event"]] else newdatas[["Event"]]   
  if (!is.null(ids)) {
    ndL <- subset_ids(object, ndL, ids)
    ndE <- subset_ids(object, ndE, ids)
  }
  id_list <- unique(ndE[[id_var]])
  if (!is.null(newdataEvent) && is.null(etimes)) {
    surv <- eval(formula(object, m = "Event")[[2L]], ndE)
    etimes  <- unclass(surv)[,"time"]
    estatus <- unclass(surv)[,"status"]    
  } else if (is.null(etimes)) {
    etimes  <- object$eventtime[as.character(id_list)]
    estatus <- object$status[as.character(id_list)]
  } else { 
    # 'etimes' are only directly specified for dynamic predictions via 
    # posterior_survfit in which case the 'etimes' correspond to the last known 
    # survival time and therefore we assume everyone has survived up to that 
    # point (ie, set estatus = 0 for all individuals), this is true even if 
    # there is an event indicated in the data supplied by the user.
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
  if (long_parts) {
    y <- lapply(1:M, function(m) eval(formula(object, m = m)[[2L]], ndL[[m]]))
    ydat <- lapply(1:M, function(m) pp_data(object, ndL[[m]], m = m))
    yX <- fetch(ydat, "x")
    yZt <- fetch(ydat, "Zt")
    yZ_names <- fetch(ydat, "Z_names")
    flist <- lapply(ndL, function(x) factor(x[[id_var]]))
    res <- c(res, nlist(y, yX, yZt, yZ_names, flist))
  }
  if (event_parts) {
    qnodes <- object$quadnodes
    qq <- get_quadpoints(qnodes)
    qtimes <- unlist(lapply(qq$points,  unstandardise_quadpoints,  0, etimes))
    qwts   <- unlist(lapply(qq$weights, unstandardise_quadweights, 0, etimes))
    starttime <- deparse(formula(object, m = "Event")[[2L]][[2L]])
    edat <- prepare_data_table(ndE, id_var, time_var = starttime)
    times <- c(etimes, qtimes) # times used to design event submodel matrices
    edat <- rolling_merge(edat, ids = rep(id_list, qnodes + 1), times = times)
    eXq  <- .pp_data_mer_x(object, newdata = edat, m = "Event")
    if (!object$coxmod_stuff$has_intercept) {
      sel <- grep("^\\(Intercept\\)$", colnames(eXq))
      if (length(sel))
        eXq <- eXq[, -sel, drop = FALSE]
    }
    assoc_parts <- lapply(1:M, function(m) {
      ymf <- prepare_data_table(ndL[[m]], id_var, time_var = time_var)
      make_assoc_parts(
        ymf, assoc = object$assoc, id_var = object$id_var, 
        time_var = object$time_var, id_list = id_list, times = times, 
        use_function = pp_data, object = object, m = m)
    })
    assoc_attr <- nlist(.Data = assoc_parts, qnodes, qtimes, qwts, etimes, estatus)
    assoc_parts <- do.call("structure", assoc_attr)
    res <- c(res, nlist(eXq, assoc_parts))
  }
  return(res)
}
