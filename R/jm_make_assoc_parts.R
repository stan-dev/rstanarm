# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2013, 2014, 2015, 2016, 2017 Trustees of Columbia University
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

# Function to construct quantities, primarily design matrices (x, Zt), that
# will be used to evaluate the longitudinal submodel contributions to the 
# association structure in the event submodel. For example, the design matrices
# evaluated at the quadpoints, quadpoints + eps, lagged quadpoints, auc quadpoints,
# and so on. Exactly what quantities are returned depends on what is specified
# in the use_function argument.
#
# @param use_function The function to call which will return the design 
#   matrices for eta, eps, lag, auc, etc. Generally either 
#   'make_assoc_parts_for_stan' or 'pp_data'.
# @param data A model frame used for constructing the design matrices
# @param assoc A named list returned by a call to validate_assoc (details
#   on the desired association structure for all longitudinal submodels)
# @param id_var The name on the ID variable
# @param time_var The name of the time variable
# @param id_list A vector of subject IDs
# @param times A vector (or possibly a list of vectors) of times at which the 
#   design matrices should be evaluated (most likely the event times and the
#   quadrature times)
# @param eps A numeric value used as the time shift for numerically evaluating
#   the slope of the longitudinal submodel using a one-sided difference
# @param auc_qnodes An integer specifying the number of quadrature nodes to
#   use when evaluating the area under the curve for the longitudinal submodel
# @param ... Additional arguments passes to use_function
# @return A named list
make_assoc_parts <- function(use_function = make_assoc_parts_for_stan, 
                             newdata, assoc, grp_stuff, ids, times, 
                             id_var, time_var, epsilon = 1E-5, 
                             auc_qnodes = 15L, ...) {
  
  if (!requireNamespace("data.table"))
    stop("the 'data.table' package must be installed to use this function")
  
  # Apply lag
  lag <- assoc[["which_lag"]]
  if (!lag == 0)
    times <- set_lag(times, lag)
  
  # Broadcast ids and times if there is lower level clustering
  if (grp_stuff$has_grp) {
    len  <- grp_stuff$qnodes + 1 # num GK quadnodes + 1
    freq <- grp_stuff$grp_freq   # num grps within each patient
    grps <- grp_stuff$grp_list   # unique grp ids
    ids   <- rep(rep(ids, freq), len)   # rep patient id sequence alongside grp ids
    times <- rep(times, rep(freq, len)) # rep time sequence alongside grp ids
    grps  <- rep(grps, len)             # rep grp ids for each quadnode
  } else grps <- NULL
  
  # Identify row in longitudinal data closest to event time or quadrature point
  #   NB if the quadrature point is earlier than the first observation time, 
  #   then covariates values are carried back to avoid missing values.
  #   In any other case, the observed covariates values from the most recent 
  #   observation time preceeding the quadrature point are carried forward to 
  #   represent the covariate value(s) at the quadrature point. (To avoid 
  #   missingness there is no limit on how far forwards or how far backwards 
  #   covariate values can be carried). If no time varying covariates are 
  #   present in the longitudinal submodel (other than the time variable) 
  #   then nothing is carried forward or backward.    
  dataQ <- rolling_merge(data = newdata, ids = ids, times = times, grps = grps)
  mod_eta <- use_function(newdata = dataQ, ...)
  
  # If association structure is based on slope, then calculate design 
  # matrices under a time shift of epsilon
  sel_slope <- grep("etaslope", names(assoc))
  if (any(unlist(assoc[sel_slope]))) {
    dataQ_pos <- dataQ_neg <- dataQ
    dataQ_neg[[time_var]] <- dataQ_neg[[time_var]] - epsilon
    dataQ_pos[[time_var]] <- dataQ_pos[[time_var]] + epsilon
    mod_neg <- use_function(newdata = dataQ_neg, ...)
    mod_pos <- use_function(newdata = dataQ_pos, ...)
    mod_eps <- mod_pos
    mod_eps$x     <- (mod_pos$x     - mod_neg$x    ) / epsilon # derivative of X
    mod_eps$xtemp <- (mod_pos$xtemp - mod_neg$xtemp) / epsilon
    mod_eps$z <- xapply(mod_pos$z, mod_neg$z,                  # derivative of z
                        FUN = function(x, y) (x - y) / epsilon)
  } else mod_eps <- NULL 
  
  # If association structure is based on area under the marker trajectory, then 
  # calculate design matrices at the subquadrature points
  sel_auc <- grep("etaauc|muauc", names(assoc))
  if (any(unlist(assoc[sel_auc]))) {
    if (grp_stuff$has_grp)
      stop2("'etaauc' and 'muauc' not yet implemented when there is a grouping ",
            "factor clustered within patients.")
    # Return a design matrix that is (qnodes * auc_qnodes * Npat) rows 
    auc_qpts <- uapply(times, function(x)
      lapply(get_quadpoints(auc_qnodes)$points, unstandardise_qpts, 0, x))
    auc_qwts <- uapply(times, function(x) 
      lapply(get_quadpoints(auc_qnodes)$weights, unstandardise_qwts, 0, x))
    ids2 <- rep(ids, each = auc_qnodes)
    dataQ_auc <- rolling_merge(data = newdata, ids = ids2, times = auc_qpts)
    mod_auc <- use_function(newdata = dataQ_auc)
  } else mod_auc <- auc_qpts <- auc_qwts <- NULL
  
  # If association structure is based on interactions with data, then calculate 
  # the design matrix which will be multiplied by etavalue, etaslope, muvalue or muslope
  sel_data <- grep("_data", names(assoc), value = TRUE)
  X_data <- xapply(sel_data, FUN = function(i) { 
    form <- assoc[["which_formulas"]][[i]]
    if (length(form)) {
      form <- as.formula(form)
      vars <- rownames(attr(terms.formula(form), "factors"))
      if (is.null(vars))
        stop2("No variables found in the formula for the '", i, "' association structure.")
      sel <- which(!vars %in% colnames(dataQ))
      if (length(sel))
        stop2("The following variables were specified in the formula for the '", i,
              "' association structure, but they cannot be found in the data: ", 
              paste0(vars[sel], collapse = ", "))
      mf <- stats::model.frame(form, data = dataQ)
      X <- stats::model.matrix(form, data = mf)
      X <- drop_intercept(X)
      if (!ncol(X))
        stop2("Bug found: A formula was specified for the '", i, "' association ", 
              "structure, but the resulting design matrix has no columns.")
    } else {
      X <- matrix(0, length(unlist(times)), 0)
    }
    X
  })
  K_data <- sapply(X_data, ncol)
  X_bind_data <- do.call(cbind, X_data)
  
  ret <- nlist(times, mod_eta, mod_eps, mod_auc, K_data, X_data, X_bind_data, grp_stuff)
  
  structure(ret, times = times, lag = lag, epsilon = epsilon, auc_qnodes = auc_qnodes,
            auc_qpts = auc_qpts, auc_qwts = auc_qwts)
}                              

# Carry out a rolling merge
#
# @param data A data.table with a set key corresponding to ids and times
# @param ids A vector of patient ids to merge against
# @param times A vector of (new) times to merge against
# @param grps A vector of cluster ids to merge against when there is clustering
#   within patient ids
# @return A data.table formed by a merge of ids, (grps), times, and the closest 
#   preceding (in terms of times) rows in data
rolling_merge <- function(data, ids, times, grps = NULL) {
  if (!requireNamespace("data.table"))
    stop("the 'data.table' package must be installed to use this function")
  key_length <- if (is.null(grps)) 2L else 3L
  if (!length(data.table::key(data)) == key_length)
    stop("Bug found: data.table key is not the same length as supplied keylist.")
  
  if (is(times, "list") && is.null(grps)) {
    return(do.call(rbind, lapply(times, FUN = function(x, ids) {
      tmp <- data.table::data.table(ids, x)
      data[tmp, roll = TRUE, rollends = c(TRUE, TRUE)]
    }, ids = ids)))     
  } else if (is(times, "list")) {
    return(do.call(rbind, lapply(times, FUN = function(x, ids, grps) {
      tmp <- data.table::data.table(ids, grps, x)
      data[tmp, roll = TRUE, rollends = c(TRUE, TRUE)]
    }, ids = ids, grps = grps)))
  } else if (is.null(grps)) {
    tmp <- data.table::data.table(ids, times)
    return(data[tmp, roll = TRUE, rollends = c(TRUE, TRUE)])       
  } else {
    tmp <- data.table::data.table(ids, grps, times)
    return(data[tmp, roll = TRUE, rollends = c(TRUE, TRUE)])       
  }
}

# Return design matrices for the longitudinal submodel. This is 
# designed to generate the design matrices evaluated at the GK
# quadrature points, because it uses a 'terms' object to generate
# the model frame, and that terms object should have been generated
# from the longitudinal submodel's model frame when it was evaluated
# at the observation times; i.e. the predvars and X_bar would have
# come from the design matrices at the observation times, not the 
# quadrature points.
#
# @param newdata A data frame; the data for the longitudinal submodel 
#   at the quadrature points.
# @param y_mod The list returned by handle_y_mod, containing info about
#   the longitudinal submodel evaluated at the observation (not quadrature)
#   times, for example, the x_bar means used for centering, the predvars 
#   attribute for the longitudinal submodel formula, and so on.
# @param include_Zt Whether to include the sparse Zt matrix in the
#   returned parts.
make_assoc_parts_for_stan <- function(newdata, y_mod, include_Zt = TRUE) {
  
  # construct model frame using predvars
  formula <- use_predvars(y_mod)
  data <- as.data.frame(newdata)
  model_frame <- stats::model.frame(lme4::subbars(formula), data)
  
  # fe design matrices
  x_form <- lme4::nobars(formula)
  x <- model.matrix(x_form, model_frame)
  xtemp <- drop_intercept(x)
  x_bar <- y_mod$x$x_bar
  xtemp <- sweep(xtemp, 2, x_bar, FUN = "-")
  
  # re design matrices
  bars <- lme4::findbars(formula)
  if (length(bars) > 2L)
    stop2("A maximum of 2 grouping factors are allowed.")
  z_parts <- lapply(bars, split_at_bars)
  z_forms <- fetch(z_parts, "re_form")
  z <- lapply(z_forms, model.matrix, model_frame)
  group_vars <- fetch(z_parts, "group_var")
  group_list <- lapply(group_vars, function(x) factor(model_frame[[x]]))
  names(z) <- names(group_list) <- group_vars
  
  ret <- nlist(x, xtemp, z, group_list, group_vars) # return list
  
  # optionally add the sparse Zt matrix
  if (include_Zt) 
    ret$Zt <- lme4::mkReTrms(bars, model_frame)$Zt
  
  ret
}
