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
# @param newdata A model frame used for constructing the design matrices
# @param assoc A list with information about the association structure for 
#   the one longitudinal submodel. 
# @param grp_stuff A list with information about any lower level grouping
#   factors that are clustered within patients and how to handle them in 
#   the association structure.
# @param ids,times The subject IDs and times vectors that correspond to the
#   event/censoring and quadrature times at which the design matrices will
#   need to be evaluated for the association structure.
# @param id_var The name on the ID variable.
# @param time_var The name of the time variable.
# @param epsilon The half-width of the central difference used for 
#   numerically calculating the derivative of the design matrix for slope
#   based association structures.
# @param auc_qnodes Integer specifying the number of GK quadrature nodes to
#   use in the integral/AUC based association structures.
# @param ... Additional arguments passes to use_function
# @return A named list
make_assoc_parts <- function(use_function = make_assoc_parts_for_stan, 
                             newdata, assoc, grp_stuff, ids, times, 
                             id_var, time_var, epsilon = 1E-5, 
                             auc_qnodes = 15L, ...) {
  
  if (!requireNamespace("data.table"))
    stop("the 'data.table' package must be installed to use this function")
  
  eps_uses_derivative_of_x <- TRUE # experimental
  
  # Apply lag
  lag <- assoc[["which_lag"]]
  if (!lag == 0)
    times <- set_lag(times, lag)
  
  # Broadcast ids and times if there is lower level clustering
  if (grp_stuff$has_grp) {
    # grps corresponding to each id
    grps <- as.vector(unlist(grp_stuff$grp_list[as.character(ids)])) 
    # freq by which to expand each ids and times element
    freq_seq <- grp_stuff$grp_freq[as.character(ids)] 
    # rep each patient id and prediction time the required num of times
    ids   <- rep(ids,   freq_seq)       
    times <- rep(times, freq_seq) 
    # indices for collapsing across clusters within patients
    grp_idx <- get_idx_array(freq_seq)  
  } else grps <- grp_idx <- NULL
  
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
    if (eps_uses_derivative_of_x) {
      # slope is evaluated by passing Stan the derivatives of the X and Z
      # design matrices directly, each evaluated using central differences 
      # with a half-width equal to epsilon
      dataQ_pos <- dataQ_neg <- dataQ
      dataQ_neg[[time_var]] <- dataQ_neg[[time_var]] - epsilon
      dataQ_pos[[time_var]] <- dataQ_pos[[time_var]] + epsilon
      mod_neg <- use_function(newdata = dataQ_neg, ...)
      mod_pos <- use_function(newdata = dataQ_pos, ...)
      mod_eps <- mod_pos
      mod_eps$x     <- (mod_pos$x     - mod_neg$x    ) / (2 * epsilon) # derivative of X
      mod_eps$xtemp <- (mod_pos$xtemp - mod_neg$xtemp) / (2 * epsilon)
      mod_eps$z <- xapply(mod_pos$z, mod_neg$z,                  # derivative of z
                          FUN = function(x, y) (x - y) / (2 * epsilon))
      if (!is.null(mod_eps$Zt))
        mod_eps$Zt <- (mod_pos$Zt - mod_neg$Zt) / (2 * epsilon)
    } else {
      # slope is evaluated by passing Stan the X and Z design matrices under
      # a time shift of epsilon and then evaluating the derivative of the
      # linear predictor in Stan using a one-sided difference
      dataQ_eps <- dataQ
      dataQ_eps[[time_var]] <- dataQ_eps[[time_var]] + epsilon
      mod_eps <- use_function(newdata = dataQ_eps, ...)
    }
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
    mod_auc <- use_function(newdata = dataQ_auc, ...)
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
      X <- matrix(0, nrow(dataQ), 0)
    }
    X
  })
  K_data <- sapply(X_data, ncol)
  X_bind_data <- do.call(cbind, X_data)
  
  ret <- nlist(times, mod_eta, mod_eps, mod_auc, K_data, X_data, X_bind_data, grp_stuff)
  
  structure(ret, times = times, lag = lag, epsilon = epsilon, grp_idx = grp_idx,
            auc_qnodes = auc_qnodes, auc_qpts = auc_qpts, auc_qwts = auc_qwts, 
            eps_uses_derivative_of_x = eps_uses_derivative_of_x)
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
#   at the event and quadrature points.
# @param y_mod The list returned by handle_y_mod, containing info about
#   the longitudinal submodel evaluated at the observation (not quadrature)
#   times, for example, the x_bar means used for centering, the predvars 
#   attribute for the longitudinal submodel formula, and so on.
# @param include_Zt Whether to include the sparse Zt matrix in the
#   returned parts.
make_assoc_parts_for_stan <- function(newdata, y_mod, include_Zt = TRUE) {
  
  # construct model frame using predvars
  formula <- use_predvars(y_mod, keep_response = FALSE)
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
