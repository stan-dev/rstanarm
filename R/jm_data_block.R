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

.datatable.aware <- TRUE # necessary for some reason when data.table is in Suggests

#--------------- Miscellaneous and helper functions

# Check input argument is a valid type, and return as a list
#
# @param arg The user input to the argument
# @param type A character vector of valid classes
# @param validate_length The required length of the returned list
# @return A list
validate_arg <- function(arg, type, validate_length = NULL) {
  nm <- deparse(substitute(arg))
  
  if (inherits(arg, type)) { 
    # input type is valid, so return as a list
    arg <- list(arg)
  } 
  else if (is(arg, "list")) { 
    # input type is a list, check each element
    check <- sapply(arg, function(x) inherits(x, type))
    if (!all(check))
      STOP_arg(nm, type)
  } 
  else {
    # input type is not valid
    STOP_arg(nm, type)
  }
  
  if (!is.null(validate_length)) {
    # return list of the specified length
    if (length(arg) == 1L)
      arg <- rep(arg, times = validate_length)
    if (!length(arg) == validate_length)
      stop2(nm, " is a list of the incorrect length.")
  }
  
  if ("data.frame" %in% type)
    arg <- lapply(arg, as.data.frame)
  if ("family" %in% type)
    arg <- lapply(arg, validate_family)
  arg
}

# Check if the user input a list of priors for the longitudinal
# submodel, and if not, then return the appropriate list
#
# @param prior The user input to the prior argument in the stan_mvmer 
#   or stan_jm call
# @param M An integer specifying the number of longitudinal submodels
broadcast_prior <- function(prior, M) {
  if (is.null(prior)) {
    return(rep(list(NULL), M))
  } 
  else if ("dist" %in% names(prior)) {
    return(rep(list(prior), M))
  } 
  else if (is.list(prior) && length(prior) == M) {
    return(prior)
  } 
  else {
    nm <- deparse(substitute(priorarg))
    stop2(nm, " appears to provide prior information separately for the ",
          "different submodels, but the list is of the incorrect length.")
  }
}

# From a vector of length M giving the number of elements (for example number
# of parameters or observations) for each submodel, create an indexing array 
# of dimension M * 2, where column 1 is the beginning index and 2 is the end index
#
# @param x A numeric vector
# @return A length(x) * 2 array
get_idx_array <- function(x) {
  as.array(do.call("rbind", lapply(1:length(x), function(i) {
    idx_beg <- ifelse(x[i] > 0L, sum(x[0:(i-1)]) + 1, 0L)
    idx_end <- ifelse(x[i] > 0L, sum(x[0:i]),         0L)
    c(idx_beg, idx_end)
  })))
}

# Function to return the range or SD of the predictors, used for scaling the priors
# This is taken from an anonymous function in stan_glm.fit
#
# @param x A vector
get_scale_value <- function(x) {
  num.categories <- n_distinct(x)
  x.scale <- 1
  if (num.categories == 2) {
    x.scale <- diff(range(x))
  } else if (num.categories > 2) {
    x.scale <- sd(x)
  }
  return(x.scale)
}

# Apply a lag to a vector of times
#
# @param x A numeric vector (e.g. observation times)
# @param lag A scalar (the lag time)
# @return A numeric vector
set_lag <- function(x, lag) {
  x <- x - lag
  x[x < 0] <- 0.0  # use baseline for lag times prior to baseline
  x
}

# Get the required number of (local) horseshoe parameters for a specified prior type
#
# @param prior_dist An integer indicating the type of prior distribution: 
#   where 1L == normal, 2L == t, 3L == hs, 4L == hs_plus
get_nvars_for_hs <- function(prior_dist) {
  if      (prior_dist <= 2L) return(0L) 
  else if (prior_dist == 3L) return(2L) 
  else if (prior_dist == 4L) return(4L)
  else return(0L)
}

# Reformulate an expression as the LHS of a model formula
# 
# @param x The expression to reformulate
# @return A model formula
reformulate_lhs <- function(x) {
  formula(substitute(LHS ~ 1, list(LHS = x)))
}

# Reformulate an expression as the RHS of a model formula
# 
# @param x The expression to reformulate
# @param subbars A logical specifying whether to call lme4::subbars
#   on the result
# @return A model formula
reformulate_rhs <- function(x, subbars = FALSE) {
  fm <- formula(substitute(~ RHS, list(RHS = x)))
  if (subbars) {
    lme4::subbars(fm)
  } else {
    fm
  }
}

#--------------- Functions related to priors

# Deal with covariance prior
#
# @param prior A list
# @param cnms A list of lists, with names of the group specific 
#   terms for each grouping factor
# @param ok_dists A list of admissible distributions
handle_cov_prior <- function(prior, cnms, ok_dists = nlist("decov", "lkj")) {
  if (!is.list(prior)) 
    stop(sQuote(deparse(substitute(prior))), " should be a named list")
  t <- length(unique(cnms)) # num grouping factors
  p <- sapply(cnms, length) # num terms for each grouping factor
  prior_dist_name <- prior$dist
  if (!prior_dist_name %in% unlist(ok_dists)) {
    stop("The prior distribution should be one of ",
         paste(names(ok_dists), collapse = ", "))
  } else if (prior_dist_name == "decov") {
    prior_shape <- as.array(maybe_broadcast(prior$shape, t))
    prior_scale <- as.array(maybe_broadcast(prior$scale, t))
    prior_concentration <- 
      as.array(maybe_broadcast(prior$concentration, sum(p[p > 1])))
    prior_regularization <- 
      as.array(maybe_broadcast(prior$regularization, sum(p > 1)))
    prior_df <- NULL
  } else if (prior_dist_name == "lkj") {
    prior_shape <- NULL
    prior_scale <- as.array(maybe_broadcast(prior$scale, sum(p)))
    prior_concentration <- NULL
    prior_regularization <- 
      as.array(maybe_broadcast(prior$regularization, sum(p > 1)))
    prior_df <- as.array(maybe_broadcast(prior$df, sum(p)))
  }
  prior_dist <- switch(prior_dist_name, decov = 1L, lkj = 2L)
  
  nlist(prior_dist_name, prior_dist, prior_shape, prior_scale, 
        prior_concentration, prior_regularization, prior_df, t, p,
        prior_autoscale = isTRUE(prior$autoscale))
}  

# Seperate the information about the covariance prior into a list
# of lists. At the top level of the returned list the elements 
# correpond to each of the grouping factors, and on the second level
# of the returned list the elements correpsond to the separate glmer
# submodels. This separation is required for autoscaling the priors 
# on the sds of group level effects, since these are autoscaled based
# on the separate Z matrices (design matrices for the random effects).
#
# @param prior_stuff The named list returned by handle_cov_prior
# @param cnms The component names for group level terms, combined across
#   all glmer submodels
# @param submodel_cnms The component names for the group level terms, 
#   separately for each glmer submodel (stored as a list of length M)
# @return A list with each element containing the covariance prior
#   information for one grouping factor
split_cov_prior <- function(prior_stuff, cnms, submodel_cnms) {
  if (!prior_stuff$prior_dist_name == "lkj") {
    return(prior_stuff) # nothing to be done for decov prior
  } else {
    M <- length(submodel_cnms) # number of submodels
    cnms_nms <- names(cnms) # names of grouping factors
    mark <- 0
    new_prior_stuff <- list()
    for (nm in cnms_nms) {
      for (m in 1:M) {
        len <- length(submodel_cnms[[m]][[nm]])
        new_prior_stuff[[nm]][[m]] <- prior_stuff 
        if (len) {
          # submodel 'm' has group level terms for group factor 'nm'
          beg <- mark + 1; end <- mark + len
          new_prior_stuff[[nm]][[m]]$prior_scale <- prior_stuff$prior_scale[beg:end]
          new_prior_stuff[[nm]][[m]]$prior_df <- prior_stuff$prior_df[beg:end]
          mark <- mark + len
        } else {
          new_prior_stuff[[nm]][[m]]$prior_scale <- NULL
          new_prior_stuff[[nm]][[m]]$prior_df <- NULL
          new_prior_stuff[[nm]][[m]]$prior_regularization <- NULL
        }
      }
    }    
  }
  new_prior_stuff
}

# Autoscaling of priors
#
# @param prior_stuff A named list returned by a call to handle_glm_prior
# @param response A vector containing the response variable, only required if
#   the priors are to be scaled by the standard deviation of the response (for
#   gaussian reponse variables only)
# @param predictors The predictor matrix, only required if the priors are to be
#   scaled by the range/sd of the predictors
# @param family A family object
# @param QR A logical specifying whether QR decomposition is used for the 
#   predictor matrix
# @param assoc A two dimensional array with information about desired association
#   structure for the joint model (returned by a call to validate_assoc). Cannot
#   be NULL if autoscaling priors for the association parameters.
# @param min_prior_scale The minimum allowed for prior scales
# @return A named list with the same structure as returned by handle_glm_prior
autoscale_prior <- function(prior_stuff, response = NULL, predictors = NULL, 
                            family = NULL, QR = FALSE, min_prior_scale = 1e-12, 
                            assoc = NULL) {
  ps <- prior_stuff
  
  if (!is.null(response) && is.gaussian(family)) { 
    # use response variable for scaling priors
    if (ps$prior_dist > 0L && ps$prior_autoscale) {
      ss <- sd(response)
      ps$prior_scale <- ss * ps$prior_scale
    }
  }
  
  if (!is.null(predictors) && !QR) {
    # use predictors for scaling priors
    if (ps$prior_dist > 0L && ps$prior_autoscale) {
      ps$prior_scale <- 
        pmax(min_prior_scale,
             ps$prior_scale / apply(predictors, 2L, get_scale_value))
    }      
  }
  
  if (!is.null(assoc)) {
    # Evaluate mean and SD of each of the association terms that will go into
    # the linear predictor for the event submodel (as implicit "covariates").
    # (NB the approximate association terms are calculated using coefs
    # from the separate longitudinal submodels estimated using glmer).
    # The mean will be used for centering each association term.
    # The SD will be used for autoscaling the prior for each association parameter.
    if (is.null(family))
      stop("'family' cannot be NULL when autoscaling association parameters.")
    assoc_terms <- make_assoc_terms(parts = mod_stuff, assoc = assoc, family = family)
    ps$a_xbar <- as.array(apply(assoc_terms, 2L, mean))
    if (ps$prior_dist > 0L && ps$prior_autoscale) {
      a_beta_scale <- apply(assoc_terms, 2L, get_scale_value)
      ps$prior_scale <- pmax(min_prior_scale, ps$prior_scale / a_beta_scale)
    }
  }
  
  ps$prior_scale <- as.array(pmin(.Machine$double.xmax, ps$prior_scale))
  ps
}

# Create "prior.info" attribute for stan_{mvmer,jm}; needed for prior_summary()
#
# @param user_* The user's priors. These should be passed in after broadcasting 
#   the df/location/scale arguments if necessary.
# @param y_has_intercept Vector of T/F, does each long submodel have an intercept?
# @param y_has_predictors Vector of T/F, does each long submodel have predictors?
# @param y_has_intercept T/F, does event submodel have an intercept?
# @param has_predictors T/F, does event submodel have predictors?
# @param adjusted_prior_*_scale adjusted scales computed if using autoscaled priors
# @param family A list of family objects.
# @return A named list with components 'prior*', 'prior*_intercept', 
#   'prior_covariance' and 'prior*_aux' each of which itself is a list
#   containing the needed values for prior_summary.
summarize_jm_prior <-
  function(user_priorLong = NULL,
           user_priorLong_intercept = NULL,
           user_priorLong_aux = NULL,
           user_priorEvent = NULL,
           user_priorEvent_intercept = NULL,
           user_priorEvent_aux = NULL,
           user_priorEvent_assoc = NULL,
           user_prior_covariance = NULL,
           y_has_intercept = NULL,
           e_has_intercept = NULL,
           y_has_predictors = NULL,
           e_has_predictors = NULL,
           has_assoc = NULL,
           adjusted_priorLong_scale = NULL,
           adjusted_priorLong_intercept_scale = NULL, 
           adjusted_priorLong_aux_scale = NULL,
           adjusted_priorEvent_scale = NULL,
           adjusted_priorEvent_intercept_scale = NULL, 
           adjusted_priorEvent_aux_scale = NULL,           
           adjusted_priorEvent_assoc_scale = NULL,
           family = NULL, 
           basehaz = NULL,
           stub_for_names = "Long") {
    if (!is.null(family) && !is(family, "list"))
      stop("'family' should be a list of family objects, one for each submodel.")
    if (!is.null(has_assoc) && !is.logical(has_assoc) && (length(has_assoc) == 1L))
      stop("'has_assoc' should be a logical vector of length 1.")
    M <- length(family)
    
    prior_list <- list()
    
    if (!is.null(user_priorLong)) {
      rescaled_coefLong <- mapply(check_if_rescaled, user_priorLong, 
                                  y_has_predictors, adjusted_priorLong_scale)
      rescaled_intLong  <- mapply(check_if_rescaled, user_priorLong_intercept, 
                                  y_has_intercept, adjusted_priorLong_intercept_scale)
      rescaled_auxLong  <- mapply(check_if_rescaled, user_priorLong_aux, 
                                  TRUE, adjusted_priorLong_aux_scale) 
      for (m in 1:M) {
        user_priorLong[[m]] <- 
          rename_t_and_cauchy(user_priorLong[[m]], y_has_predictors[m])
        user_priorLong_intercept[[m]] <-
          rename_t_and_cauchy(user_priorLong_intercept[[m]], y_has_intercept[m])
        user_priorLong_aux[[m]] <-
          rename_t_and_cauchy(user_priorLong_aux[[m]], TRUE)
      }
      prior_list$priorLong <- list_nms(lapply(1:M, function(m) {
        if (!y_has_predictors[m]) NULL else with(user_priorLong[[m]], list(
          dist = prior_dist_name,
          location = prior_mean,
          scale = prior_scale,
          adjusted_scale = if (rescaled_coefLong[m])
            adjusted_priorLong_scale[[m]] else NULL,
          df = if (prior_dist_name %in% c
                   ("student_t", "hs", "hs_plus", "lasso", "product_normal"))
            prior_df else NULL
        ))        
      }), M, stub = stub_for_names)
      prior_list$priorLong_intercept <- list_nms(lapply(1:M, function(m) {
        if (!y_has_intercept[m]) NULL else with(user_priorLong_intercept[[m]], list(
          dist = prior_dist_name,
          location = prior_mean,
          scale = prior_scale,
          adjusted_scale = if (rescaled_intLong[m]) 
            adjusted_priorLong_intercept_scale[[m]] else NULL,
          df = if (prior_dist_name %in% "student_t") 
            prior_df else NULL
        ))
      }), M, stub = stub_for_names)      
      aux_name <- lapply(family, .rename_aux)
      prior_list$priorLong_aux <- list_nms(lapply(1:M, function(m) {
        if (is.na(aux_name[[m]])) NULL else with(user_priorLong_aux[[m]], list(
          dist = prior_dist_name,
          location = if (!is.na(prior_dist_name) && 
                         prior_dist_name != "exponential")
            prior_mean else NULL,
          scale = if (!is.na(prior_dist_name) && 
                      prior_dist_name != "exponential")
            prior_scale else NULL,
          adjusted_scale = if (rescaled_auxLong[m])
            adjusted_priorLong_aux_scale[[m]] else NULL,
          df = if (!is.na(prior_dist_name) && 
                   prior_dist_name %in% "student_t")
            prior_df else NULL, 
          rate = if (!is.na(prior_dist_name) && 
                     prior_dist_name %in% "exponential")
            1 / prior_scale else NULL,
          aux_name = aux_name[[m]]
        ))
      }), M, stub = stub_for_names)     
    }
    
    if (!is.null(user_priorEvent)) {
      rescaled_coefEvent <- check_if_rescaled(user_priorEvent, e_has_predictors,
                                              adjusted_priorEvent_scale)
      rescaled_intEvent  <- check_if_rescaled(user_priorEvent_intercept, e_has_intercept, 
                                              adjusted_priorEvent_intercept_scale)
      rescaled_auxEvent  <- check_if_rescaled(user_priorEvent_aux, TRUE, 
                                              adjusted_priorEvent_aux_scale)
      user_priorEvent <- 
        rename_t_and_cauchy(user_priorEvent, e_has_predictors)  
      user_priorEvent_intercept <- 
        rename_t_and_cauchy(user_priorEvent_intercept, e_has_intercept)  
      user_priorEvent_aux <- 
        rename_t_and_cauchy(user_priorEvent_aux, TRUE)     
      prior_list$priorEvent <-
        if (!e_has_predictors) NULL else with(user_priorEvent, list(
          dist = prior_dist_name,
          location = prior_mean,
          scale = prior_scale,
          adjusted_scale = if (rescaled_coefEvent)
            adjusted_priorEvent_scale else NULL,
          df = if (prior_dist_name %in% c
                   ("student_t", "hs", "hs_plus", "lasso", "product_normal"))
            prior_df else NULL
        ))
      prior_list$priorEvent_intercept <-
        if (!e_has_intercept) NULL else with(user_priorEvent_intercept, list(
          dist = prior_dist_name,
          location = prior_mean,
          scale = prior_scale,
          adjusted_scale = if (rescaled_intEvent)
            adjusted_priorEvent_intercept_scale else NULL,
          df = if (prior_dist_name %in% "student_t")
            prior_df else NULL
        ))
      e_aux_name <- .rename_e_aux(basehaz) 
      prior_list$priorEvent_aux <-
        with(user_priorEvent_aux, list(
          dist = prior_dist_name,
          location = prior_mean,
          scale = prior_scale,
          adjusted_scale = if (rescaled_auxEvent)
            adjusted_priorEvent_aux_scale else NULL,
          df = if (!is.na(prior_dist_name) && 
                   prior_dist_name %in% "student_t")
            prior_df else NULL, 
          aux_name = e_aux_name
        ))      
    }
    
    if (!is.null(user_priorEvent_assoc)) {
      rescaled_coefAssoc <- check_if_rescaled(user_priorEvent_assoc, has_assoc, 
                                              adjusted_priorEvent_assoc_scale)
      user_priorEvent_assoc <- rename_t_and_cauchy(user_priorEvent_assoc, has_assoc)        
      prior_list$priorEvent_assoc <-
        if (!has_assoc) NULL else with(user_priorEvent_assoc, list(
          dist = prior_dist_name,
          location = prior_mean,
          scale = prior_scale,
          adjusted_scale = if (rescaled_coefAssoc)
            adjusted_priorEvent_assoc_scale else NULL,
          df = if (prior_dist_name %in% c
                   ("student_t", "hs", "hs_plus", "lasso", "product_normal"))
            prior_df else NULL
        ))
    }
    
    if (length(user_prior_covariance))
      prior_list$prior_covariance <- user_prior_covariance
    
    return(prior_list)
  }

# Get name of auxiliary parameters for event submodel
#
# @param basehaz A list with information about the baseline hazard
.rename_e_aux <- function(basehaz) {
  nm <- basehaz$type_name
  if (nm == "weibull") "weibull-shape" else
    if (nm == "bs") "spline-coefficients" else
      if (nm == "piecewise") "piecewise-coefficients" else NA
}

# Check if priors were autoscaled
#
# @param prior_stuff A list with prior info returned by handle_glm_prior
# @param has A logical checking, for example, whether the model has_predictors, 
#   has_intercept, has_assoc, etc
# @param adjusted_prior_scale The prior scale after any autoscaling
check_if_rescaled <- function(prior_stuff, has, adjusted_prior_scale) {
  prior_stuff$prior_autoscale && has &&
    !is.na(prior_stuff$prior_dist_name) &&
    !all(prior_stuff$prior_scale == adjusted_prior_scale)      
}

# Rename the t prior as being student-t or cauchy
#
# @param prior_stuff A list with prior info returned by handle_glm_prior
# @param has A logical checking, for example, whether the model has_predictors, 
#   has_intercept, has_assoc, etc
rename_t_and_cauchy <- function(prior_stuff, has) {
  if (has && prior_stuff$prior_dist_name %in% "t") {
    if (all(prior_stuff$prior_df == 1)) {
      prior_stuff$prior_dist_name <- "cauchy"
    } else {
      prior_stuff$prior_dist_name <- "student_t"
    }
  }
  return(prior_stuff)
}

#--------------- Functions related to longitudinal submodel

# Construct a list with information on the glmer submodel
#
# @param formula The model formula for the glmer submodel
# @param data The data for the glmer submodel
# @param family The family object for the glmer submodel
# @return A named list with the following elements:
#   Y: named list with the reponse vector and related info
#   X: named list with the fe design matrix and related info
#   Z: named list with the re design matrices and related info
#   terms: the model.frame terms object with bars "|" replaced by "+".
#   family: the modified family object for the glmer submodel
#   intercept_type: named list with info about the type of 
#     intercept required for the glmer submodel
#   has_aux: logical specifying whether the glmer submodel 
#     requires an auxiliary parameter
handle_y_mod <- function(formula, data, family) {
  mf <- stats::model.frame(lme4::subbars(formula), data)
  if (!length(formula) == 3L)
    stop2("An outcome variable must be specified.")
  
  # Response vector, design matrices
  Y <- make_Y(formula, mf, family) 
  X <- make_X(formula, mf, drop_intercept = TRUE, centre = TRUE)
  Z <- make_Z(formula, mf) 
  
  # Terms
  terms <- attr(mf, "terms")
  terms <- append_predvars_attribute(terms, formula, data)
  
  # Binomial with >1 trials not allowed by stan_{mvmver,jm}
  is_binomial <- is.binomial(family$family)
  is_bernoulli <- is_binomial && NCOL(Y) == 1L && all(Y %in% 0:1)
  if (is_binomial && !is_bernoulli)
    STOP_binomial()
  
  # Various flags
  intercept_type <- check_intercept_type(X, family)
  has_aux <- check_for_aux(family)
  family <- append_mvmer_famlink(family, is_bernoulli)
  
  nlist(Y, X, Z, terms, family, intercept_type, has_aux)
}

# Return the response vector
#
# @param formula The model formula
# @param model_frame The model frame
# @param family A family object
# @return A named list with the following elements:
#   Y: the response vector
#   real: the response vector if real, else numeric(0) 
#   integer: the response vector if integer, else integer(0)
#   is_real: a logical indicating whether the response is real
make_Y <- function(formula, model_frame, family) {
  Y <- as.vector(model.response(model_frame))
  Y <- validate_glm_outcome_support(Y, family)
  is_real <- check_response_real(family)
  real <- if (is_real) Y else numeric(0) 
  integer <- if (!is_real) Y else integer(0) 
  nlist(Y, real, integer, is_real)
}

# Return the design matrix, possibly centred
#
# @param formula The model formula
# @param model_frame The model frame
# @param drop_intercept Logical specifying whether to drop the intercept
#   from the returned model matrix
# @param centre Logical specifying whether to centre the predictors
# @return A named list with the following elements:
#   X: the fixed effects model matrix, possibly centred
#   X_bar: The column means of the model matrix
#   has_intercept: Logical flag for whether the submodel has an intercept
#   N,K: number of rows (observations) and columns (predictors) in the
#     fixed effects model matrix
make_X <- function(formula, model_frame = NULL, drop_intercept = TRUE, 
                   centre = TRUE) {
  X_formula <- lme4::nobars(formula)
  X <- model.matrix(X_formula, model_frame)
  has_intercept <- check_for_intercept(X, logical = TRUE)
  if (drop_intercept)
    X <- drop_intercept(X)
  if (centre) {
    if (!drop_intercept)
      stop2("Cannot centre 'x' without dropping the intercept.")
    X_bar <- colMeans(X)
    X <- sweep(X, 2, X_bar, FUN = "-")
  } else {
    X_bar <- rep(0, ncol(X))
  }
  # drop any column of x with < 2 unique values (empty interaction levels)
  sel <- (2 > apply(X, 2L, function(x) length(unique(x))))
  if (any(sel)) {
    warning("Dropped empty interaction levels: ",
            paste(colnames(X)[sel], collapse = ", "))
    X <- X[, !sel, drop = FALSE]
    X_bar <- X_bar[!sel]
  }    
  nlist(X, X_bar, has_intercept, N = NROW(X), K = NCOL(X))
}

# Return the design matrices for the group level terms
#
# @param formula The model formula
# @param model_frame The model frame
# @return A named list with the following elements:
#   Z: a list with each element containing the random effects model 
#     matrix for one grouping factor
#   group_vars: a character vector with the name of each of the
#     grouping factors
#   group_cnms: a list with each element containing the names of the
#     group level parameters for one grouping factor
#   group_list: a list with each element containing the vector of group 
#     IDs for the rows of Z
#   nvars: a vector with the number of group level parameters for each
#     grouping factor
#   ngrps: a vector with the number of groups for each grouping factor 
make_Z <- function(formula, model_frame = NULL) {
  bars <- lme4::findbars(formula)
  if (length(bars) > 2L)
    stop2("A maximum of 2 grouping factors are allowed.")
  Z_parts <- lapply(bars, split_at_bars)
  Z_forms <- fetch(Z_parts, "re_form")
  Z <- lapply(Z_forms, model.matrix, model_frame)
  group_cnms <- lapply(Z, colnames)
  group_vars <- fetch(Z_parts, "group_var")
  group_list <- lapply(group_vars, function(x) model_frame[[x]])
  nvars <- sapply(group_cnms, length)
  ngrps <- sapply(group_list, n_distinct)
  names(Z) <- names(Z_forms) <- names(group_cnms) <- 
    names(group_list) <- names(nvars) <- names(ngrps) <- group_vars
  nlist(Z, Z_forms, group_vars, group_cnms, group_list, nvars, ngrps)
}

# Return info on the required type of intercept
#
# @param X The model matrix
# @param family A family object
# @return A character string specifying the type of bounds 
#   required for the intercept term
check_intercept_type <- function(X, family) {
  fam <- family$family
  link <- family$link
  if (!X$has_intercept) { # no intercept
    type <- "none"
    needs_intercept <- 
      (!is.gaussian(fam) && link == "identity") ||
      (is.gamma(fam) && link == "inverse") ||
      (is.binomial(fam) && link == "log")
    if (needs_intercept)
      stop2("To use the specified combination of family and link (", fam, 
            ", ", link, ") the model must have an intercept.")
  } else if (fam == "binomial" && link == "log") { # binomial, log
    type <- "upper_bound" 
  } else if (fam == "binomial") { # binomial, !log
    type <- "no_bound"
  } else if (link == "log") { # gamma/inv-gaus/poisson/nb, log
    type <- "no_bound"  
  } else if (fam == "gaussian") { # gaussian, !log
    type <- "no_bound"  
  } else { # gamma/inv-gaus/poisson/nb, !log 
    type <- "lower_bound"  
  }
  number <- switch(type, none = 0L, no_bound = 1L,
                   lower_bound = 2L, upper_bound = 3L)
  nlist(type, number) 
}

# Check the family and link function are supported by stan_{mvmer,jm}
#
# @param family A family object
# @param supported_families A character vector of supported family names
# @return A family object
validate_famlink <- function(family, supported_families) {
  famname <- family$family
  fam <- which(supported_families == famname)
  if (!length(fam)) 
    stop2("'family' must be one of ", paste(supported_families, collapse = ", "))
  supported_links <- supported_glm_links(famname)
  link <- which(supported_links == family$link)
  if (!length(link)) 
    stop("'link' must be one of ", paste(supported_links, collapse = ", "))
  return(family)
}

# Append a family object with numeric family and link information used by Stan
#
# @param family The existing family object
# @param is_bernoulli Logical specifying whether the family should be bernoulli
# @return A family object with two appended elements: 
#   mvmer_family: an integer telling Stan which family
#   mvmer_link: an integer telling Stan which link function (varies by family!)
append_mvmer_famlink <- function(family, is_bernoulli = FALSE) {
  famname <- family$family
  family$mvmer_family <- switch(
    famname, 
    gaussian = 1L, 
    Gamma = 2L,
    inverse.gaussian = 3L,
    binomial = 5L, # bernoulli = 4L changed later
    poisson = 6L,
    "neg_binomial_2" = 7L)
  if (is_bernoulli)
    family$mvmer_family <- 4L
  supported_links <- supported_glm_links(famname)
  link <- which(supported_links == family$link)
  family$mvmer_link <- link
  return(family)
}

# Split the random effects part of a model formula into
#   - the formula part (ie. the formula on the LHS of "|"), and 
#   - the name of the grouping factor (ie. the variable on the RHS of "|")
#
# @param x Random effects part of a model formula, as returned by lme4::findbars
# @return A named list with the following elements:
#   re_form: a formula specifying the random effects structure
#   group_var: the name of the grouping factor
split_at_bars <- function(x) {
  terms <- strsplit(deparse(x), "\\s\\|\\s")[[1L]]
  if (!length(terms) == 2L)
    stop2("Could not parse the random effects formula.")
  re_form <- formula(paste("~", terms[[1L]]))
  group_var <- terms[[2L]]
  nlist(re_form, group_var)
}

# Function to check if the response vector is real or integer
#
# @param family A family object
# @return A logical specify whether the response is real (TRUE) or integer (FALSE)
check_response_real <- function(family) {
  !(family$family %in% c("binomial", "poisson", "neg_binomial_2"))
}

# Function to check if the submodel should include a auxiliary term
#
# @param family A family object
# @return A logical specify whether the submodel includes a auxiliary term
check_for_aux <- function(family) {
  !(family$family %in% c("binomial", "poisson"))
}

# Function to return a single cnms object for all longitudinal submodels
#
# @param x A list, with each element being a cnms object returned by (g)lmer
get_common_cnms <- function(x, stub = "Long") {
  nms <- lapply(x, names)
  unique_nms <- unique(unlist(nms))
  cnms <- lapply(seq_along(unique_nms), function(i) {
    nm <- unique_nms[i]
    unlist(lapply(1:length(x), function(m) 
      if (nm %in% nms[[m]]) paste0(stub, m, "|", x[[m]][[nm]])))
  })
  names(cnms) <- unique_nms
  cnms
}

# Take a list of cnms objects (each element containing the cnms for one 
# submodel) and assess whether the specified variable is included as a 
# grouping factor in all of the submodels
#
# @param y_cnms A list with each element containing the cnms object for
#   one submodel.
# @param group_var The name of the grouping factor variable.
# @return The name of the grouping factor, or an error if it doesn't 
#   appear in every submodel.
validate_grouping_factor <- function(y_cnms, group_var) {
  check <- sapply(y_cnms, function(x) group_factor %in% names(x))
  if (!all(check)) {
    nm <- deparse(substitute(group_var))
    stop2(nm, " must be a grouping factor in all longitudinal submodels.")
  }
  group_factor
}

# Check the factor list corresponding to subject ID is the same in each 
# of the longitudinal submodels
#
# @param id_var The name of the ID variable
# @param y_flist A list containing the flist objects returned for each 
#   separate longitudinal submodel
# @return A vector of factor levels corresponding to the IDs appearing
#   in the longitudinal submodels
check_id_list <- function(id_var, y_flist) {
  id_list <- unique(lapply(y_flist, function(x) levels(x[[id_var]])))
  if (length(id_list) > 1L)
    stop2("The subject IDs are not the same in all longitudinal submodels.")
  unlist(id_list)  
}

# Take the model frame terms object and append with attributes
# that provide the predvars for the fixed and random effects 
# parts, based on the model formula and data
#
# @param terms The existing model frame terms object
# @param formula The formula that was used to build the model frame
#   (but prior to having called lme4::subbars on it!)
# @param data The data frame that was used to build the model frame
# @return A terms object with predvars.fixed and predvars.random as
#   additional attributes
append_predvars_attribute <- function(terms, formula, data) {
  fe_form <- lme4::nobars(formula)
  re_form <- lme4::subbars(justRE(formula, response = TRUE))
  fe_frame <- stats::model.frame(fe_form, data)
  re_frame <- stats::model.frame(re_form, data)
  fe_terms <- attr(fe_frame, "terms")
  re_terms <- attr(re_frame, "terms")
  fe_predvars <- attr(fe_terms, "predvars")
  re_predvars <- attr(re_terms, "predvars")
  attr(terms, "predvars.fixed")  <- attr(fe_terms, "predvars")
  attr(terms, "predvars.random") <- attr(re_terms, "predvars")
  terms
}

# Function to substitute variables in the formula of a fitted model
# with the corresponding predvars based on the terms object for the model.
# (This is useful since lme4::glFormula doesn't allow a terms object to be 
# passed as the first argument instead of a model formula).
#
# @param mod A (g)lmer model object from which to extract the formula and terms
# @return A reformulated model formula with variables replaced by predvars
use_predvars <- function(mod) {
  fm <- formula(mod)
  ff <- grep("", attr(terms(mod, fixed.only = TRUE), "variables"), value = TRUE)[-1]
  fr <- grep("", attr(terms(mod, random.only = TRUE), "variables"), value = TRUE)[-1]
  pf <- grep("", attr(terms(mod, fixed.only = TRUE), "predvars"), value = TRUE)[-1]
  pr <- grep("", attr(terms(mod, random.only = TRUE), "predvars"), value = TRUE)[-1]
  if (!identical(c(ff, fr), c(pf, pr))) {
    for (j in 1:length(ff))
      fm <- gsub(ff[[j]], pf[[j]], fm, fixed = TRUE)    
    for (j in 1:length(fr))
      fm <- gsub(fr[[j]], pr[[j]], fm, fixed = TRUE)    
    fm <- reformulate(fm[[3L]], response = formula(mod)[[2L]])
  }
  fm
}

#--------------- Functions related to event submodel

# Construct a list with information on the event submodel
#
# @param formula The model formula for the event submodel
# @param data The data for the event submodel
# @param qnodes An integer specifying the number of GK quadrature nodes
# @param id_var The name of the ID variable
# @param y_id_list A character vector with a unique list of subject IDs 
#   (factor levels) that appeared in the longitudinal submodels
# @return A named list with the following elements:
#
handle_e_mod <- function(formula, data, qnodes, id_var, y_id_list) {
  if (!requireNamespace("survival"))
    stop("the 'survival' package must be installed to use this function")
  if (!requireNamespace("data.table"))
    stop("the 'data.table' package must be installed to use this function")
  
  mod <- survival::coxph(formula, data = data, x = TRUE)
  mf1 <- expand.model.frame(mod, id_var, na.expand = TRUE)
  #mf1[[id_var]] <- promote_to_factor(mf1[[id_var]]) # same as lme4
  mf2 <- unclass_Surv_column(mf1) 
  if (attr(mod$y, "type") == "counting") {
    tvc <- TRUE; t0_var <- "start"; t1_var <- "stop"
  } else if (attr(mod$y, "type") == "right") {
    tvc <- FALSE; t0_var <- "time"; t1_var <- "time"
  } else {
    stop2("Only 'right' or 'counting' type Surv objects are allowed ", 
          "on the LHS of 'formulaEvent'.")
  }
  
  # Split model frame and find event time and status
  mf_by_id <- split(mf2, mf2[, id_var])
  mf_entry <- do.call(rbind, lapply(
    mf_by_id, FUN = function(x) x[which.min(x[, t0_var]), ]))
  mf_event <- do.call(rbind, lapply(
    mf_by_id, FUN = function(x) x[which.max(x[, t1_var]), ]))
  entrytime <- mf_entry[[t0_var]]
  if (tvc && (any(entrytime) > 0))
    warning("Note that delayed entry is not yet implemented. It will ",
            "be assumed that all individuals were at risk from time 0.")
  entrytime <- rep(0, length(entrytime)) # no delayed entry
  eventtime <- mf_event[[t1_var]]
  status    <- mf_event[["status"]]  
  id_list   <- mf_event[[id_var]]
  names(entrytime) <- names(eventtime) <- names(status) <- ids
  
  # Mean log incidence rate - used for shifting log baseline hazard
  norm_const <- log(sum(status) / sum(eventtime))
  
  # Error checks for the ID variable
  if (!identical(y_id_list, levels(factor(id_list))))
    stop2("The patient IDs (levels of the grouping factor) included ",
          "in the longitudinal and event submodels do not match")
  if (is.unsorted(factor(id_list)))
    stop2("'dataEvent' needs to be sorted by the subject ",
          "ID/grouping variable")
  if (!identical(length(y_id_list), length(id_list)))
    stop2("The number of patients differs between the longitudinal and ",
          "event submodels. Perhaps you intended to use 'start/stop' notation ",
          "for the Surv() object.")
  
  # Quadrature weights/times/ids
  qq <- get_quadpoints(qnodes)
  qwts <- uapply(qq$weights, unstandardise_qwts, entrytime, eventtime)
  qpts <- uapply(qq$points, unstandardise_qpts, entrytime, eventtime)
  qids <- rep(id_list, qnodes)
  
  # Event times/ids (for failures only)
  epts <- eventtime[status == 1] # event times (for failures only)
  eids <- id_list[status == 1]   # subject ids (for failures only)
  
  # Both event times/ids and quadrature times/ids
  cpts <- c(epts, qpts)
  cids <- unlist(list(eids, qids)) # NB using c(.) demotes factors to integers
  
  # Evaluate design matrix at event and quadrature times
  if (ncol(mod$x)) {
    # Convert model frame from Cox model into a data.table
    mf2 <- data.table::data.table(mf2, key = c(id_var, t0_var))
    mf2[[t0_var]] <- as.numeric(mf2[[t0_var]])
    # Obtain rows of the model frame that are as close as possible to 
    # the event times (failures only) and quadrature times                      
    mf2_q <- mf2[data.table::data.table(cids, cpts), 
                 roll = TRUE, rollends = c(TRUE, TRUE)]
    # Construct design matrix evaluated at event and quadrature times
    fm_RHS <- reformulate(attr(terms(mod), "term.labels"))
    Xq <- model.matrix(fm_RHS, data = mf2_q)
    Xq <- Xq[, -1L, drop = FALSE] # drop intercept
    # Centre the design matrix
    Xbar <- colMeans(Xq)
    Xq <- sweep(Xq, 2, Xbar, FUN = "-")
    sel <- (2 > apply(Xq, 2L, function(x) length(unique(x))))
    if (any(sel)) {
      # drop any column of x with < 2 unique values (empty interaction levels)
      warning("Dropped empty interaction levels: ",
              paste(colnames(Xq)[sel], collapse = ", "))
      Xq <- Xq[, !sel, drop = FALSE]
      Xbar <- Xbar[!sel]
    }
  } else {
    Xq <- matrix(0,0L,0L)
    Xbar <- rep(0,0L)
  }
  
  nlist(mod, entrytime, eventtime, status, Npat = length(eventtime), 
        Nevents = sum(status), id_list, qnodes, qwts, qpts, qids, 
        epts, eids, cpts, cids, Xq, Xbar, K = ncol(Xq), norm_const, 
        model_frame = mf1, tvc)
}

# Deal with the baseline hazard
#
# @param basehaz A string specifying the type of baseline hazard
# @param basehaz_ops A named list with elements df, knots 
# @param ok_basehaz A list of admissible baseline hazards
# @param eventtime A numeric vector with eventtimes for each individual
# @param status A numeric vector with event indicators for each individual
# @return A named list with the following elements:
#
handle_basehaz <- function(basehaz, basehaz_ops, 
                           ok_basehaz = nlist("weibull", "bs", "piecewise"),
                           ok_basehaz_ops = nlist("df", "knots"),
                           eventtime, status) {
  
  if (!basehaz %in% unlist(ok_basehaz))
    stop("The baseline hazard should be one of ", paste(names(ok_basehaz), collapse = ", "))
  if (!all(names(basehaz_ops) %in% unlist(ok_basehaz_ops)))
    stop("The baseline hazard options list can only include ", paste(names(ok_basehaz_ops), collapse = ", "))
  
  type <- switch(basehaz, weibull = 1L, bs = 2L, piecewise = 3L)
  type_name <- basehaz
  user_df   <- basehaz_ops$df
  df        <- basehaz_ops$df
  knots     <- basehaz_ops$knots
  bs_basis  <- NULL
  
  if (type_name == "weibull") {
    # handle df and knots
    if (!is.null(df))
      warning("'df' will be ignored since baseline hazard was set to weibull.", 
              immediate. = TRUE, call. = FALSE)
    if (!is.null(knots))
      warning("'knots' will be ignored since baseline hazard was set to weibull.", 
              immediate. = TRUE, call. = FALSE) 
    user_df <- NULL
    df      <- 1L
    knots   <- NULL
  } else if (type_name %in% c("bs", "piecewise")) {
    # handle df and knots
    if (!any(is.null(df), is.null(knots))) { 
      # both specified
      stop("Cannot specify both 'df' and 'knots' for the baseline hazard.", call. = FALSE)
    } else if (all(is.null(df), is.null(knots))) { 
      # both null -- use default df
      user_df <- df <- 6L
      knots <- NULL
    } else if (!is.null(df)) { 
      # only df specified
      if (type == 2L) {
        if (df < 3) stop("'df' must be at least 3 for B-splines baseline hazard.")
        user_df <- df <- df + 1
      }
    } else if (!is.null(knots)) {          
      # only knots specified
      if (!is.numeric(knots)) stop("'knots' vector must be numeric", call. = FALSE)
      if (any(knots < 0)) stop("'knots' must be non-negative", call. = FALSE)      
      if (type == 2L) df <- length(knots) + 4
      else if (type == 3L) df <- length(knots) + 1
    } else {
      stop("Bug found: unable to reconcile 'df' and 'knots' arguments.", call. = FALSE) 
    }
  }  
  
  # Evaluate spline basis (knots, df, etc) based on distribution of observed event times
  # or evaluate cut points for piecewise constant baseline hazard
  if (type == 2L) {
    bs_basis <- splines::bs(eventtime[(status > 0)], df = user_df, knots = knots, 
                            Boundary.knots = c(0, max(eventtime)), intercept = TRUE)
  } else if (type == 3L) {
    if (is.null(knots)) {
      knots <- quantile(eventtime[(status > 0)], probs = seq(0, 1, 1 / df))
      knots[[1]] <- 0
      knots[[length(knots)]] <- max(eventtime)
    } else {
      if (any(knots > max(eventtime)))
        stop("'knots' for the baseline hazard cannot be greater than the ",
             "largest event time.", call. = FALSE)
      knots <- c(0, knots, max(eventtime))
    }
  }  
  
  nlist(type, type_name, user_df, df, knots, bs_basis)   
}

# Return the design matrix for the baseline hazard
#
# @param times A vector of times at which to evaluate the baseline hazard
# @param basehaz A named list with info about the baseline hazard,
#   returned by a call to handle_basehaz
# @return A matrix
make_basehaz_X <- function(times, basehaz_info) {
  if (basehaz$type_name == "weibull") {
    X <- matrix(log(times), nrow = length(times), ncol = 1) 
  } else if (basehaz$type_name == "bs") {
    basis <- basehaz$bs_basis
    if (is.null(basis))
      stop2("Bug found: could not find info on B-splines basis terms.")
    X <- as.array(predict(basis, times)) 
  } else if (basehaz$type_name == "piecewise") {
    knots <- basehaz$knots
    df <- basehaz$df
    if (is.null(knots) || is.null(df))
      stop2("Bug found: could not find info on basehaz df and knot locations.")
    times_quantiles <- cut(times, knots, include.lowest = TRUE, labels = FALSE)
    X <- matrix(NA, length(times_quantiles), df)
    for (i in 1:df) 
      X[, i] <- ifelse(times_quantiles == i, 1, 0)
    X <- as.array(X)
  } else {
    stop2("Bug found: type of baseline hazard unknown.") 
  }
  X
}

# Function to return standardised GK quadrature points and weights
#
# @param nodes The required number of quadrature nodes
# @return A list with two named elements (points and weights) each
#   of which is a numeric vector with length equal to the number of
#   quadrature nodes
get_quadpoints <- function(nodes = 15) {
  if (!is.numeric(nodes) || (length(nodes) > 1L)) {
    stop("'qnodes' should be a numeric vector of length 1.")
  } else if (nodes == 15) {
    list(
      points = c(
        -0.991455371120812639207,
        -0.949107912342758524526,
        -0.86486442335976907279,
        -0.7415311855993944398639,
        -0.5860872354676911302941,
        -0.4058451513773971669066,
        -0.2077849550078984676007,
        0,
        0.2077849550078984676007,
        0.405845151377397166907,
        0.5860872354676911302941,
        0.741531185599394439864,
        0.86486442335976907279,
        0.9491079123427585245262,
        0.991455371120812639207),
      weights = c(
        0.0229353220105292249637,
        0.063092092629978553291,
        0.10479001032225018384,
        0.140653259715525918745,
        0.1690047266392679028266,
        0.1903505780647854099133,
        0.204432940075298892414,
        0.209482141084727828013,
        0.204432940075298892414,
        0.1903505780647854099133,
        0.169004726639267902827,
        0.140653259715525918745,
        0.1047900103222501838399,
        0.063092092629978553291,
        0.0229353220105292249637))      
  } else if (nodes == 11) {
    list(
      points = c(
        -0.984085360094842464496,
        -0.906179845938663992798,
        -0.754166726570849220441,
        -0.5384693101056830910363,
        -0.2796304131617831934135,
        0,
        0.2796304131617831934135,
        0.5384693101056830910363,
        0.754166726570849220441,
        0.906179845938663992798,
        0.984085360094842464496),
      weights = c(
        0.042582036751081832865,
        0.1152333166224733940246,
        0.186800796556492657468,
        0.2410403392286475866999,
        0.272849801912558922341,
        0.2829874178574912132043,
        0.272849801912558922341,
        0.241040339228647586701,
        0.186800796556492657467,
        0.115233316622473394025,
        0.042582036751081832865))     
  } else if (nodes == 7) {
    list(
      points = c(
        -0.9604912687080202834235,
        -0.7745966692414833770359,
        -0.4342437493468025580021,
        0,
        0.4342437493468025580021,
        0.7745966692414833770359,
        0.9604912687080202834235),
      weights = c(
        0.1046562260264672651938,
        0.268488089868333440729,
        0.401397414775962222905,
        0.450916538658474142345,
        0.401397414775962222905,
        0.268488089868333440729,
        0.104656226026467265194))      
  } else stop("'qnodes' must be either 7, 11 or 15.")  
}

# Remove the "Surv" class attribute from the first column 
# of the model frame after a survival::coxph call
#
# @param data A model frame with the first column being the Surv() response
unclass_Surv_column <- function(data) {
  cbind(unclass(mf1[,1]), mf1[, -1, drop = FALSE], stringsAsFactors = FALSE)
}

#--------------- Functions related to association structure

# Return a named list with information about the specified association structure 
# 
# @param user_x A character vector or NULL, being the user input to the
#   assoc argument (for one submodel) in the stan_jm call
# @param y_mod_stuff A list returned by a call to handle_glmod
# @param id_var The name of the ID variable 
# @param M Integer specifying the total number of longitudinal submodels
# @return A list with information about the desired association structure
validate_assoc <- function(user_x, y_mod_stuff, ok_assoc, ok_assoc_data,
                           ok_assoc_interactions, lag, id_var, M) {
  
  ok_inputs <- c(ok_assoc, paste0(ok_assoc_data, "_data"),
                 unlist(lapply(ok_assoc_interactions, paste0, "_", ok_assoc_interactions))) 
  
  # Check user input to assoc argument
  trimmed_x <- trim_assoc(user_x, ok_assoc_data, ok_assoc_interactions)
  if (is.null(user_x) || all(trimmed_x %in% ok_inputs)) {
    assoc <- sapply(ok_inputs, `%in%`, trimmed_x, simplify = FALSE)
    if (is.null(user_x)) {
      assoc$null <- TRUE
    } else if (is.vector(user_x) && is.character(user_x)) {
      if ((assoc$null) && (length(user_x) > 1L))
        stop("In assoc, 'null' cannot be specified in conjuction ",
             "with another association type", call. = FALSE)
      STOP_combination_not_allowed(assoc, "etavalue", "muvalue")
      STOP_combination_not_allowed(assoc, "etaslope", "muslope")
      STOP_combination_not_allowed(assoc, "etaauc",   "muauc")
    } else {
      stop("'assoc' argument should be a character vector or, for a multivariate ",
           "joint model, possibly a list of character vectors.", call. = FALSE)    
    }    
  } else {
    stop("An unsupported association type has been specified. The ",
         "'assoc' argument can only include the following association ", 
         "types: ", paste(ok_assoc, collapse = ", "), ", as well as ",
         "possible interactions either between association terms or ",
         "with observed data.", call. = FALSE)  
  }
  
  # Parse suffix specifying indices for shared random effects
  cnms <- y_mod_stuff$cnms[[id_var]] # names of random effect terms
  assoc$which_b_zindex    <- parse_assoc_sharedRE("shared_b",    user_x, max_index = length(cnms), cnms)
  assoc$which_coef_zindex <- parse_assoc_sharedRE("shared_coef", user_x, max_index = length(cnms), cnms)
  
  if (length(intersect(assoc$which_b_zindex, assoc$which_coef_zindex)))
    stop("The same random effects indices should not be specified in both ",
         "'shared_b' and 'shared_coef'. Specifying indices in 'shared_coef' ",
         "will include both the fixed and random components.", call. = FALSE)
  
  if (length(assoc$which_coef_zindex)) {
    if (length(y_mod_stuff$cnms) > 1L)
      stop("'shared_coef' association structure cannot be used when there is ",
           "clustering at levels other than the individual-level.", call. = FALSE)
    b_nms <- names(assoc$which_coef_zindex)
    assoc$which_coef_xindex <- sapply(b_nms, function(y, beta_nms) {
      beta_match <- grep(y, beta_nms, fixed = TRUE)
      if (!length(beta_match)) {
        stop("In association structure 'shared_coef', no matching fixed effect ",
             "component could be found for the following random effect: ", y, 
             ". Perhaps consider using 'shared_b' association structure instead.")
      } else if (length(beta_match) > 1L) {
        stop("Bug found: In association structure 'shared_coef', multiple ",
             "fixed effect components have been found to match the following ",
             "random effect: ", y)
      }  
      beta_match
    }, beta_nms = colnames(y_mod_stuff$x))
  } else assoc$which_coef_xindex <- numeric(0)
  
  if (!identical(length(assoc$which_coef_zindex), length(assoc$which_coef_xindex)))
    stop("Bug found: the lengths of the fixed and random components of the ",
         "'shared_coef' association structure are not the same.")
  
  # Parse suffix specifying formula for interactions with data
  ok_inputs_data <- paste0(ok_assoc_data, "_data")
  assoc$which_formulas <- sapply(ok_inputs_data, parse_assoc_data, user_x, simplify = FALSE) 
  
  # Parse suffix specifying indices for interactions between association terms
  ok_inputs_interactions <- unlist(lapply(ok_assoc_interactions, paste0, "_", ok_assoc_interactions))
  assoc$which_interactions <- sapply(ok_inputs_interactions, parse_assoc_interactions, 
                                     user_x, max_index = M, simplify = FALSE)
  
  # Lag for association structure
  assoc$which_lag <- lag
  
  assoc
}

# Validate the user input to the lag_assoc argument of stan_jm
#
# @param lag_assoc The user input to the lag_assoc argument
# @param M Integer specifying the number of longitudinal submodels
# @return Nothing
validate_lag_assoc <- function(lag_assoc, M) {
  if (length(lag_assoc) == 1L)
    lag_assoc <- rep(lag_assoc, M)
  if (!length(lag_assoc) == M)
    stop2("'lag_assoc' should length 1 or length equal to the ",
          "number of markers (", M, ").")
  if (!is.numeric(lag_assoc))
    stop2("'lag_assoc' must be numeric.")
  if (any(lag_assoc < 0))
    stop2("'lag_assoc' must be non-negative.")
}

# Remove suffixes from the user inputted assoc argument
#
# @param x A character vector, being the user input to the 
#   assoc argument in the stan_jm call
# @param ok_assoc_data A character vector specifying which types
#   of association terms are allowed to be interacted with data
# @param ok_assoc_interactions A character vector specifying which types
#   of association terms are allowed to be interacted with other 
#   association terms
trim_assoc <- function(x, ok_assoc_data, ok_assoc_interactions) {
  x <- gsub("^shared_b\\(.*",    "shared_b",    x) 
  x <- gsub("^shared_coef\\(.*", "shared_coef", x) 
  for (i in ok_assoc_data)
    x <- gsub(paste0("^", i, "_data\\(.*"),    paste0(i, "_data"), x)
  for (i in ok_assoc_interactions) for (j in ok_assoc_interactions)
    x <- gsub(paste0("^", i, "_", j, "\\(.*"), paste0(i, "_", j),  x) 
  x     
}

# Parse the formula for specifying a data interaction with an association term
#
# @param x A character string corresponding to one of the allowed
#   association structures for interactions with data, for example, 
#   "etavalue_data" or "etaslope_data"
# @param user_x A character vector, being the user input to the assoc
#   argument in the stan_jm call
# @return The parsed formula (which can be used for constructing a 
#   design matrix for interacting data with association type x) or NULL
parse_assoc_data <- function(x, user_x) {
  val <- grep(paste0("^", x, ".*"), user_x, value = TRUE)
  if (length(val)) {
    val2 <- unlist(strsplit(val, x))[-1]
    fm <- tryCatch(eval(parse(text = val2)), error = function(e) 
      stop(paste0("Incorrect specification of the formula in the '", x,
                  "' association structure. See Examples in the help file."), call. = FALSE))
    if (!is(fm, "formula"))
      stop(paste0("Suffix to '", x, "' association structure should include ",
                  "a formula within parentheses."), call. = FALSE)
    if (identical(length(fm), 3L))
      stop(paste0("Formula specified for '", x, "' association structure should not ",
                  "include a response."), call. = FALSE)
    if (length(lme4::findbars(fm)))
      stop(paste0("Formula specified for '", x, "' association structure should only ",
                  "include fixed effects."), call. = FALSE)
    if (fm[[2L]] == 1)
      stop(paste0("Formula specified for '", x, "' association structure cannot ",
                  "be an intercept only."), call. = FALSE)
    return(fm)
  } else numeric(0)
}

# Parse the indices specified for shared random effects
#
# @param x A character string corresponding to one of the allowed
#   association structures for shared random effects
# @param user_x A character vector, being the user input to the assoc
#   argument in the stan_jm call
# @param max_index An integer specifying the total number of random effects
#   in the longitudinal submodel, and therefore the maximum allowed index for
#   the shared random effects
# @param cnms The names of the random effects corresponding to the 
#   individual-level (id_var) of clustering
# @return A numeric vector specifying indices for the shared random effects
parse_assoc_sharedRE <- function(x, user_x, max_index, cnms) {
  val <- grep(paste0("^", x, ".*"), user_x, value = TRUE)
  if (length(val)) {
    val2 <- unlist(strsplit(val, x))[-1]
    if (length(val2)) {
      index <- tryCatch(eval(parse(text = paste0("c", val2))), error = function(e) 
        stop("Incorrect specification of the '", x, "' association structure. ",
             "See Examples in help file.", call. = FALSE))
      if (any(index > max_index))
        stop(paste0("The indices specified for the '", x, "' association structure are ",
                    "greater than the number of subject-specific random effects."), call. = FALSE)
    } else index <- seq_len(max_index)
    names(index) <- cnms[index]
    return(index)   
  } else numeric(0)
}

# Parse the indices specified for interactions between association terms
#
# @param x A character string corresponding to one of the allowed
#   association structures
# @param user_x A character vector, being the user input to the assoc
#   argument in the stan_jm call
# @param max_index An integer specifying the maximum allowed index
# @return A numeric vector specifying indices
parse_assoc_interactions <- function(x, user_x, max_index) {
  val <- grep(paste0("^", x, ".*"), user_x, value = TRUE)
  if (length(val)) {
    val2 <- unlist(strsplit(val, x))[-1]
    if (length(val2)) {
      index <- tryCatch(eval(parse(text = paste0("c", val2))), error = function(e) 
        stop("Incorrect specification of the '", x, "' association structure. It should ",
             "include a suffix with parentheses specifying the indices of the association ",
             "terms you want to include in the interaction. See Examples in the help file.", call. = FALSE))
      if (any(index > max_index))
        stop("The indices specified for the '", x, "' association structure ",
             "cannot be greater than the number of longitudinal submodels.", call. = FALSE)     
      return(index)
    } else
      stop("Incorrect specification of the '", x, "' association structure. It should ",
           "include a suffix with parentheses specifying the indices of the association ",
           "terms you want to include in the interaction. See Examples in the help file.", call. = FALSE)
  } else numeric(0)      
}

# Make sure that interactions between association terms (for example
# etavalue_etaslope or mu_value_muvalue etc) are always ordered so that
# the first listed association term is for the submodel with the smallest
# index. For example, etavalue1_etavalue2 NOT etavalue2_etavalue1. This
# is to ensure there is no replication such as including both 
# etavalue1_etavalue2 AND etavalue2_etavalue1 when passing to Stan.
#
# @param assoc A named list of named lists, returned by a call to validate_assoc
# @param ok_assoc_interactions A character vector, specifying which association
#   structures are allowed to be used in interactions
check_order_of_assoc_interactions <- function(assoc, ok_assoc_interactions) {
  M <- ncol(assoc)
  for (i in ok_assoc_interactions) {
    for (j in ok_assoc_interactions) {
      header <- paste0(i, "_", j)
      header_reversed <- paste0(j, "_", i)
      for (m in 1:M) {
        if (assoc[header,][[m]]) {
          indices <- assoc["which_interactions",][[m]][[header]]
          sel <- which(indices < m)
          if (length(sel)) {
            # Remove indices for submodels before the current submodel m
            new_indices <- indices[-sel]
            assoc["which_interactions", ][[m]][[header]] <- new_indices
            assoc[header,][[m]] <- (length(new_indices) > 0L)
            # Replace those indices by reversing the order of association terms
            for (k in indices[sel]) {
              assoc["which_interactions",][[k]][[header_reversed]] <- 
                unique(c(assoc["which_interactions",][[k]][[header_reversed]], m))
              assoc[header_reversed,][[k]] <- 
                (length(assoc["which_interactions",][[k]][[header_reversed]]) > 0L)
            }
          }
        }
      }       
    }
  }
  assoc
}

# Return design matrices for evaluating longitudinal submodel quantities 
# at specified quadrature points/times
#
# @param m Integer specifying the longitudinal submodel that the association
#   structure is related to
# @param mc The matched call for the longitudinal submodel
# @param y_mod_stuff A named list returned by a call to handle_glmod (the
#   fit for a single longitudinal submodel)
# @param assoc A named list returned by a call to validate_assoc (details
#   on the desired association structure for all longitudinal submodels)
# @param id_var The name on the ID variable
# @param time_var The name of the time variable
# @param eps The time shift used for the numerical derivative calculation 
#   based on a one-sided different
# @param dataAssoc An optional data frame containing data for interactions within
#   the association terms
handle_assocmod <- function(m, data, assoc, ids, times, id_var, time_var,  
                            clust_stuff, epsilon, auc_qnodes) {
  if (!requireNamespace("data.table"))
    stop2("the 'data.table' package must be installed to use this function.")
  
  # Declare data as a data.table for merging with quadrature points
  if (clust_stuff$has_clust) {
    key_vars <- c(id_var, clust_stuff$clust_var, time_var)
  } else {
    key_vars <- c(id_var, time_var)
  }
  df <- data.table::data.table(data, key = key_vars)
  df[[time_var]] <- as.numeric(df[[time_var]]) # ensures no rounding on merge
  
  # Design matrices for calculating association structure based on 
  # (possibly lagged) eta, slope, auc and any interactions with data
  parts <- make_assoc_parts(newdata = df, assoc = assoc, ids = ids, 
                            times = times, id_var = id_var, 
                            time_var = time_var, clust_stuff = clust_stuff, 
                            epsilon = epsilon, auc_qnodes = auc_qnodes)
  
  # If association structure is based on shared random effects or shared 
  # coefficients then construct a matrix with the estimated b parameters
  # from the separate glmod (for the id_var grouping factor only). Note this
  # matrix is not passed to standata, but just used for autoscaling the 
  # priors for association parameters.
  sel_shared <- grep("^shared", rownames(assoc))
  if (any(unlist(assoc[sel_shared,]))) {
    # flist for long submodel
    flist_tmp <- lme4::getME(y_mod_stuff$mod, "flist")
    # which grouping factor is id_var
    Gp_sel <- which(names(flist_tmp) == id_var) 
    # grouping factor indices
    Gp <- lme4::getME(y_mod_stuff$mod, "Gp")  
    b_beg <- Gp[[Gp_sel]] + 1
    b_end <- Gp[[Gp_sel + 1]]
    # b vector for grouping factor = id_var
    b_vec <- lme4::getME(y_mod_stuff$mod, "b")[b_beg:b_end]
    # convert to Npat * n_re matrix
    b_mat <- matrix(b_vec, nrow = length(levels(flist_tmp[[Gp_sel]])), byrow = TRUE)
  } else b_mat <- NULL
  
  parts$b_mat <- b_mat
  return(parts)
}

# Function to construct quantities, primarily design matrices (X, Zt), that
# will be used to evaluate the longitudinal submodel contributions to the 
# association structure in the event submodel. For example, the design matrices
# evaluated at the quadpoints, the quadpoints, lagged quadpoints, auc quadpoints,
# and so on. Exactly what quantities are returned depends on what is specified
# in the use_function argument.
#
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
# @param use_function The function to call which will return the design 
#   matrices for eta, eps, lag, auc, etc.
# @param ... Additional arguments passes to use_function
# @return A named list
make_assoc_parts <- function(newdata, assoc, ids, times, id_var, time_var, 
                             clust_stuff, epsilon = 1E-5, auc_qnodes = 15L) {
  if (!requireNamespace("data.table"))
    stop("the 'data.table' package must be installed to use this function")
  
  dots <- list(...)
  m <- dots$m 
  if (is.null(m)) stop("Argument m must be specified in dots.")
  
  # Apply lag
  lag <- assoc["which_lag",][[m]]
  if (!lag == 0)
    times <- set_lag(times, lag)

  # Broadcast ids and times if there is lower level clustering
  if (clust_stuff$has_clust) {
    clust_var <- clust_stuff$clust_var
    if (is(times, "list")) {
      ids <- rep(ids, clust_stuff$clust_freq)
      times <- lapply(times, rep, clust_stuff$clust_freq)
      clust <- clust_stuff$clust_list
    } else {
      ids <- rep(rep(ids, clust_stuff$clust_freq), clust_stuff$qnodes + 1)
      times <- rep(times, rep(clust_stuff$clust_freq, clust_stuff$qnodes + 1))
      clust <- rep(clust_stuff$clust_list, clust_stuff$qnodes + 1)
    }
  } else clust <- NULL
  
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
  dataQ <- rolling_merge(data = newdata, ids = ids, times = times, clust = clust)
  mod_eta <- make_X_and_Z_fast(terms, data = dataQ, X_form, Z_forms, X_bar)
  
  # If association structure is based on slope, then calculate design 
  # matrices under a time shift of epsilon
  sel_slope <- grep("etaslope", rownames(assoc))
  if (any(unlist(assoc[sel_slope,]))) {
    dataQ_pos <- dataQ_neg <- dataQ
    dataQ_neg[[time_var]] <- dataQ_neg[[time_var]] - epsilon
    dataQ_pos[[time_var]] <- dataQ_pos[[time_var]] + epsilon
    mod_neg <- make_X_and_Z_fast(terms, data = dataQ_neg, X_form, Z_forms)
    mod_pos <- make_X_and_Z_fast(terms, data = dataQ_pos, X_form, Z_forms)
    mod_eps <- mod_pos
    mod_eps$X <- (mod_pos$X - mod_neg$X) / epsilon # derivative of X
    mod_eps$Z <- xapply(mod_pos$Z, mod_neg$Z,      # derivative of Z
                        FUN = function(x, y) (x - y) / epsilon)
  } else mod_eps <- NULL 
  
  # If association structure is based on area under the marker trajectory, then 
  # calculate design matrices at the subquadrature points
  sel_auc <- grep("etaauc|muauc", rownames(assoc))
  if (any(unlist(assoc[sel_auc,]))) {
    if (clust_stuff$has_clust)
      stop2("'etaauc' and 'muauc' not yet implemented when there is a grouping ",
            "factor clustered within patients.")
    # Return a design matrix that is (qnodes * auc_qnodes * Npat) rows 
    auc_qtimes <- 
      lapply(times, function(x) unlist(
        lapply(x, function(y) 
          lapply(get_quadpoints(auc_qnodes)$points, unstandardise_qpts, 0, y))))
    auc_qwts <- 
      lapply(times, function(x) unlist(
        lapply(x, function(y) 
          lapply(get_quadpoints(auc_qnodes)$weights, unstandardise_qwts, 0, y))))
    ids2 <- rep(ids, each = auc_qnodes)
    dataQ_auc <- rolling_merge(data = newdata, ids = ids2, times = auc_qtimes)
    mod_auc <- use_function(newdata = dataQ_auc, ...)
  } else mod_auc <- auc_qtimes <- auc_qwts <- NULL
  
  # If association structure is based on interactions with data, then calculate 
  # the design matrix which will be multiplied by etavalue, etaslope, muvalue or muslope
  sel_data <- grep("_data", rownames(assoc), value = TRUE)
  X_data <- xapply(sel_data, FUN = function(i) { 
    form <- assoc["which_formulas",][[m]][[i]]
    if (length(form)) {
      form <- as.formula(form)
      vars <- rownames(attr(terms.formula(form), "factors"))
      if (is.null(vars))
        stop2("No variables found in the formula for the '", i, "' association structure.")
      sel <- which(!vars %in% colnames(dataQ))
      if (length(sel))
        stop2("The following variables were specified in the formula for the '", i,
              "' association structure, but they cannot be found in the data: ", 
              paste0(sel, collapse = ", "))
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
  
  ret <- nlist(times, mod_eta, mod_eps, mod_auc, X_data, K_data, X_bind_data, clust_stuff)
  
  structure(ret, times = times, lag = lag, epsilon = epsilon, auc_qnodes = auc_qnodes,
            auc_qtimes = auc_qtimes, auc_qwts = auc_qwts)
}                              

# Carry out a rolling merge
#
# @param data A data.table with a set key corresponding to ids and times
# @param ids A vector of patient ids to merge against
# @param times A vector of (new) times to merge against
# @param clust A vector of cluster ids to merge against when there is clustering
#   within patient ids
# @return A data.table formed by a merge of ids, (clust), times, and the closest 
#   preceding (in terms of times) rows in data
rolling_merge <- function(data, ids, times, clust = NULL) {
  if (!requireNamespace("data.table"))
    stop("the 'data.table' package must be installed to use this function")
  key_length <- if (is.null(clust)) 2L else 3L
  if (!length(data.table::key(data)) == key_length)
    stop("Bug found: data.table key is not the same length as supplied keylist.")
  
  if (is(times, "list") && is.null(clust)) {
    return(do.call(rbind, lapply(times, FUN = function(x, ids) {
      tmp <- data.table::data.table(ids, x)
      data[tmp, roll = TRUE, rollends = c(TRUE, TRUE)]
    }, ids = ids)))     
  } else if (is(times, "list")) {
    return(do.call(rbind, lapply(times, FUN = function(x, ids, clust) {
      tmp <- data.table::data.table(ids, clust, x)
      data[tmp, roll = TRUE, rollends = c(TRUE, TRUE)]
    }, ids = ids, clust = clust)))
  } else if (is.null(clust)) {
    tmp <- data.table::data.table(ids, times)
    return(data[tmp, roll = TRUE, rollends = c(TRUE, TRUE)])       
  } else {
    tmp <- data.table::data.table(ids, clust, times)
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
# @param terms A terms object corresponding to the model frame for the
#   longitudinal submodel evaluated at the longitudinal observation times.
# @param data A data frame; the data for the longitudinal submodel 
#   at the quadrature points.
# @param X_form The formula for the fe design matrix.
# @param Z_forms A list of formulas for the re design matrices for each
#   grouping factor
# @param X_bar An optional vector of column means used for centering X.
#   Only relevant if the design matrices are going to be used in the stan code.
make_X_and_Z_fast <- function(terms, data, X_form, Z_forms, X_bar = NULL) {
  model_frame <- stats::model.frame(terms, data)
  X <- model.matrix(X_form, model_frame)
  if (!is.null(X_bar)) { # assoc parts are for passing to jm.stan
    X <- drop_intercept(X)
    X <- sweep(X, 2, X_bar, FUN = "-")
  }
  group_vars <- names(Z_forms)
  group_list <- lapply(group_vars, function(x) model_frame[[x]])
  Z <- lapply(Z_forms, model.matrix, model_frame)
  names(Z) <- names(group_list) <- group_vars
  nlist(X, Z, group_list, group_vars)
}

# Get the information need for combining the information in lower-level units
# clustered within an individual, when the patient-level is not the only 
# clustering level in the longitudinal submodel
#
# @param cnms The component names for a single longitudinal submodel
# @param flist The flist for a single longitudinal submodel
# @param id_var The name of the ID variable
# @param qnodes Integer specifying the number of qnodes being used for 
#   the GK quadrature in the stan_jm call
# @param grp_assoc Character string specifying the association structure used
#   for combining information in the lower level units clustered within an
#   individual
# @return A named list with the following elements:
#   has_clust: logical specifying whether the submodel has a grouping factor
#     that is clustered with patients
#   clust_var: the name of any grouping factor that is clustered with patients
#   clust_list: a vector of unique group ids for the clust_var grouping 
#     factor; this is in the same order as the original data.
#   clust_ids: a vector of patient ids that provides the patient that each
#     unique level of clust_list is clustered within.
#   clust_ids_q: the clust_ids vector repeated qnodes+1 times. This is used
#     in the stan code, to identify which elements of the longitudinal 
#     linear predictor need to be collapsed across for each patient.
#   qnodes: integer specifying the number of GK quadrature nodes.
get_basic_clust <- function(cnms, flist, id_var) {
  cnms_nms <- names(cnms)
  tally <- xapply(cnms_nms, FUN = function(x) 
    # within each ID, count the number of levels for the grouping factor x
    tapply(flist[[x]], flist[[id_var]], FUN = n_distinct))
  sel <- which(sapply(tally, function(x) !all(x == 1L)) == TRUE)
  has_clust <- as.logical(length(sel))
  if (!has_clust) {
    return(nlist(has_clust))
  } else {
    if (length(sel) > 1L)
      stop("There can only be one grouping factor clustered within 'id_var'.")
    clust_var <- cnms_nms[sel] 
    return(nlist(has_clust, clust_var))
  }
}

get_extra_clust <- function(basic_info, flist, id_var, qnodes, grp_assoc, 
                            ok_grp_assocs = c("sum", "mean", "min", "max")) {
  has_clust <- basic_info$has_clust
  clust_var <- basic_info$clust_var
  if (!has_clust) { # no grouping factor clustered within patients
    return(basic_info)
  } else { # submodel has a grouping factor clustered within patients
    if (is.null(clust_var))
      stop2("Bug found: could not find 'clust_var' in basic_info.")
    if (is.null(grp_assoc))
      stop2("'grp_assoc' cannot be NULL when there is a grouping factor ",
            "clustered within patients.")       
    if (!grp_assoc %in% ok_grp_assocs)
      stop2("'grp_assoc' must be one of: ", paste(ok_grp_assocs, collapse = ", "))
    
    # cluster and patient ids for each row of the Z matrix
    factor_clust <- factor(flist[[clust_var]]) 
    factor_ids <- factor(flist[[id_var]])
    
    # num clusters within each patient
    clust_freq <- tapply(factor_clust, factor_ids, FUN = n_distinct)
    
    # a vector of unique cluster ids
    clust_list <- unique(factor_clust)
    
    # a vector of patient ids indexed alongside clust_list
    clust_ids <- rep(unique(factor_ids), clust_freq) 
    
    # broadcast clust_ids for each quadnode
    clust_ids_q <- rep(clust_ids, qnodes + 1)
    
    basic_info <- nlist(has_clust, clust_var)
    extra_info <- nlist(clust_freq, clust_list, clust_ids, clust_ids_q, qnodes)
    return(c(basic_info, extra_info))
  }
}

# Function to calculate the number of association parameters in the model
#
# @param assoc A list of length M with information about the association structure
#   type for each submodel, returned by an mapply call to validate_assoc
# @param a_mod_stuff A list of length M with the design matrices related to
#   the longitudinal submodels in the GK quadrature, returned by an mapply 
#   call to handle_assocmod
# @return Integer indicating the number of association parameters in the model 
get_num_assoc_pars <- function(assoc, a_mod_stuff) {
  sel1 <- c("etavalue", "etaslope", "etaauc", 
            "muvalue", "muslope", "muauc")
  sel2 <- c("which_b_zindex", "which_coef_zindex")
  sel3 <- c("which_interactions")
  K1 <- sum(as.integer(assoc[sel1,]))
  K2 <- length(unlist(assoc[sel2,]))
  K3 <- length(unlist(assoc[sel3,]))
  K4 <- sum(fetch_(a_mod_stuff, "K_data"))
  K1 + K2 + K3 + K4
}



# Function to construct a design matrix for the association structure in
# the event submodel, to be multiplied by a vector of association parameters
#
# @param assoc An array with information about the desired association 
#   structure, returned by a call to validate_assoc
# @param parts A list equal in length to the number of markers. Each element
#   parts[[m]] should contain a named list with components $mod_eta, $mod_eps,
#   $mod_auc, which each contain either the linear predictor at quadtimes, 
#   quadtimes + eps, and auc quadtimes, or the design matrices
#   used for constructing the linear predictor. Each element parts[[m]] should 
#   also contain $xmat_data and $K_data.
# @param family A list of family objects, equal in length to the number of 
#   longitudinal submodels
# @param ... If parts does not contain the linear predictors, then this should
#   include elements beta and b, each being a length M list of parameters for the
#   longitudinal submodels
# @return A design matrix containing the association terms to be multiplied by
#   the association paramters.
make_assoc_terms <- function(parts, assoc, family, ...) {
  M <- length(parts)
  a_X <- list()
  mark <- 1
  for (m in 1:M) {
    times  <- attr(parts[[m]], "times")
    eps    <- attr(parts[[m]], "eps")  
    qnodes <- attr(parts[[m]], "auc_qnodes")
    qwts   <- unlist(attr(parts[[m]], "auc_qwts"))
    
    if (!assoc["null",][[m]]) {
      invlink_m <- family[[m]]$linkinv    
      eta_m <- get_element(parts, m = m, "eta", ...)
      eps_m <- get_element(parts, m = m, "eps", ...)
      auc_m <- get_element(parts, m = m, "auc", ...)
      data_m <- get_element(parts, m = m, "xmat_data", ...)
      K_data_m <- get_element(parts, m = m, "K_data", ...)
      has_clust_m <- parts[[m]]$clust_stuff$has_clust
      if (has_clust_m)
        clust_mat_m <- parts[[m]]$clust_stuff$clust_mat
      
      # etavalue and any interactions      
      if (assoc["etavalue",][[m]]) { # etavalue
        if (has_clust_m) {
          a_X[[mark]] <- clust_mat_m %*% eta_m
        } else {
          a_X[[mark]] <- eta_m
        }
        mark <- mark + 1
      }
      if (assoc["etavalue_data",][[m]]) { # etavalue*data
        idx_ev <- which(names(K_data_m) == "etavalue_data")
        cbeg  <- sum(K_data_m[0:(idx_ev-1)]) + 1
        cend  <- sum(K_data_m[0: idx_ev   ])
        val <- as.vector(eta_m) * data_m[, cbeg:cend, drop = FALSE]
        if (has_clust_m) {
          a_X[[mark]] <- do.call(cbind, lapply(
            1:ncol(val), function(i) clust_mat_m %*% as.vector(val[,i])))
        } else {
          a_X[[mark]] <- val
        }
        mark <- mark + 1
      }
      if (assoc["etavalue_etavalue",][[m]]) { # etavalue*etavalue
        sel <- assoc["which_interactions",][[m]][["etavalue_etavalue"]]
        for (j in sel) {
          eta_j <- get_element(parts, m = j, "eta", ...)
          val   <- eta_m * eta_j 
          a_X[[mark]] <- val
          mark <- mark + 1
        }
      } 
      if (assoc["etavalue_muvalue",][[m]]) { # etavalue*muvalue
        sel <- assoc["which_interactions",][[m]][["etavalue_muvalue"]]
        for (j in sel) {
          eta_j <- get_element(parts, m = j, "eta", ...)
          invlink_j <- family[[j]]$linkinv
          val <- eta_m * invlink_j(eta_j) 
          a_X[[mark]] <- val
          mark <- mark + 1             
        }
      }       
      # etaslope and any interactions
      if (assoc["etaslope",][[m]]) { # etaslope
        dydt_m <- (eps_m - eta_m) / eps
        if (has_clust_m) {
          a_X[[mark]] <- clust_mat_m %*% dydt_m
        } else {
          a_X[[mark]] <- dydt_m
        } 
        mark <- mark + 1             
      }
      if (assoc["etaslope_data",][[m]]) { # etaslope*data
        dydt_m <- (eps_m - eta_m) / eps
        idx_es <- which(names(K_data_m) == "etaslope_data")
        cbeg  <- sum(K_data_m[0:(idx_es-1)]) + 1
        cend  <- sum(K_data_m[0: idx_es   ])
        val <- as.vector(dydt_m) * data_m[, cbeg:cend, drop = FALSE]
        if (has_clust_m) {
          a_X[[mark]] <- do.call(cbind, lapply(
            1:ncol(val), function(i) clust_mat_m %*% as.vector(val[,i])))
        } else {
          a_X[[mark]] <- val
        }
        mark <- mark + 1            
      }
      # etaauc
      if (assoc["etaauc",][[m]]) { # etaauc
        val   <- c()
        for (j in 1:length(eta_m)) {
          wgt_j <- qwts[((j-1) * qnodes + 1):(j * qnodes)]
          auc_j <- auc_m[((j-1) * qnodes + 1):(j * qnodes)]
          val[j] <- sum(wgt_j * auc_j)
        }
        a_X[[mark]] <- val
        mark <- mark + 1            
      }        
      # muvalue and any interactions
      if (assoc["muvalue",][[m]]) { # muvalue
        val <- invlink_m(eta_m) 
        a_X[[mark]] <- val
        mark <- mark + 1            
      }
      if (assoc["muvalue_data",][[m]]) { # muvalue*data
        idx_mv <- which(names(K_data_m) == "muvalue_data")
        cbeg  <- sum(K_data_m[0:(idx_mv-1)]) + 1
        cend  <- sum(K_data_m[0: idx_mv   ])
        val   <- as.vector(invlink_m(eta_m)) * data_m[, cbeg:cend, drop = FALSE] 
        a_X[[mark]] <- val
        mark <- mark + 1           
      }
      if (assoc["muvalue_etavalue",][[m]]) { # muvalue*etavalue
        sel <- assoc["which_interactions",][[m]][["muvalue_etavalue"]]
        for (j in sel) {
          eta_j <- get_element(parts, m = j, "eta", ...)
          val   <- invlink_m(eta_m) * eta_j 
          a_X[[mark]] <- val
          mark <- mark + 1           
        }
      } 
      if (assoc["muvalue_muvalue",][[m]]) { # muvalue*muvalue
        sel <- assoc["which_interactions",][[m]][["muvalue_muvalue"]]
        for (j in sel) {
          eta_j <- get_element(parts, m = j, "eta", ...)
          invlink_j <- family[[j]]$linkinv
          val <- invlink_m(eta_m) * invlink_j(eta_j) 
          a_X[[mark]] <- val
          mark <- mark + 1                   
        }
      }       
      # muslope and any interactions
      if (assoc["muslope",][[m]]) { # muslope
        val <- (invlink_m(eps_m) - invlink_m(eta_m)) / eps
        a_X[[mark]] <- val
        mark <- mark + 1                   
      }
      if (assoc["muslope_data",][[m]]) { # muslope*data
        dydt_m <- (invlink_m(eps_m) - invlink_m(eta_m)) / eps
        idx_ms <- which(names(K_data_m) == "muslope_data")
        cbeg  <- sum(K_data_m[0:(idx_ms-1)]) + 1
        cend  <- sum(K_data_m[0: idx_ms   ])
        val   <- as.vector(dydt_m) * data_m[, cbeg:cend, drop = FALSE] 
        a_X[[mark]] <- val
        mark <- mark + 1              
      }    
      # muauc
      if (assoc["muauc",][[m]]) { # muauc
        val   <- c()
        for (j in 1:length(eta_m)) {
          wgt_j <- qwts[((j-1) * qnodes + 1):(j * qnodes)]
          auc_j <- invlink_m(auc_m[((j-1) * qnodes + 1):(j * qnodes)])
          val[j] <- sum(wgt_j * auc_j)
        }
        a_X[[mark]] <- val
        mark <- mark + 1 
      }
    }
  }
  for (m in 1:M) {
    # shared_b
    if (assoc["shared_b",][[m]]) {
      sel <- assoc["which_b_zindex",][[m]]
      val <- get_element(parts, m = m, "b_mat", ...)[,sel]
      a_X[[mark]] <- val
      mark <- mark + 1                   
    }
  }    
  for (m in 1:M) {
    # shared_coef
    if (assoc["shared_coef",][[m]]) {
      sel <- assoc["which_coef_zindex",][[m]]
      val <- get_element(parts, m = m, "b_mat", ...)[,sel]
      a_X[[mark]] <- val
      mark <- mark + 1                   
    }
  }
  return(do.call("cbind", a_X))
}

# Function to get linear predictor and other bits
#
# @param parts A named list containing the parts for constructing the association 
#   structure. It may contain elements $mod_eta, $mod_eps, $mod_auc, etc. as 
#   well as $xmat_data, $K_data
# @param which
get_element <- function(parts, m = 1, which = "eta", ...) {
  dots <- list(...)
  ok_which_args <- c("eta", "eps", "auc", "xmat_data", "K_data", "b_mat")
  if (!which %in% ok_which_args)
    stop("'which' must be one of: ", paste(ok_which_args, collapse = ", "))
  if (which %in% c("eta", "eps", "auc")) {
    part <- parts[[m]][[paste0("mod_", which)]]
    if (is.null(part)) { # model doesn't include an assoc related to 'part'
      return(NULL)
    } else if (!is.null(part$linpred)) { # linpred already provided in object
      return(part$linpred)
    } else { # need to construct linpred
      x <- part$x
      Zt <- part$Zt
      Znames  <- part$Z_names
      if (is.null(x) || is.null(Zt))
        stop(paste0("Bug found: cannot find x and Zt in object. They are ",
                    "required to build the linear predictor for '", which, "'."))          
      beta <- dots$beta[[m]]
      b <- dots$b[[m]]
      if (is.null(beta) || is.null(b))
        stop("Bug found: beta and b must be provided to construct linpred.")
      return(linear_predictor.default(beta, x) + as.vector(b %*% Zt))
    }
  } else if (which %in% c("xmat_data", "K_data", "b_mat")) {
    return(parts[[m]][[which]])
  } else {
    stop("'which' argument doesn't include a valid entry.")
  }
}

#--------------- Functions related to generating initial values

# Create a function that can be used to generate the model-based initial values for Stan
#
# @param y_mod_stuff A list, with each element containing the list object returned by
#   a call to the handle_glmod function
# @param e_mod_stuff A list object returned by a call to the handle_coxmod function
# @param standata The data list that will be passed to Stan
generate_init_function <- function(y_mod_stuff, e_mod_stuff, standata) {
  
  # Initial values for intercepts, coefficients and aux parameters
  est <- lapply(y_mod_stuff, function(x) summary(x$mod)$coefficients[, "Estimate"])
  xbar      <- fetch(y_mod_stuff, "xbar")
  gamma     <- lapply(seq_along(y_mod_stuff), function(m) 
    return_intercept(est[[m]]) - xbar[[m]] %*% drop_intercept(est[[m]]))
  gamma_nob <- gamma[as.logical(standata$has_intercept_nob)]
  gamma_lob <- gamma[as.logical(standata$has_intercept_lob)]
  gamma_upb <- gamma[as.logical(standata$has_intercept_upb)]
  beta      <- lapply(est, drop_intercept)
  aux       <- lapply(y_mod_stuff, function(x) sigma(x$mod))
  e_beta    <- e_mod_stuff$mod$coef
  e_aux     <- if (standata$basehaz_type == 1L) runif(1, 0.5, 3) else rep(0, standata$basehaz_df)
  z_beta        <- standardise_coef(unlist(beta), standata$prior_mean,           standata$prior_scale)
  aux_unscaled  <- standardise_coef(unlist(aux),  standata$prior_mean_for_aux,   standata$prior_scale_for_aux)
  aux_unscaled  <- aux_unscaled[as.logical(standata$has_aux)] # only keep aux where relevant
  e_z_beta      <- standardise_coef(e_beta,       standata$e_prior_mean,         standata$e_prior_scale) 
  e_aux_unscaled<- standardise_coef(e_aux,        standata$e_prior_mean_for_aux, standata$e_prior_scale_for_aux)
  b_Cov         <- lapply(y_mod_stuff, function(x) lme4::VarCorr(x$mod)[[1L]])
  sel           <- sapply(y_mod_stuff, function(x) length(lme4::VarCorr(x$mod)) > 1L)
  #if (any(sel)) stop("Model-based initial values cannot yet be used with more ",
  #                   "than one clustering level.", call. = FALSE)
  
  # Initial values for random effects distribution
  scale <- standata$scale
  t     <- standata$t
  p     <- standata$p
  
  # Cholesky decomp of b_Cov combined across all submodel
  L_b_Cov         <- t(suppressWarnings(chol(as.matrix(Matrix::bdiag(b_Cov)), pivot = TRUE))) 
  diag_L_b_Cov    <- diag(L_b_Cov)
  
  # Dimensions
  len_z_T <- 0
  for (i in 1:t) {
    if (p[i] > 2) 
      for (j in 3:p[i]) len_z_T <- len_z_T + p[i] - 1;
  }
  len_rho <- ifelse((sum(p) - t) > 0, sum(p) - t, 0)
  
  # Construct initial values for theta_L matrix
  # ** Much room for improvement here! **
  tau              <- c()
  rho              <- c()
  normalised_zetas <- c()
  rho_mark         <- 1
  for (i in 1:t) {
    trace <- sum(diag_L_b_Cov)  # equal to variance of RE if only one RE
    normalised_zetas_tmp <- diag_L_b_Cov / trace  # equal to 1 if only one random effect
    tau[i] <- (sqrt(trace / p[i])) / scale[i]
    std_dev1 <- sqrt(L_b_Cov[1,1])
    if (p[i] > 1) {
      std_dev2 <- sqrt(L_b_Cov[2,2])
      T21 <- L_b_Cov[2,1] / std_dev2
      rho[rho_mark] <- (T21 + 1) / 2
      rho_mark <- rho_mark + 1
      normalised_zetas <- c(normalised_zetas, normalised_zetas_tmp)
    }
  }
  
  # Parameters related to priors
  len_global <- sum((2 * (standata$prior_dist == 3)) + (4 * (standata$prior_dist == 4)))
  len_local2 <- sum((standata$prior_dist == 3) * standata$KM) 
  len_local4 <- sum((standata$prior_dist == 4) * standata$KM)
  len_mix   <- sum((standata$prior_dist %in% c(5,6)) * standata$KM)
  len_ool <- sum(standata$prior_dist == 6)
  len_noise <- sum((standata$family == 8) * standata$NM)
  
  # Function to generate model based initial values
  model_based_inits <- Filter(function(x) (!is.null(x)), list(
    gamma_nob      = array_else_double(gamma_nob),
    gamma_lob      = array_else_double(gamma_lob),
    gamma_upb      = array_else_double(gamma_upb),
    z_beta         = array_else_double(z_beta),
    aux_unscaled   = array_else_double(aux_unscaled),
    e_z_beta       = array_else_double(e_z_beta),
    e_aux_unscaled = array_else_double(e_aux_unscaled),
    e_gamma  = array_else_double(rep(0, standata$e_has_intercept)),
    a_z_beta = array_else_double(rep(0, standata$a_K)),
    z_b      = array_else_double(runif(standata$q, -0.5, 0.5)),
    global   = array_else_double(runif(len_global)),
    e_global = array_else_double(runif(get_nvars_for_hs(standata$e_prior_dist))),
    a_global = array_else_double(runif(get_nvars_for_hs(standata$a_prior_dist))),
    local2   = matrix_of_uniforms(nrow = 2, ncol = len_local2),
    local4   = matrix_of_uniforms(nrow = 4, ncol = len_local4),
    e_local  = matrix_of_uniforms(nrow = get_nvars_for_hs(standata$e_prior_dist), ncol = standata$e_K),
    a_local  = matrix_of_uniforms(nrow = get_nvars_for_hs(standata$a_prior_dist), ncol = standata$a_K),
    mix   = if (len_mix > 0) matrix(rep(1, len_mix), 1, len_mix) else matrix(0,0,0),
    e_mix = if (standata$e_prior_dist %in% c(5,6)) matrix(rep(1, standata$e_K), 1, standata$e_K) else matrix(0,0,standata$e_K),
    a_mix = if (standata$a_prior_dist %in% c(5,6)) matrix(rep(1, standata$a_K), 1, standata$a_K) else matrix(0,0,standata$a_K),
    ool   = if (len_ool > 0) as.array(len_ool) else as.array(double(0)), 
    e_ool = if (standata$e_prior_dist == 6) as.array(1) else as.array(double(0)), 
    a_ool = if (standata$a_prior_dist == 6) as.array(1) else as.array(double(0)),
    noise = if (len_noise > 0) matrix(runif(len_noise), 1, len_noise) else matrix(0,0,0),
    z_T   = array_else_double(rep(sqrt(1 / len_z_T), len_z_T)),
    rho   = array_else_double(rep(1 / (len_rho + 1), len_rho)),
    zeta  = array_else_double(normalised_zetas),
    tau   = array_else_double(tau)))
  
  return(function() model_based_inits)
}


#--------------- Functions related to standata and sampling

# Set arguments for sampling for stan_jm
#
# Prepare a list of arguments to use with \code{rstan::sampling} via
# \code{do.call}.
#
# *Note that this differs from the set_sampling_args function in that
# it uses different default numbers of iterations and chains, and 
# different default adapt_delta and max_treedepth. Using the shorter 
# treedepth is to stop the sampler trailing off during early iterations. 
# This can drastically reduce the model estimation time, and in most
# examples using a shorter treedepth hasn't compromised the sampler
# at later interations (ie, at later iterations the sampler doesn't
# hit the maximum treedepth).
#
# @param object The stanfit object to use for sampling.
# @param cnms The component names for the group level parameters combined
#   across all glmer submodels. This is used to determine the maximum number
#   of parameters for any one grouping factor in the model, which in turn is
#   used to determine the default adapt_delta.
# @param user_dots The contents of \code{...} from the user's call to
#   the \code{stan_jm} modeling function.
# @param user_adapt_delta The value for \code{adapt_delta} specified by the
#   user.
# @param user_max_treedepth The value for \code{max_treedepth} specified by the
#   user.
# @param ... Other arguments to \code{\link[rstan]{sampling}} not coming from
#   \code{user_dots} (e.g. \code{pars}, \code{init}, etc.)
# @return A list of arguments to use for the \code{args} argument for 
#   \code{do.call(sampling, args)}.
set_jm_sampling_args <- function(object, cnms, user_dots = list(), 
                                 user_adapt_delta = NULL, 
                                 user_max_treedepth = NULL, 
                                 ...) {
  args <- list(object = object, ...)
  unms <- names(user_dots)
  for (j in seq_along(user_dots)) {
    args[[unms[j]]] <- user_dots[[j]]
  }
  
  max_p <- max(sapply(cnms, length))
  
  default_adapt_delta <- if (max_p > 2) 0.99 else 0.95
  default_max_treedepth <- 11L
  
  if (!is.null(user_adapt_delta))
    args$control$adapt_delta <- user_adapt_delta else 
      if (is.null(args$control$adapt_delta))
        args$control$adapt_delta <- default_adapt_delta
  
  if (!is.null(user_max_treedepth))
    args$control$max_treedepth <- user_max_treedepth else
      if (is.null(args$control$max_treedepth))
        args$control$max_treedepth <- default_max_treedepth
  
  if (!"save_warmup" %in% unms) 
    args$save_warmup <- TRUE  
  
  return(args)
}  

#--------------- Functions related to observation weights

# Check the weights argument for stan_jm
#
# @param weights The data frame passed via the weights argument
# @param id_var The name of the ID variable
check_weights <- function(weights, id_var) {
  
  # Check weights are an appropriate data frame
  if ((!is.data.frame(weights)) || (!ncol(weights) == 2))
    stop("'weights' argument should be a data frame with two columns: the first ",
         "containing patient IDs, the second containing their corresponding ",
         "weights.", call. = FALSE)
  if (!id_var %in% colnames(weights))
    stop("The data frame supplied in the 'weights' argument should have a ",
         "column named ", id_var, call. = FALSE)
  weight_var <- setdiff(colnames(weights), id_var)
  
  # Check weights are positive and numeric
  wts <- weights[[weight_var]]
  if (!is.numeric(wts)) 
    stop("The weights supplied must be numeric.", call. = FALSE)
  if (any(wts < 0)) 
    stop("Negative weights are not allowed.", call. = FALSE)
  
  # Check only one weight per ID
  n_weights_per_id <- tapply(weights[[weight_var]], weights[[id_var]], length)
  if (!all(n_weights_per_id == 1L))
    stop("The data frame supplied in the 'weights' argument should only have ",
         "one row (ie, one weight) per patient ID.", call. = FALSE)
}

# Return the vector of prior weights for one of the submodels
#
# @param mod_stuff A named list with elements: y, flist, ord
# @param weights The data frame passed via the weights argument
# @param id_var The name of the ID variable
handle_weights <- function(mod_stuff, weights, id_var) {
  
  is_glmod <- (is.null(mod_stuff$eventtime))
  
  # No weights provided by user
  if (is.null(weights)) {
    len <- if (is_glmod) length(mod_stuff$y) else length(mod_stuff$eventtime)
    return(rep(0.0, len)) 
  }
  
  # Check for IDs with no weight supplied
  weights[[id_var]] <- factor(weights[[id_var]])
  ids <- if (is_glmod) mod_stuff$flist[[id_var]] else factor(mod_stuff$flist)
  sel <- which(!ids %in% weights[[id_var]])
  if (length(sel)) {
    if (length(sel) > 30L) sel <- sel[1:30]
    stop(paste0("The following patient IDs are used in fitting the model, but ",
                "do not have weights supplied via the 'weights' argument: ",
                paste(ids[sel], collapse = ", ")), call. = FALSE)
  }
  
  # Obtain length and ordering of weights vector using flist
  wts_df  <- merge(data.frame(id = ids), weights, by.x = "id", by.y = id_var, sort = FALSE)
  wts_var <- setdiff(colnames(weights), id_var)
  wts     <- wts_df[[wts_var]]
  
  wts
}



