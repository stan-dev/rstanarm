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

#' @importFrom survival Surv
#' @export
survival::Surv

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
# @param min_prior_scale The minimum allowed for prior scales
# @param assoc A two dimensional array with information about desired association
#   structure for the joint model (returned by a call to validate_assoc). Cannot
#   be NULL if autoscaling priors for the association parameters.
# @param ... Other arguments passed to make_assoc_terms. If autoscaling priors 
#   for the association parameters then this should include 'parts' which 
#   is a list containing the design matrices for the longitudinal submodel 
#   evaluated at the quadrature points, as well as 'beta' and 'b' which are
#   the parameter values to use when constructing the linear predictor(s) in
#   make_assoc_terms.
# @return A named list with the same structure as returned by handle_glm_prior
autoscale_prior <- function(prior_stuff, response = NULL, predictors = NULL, 
                            family = NULL, QR = FALSE, min_prior_scale = 1e-12, 
                            assoc = NULL, ...) {
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
    assoc_terms <- make_assoc_terms(family = family, assoc = assoc, ...)
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
# @param e_has_intercept T/F, does event submodel have an intercept?
# @param e_has_predictors T/F, does event submodel have predictors?
# @param has_assoc Logical specifying whether the model has an association 
#   structure. Can be NULL if the prior summary is not for a joint model.
# @param adjusted_prior_*_scale Adjusted scales computed if using autoscaled priors
# @param family A list of family objects.
# @param basehaz A list with information about the baseline hazard.
# @param stub_for_names Character string with the text stub to use in the 
#   names identifying the glmer or longitudinal submodels. 
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
           b_user_prior_stuff = NULL,
           b_prior_stuff = NULL,
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
      if (!is.na(e_aux_name)) {
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
    
    if (length(user_prior_covariance)) {
      if (user_prior_covariance$dist == "decov") {
        prior_list$prior_covariance <- user_prior_covariance
      } else if (user_prior_covariance$dist == "lkj") {
        # lkj prior for correlation matrix
        prior_list$prior_covariance <- user_prior_covariance
        # half-student_t prior on SD for each ranef (possibly autoscaled)
        prior_list$prior_covariance$df <- b_user_prior_stuff$prior_df
        prior_list$prior_covariance$scale <- b_user_prior_stuff$prior_scale
        adj_scales <- uapply(b_prior_stuff, FUN = uapply, '[[', "prior_scale")
        if (!all(b_user_prior_stuff$prior_scale == adj_scales)) {
          prior_list$prior_covariance$adjusted_scale <- adj_scales
        } else {
          prior_list$prior_covariance$adjusted_scale <- NULL
        }
      } else {
        prior_list$prior_covariance <- NULL
      }
    }
    
    if (!stub_for_names == "Long") {
      nms <- names(prior_list)
      new_nms <- gsub("Long", "", nms)
      names(prior_list) <- new_nms
    }
    
    return(prior_list)
  }

# Get name of auxiliary parameters for event submodel
#
# @param basehaz A list with information about the baseline hazard
.rename_e_aux <- function(basehaz) {
  nm <- basehaz$type_name
  switch(nm,
         weibull   = "weibull-shape",
         gompertz  = "gompertz-scale",
         bs        = "B-spline-coefficients",
         ms        = "M-spline-coefficients",
         piecewise = "piecewise-coefficients",
         NA)
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

# Return the response vector for passing to Stan
#
# @param formula The model formula
# @param model_frame The model frame
# @param family A family object
# @return A named list with the following elements:
#   y: the response vector
#   real: the response vector if real, else numeric(0)
#   integer: the response vector if integer, else integer(0)
#   resp_type: 1L if response is real, 2L is response is integer
make_y_for_stan <- function(formula, model_frame, family) {
  y <- as.vector(model.response(model_frame))
  y <- validate_glm_outcome_support(y, family)
  resp_type <- if (check_response_real(family)) 1L else 2L
  real    <- if (resp_type == 1L) y else numeric(0) 
  integer <- if (resp_type == 2L) y else integer(0) 
  nlist(y, real, integer, resp_type)
}

# Return the design matrix for passing to Stan
#
# @param formula The model formula.
# @param model_frame The model frame.
# @return A named list with the following elements:
#   x: the fe model matrix, not centred and may have intercept.
#   xtemp: fe model matrix, centred and no intercept.
#   x_form: the formula for the fe model matrix.
#   x_bar: the column means of the model matrix.
#   has_intercept: logical for whether the submodel has an intercept
#   N,K: number of rows (observations) and columns (predictors) in the
#     fixed effects model matrix
make_x_for_stan <- function(formula, model_frame) {
  x_form <- lme4::nobars(formula)
  x <- model.matrix(x_form, model_frame)
  has_intercept <- check_for_intercept(x, logical = TRUE)
  xtemp <- drop_intercept(x)
  x_bar <- colMeans(xtemp)
  xtemp <- sweep(xtemp, 2, x_bar, FUN = "-")
  # identify any column of x with < 2 unique values (empty interaction levels)
  sel <- (2 > apply(xtemp, 2L, function(x) length(unique(x))))
  if (any(sel))
    stop2("Cannot deal with empty interaction levels found in columns: ",
            paste(colnames(xtemp)[sel], collapse = ", "))
  nlist(x, xtemp, x_form, x_bar, has_intercept, N = NROW(xtemp), K = NCOL(xtemp))
}

# Return design matrices for the group level terms for passing to Stan
#
# @param formula The model formula
# @param model_frame The model frame
# @return A named list with the following elements:
#   z: a list with each element containing the random effects model 
#     matrix for one grouping factor.
#   z_forms: a list with each element containing the model formula for 
#     one grouping factor.
#   group_vars: a character vector with the name of each of the
#     grouping factors
#   group_cnms: a list with each element containing the names of the
#     group level parameters for one grouping factor
#   group_list: a list with each element containing the vector of group 
#     IDs for the rows of z
#   nvars: a vector with the number of group level parameters for each
#     grouping factor
#   ngrps: a vector with the number of groups for each grouping factor 
make_z_for_stan <- function(formula, model_frame) {
  bars <- lme4::findbars(formula)
  if (length(bars) > 2L)
    stop2("A maximum of 2 grouping factors are allowed.")
  z_parts <- lapply(bars, split_at_bars)
  z_forms <- fetch(z_parts, "re_form")
  z <- lapply(z_forms, model.matrix, model_frame)
  group_cnms <- lapply(z, colnames)
  group_vars <- fetch(z_parts, "group_var")
  group_list <- lapply(group_vars, function(x) factor(model_frame[[x]]))
  nvars <- lapply(group_cnms, length)
  ngrps <- lapply(group_list, n_distinct)
  names(z) <- names(z_forms) <- names(group_cnms) <- 
    names(group_list) <- names(nvars) <- names(ngrps) <- group_vars
  nlist(z, z_forms, group_vars, group_cnms, group_list, nvars, ngrps)
}

# Return info on the required type of intercept
#
# @param X The model matrix
# @param family A family object
# @return A named list with the following elements:
#   type: character string specifying the type of bounds to use
#     for the intercept.
#   number: an integer specifying the type of bounds to use
#     for the intercept where 0L = no intercept, 1L = no bounds 
#     on intercept, 2L = lower bound, 3L = upper bound.
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

# Check the id_var argument is valid and is included appropriately in the
# formulas for each of the longitudinal submodels
#
# @param id_var The character string that the user specified for the id_var
#   argument -- will have been set to NULL if the argument was missing.
# @param y_cnms A list of length M with the cnms for each longitudinal submodel
# @param y_flist A list of length M with the flist for each longitudinal submodel
# @return Returns the character string corresponding to the appropriate id_var.
#   This will either be the user specified id_var argument or the only grouping
#   factor.
check_id_var <- function(y_mod, id_var) {
  
  y_cnms <- fetch(y_mod, "z", "group_cnms")

  len_cnms <- sapply(y_cnms, length) # num grouping factors in each submodel
  
  if (any(len_cnms > 1L)) { # multiple grouping factors
    
    if (is.null(id_var))
      stop2("With more than one grouping factor 'id_var' must be specified.")
    lapply(y_cnms, function(x)  if (!(id_var %in% names(x)))
      stop2("'id_var' must be a grouping factor in each longitudinal submodel.")) 
    return(id_var)
    
  } else { # only one grouping factor (assumed to be subject ID)
    
    only_cnm <- unique(sapply(y_cnms, names))
    if (length(only_cnm) > 1L)
      stop2("The grouping factor (ie, subject ID variable) is not the ",
           "same in all longitudinal submodels.")
    if (not.null(id_var) && !identical(id_var, only_cnm))
      warning2("The user specified 'id_var' (", paste(id_var), 
              ") and the assumed ID variable based on the single ",
              "grouping factor (", paste(only_cnm), ") are not the same; ", 
              "'id_var' will be ignored")
    return(only_cnm)
  }
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
  terms <- strsplit(deparse(x, 500), "\\s\\|\\s")[[1L]]
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
get_common_cnms <- function(y_mod, stub = "Long") {
  y_cnms <- fetch(y_mod, "z", "group_cnms")
  nms <- lapply(y_cnms, names)
  unique_nms <- unique(unlist(nms))
  cnms <- lapply(seq_along(unique_nms), function(i) {
    nm <- unique_nms[i]
    unlist(lapply(1:length(y_cnms), function(m) 
      if (nm %in% nms[[m]]) paste0(stub, m, "|", y_cnms[[m]][[nm]])))
  })
  names(cnms) <- unique_nms
  if (length(cnms) > 2L)
    stop("A maximum of 2 grouping factors are allowed.")
  cnms
}

# Function to return a single list with the factor levels for each
# grouping factor, but collapsed across all longitudinal submodels
# 
# @param x A list containing the flist object for each of the submodels
get_common_flevels <- function(y_mod) {
  y_flevels <- fetch(y_mod, "z", "group_list")
  nms <- lapply(y_flevels, names)
  unique_nms <- unique(unlist(nms))
  flevels <- lapply(seq_along(unique_nms), function(i) {
    nm <- unique_nms[i]
    flevels_nm <- lapply(1:length(y_flevels), function(m) 
      if (nm %in% nms[[m]]) levels(y_flevels[[m]][[nm]]))
    flevels_nm <- rm_null(unique(flevels_nm))
    if (length(flevels_nm) > 1L)
      stop2("The group factor levels must be the same for all submodels.")
    flevels_nm[[1L]]
  })
  names(flevels) <- unique_nms
  flevels
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
  check <- sapply(y_cnms, function(x) group_var %in% names(x))
  if (!all(check)) {
    nm <- deparse(substitute(group_var))
    stop2(nm, " must be a grouping factor in all longitudinal submodels.")
  }
  group_var
}

# Check the factor list corresponding to subject ID is the same in each 
# of the longitudinal submodels
#
# @param id_var The name of the ID variable
# @param y_flist A list containing the flist objects returned for each 
#   separate longitudinal submodel
# @return A vector of factor levels corresponding to the IDs appearing
#   in the longitudinal submodels
check_id_list <- function(y_mod, id_var) {
  y_flist <- fetch(y_mod, "z", "group_list")
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
use_predvars <- function(mod, keep_response = TRUE) {
  fm <- formula(mod)
  ff <- lapply(attr(terms(mod, fixed.only  = TRUE), "variables"), deparse, 500)[-1]
  fr <- lapply(attr(terms(mod, random.only = TRUE), "variables"), deparse, 500)[-1]
  pf <- lapply(attr(terms(mod, fixed.only  = TRUE), "predvars"),  deparse, 500)[-1]
  pr <- lapply(attr(terms(mod, random.only = TRUE), "predvars"),  deparse, 500)[-1]
  if (!identical(c(ff, fr), c(pf, pr))) {
    for (j in 1:length(ff))
      fm <- gsub(ff[[j]], pf[[j]], fm, fixed = TRUE)    
    for (j in 1:length(fr))
      fm <- gsub(fr[[j]], pr[[j]], fm, fixed = TRUE)    
  }
  rhs <- fm[[length(fm)]]
  if (is(rhs, "call")) 
    rhs <- deparse(rhs, 500L)
  if (keep_response && length(fm) == 3L) {
    fm <- reformulate(rhs, response = formula(mod)[[2L]])
  } else if (keep_response && length(fm) == 2L) {
    warning2("No response variable found, reformulating RHS only.")
    fm <- reformulate(rhs, response = NULL)
  } else {
    fm <- reformulate(rhs, response = NULL)
  }
  fm
}

# Check that the observation times for the longitudinal submodel are all
# positive and not observed after the individual's event time
#
# @param data A data frame (data for one longitudinal submodel)
# @param eventtimes A named numeric vector with the event time for each
#   individual. The vector names should be the individual ids.
# @param id_var,time_var The ID and time variable in the longitudinal data.
# @return Nothing.
validate_observation_times <-function(data, exittime, id_var, time_var) {
  if (!time_var %in% colnames(data)) 
    STOP_no_var(time_var)
  if (!id_var %in% colnames(data)) 
    STOP_no_var(id_var)
  if (any(data[[time_var]] < 0))
    stop2("Values for the time variable (", time_var, ") should not be negative.")
  mt  <- tapply(data[[time_var]], factor(data[[id_var]]), max) # max observation time
  nms <- names(exittime)                                       # patient IDs
  if (is.null(nms))
    stop2("Bug found: cannot find names in the vector of exit times.")
  sel <- which(sapply(nms, FUN = function(i) mt[i] > exittime[i]))
  if (length(sel))
    stop2("The following individuals have observation times in the longitudinal ",
          "data that are later than their event time: ", comma(nms[sel]))     
}


#--------------- Functions related to event submodel

# Construct a list with information about the baseline hazard
#
# @param basehaz A string specifying the type of baseline hazard
# @param basehaz_ops A named list with elements df, knots 
# @param ok_basehaz A list of admissible baseline hazards
# @param times A numeric vector with eventtimes for each individual
# @param status A numeric vector with event indicators for each individual
# @param upper_times A numeric vector (or NULL) with the upper limit for any
#   observations with interval censoring.
# @return A named list with the following elements:
#   type: integer specifying the type of baseline hazard, 1L = weibull,
#     2L = b-splines, 3L = piecewise.
#   type_name: character string specifying the type of baseline hazard.
#   user_df: integer specifying the input to the df argument
#   df: integer specifying the number of parameters to use for the 
#     baseline hazard.
#   knots: the knot locations for the baseline hazard.
#   bs_basis: The basis terms for the B-splines. This is passed to Stan
#     as the "model matrix" for the baseline hazard. It is also used in
#     post-estimation when evaluating the baseline hazard for posterior
#     predictions since it contains information about the knot locations
#     for the baseline hazard (this is implemented via splines::predict.bs). 
handle_basehaz <- function(basehaz, 
                           basehaz_ops, 
                           ok_basehaz     = c("weibull", "bs", "piecewise"),
                           ok_basehaz_ops = c("df", "knots"),
                           times, 
                           status,
                           upper_times) {
  
  if (!basehaz %in% ok_basehaz)
    stop2("'basehaz' should be one of: ", comma(ok_basehaz))
  
  if (!all(names(basehaz_ops) %in% ok_basehaz_ops))
    stop2("'basehaz_ops' can only include: ", comma(ok_basehaz_ops))
  
  tt <- times[(status == 1)] # uncensored event times
  
  if (basehaz == "exp") {
    
    bknots <- NULL # boundary knot locations
    iknots <- NULL # internal knot locations
    basis  <- NULL # spline basis
    nvars  <- 0L   # number of aux parameters, none
    
  } else if (basehaz == "gompertz") {
    
    bknots <- NULL # boundary knot locations
    iknots <- NULL # internal knot locations
    basis  <- NULL # spline basis
    nvars  <- 1L   # number of aux parameters, Gompertz scale
    
  } else if (basehaz == "weibull") {
    
    bknots <- NULL # boundary knot locations
    iknots <- NULL # internal knot locations
    basis  <- NULL # spline basis
    nvars  <- 1L   # number of aux parameters, Weibull shape
    
  } else if (basehaz == "bs") {
    
    df    <- basehaz_ops$df
    knots <- basehaz_ops$knots
    
    if (!is.null(df) && !is.null(knots))
      stop2("Cannot specify both 'df' and 'knots' for the baseline hazard.")
    
    if (is.null(df))
      df <- 5L # default df for B-splines, assuming no intercept
    # NB this is ignored if the user specified knots
    
    bknots <- get_bknots(c(times, upper_times))
    iknots <- get_iknots(tt, df = df, iknots = knots)
    basis  <- get_basis(tt, iknots = iknots, bknots = bknots, type = "bs")      
    nvars  <- ncol(basis)  # number of aux parameters, basis terms
    
  } else if (basehaz == "ms") {
    
    df    <- basehaz_ops$df
    knots <- basehaz_ops$knots
    
    if (!is.null(df) && !is.null(knots)) {
      stop2("Cannot specify both 'df' and 'knots' for the baseline hazard.")
    }
    
    if (is.null(df)) {
      df <- 5L # default df for B-splines, assuming no intercept
      # NB this is ignored if the user specified knots
    }
    
    bknots <- get_bknots(c(times, upper_times))
    iknots <- get_iknots(tt, df = df, iknots = knots)
    basis  <- get_basis(tt, iknots = iknots, bknots = bknots, type = "ms")      
    nvars  <- ncol(basis)  # number of aux parameters, basis terms
    
  } else if (basehaz == "piecewise") {
    
    df    <- basehaz_ops$df
    knots <- basehaz_ops$knots
    
    if (!is.null(df) && !is.null(knots)) {
      stop2("Cannot specify both 'df' and 'knots' for the baseline hazard.")
    }
    
    if (is.null(df)) {
      df <- 6L # default number of segments for piecewise constant
      # NB this is ignored if the user specified knots
    }
    
    bknots <- get_bknots(c(times, upper_times))
    iknots <- get_iknots(tt, df = df, iknots = knots)
    basis  <- NULL               # spline basis
    nvars  <- length(iknots) + 1 # number of aux parameters, dummy indicators
    
  }  
  
  nlist(type_name = basehaz, 
        type = basehaz_for_stan(basehaz), 
        nvars, 
        iknots, 
        bknots, 
        basis,
        df = nvars,
        user_df = nvars,
        knots = if (basehaz == "bs") iknots else c(bknots[1], iknots, bknots[2]),
        bs_basis = basis)
}

# Check that the ids in the longitudinal and survival models match
validate_jm_ids <- function(y_ids, e_ids) {
  if (!identical(y_ids, levels(factor(e_ids))))
    stop2("The patient IDs (levels of the grouping factor) included ",
          "in the longitudinal and event submodels do not match")
  if (is.unsorted(factor(e_ids)))
    stop2("'dataEvent' needs to be sorted by the subject ID/grouping variable.")
  if (!identical(length(y_ids), length(e_ids)))
    stop2("The number of patients differs between the longitudinal and ",
          "event submodels. Perhaps you intended to use 'start/stop' notation ",
          "for the Surv() object.")
}

# Return the integer respresentation for the baseline hazard, used by Stan
#
# @param basehaz_name A character string, the type of baseline hazard.
# @return An integer, or NA if unmatched.
basehaz_for_stan <- function(basehaz_name) {
  switch(basehaz_name, 
         weibull   = 1L, 
         bs        = 2L,
         piecewise = 3L,
         ms        = 4L,
         exp       = 5L,
         gompertz  = 6L,
         NA)
}

# Return a vector with boundary knots for 'x'
#
# @param x A numeric vector.
# @return A numeric vector of length 2, the boundary knot locations.
get_bknots <- function(x) {
  c(0, max(x))
}

# Return a vector with internal knots for 'x', based on evenly spaced quantiles
#
# @param x A numeric vector.
# @param df The degrees of freedom. If specified, then 'df - degree - intercept'.
#   knots are placed at evenly spaced percentiles of 'x'. If 'iknots' is 
#   specified then 'df' is ignored.
# @return A numeric vector of internal knot locations, or NULL if there are
#   no internal knots.
get_iknots <- function(x, df = 6L, degree = 3L, iknots = NULL, intercept = TRUE) {
  
  # obtain number of internal knots
  if (is.null(iknots)) {
    nk <- df - degree - intercept
  } else {
    nk <- length(iknots)
  }
  
  # validate number of internal knots
  if (nk < 0) {
    stop2("Number of internal knots cannot be negative.")
  }
  
  # obtain default knot locations if necessary
  if (is.null(iknots)) {
    iknots <- qtile(x, nq = nk + 1)  # evenly spaced percentiles
  }
  
  # return internal knot locations, ensuring they are positive
  validate_positive_scalar(iknots)
  
  return(iknots)
}

# Return a vector with valid names for elements in the list passed to the
# 'basehaz_ops' argument of a 'stan_jm' or 'stan_surv' call
#
# @param basehaz_name A character string, the type of baseline hazard.
# @return A character vector, or NA if unmatched.
get_ok_basehaz_ops <- function(basehaz_name) {
  switch(basehaz_name,
         weibull   = c(),
         bs        = c("df", "knots"),
         piecewise = c("df", "knots"),
         ms        = c("df", "knots"),
         NA)
}

# Return the design matrix for the baseline hazard
#
# @param times A vector of times at which to evaluate the baseline hazard
# @param basehaz A named list with info about the baseline hazard,
#   returned by a call to handle_basehaz
# @return A matrix
make_basehaz_X <- function(times, basehaz, deriv = FALSE) {
  
	if (basehaz$type_name == "exponential") {
  
    X <- matrix(0, nrow = length(times), ncol = 0L) # dud matrix for Stan
  
	} else if (basehaz$type_name == "weibull") {
  
    X <- matrix(log(times), nrow = length(times), ncol = 1) 
  
	} else if (basehaz$type_name == "bs") {
  
    if (is.null(basehaz$basis))
      stop2("Bug found: could not find info on spline basis.")
	  
	  X <- as.array(predict(basehaz$basis, times))
	  
  } else if (basehaz$type_name == "piecewise") {
	
    if (is.null(basehaz$knots))
      stop2("Bug found: could not find info on knot locations.")
    
    X <- dummy_matrix(times, knots = basehaz$knots)
    
  } else if (basehaz$type_name == "ms") {
    
    timescale <- basehaz$timescale
    if (is.null(timescale)) {
      tt <- times
    } else if (timescale == "log") {
      tt <- log(times)
    }
    
    if (is.null(basehaz$basis))
      stop2("Bug found: could not find info on spline basis.")
    
    if (deriv) { # M-splines, i.e. derivative of I-spline basis
      X <- as.array(deriv(predict(basehaz$basis, tt)))
    } else { # I-splines
      X <- as.array(predict(basehaz$basis, tt))
    }
    
  } else {
    
    stop2("Bug found: type of baseline hazard unknown.") 
  
  }
  X
}

# Create a dummy indicator matrix for time intervals defined by 'knots'
#
# @param x A numeric vector with the original data.
# @param knots The cutpoints defining the desired categories of 'x'.
# @return A dummy matrix.
dummy_matrix <- function(x, knots) {
  n_intervals <- length(knots) - 1
  interval <- cut(x, knots, include.lowest = TRUE, labels = FALSE)
  out <- matrix(NA, length(interval), n_intervals)
  for (i in 1:nvars) 
    out[, i] <- ifelse(interval == i, 1, 0)
  as.matrix(out)
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
  cbind(unclass(data[,1]), data[, -1, drop = FALSE], stringsAsFactors = FALSE)
}

#--------------- Functions related to association structure

# Return a named list summarising the association for one submodel
# 
# @param user_x A character vector or NULL, being the user input to the
#   assoc argument (for one submodel) in the stan_jm call
# @param y_mod A list returned by a call to handle_glmod
# @param id_var The name of the ID variable 
# @param metastuff A list with meta information about the joint model.
# @return A list with information about the desired association structure.
handle_assoc <- function(user_assoc, user_lag, y_mod, meta) {
  
  M              <- meta$M
  id_var         <- meta$id_var
  has_icens      <- meta$has_icens
  ok_assoc       <- meta$ok_assoc
  ok_assoc_data  <- meta$ok_assoc_data
  ok_assoc_int   <- meta$ok_assoc_int
  ok_assoc_icens <- meta$ok_assoc_icens
  
  ok_inputs <- c(ok_assoc, paste0(ok_assoc_data, "_data"),
                 uapply(ok_assoc_int, paste0, "_", ok_assoc_int))
  
  disallowed <- c("muslope", "shared_b", "shared_coef")

  trimmed_assoc <- trim_assoc(user_assoc, ok_assoc_data, ok_assoc_int)
  
  if (not.null(user_assoc) && !all(trimmed_assoc %in% ok_inputs))
    stop2("An unsupported association type has been specified. The ",
         "'assoc' argument can only include the following association ", 
         "types: ", comma(ok_assoc), ", as well as possible interactions ",
         "either between association terms or with observed data.")      
  
  if (any(trimmed_assoc %in% disallowed))
    stop2("The following association structures are temporarily disallowed: ",
          "and may be reinstated in a future release: ", comma(disallowed))
    
  if (has_icens && !all(trimmed_assoc %in% ok_assoc_icens))
    stop2("When interval censoring is present, only the following association ",
          "structures are allowed:", comma(ok_assoc_icens))
  
  if (not.null(user_assoc) && !is.character(user_assoc))
    stop2("The 'assoc' argument should be a character vector or, for a ",
          "multivariate joint model, possibly a list of character vectors.") 
  
  assoc <- sapply(ok_inputs, `%in%`, trimmed_assoc, simplify = FALSE)
  if (is.null(user_assoc)) {
    assoc$null <- TRUE
  } else {
    if (assoc$null && (length(user_assoc) > 1L))
      stop2("In assoc, 'null' cannot be specified in conjuction ",
           "with another association type.")
    STOP_combination_not_allowed(assoc, "etavalue", "muvalue")
    STOP_combination_not_allowed(assoc, "etaslope", "muslope")
    STOP_combination_not_allowed(assoc, "etaauc",   "muauc")
  }
  
  # Parse suffix specifying indices for shared random effects
  cnms <- y_mod$z$group_cnms
  cnms_id <- cnms[[id_var]] # names of random effect terms
  assoc$which_b_zindex <- parse_assoc_sharedRE("shared_b",    user_assoc, 
                                                  max_index = length(cnms_id), cnms_id)
  assoc$which_coef_zindex <- parse_assoc_sharedRE("shared_coef", user_assoc, 
                                                  max_index = length(cnms_id), cnms_id)
  
  if (length(intersect(assoc$which_b_zindex, assoc$which_coef_zindex)))
    stop("The same random effects indices should not be specified in both ",
         "'shared_b' and 'shared_coef'. Specifying indices in 'shared_coef' ",
         "will include both the fixed and random components.", call. = FALSE)
  
  if (length(assoc$which_coef_zindex)) {
    if (length(cnms) > 1L)
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
    }, beta_nms = colnames(y_mod$X$X))
  } else assoc$which_coef_xindex <- numeric(0)
  
  if (!identical(length(assoc$which_coef_zindex), length(assoc$which_coef_xindex)))
    stop("Bug found: the lengths of the fixed and random components of the ",
         "'shared_coef' association structure are not the same.")
  
  # Parse suffix specifying formula for interactions with data
  ok_inputs_data <- paste0(ok_assoc_data, "_data")
  assoc$which_formulas <- sapply(ok_inputs_data, parse_assoc_data, user_assoc, simplify = FALSE) 
  
  # Parse suffix specifying indices for interactions between association terms
  ok_inputs_interactions <- unlist(lapply(ok_assoc_int, paste0, "_", ok_assoc_int))
  assoc$which_interactions <- sapply(ok_inputs_interactions, parse_assoc_interactions, 
                                     user_assoc, max_index = M, simplify = FALSE)
  
  # Lag for association structure
  assoc$which_lag <- user_lag
  
  assoc
}

# Check whether an association structure was specified that is not allowed
# when there is an additional grouping factor clustered within patients
#
# @param has_grp Logical vector specifying where each of the 1:M submodels
#   has a grp factor clustered within patients or not.
# @param assoc A two dimensional array with information about desired association
#   structure for the joint model (returned by a call to validate_assoc). 
# @param ok_assocs_with_grp A character vector with the rownames in assoc
#   that are allowed association structures when there is a grp factor 
#   clustered within patients.
validate_assoc_with_grp <- function(grp_stuff, assoc, ok_assoc_with_grp) {
  
  has_grp       <- fetch_(grp_stuff, "has_grp")
  grp_structure <- fetch (grp_stuff, "grp_list")[has_grp]
  
  if (n_distinct(grp_structure) > 1L)
    stop2("Any longitudinal submodels with a grouping factor clustered within ",
          "patients must use the same clustering structure; that is, the same ",
          "clustering variable and the same number of units clustered within a ",
          "given patient.")
  
  all_nms <- grep("which|null", rownames(assoc), invert = TRUE, value = TRUE)
  disallowed_nms <- setdiff(all_nms, ok_assoc_with_grp)
  sel <- unlist(assoc[disallowed_nms, which(has_grp)])
  if (any(sel))
    stop2("Only the following association structures are allowed when there is a ",
          "grouping factor clustered within individuals: ", comma(ok_assoc_with_grp))
}

# Validate the user input to the lag_assoc argument of stan_jm
#
# @param lag_assoc The user input to the lag_assoc argument
# @param M Integer specifying the number of longitudinal submodels
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
  lag_assoc
}

# Remove suffixes from the user inputted assoc argument
#
# @param x A character vector, being the user input to the 
#   assoc argument in the stan_jm call
# @param ok_assoc_data A character vector specifying which types
#   of association terms are allowed to be interacted with data
# @param ok_assoc_int A character vector specifying which types
#   of association terms are allowed to be interacted with other 
#   association terms
trim_assoc <- function(x, ok_assoc_data, ok_assoc_int) {
  x <- gsub("^shared_b\\(.*",    "shared_b",    x) 
  x <- gsub("^shared_coef\\(.*", "shared_coef", x) 
  for (i in ok_assoc_data)
    x <- gsub(paste0("^", i, "_data\\(.*"),    paste0(i, "_data"), x)
  for (i in ok_assoc_int) for (j in ok_assoc_int)
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
# @param assoc A two dimensional array with information about desired association
#   structure for the joint model (returned by a call to validate_assoc). 
# @param ok_assoc_int A character vector, specifying which association
#   structures are allowed to be used in interactions
check_order_of_assoc_interactions <- function(assoc, ok_assoc_int) {
  M <- ncol(assoc)
  for (i in ok_assoc_int) {
    for (j in ok_assoc_int) {
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

# Return a matrix with information for stan about which components are
# required for building the association structure of the joint model
#
# @param assoc An array with information on the joint model association struct.
# @return An indicator matrix, with a row for each component type and a 
#   column for each longitudinal submodel.
make_assoc_component_flags <- function(assoc) {
  component_types <- c("etavalue", 
                       "etaslope", 
                       "etaauc", 
                       "muvalue", 
                       "muslope", 
                       "muauc")
  assoc_uses <- sapply(component_types,
    function(x, assoc) {
      nm_check <- switch(x,
                         etavalue = "^eta|^mu",
                         etaslope = "etaslope|muslope",
                         etaauc   = "etaauc|muauc",
                         muvalue  = "muvalue|muslope",
                         muslope  = "muslope",
                         muauc    = "muauc")
      sel <- grep(nm_check, rownames(assoc))
      tmp <- assoc[sel, , drop = FALSE]
      tmp <- pad_matrix(tmp, cols = 3L, value = FALSE)
      ai(as.logical(colSums(tmp > 0)))
    }, assoc = assoc)
  t(assoc_uses)
}

# Return a matrix with information for stan about which association structures
# are to be used in the association structure
#
# @param assoc An array with information on the joint model association struct.
# @return An indicator matrix, with a row for each possible type of association
#   structure and a column for each longitudinal submodel.
make_assoc_type_flags <- function(assoc) {
  nms <- rownames(assoc)
  sel <- grep("which|null", nms, invert = TRUE)
  matrix(ai(assoc[sel,]), ncol = ncol(assoc))
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
#   has_grp: logical specifying whether the submodel has a grouping factor
#     that is clustered with patients.
#   grp_var: the name of any grouping factor that is clustered with patients.
#   grp_assoc: the user input to the grp_assoc argument in the stan_jm call.
#   grp_freq: a named vector with the number of lower level units clustered
#     within each individual.
#   grp_list: a named list containing the unique names for the lower level 
#     units clustered within each individual.
get_basic_grp_info <- function(y_mod, id_var) {
  cnms  <- y_mod[["z"]][["group_cnms"]]
  flist <- y_mod[["z"]][["group_list"]]
  cnms_nms <- names(cnms)
  tally <- xapply(cnms_nms, FUN = function(x) 
    # within each ID, count the number of levels for the grouping factor x
    tapply(flist[[x]], flist[[id_var]], FUN = n_distinct))
  sel <- which(sapply(tally, function(x) !all(x == 1L)) == TRUE)
  has_grp <- as.logical(length(sel))
  if (!has_grp) {
    return(nlist(has_grp))
  } else {
    if (length(sel) > 1L)
      stop("There can only be one grouping factor clustered within 'id_var'.")
    grp_var <- cnms_nms[sel] 
    return(nlist(has_grp, grp_var))
  }
}

get_extra_grp_info <- function(basic_info, flist, id_var, grp_assoc,
                               ok_assoc_grp = c("sum", "mean", "min", "max")) {
  has_grp <- basic_info$has_grp
  grp_var <- basic_info$grp_var
  if (!has_grp) { # no grouping factor clustered within patients
    return(basic_info)
  } else { # submodel has a grouping factor clustered within patients
    if (is.null(grp_var))
      stop2("Bug found: could not find 'grp_var' in basic_info.")
    if (is.null(grp_assoc))
      stop2("'grp_assoc' cannot be NULL when there is a grouping factor ",
            "clustered within patients.")       
    if (!grp_assoc %in% ok_assoc_grp)
      stop2("'grp_assoc' must be one of: ", comma(ok_assoc_grp))
    
    # cluster and patient ids for each row of the z matrix
    factor_grp <- factor(flist[[grp_var]]) 
    factor_ids <- factor(flist[[id_var]])
    
    # num clusters within each patient
    grp_freq <- tapply(factor_grp, factor_ids, FUN = n_distinct, simplify = FALSE)
    grp_freq <- unlist(grp_freq)
    
    # unique cluster ids for each patient id
    grp_list <- tapply(factor_grp, factor_ids, FUN = unique, simplify = FALSE)

    basic_info <- nlist(has_grp, grp_var)
    extra_info <- nlist(grp_assoc, grp_freq, grp_list)
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


#--------------- Functions related to generating initial values

# Create a function that can be used to generate the model-based initial values for Stan
#
# @param e_mod_stuff A list object returned by a call to the handle_coxmod function
# @param standata The data list that will be passed to Stan
get_prefit_inits <- function(init_fit, e_mod, standata) {
  
  init_means <- rstan::get_posterior_mean(init_fit)
  init_nms   <- rownames(init_means)
  inits <- generate_init_function(e_mod, standata)()
  
  sel_b1 <- grep(paste0("^z_bMat1\\."), init_nms)
  if (length(sel_b1))
    inits[["z_bMat1"]] <- matrix(init_means[sel_b1,], nrow = standata$bK1)
  
  sel_b2 <- grep(paste0("^z_bMat2\\."), init_nms)
  if (length(sel_b2))
    inits[["z_bMat2"]] <- matrix(init_means[sel_b2,], nrow = standata$bK2)
  
  sel_bC1 <- grep(paste0("^bCholesky1\\."), init_nms)
  if (length(sel_bC1) > 1) {
    inits[["bCholesky1"]] <- matrix(init_means[sel_bC1,], nrow = standata$bK1)
  } else if (length(sel_bC1) == 1) {
    inits[["bCholesky1"]] <- aa(init_means[sel_bC1,])
  }
  
  sel_bC2 <- grep(paste0("^bCholesky2\\."), init_nms)
  if (length(sel_bC2) > 1) {
    inits[["bCholesky2"]] <- matrix(init_means[sel_bC2,], nrow = standata$bK2)
  } else if (length(sel_bC1) == 1) {
    inits[["bCholesky2"]] <- aa(init_means[sel_bC2,])
  }      
  
  sel <- c("yGamma1", "yGamma2", "yGamma3", 
           "z_yBeta1", "z_yBeta2", "z_yBeta3",
           "yAux1_unscaled", "yAux2_unscaled", "yAux3_unscaled", 
           "bSd1", "bSd2", "z_b", "z_T", "rho", "zeta", "tau", 
           "yGlobal1", "yGlobal2", "yGlobal3", 
           "yLocal1", "yLocal2", "yLocal3", 
           "yMix1", "yMix2", "yMix3", 
           "yOol1", "yOol2", "yOol3")
  for (i in sel) {
    sel_i <- grep(paste0("^", i, "\\."), init_nms)
    if (length(sel_i))
      inits[[i]] <- aa(init_means[sel_i,])
  }
  return(function() inits)
}

generate_init_function <- function(e_mod_stuff, standata) {
  
  # Initial values for intercepts, coefficients and aux parameters
  if (e_mod_stuff$surv_type %in% c("right", "counting")) {
    e_beta <- e_mod_stuff$mod$coef
  } else if (e_mod_stuff$surv_type %in% c("interval", "interval2")) {
    e_beta <- -drop_intercept(e_mod_stuff$mod$coef) * e_mod_stuff$mod$scale
  } else {
    stop("Bug found: Invalid Surv type.")
  }
  e_aux <- if (standata$basehaz_type == 1L) runif(1, 0.5, 3) else rep(0, standata$basehaz_nvars)
  e_z_beta       <- standardise_coef(x        = e_beta, 
                                     location = standata$e_prior_mean, 
                                     scale    = standata$e_prior_scale) 
  e_aux_unscaled <- standardise_coef(x        = e_aux, 
                                     location = standata$e_prior_mean_for_aux, 
                                     scale    = standata$e_prior_scale_for_aux)

  # Function to generate model based initial values
  model_based_inits <- rm_null(list(
    e_z_beta       = array_else_double(e_z_beta),
    e_aux_unscaled = array_else_double(e_aux_unscaled),
    e_gamma        = array_else_double(rep(0, standata$e_has_intercept))))
  
  return(function() model_based_inits)
}


#--------------- Functions related to standata and sampling

# Set arguments for sampling for stan_jm
#
# Prepare a list of arguments to use with \code{rstan::sampling} via
# \code{do.call}.
#
# *Note that this differs from the set_sampling_args function in that
# it uses a different default adapt_delta and max_treedepth. Using a 
# shorter treedepth seems to stop the sampler trailing off during early 
# iterations and can drastically reduce the model estimation time, and 
# in most examples using a shorter treedepth hasn't compromised the sampler
# at later interations (ie, at later iterations the sampler doesn't
# hit the maximum treedepth). The default adapt_delta depends on the 
# largest number of group-specific parameters for any single grouping
# factor in the model.
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
  
  default_adapt_delta <- if (max_p > 2) 0.85 else 0.80
  default_max_treedepth <- 10L
  
  if (!is.null(user_adapt_delta))
    args$control$adapt_delta <- user_adapt_delta else 
      if (is.null(args$control$adapt_delta))
        args$control$adapt_delta <- default_adapt_delta
  
  if (!is.null(user_max_treedepth))
    args$control$max_treedepth <- user_max_treedepth else
      if (is.null(args$control$max_treedepth))
        args$control$max_treedepth <- default_max_treedepth
  
  if (!"save_warmup" %in% unms) 
    args$save_warmup <- FALSE  
  
  return(args)
}  

# Return the list of pars for Stan to monitor
# 
# @param standata The list of data to pass to Stan
# @param is_jm A logical
# @return A character vector
pars_to_monitor <- function(standata, is_jm = FALSE) {
  c(if (standata$M > 0 && standata$intercept_type[1]) "yAlpha1", 
    if (standata$M > 1 && standata$intercept_type[2]) "yAlpha2", 
    if (standata$M > 2 && standata$intercept_type[3]) "yAlpha3", 
    if (standata$M > 0 && standata$yK[1]) "yBeta1",
    if (standata$M > 1 && standata$yK[2]) "yBeta2",
    if (standata$M > 2 && standata$yK[3]) "yBeta3",
    if (is_jm) "e_alpha",
    if (is_jm && standata$e_K) "e_beta",
    if (is_jm && standata$a_K) "a_beta",
    if (standata$bK1 > 0) "b1",
    if (standata$bK2 > 0) "b2",
    if (standata$M > 0 && standata$has_aux[1]) "yAux1",
    if (standata$M > 1 && standata$has_aux[2]) "yAux2",
    if (standata$M > 2 && standata$has_aux[3]) "yAux3",
    if (is_jm && length(standata$basehaz_nvars)) "e_aux",
    if (standata$prior_dist_for_cov == 2 && standata$bK1 > 0) "bCov1",
    if (standata$prior_dist_for_cov == 2 && standata$bK2 > 0) "bCov2",
    if (standata$prior_dist_for_cov == 1 && standata$len_theta_L) "theta_L",
    "mean_PPD")
}

# Change the MCMC samples for theta_L to Sigma
#
# @param stanfit The stanfit object from the fitted model
# @param cnms The component names for the group level terms, combined
#   across all glmer submodels
# @return A stanfit object
evaluate_Sigma <- function(stanfit, cnms) {
  nc <- sapply(cnms, FUN = length)
  nms <- names(cnms) 
  thetas <- extract(stanfit, pars = "theta_L", inc_warmup = TRUE, 
                    permuted = FALSE)
  Sigma <- apply(thetas, 1:2, FUN = function(theta) {
    Sigma <- mkVarCorr(sc = 1, cnms, nc, theta, nms)
    unlist(sapply(Sigma, simplify = FALSE, 
                  FUN = function(x) x[lower.tri(x, TRUE)]))
  })
  l <- length(dim(Sigma))
  end <- tail(dim(Sigma), 1L)
  shift <- grep("^theta_L", names(stanfit@sim$samples[[1]]))[1] - 1L
  if (l == 3) for (chain in 1:end) for (param in 1:nrow(Sigma)) {
    stanfit@sim$samples[[chain]][[shift + param]] <- Sigma[param, , chain] 
  }
  else for (chain in 1:end) {
    stanfit@sim$samples[[chain]][[shift + 1]] <- Sigma[, chain]
  }
  stanfit
}

# Get the names for the Sigma var-cov matrix
#
# @param cnms The component names for the group level terms, combined
#   across all glmer submodels
# @return A character vector
get_Sigma_nms <- function(cnms) {
  nms <- names(cnms) 
  Sigma_nms <- lapply(cnms, FUN = function(grp) {
    nm <- outer(grp, grp, FUN = paste, sep = ",")
    nm[lower.tri(nm, diag = TRUE)]
  })
  for (j in seq_along(Sigma_nms)) {
    Sigma_nms[[j]] <- paste0(nms[j], ":", Sigma_nms[[j]])
  }
  paste0("Sigma[", unlist(Sigma_nms), "]")
}


#--------------- Functions related to observation weights

# Check the weights argument for stan_jm
#
# @param weights The data frame passed via the weights argument
# @param id_var The name of the ID variable
check_weights <- function(weights, id_var) {
  
  if (is.null(weights))
    return(weights)
  
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
  
  is_glmod <- (is.null(mod_stuff$exittime))
  
  # No weights provided by user
  if (is.null(weights)) {
    len <- if (is_glmod) length(mod_stuff$Y$Y) else 0
    return(rep(0.0, len)) 
  }
  
  # Check for IDs with no weight supplied
  weights[[id_var]] <- factor(weights[[id_var]])
  ids <- if (is_glmod) mod_stuff$Z$group_list[[id_var]] else factor(mod_stuff$id_list)
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

# Parse the model formula
#
# @param formula The user input to the formula argument.
# @param data The user input to the data argument (i.e. a data frame).
parse_formula <- function(formula, data) {
  
  formula <- validate_formula(formula, needs_response = TRUE)
  
  lhs        <- lhs(formula) # full LHS of formula
  lhs_form   <- reformulate_lhs(lhs)
  
  rhs        <- rhs(formula)         # RHS as expression
  rhs_form   <- reformulate_rhs(rhs) # RHS as formula
  rhs_terms  <- terms(rhs_form, specials = "tde")
  rhs_vars   <- rownames(attr(rhs_terms, "factors"))
  
  allvars      <- all.vars(formula)
  allvars_form <- reformulate(allvars)

  surv <- eval(lhs, envir = data) # Surv object
  surv <- validate_surv(surv)
  type <- attr(surv, "type")
  
  if (type == "right") {
    tvar_beg <- NULL
    tvar_end <- as.character(lhs[[2L]])
    dvar     <- as.character(lhs[[3L]])
    min_t    <- 0
    max_t    <- max(surv[, "time"])
  } else if (type == "counting") {
    tvar_beg <- as.character(lhs[[2L]])
    tvar_end <- as.character(lhs[[3L]])
    dvar     <- as.character(lhs[[4L]])
    min_t    <- min(surv[, "start"])
    max_t    <- max(surv[, "stop"])
  } else if (type == "interval") {
    tvar_beg <- NULL
    tvar_end <- as.character(lhs[[2L]])
    dvar     <- as.character(lhs[[4L]])
    min_t    <- 0
    max_t    <- max(surv[, c("time1", "time2")])
  } else if (type == "interval2") {
    tvar_beg <- NULL
    tvar_end <- as.character(lhs[[2L]])
    dvar     <- as.character(lhs[[4L]])
    min_t    <- 0
    max_t    <- max(surv[, c("time1", "time2")])
  }
  
  sel <- attr(rhs_terms, "specials")$tde

  if (!is.null(sel)) { # model has tde
    
    # replace 'tde(x, ...)' in formula with 'x'
    tde_oldvars <- rhs_vars
    tde_newvars <- sapply(tde_oldvars, function(oldvar) {
      if (oldvar %in% rhs_vars[sel]) {
        tde <- function(newvar, ...) { # define tde function locally
          safe_deparse(substitute(newvar)) 
        }
        eval(parse(text = oldvar))
      } else oldvar
    }, USE.NAMES = FALSE)
    for (i in sel) {
      sel_terms <- which(attr(rhs_terms, "factors")[i, ] > 0)
      for (j in sel_terms) {
        term_labels[j] <- gsub(tde_oldvars[i], 
                               tde_newvars[i], 
                               term_labels[j], 
                               fixed = TRUE)
      }
    }
    tf_form <- reformulate(term_labels, response = lhs)
    
    # extract 'tde(x, ...)' from formula and construct 'bs(times, ...)'
    tde_terms <- lapply(rhs_vars[sel], function(x) {
      tde <- function(vn, ...) { # define tde function locally
        dots <- list(...)
        ok_args <- c("df")
        if (!isTRUE(all(names(dots) %in% ok_args)))
          stop2("Invalid argument to 'tde' function. ",
                "Valid arguments are: ", comma(ok_args))
        df <- if (is.null(dots$df)) 3 else dots$df
        degree <- 3
        if (df == 3) {
          dots[["knots"]] <- numeric(0)
        } else {
          dx <- (max_t - min_t) / (df - degree + 1)
          dots[["knots"]] <- seq(min_t + dx, max_t - dx, dx)
        }
        dots[["Boundary.knots"]] <- c(min_t, max_t) 
        sub("^list\\(", "bs\\(times__, ", deparse(dots))
      }
      tde_calls <- eval(parse(text = x))
      sel_terms <- which(attr(rhs_terms, "factors")[x, ] > 0)
      new_calls <- sapply(seq_along(sel_terms), function(j) {
        paste0(term_labels[sel_terms[j]], ":", tde_calls)
      })
      nlist(tde_calls, new_calls)
    })
    td_basis <- fetch(tde_terms, "tde_calls")
    new_calls <- fetch_(tde_terms, "new_calls")
    td_form <- reformulate(new_calls, response = NULL, intercept = FALSE)
    
  } else { # model doesn't have tde
    tf_form  <- formula
    td_form  <- NULL
    td_basis <- NULL
  }
  
  nlist(formula,
        lhs,
        rhs,
        lhs_form,
        rhs_form,
        tf_form,
        td_form,
        td_basis,
        fe_form = rhs_form, # no re terms accommodated yet
        re_form = NULL,     # no re terms accommodated yet
        allvars,
        allvars_form,
        tvar_beg,
        tvar_end,
        dvar,
        surv_type = attr(surv, "type"))
}

# Data for handling lower-level clustering in association structure
standata_add_assoc_grp <- function(standata, a_mod, grp_stuff) {
  has_grp <- fetch_(grp_stuff, "has_grp")
  if (any(has_grp)) {
    sel <- which(has_grp)[[1L]]
    standata$idx_grp <- attr(a_mod[[sel]], "grp_idx")
    standata$grp_assoc <- switch(grp_assoc, 
                                 sum  = 1L,
                                 mean = 2L,
                                 min  = 3L,
                                 max  = 4L,
                                 0L)
  } else { # no lower level clustering
    standata$idx_grp   <- matrix(0L, standata$len_cpts, 2L)
    standata$grp_assoc <- 0L
  } 
  standata$has_grp <- aa(ai(has_grp))
  standata
}


# Data (design matrices) for the main types of association structure
standata_add_assoc_xz <- function(standata, a_mod, meta, assoc) {
  cnms_nms <- names(meta$cnms)
  N_tmp <- sapply(a_mod, function(x) NROW(x$mod_eta$xtemp))
  N_tmp <- c(N_tmp, rep(0, 3 - length(N_tmp)))
  standata$y_qrows <- as.array(as.integer(N_tmp))
  for (m in 1:3) {
    for (i in c("eta", "eps", "auc")) {
      nm_check <- switch(i,
                         eta = "^eta|^mu",
                         eps = "slope",
                         auc = "auc")
      sel <- grep(nm_check, rownames(assoc))
      if (m <= meta$M && any(unlist(assoc[sel,m]))) {
        tmp_stuff <- a_mod[[m]][[paste0("mod_", i)]]
        
        # fe design matrix 
        X_tmp <- tmp_stuff$xtemp
        
        # re design matrix, group factor 1
        Z1_tmp    <- tmp_stuff$z[[cnms_nms[1L]]]
        Z1_tmp    <- transpose(Z1_tmp)
        Z1_tmp    <- convert_null(Z1_tmp, "matrix")
        Z1_tmp_id <- tmp_stuff$group_list[[cnms_nms[1L]]]
        Z1_tmp_id <- groups(Z1_tmp_id)
        Z1_tmp_id <- convert_null(Z1_tmp_id, "arrayinteger")
        
        # re design matrix, group factor 2
        if (length(cnms_nms) > 1L) {
          Z2_tmp    <- tmp_stuff$z[[cnms_nms[2L]]]
          Z2_tmp    <- transpose(Z2_tmp)
          Z2_tmp    <- convert_null(Z2_tmp, "matrix")
          Z2_tmp_id <- tmp_stuff$group_list[[cnms_nms[2L]]]
          Z2_tmp_id <- groups(Z2_tmp_id)
          Z2_tmp_id <- convert_null(Z2_tmp_id, "arrayinteger")
        } else {
          Z2_tmp    <- matrix(0,standata$bK2_len[m],0) 
          Z2_tmp_id <- aa(integer(0))
        }
      } else {
        X_tmp     <- matrix(0,0,standata$yK[m])
        Z1_tmp    <- matrix(0,standata$bK1_len[m],0) 
        Z2_tmp    <- matrix(0,standata$bK2_len[m],0) 
        Z1_tmp_id <- aa(integer(0)) 
        Z2_tmp_id <- aa(integer(0)) 
      }
      standata[[paste0("y", m, "_x_",     i)]] <- X_tmp
      standata[[paste0("y", m, "_z1_",    i)]] <- Z1_tmp
      standata[[paste0("y", m, "_z2_",    i)]] <- Z2_tmp
      standata[[paste0("y", m, "_z1_id_", i)]] <- Z1_tmp_id
      standata[[paste0("y", m, "_z2_id_", i)]] <- Z2_tmp_id              
    }             
  }
  standata
}

# Dimensions for auc association structure
standata_add_assoc_auc <- function(standata, a_mod, meta) {
  
  count_rows <- function(i) { 
    nr <- NROW(i$mod_auc$x)
    if (nr > 0) nr else NULL 
  }
  
  nrow_y_x_auc <- unique(uapply(a_mod, count_rows))
  if (length(nrow_y_x_auc) > 1L)
    stop2("Bug found: nrows for auc should be the same for all submodels.")
  qnodes <- meta$auc_qnodes
  uses_auc_assoc <- any(standata$assoc_uses[3,] > 0)
  if (uses_auc_assoc) {
    qw   <- get_quadpoints(qnodes)$weights
    qfun <- function(i) lapply(qw, unstandardise_qwts, 0, i)
    standata$auc_qwts <- aa(uapply(standata$cpts, qfun))
  } else {
    standata$auc_qwts <- double(0)
  }
  standata$auc_qnodes <- ai(qnodes)
  standata$y_qrows_for_auc <- ai(nrow_y_x_auc %ORifNULL% 0)
  standata
}

# Data for interaction-based and shared parameter association structures
standata_add_assoc_extras <- function(standata, a_mod, assoc) {
  
  # interactions between association terms and data, with the following objects:
  #   a_K_data: number of columns in y_Xq_data corresponding to each interaction 
  #     type (ie, etavalue, etaslope, muvalue, muslope) for each submodel
  #   idx_q: indexing for the rows of Xq_data that correspond to each submodel, 
  #     since it is formed as a block diagonal matrix
  xq_data <- fetch(a_mod, "X_bind_data") # design mat for the interactions
  standata$y_x_data <- aa(am(Matrix::bdiag(xq_data)))
  standata$a_K_data <- fetch_array(a_mod, "K_data")
  standata$idx_data <- get_idx_array(standata$y_qrows)
  
  # interactions between association terms
  wi <- "which_interactions"
  standata$which_interactions      <- aa(unlist(assoc[wi,]))
  standata$size_which_interactions <-  c(sapply(assoc[wi,], sapply, length))
  
  # shared random effects
  wbz <- "which_b_zindex"
  wcz <- "which_coef_zindex"
  wcx <- "which_coef_zindex"
  standata$which_b_zindex    <- aa(unlist(assoc[wbz,]))
  standata$which_coef_zindex <- aa(unlist(assoc[wcz,]))
  standata$which_coef_xindex <- aa(unlist(assoc[wcx,]))
  standata$size_which_b      <- aa(sapply(assoc[wbz,], length))
  standata$size_which_coef   <- aa(sapply(assoc[wcz,], length))
  
  # sum dimensions
  for (i in c("a_K_data", paste0("size_which_", c("b", "coef", "interactions")))) {
    standata[[paste0("sum_", i)]] <- ai(sum(standata[[i]]))
  }
  standata
}
