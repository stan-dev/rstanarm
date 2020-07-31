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
  
  if (!identical(NULL, response) && is.gaussian(family$family)) { 
    # use response variable for scaling priors
    if (ps$prior_dist > 0L && ps$prior_autoscale) {
      ss <- sd(response)
      ps$prior_scale <- ss * ps$prior_scale
    }
  }
  
  if (!identical(NULL, predictors) && !QR) {
    # use predictors for scaling priors
    if (ps$prior_dist > 0L && ps$prior_autoscale) {
      ps$prior_scale <- 
        pmax(min_prior_scale,
             ps$prior_scale / apply(predictors, 2L, get_scale_value))
    }      
  }
  
  if (!identical(NULL, assoc)) {
    # Evaluate mean and SD of each of the association terms that will go into
    # the linear predictor for the event submodel (as implicit "covariates").
    # (NB the approximate association terms are calculated using coefs
    # from the separate longitudinal submodels estimated using glmer).
    # The mean will be used for centering each association term.
    # The SD will be used for autoscaling the prior for each association parameter.
    if (identical(NULL, family))
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
# @param formula The model formula for the glmer submodel.
# @param data The data for the glmer submodel.
# @param family The family object for the glmer submodel.
# @return A named list with the following elements:
#   y: named list with the reponse vector and related info.
#   x: named list with the fe design matrix and related info.
#   z: named list with the re design matrices and related info.
#   terms: the model.frame terms object with bars "|" replaced by "+".
#   model_frame: The model frame with all variables used in the 
#     model formula.
#   formula: The model formula.
#   reTrms: returned by lme4::glFormula$reTrms.
#   family: the (modified) family object for the glmer submodel.
#   intercept_type: named list with info about the type of 
#     intercept required for the glmer submodel.
#   has_aux: logical specifying whether the glmer submodel 
#     requires an auxiliary parameter.
handle_y_mod <- function(formula, data, family) {
  mf <- stats::model.frame(lme4::subbars(formula), data)
  if (!length(formula) == 3L)
    stop2("An outcome variable must be specified.")
  
  # lme4 parts
  lme4_parts <- lme4::glFormula(formula, data)
  reTrms <- lme4_parts$reTrms
  
  # Response vector, design matrices
  y <- make_y_for_stan(formula, mf, family) 
  x <- make_x_for_stan(formula, mf)
  z <- make_z_for_stan(formula, mf) 
  
  # Terms
  terms <- attr(mf, "terms")
  terms <- append_predvars_attribute(terms, formula, data)
  
  # Binomial with >1 trials not allowed by stan_{mvmver,jm}
  is_binomial <- is.binomial(family$family)
  is_bernoulli <- is_binomial && NCOL(y$y) == 1L && all(y$y %in% 0:1)
  if (is_binomial && !is_bernoulli)
    STOP_binomial()
  
  # Various flags
  intercept_type <- check_intercept_type(x, family)
  has_aux <- check_for_aux(family)
  family <- append_mvmer_famlink(family, is_bernoulli)
  
  nlist(y, x, z, reTrms, model_frame = mf, formula, terms, 
        family, intercept_type, has_aux)
}

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
check_id_var <- function(id_var, y_cnms, y_flist) {
  len_cnms <- sapply(y_cnms, length)
  if (any(len_cnms > 1L)) {  # more than one grouping factor
    if (is.null(id_var)) {
      stop("'id_var' must be specified when using more than one grouping factor",
           call. = FALSE)
    } else {
      lapply(y_cnms, function(x)  if (!(id_var %in% names(x)))
        stop("'id_var' must be included as a grouping factor in each ",
             "of the longitudinal submodels", call. = FALSE)) 
    }
    return(id_var)
  } else {  # only one grouping factor (assumed to be subject ID)
    only_cnm <- unique(sapply(y_cnms, names))
    if (length(only_cnm) > 1L)
      stop("The grouping factor (ie, subject ID variable) is not the ",
           "same in all longitudinal submodels", call. = FALSE)
    if ((!is.null(id_var)) && (!identical(id_var, only_cnm)))
      warning("The user specified 'id_var' (", paste(id_var), 
              ") and the assumed ID variable based on the single ",
              "grouping factor (", paste(only_cnm), ") are not the same; ", 
              "'id_var' will be ignored", call. = FALSE, immediate. = TRUE)
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

# Function to return a single list with the factor levels for each
# grouping factor, but collapsed across all longitudinal submodels
# 
# @param x A list containing the flist object for each of the submodels
get_common_flevels <- function(x) {
  nms <- lapply(x, names)
  unique_nms <- unique(unlist(nms))
  flevels <- lapply(seq_along(unique_nms), function(i) {
    nm <- unique_nms[i]
    flevels_nm <- lapply(1:length(x), function(m) 
      if (nm %in% nms[[m]]) levels(x[[m]][[nm]]))
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
validate_observation_times <-function(data, eventtimes, id_var, time_var) {
  if (!time_var %in% colnames(data)) 
    STOP_no_var(time_var)
  if (!id_var %in% colnames(data)) 
    STOP_no_var(id_var)
  if (any(data[[time_var]] < 0))
    stop2("Values for the time variable (", time_var, ") should not be negative.")
  mt <- tapply(data[[time_var]], factor(data[[id_var]]), max)  # max observation time
  nms <- names(eventtimes)                                     # patient IDs
  if (is.null(nms))
    stop2("Bug found: cannot find names in the vector of event times.")
  sel <- which(sapply(nms, FUN = function(i) mt[i] > eventtimes[i]))
  if (length(sel))
    stop2("The following individuals have observation times in the longitudinal data ",
          "are later than their event time: ", paste(nms[sel], collapse = ", "))      
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
#   mod: The fitted Cox model.
#   entrytime: Named vector of numeric entry times.
#   eventtime: Named vector of numeric event times.
#   status: Named vector of event/failure indicators.
#   Npat: Number of individuals.
#   Nevents: Total number of events/failures.
#   id_list: A vector of unique subject IDs, as a factor.
#   qnodes: The number of GK quadrature nodes.
#   qwts,qpts: Vector of unstandardised quadrature weights and points.
#     The vector is ordered such that the first Npat items are the
#     weights/locations of the first quadrature point, then the second
#     Npat items are the weights/locations for the second quadrature
#     point, and so on. 
#   qids: The subject IDs corresponding to each element of qwts/qpts.
#   epts: The event times, but only for individuals who were NOT censored
#     (i.e. those individual who had an event).
#   eids: The subject IDs corresponding to each element of epts.
#   cpts: Combined vector of failure and quadrature times: c(epts, qpts).
#   cids: Combined vector subject IDs: c(eids, qids).
#   Xq: The model matrix for the event submodel, centred and no intercept.
#   Xbar: Vector of column means for the event submodel model matrix.
#   K: Number of predictors for the event submodel.
#   norm_const: Scalar, the constant used to shift the event submodel
#     linear predictor (equal to the log of the mean incidence rate). 
#   model_frame: The model frame for the fitted Cox model, but with the
#     subject ID variable also included.
#   tvc: Logical, if TRUE then a counting type Surv() object was used
#     in the fitted Cox model (ie. time varying covariates). 
handle_e_mod <- function(formula, data, qnodes, id_var, y_id_list) {
  if (!requireNamespace("survival"))
    stop("the 'survival' package must be installed to use this function")
  if (!requireNamespace("data.table"))
    stop("the 'data.table' package must be installed to use this function")
  
  mod <- survival::coxph(formula, data = data, x = TRUE)
  RHS_with_id <- paste(deparse(formula[[3L]]), "+", id_var)
  formula_with_id <- reformulate(RHS_with_id, response = formula[[2L]])
  mf1 <- model.frame(formula_with_id, data = data)
  mf1[[id_var]] <- promote_to_factor(mf1[[id_var]]) # same as lme4
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
  id_list   <- factor(mf_event[[id_var]])
  names(entrytime) <- names(eventtime) <- names(status) <- id_list
  
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
    dt <- prepare_data_table(mf2, id_var = id_var, time_var = t0_var)
    # Obtain rows of the model frame that are as close as possible to 
    # the event times (failures only) and quadrature times                      
    mf2 <- rolling_merge(dt, ids = cids, times = cpts)
    # Construct design matrix evaluated at event and quadrature times
    fm_RHS <- reformulate(attr(terms(mod), "term.labels"))
    Xq <- model.matrix(fm_RHS, data = mf2)
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
make_basehaz_X <- function(times, basehaz) {
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
  cbind(unclass(data[,1]), data[, -1, drop = FALSE], stringsAsFactors = FALSE)
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
    
    temporarily_disallowed <- c("muslope", "shared_b", "shared_coef")
    if (any(trimmed_x %in% temporarily_disallowed))
      stop2("The following association structures have been temporarily disallowed ",
            "and will be reinstated in a future release: ", 
            paste(temporarily_disallowed, collapse = ", "))
    
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
  cnms <- y_mod_stuff$z$group_cnms
  cnms_id <- cnms[[id_var]] # names of random effect terms
  assoc$which_b_zindex <- parse_assoc_sharedRE("shared_b",    user_x, 
                                                  max_index = length(cnms_id), cnms_id)
  assoc$which_coef_zindex <- parse_assoc_sharedRE("shared_coef", user_x, 
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
    }, beta_nms = colnames(y_mod_stuff$X$X))
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
validate_assoc_with_grp <- function(has_grp, assoc, ok_assocs_with_grp) {
  all_rownames <- grep("which|null", rownames(assoc), 
                       invert = TRUE, value = TRUE)
  disallowed_rows <- setdiff(all_rownames, ok_assocs_with_grp)
  sel <- which(has_grp)
  check <- unlist(assoc[disallowed_rows, sel])
  if (any(check))
    stop2("Only the following association structures are allowed when ",
          "there is a grouping factor clustered within individuals: ",
          paste(ok_assocs_with_grp, collapse = ", "))
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
# @param assoc A two dimensional array with information about desired association
#   structure for the joint model (returned by a call to validate_assoc). 
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
# @param data A data frame, the data for the longitudinal submodel.
# @param assoc A list with information about the association structure for 
#   the one longitudinal submodel. 
# @param y_mod A named list returned by a call to handle_y_mod (the
#   fit for a single longitudinal submodel)
# @param grp_stuff A list with information about any lower level grouping
#   factors that are clustered within patients and how to handle them in 
#   the association structure.
# @param ids,times The subject IDs and times vectors that correspond to the
#   event and quadrature times at which the design matrices will
#   need to be evaluated for the association structure.
# @param id_var The name on the ID variable.
# @param time_var The name of the time variable.
# @param epsilon The half-width of the central difference used for 
#   numerically calculating the derivative of the design matrix for slope
#   based association structures.
# @param auc_qnodes Integer specifying the number of GK quadrature nodes to
#   use in the integral/AUC based association structures.
# @return The list returned by make_assoc_parts.
handle_assocmod <- function(data, assoc, y_mod, grp_stuff, ids, times,  
                            id_var, time_var, epsilon, auc_qnodes) {
  
  if (!requireNamespace("data.table"))
    stop2("the 'data.table' package must be installed to use this function.")
  
  # Before turning data into a data.table (for a rolling merge
  # against the quadrature points) we want to make sure that the 
  # data does not include any NAs for the predictors or assoc formula variables
  tt <- y_mod$terms
  assoc_interaction_forms <- assoc[["which_formulas"]]
  extra_vars <- uapply(assoc_interaction_forms, function(i) {
    # loop over the four possible assoc interaction formulas and 
    # collect any variables used
    if (length(i)) {
      rownames(attr(terms.formula(i), "factors")) 
    } else NULL
  })
  rhs <- deparse(tt[[3L]], 500L)
  if (!is.null(extra_vars))
    rhs <- c(rhs, extra_vars)
  form_new <- reformulate(rhs, response = NULL)
  df <- get_all_vars(form_new, data)
  df <- df[complete.cases(df), , drop = FALSE]

  # Declare df as a data.table for merging with quadrature points
  dt <- prepare_data_table(df, id_var = id_var, time_var = time_var, 
                           grp_var = grp_stuff$grp_var) # NB grp_var may be NULL

  # Design matrices for calculating association structure based on 
  # (possibly lagged) eta, slope, auc and any interactions with data
  parts <- make_assoc_parts(use_function = make_assoc_parts_for_stan,
                            newdata = dt, assoc = assoc, id_var = id_var, 
                            time_var = time_var, grp_stuff = grp_stuff, 
                            ids = ids, times = times, epsilon = epsilon, 
                            auc_qnodes = auc_qnodes, y_mod = y_mod)
  
  # If association structure is based on shared random effects or shared 
  # coefficients then construct a matrix with the estimated b parameters
  # from the separate glmod (for the id_var grouping factor only). Note this
  # matrix is not passed to standata, but just used for autoscaling the 
  # priors for association parameters.
  sel_shared <- grep("^shared", rownames(assoc))
  if (any(unlist(assoc[sel_shared]))) {
    # flist for long submodel
    flist_tmp <- lme4::getME(y_mod$mod, "flist")
    # which grouping factor is id_var
    Gp_sel <- which(names(flist_tmp) == id_var) 
    # grouping factor indices
    Gp <- lme4::getME(y_mod$mod, "Gp")  
    b_beg <- Gp[[Gp_sel]] + 1
    b_end <- Gp[[Gp_sel + 1]]
    # b vector for grouping factor = id_var
    b_vec <- lme4::getME(y_mod$mod, "b")[b_beg:b_end]
    # convert to Npat * n_re matrix
    b_mat <- matrix(b_vec, nrow = length(levels(flist_tmp[[Gp_sel]])), byrow = TRUE)
  } else b_mat <- NULL
  
  parts$b_mat <- b_mat
  return(parts)
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
get_basic_grp_info <- function(cnms, flist, id_var) {
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
                               ok_grp_assocs = c("sum", "mean", "min", "max")) {
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
    if (!grp_assoc %in% ok_grp_assocs)
      stop2("'grp_assoc' must be one of: ", paste(ok_grp_assocs, collapse = ", "))
    
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
generate_init_function <- function(e_mod_stuff, standata) {
  
  # Initial values for intercepts, coefficients and aux parameters
  e_beta    <- e_mod_stuff$mod$coef
  e_aux     <- if (standata$basehaz_type == 1L) runif(1, 0.5, 3) else rep(0, standata$basehaz_df)
  e_z_beta      <- standardise_coef(e_beta, standata$e_prior_mean, standata$e_prior_scale) 
  e_aux_unscaled<- standardise_coef(e_aux, standata$e_prior_mean_for_aux, standata$e_prior_scale_for_aux)

  # Function to generate model based initial values
  model_based_inits <- rm_null(list(
    e_z_beta       = array_else_double(e_z_beta),
    e_aux_unscaled = array_else_double(e_aux_unscaled),
    e_gamma  = array_else_double(rep(0, standata$e_has_intercept))))
  
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
    if (is_jm && length(standata$basehaz_X)) "e_aux",
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
  unlist(Sigma_nms)
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
    len <- if (is_glmod) length(mod_stuff$Y$Y) else length(mod_stuff$eventtime)
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

