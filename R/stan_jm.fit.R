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

# Internal model fitting function for models estimated using 
# \code{stan_mvmer} or \code{stan_jm}.
# 
# See \code{stan_jm} for a description of the arguments to the 
# \code{stan_jm.fit} function call.
#
stan_jm.fit <- function(formulaLong          = NULL, 
                        dataLong             = NULL, 
                        formulaEvent         = NULL, 
                        dataEvent            = NULL, 
                        time_var, 
                        id_var,  
                        family               = gaussian, 
                        assoc                = "etavalue", 
                        lag_assoc            = 0, 
                        grp_assoc, 
                        epsilon              = 1E-5, 
                        basehaz              = c("bs", "weibull", "piecewise"), 
                        basehaz_ops, 
                        qnodes               = 15, 
                        init                 = "prefit", 
                        weights,					          
                        priorLong            = normal(),
                        priorLong_intercept  = normal(), 
                        priorLong_aux        = cauchy(0, 5), 
                        priorEvent           = normal(), 
                        priorEvent_intercept = normal(), 
                        priorEvent_aux       = cauchy(),
                        priorEvent_assoc     = normal(), 
                        prior_covariance     = lkj(), 
                        prior_PD             = FALSE, 
                        algorithm            = c("sampling", "meanfield", "fullrank"), 
                        adapt_delta          = NULL, 
                        max_treedepth        = 10L, 
                        QR                   = FALSE, 
                        sparse               = FALSE, 
                        ...) {
  
  #-----------------------------
  # Pre-processing of arguments
  #-----------------------------  
  
  if (!requireNamespace("survival"))
    stop("the 'survival' package must be installed to use this function.")
  
  # Set seed if specified
  dots <- list(...)
  if ("seed" %in% names(dots))
    set.seed(dots$seed)
  
  algorithm <- match.arg(algorithm)
  basehaz   <- match.arg(basehaz)
  
  if (missing(basehaz_ops)) basehaz_ops <- NULL
  if (missing(weights))     weights     <- NULL
  if (missing(id_var))      id_var      <- NULL
  if (missing(time_var))    time_var    <- NULL
  if (missing(grp_assoc))   grp_assoc   <- NULL

  if (!is.null(weights)) 
    stop("'weights' are not yet implemented.")
  if (QR)               
    stop("'QR' decomposition is not yet implemented.")
  if (sparse)
    stop("'sparse' option is not yet implemented.")

  # Error if args not supplied together
  supplied_together(formulaLong,  dataLong,  error = TRUE)
  supplied_together(formulaEvent, dataEvent, error = TRUE)
  
  # Determine whether a joint longitudinal-survival model was specified
  is_jm <- supplied_together(formulaLong, formulaEvent)
  stub  <- ifelse(is_jm, "Long", "y")

  if (is_jm && is.null(time_var))
    stop("'time_var' must be specified.")

  # Formula
  formulaLong <- validate_arg(formulaLong, "formula"); M <- length(formulaLong)
  
  # Data
  dataLong <- validate_arg(dataLong, "data.frame", validate_length = M)  
  if (is_jm)
    dataEvent <- as.data.frame(dataEvent)
  
  # Family
  ok_classes <- c("function", 
                  "family", 
                  "character")
  ok_families <- c("binomial", 
                   "gaussian", 
                   "Gamma", 
                   "inverse.gaussian", 
                   "poisson", 
                   "neg_binomial_2")
  family <- validate_arg(family, ok_classes, validate_length = M)
  family <- lapply(family, validate_famlink, ok_families)
  family <- lapply(family, append_mvmer_famlink)
  
  # Observation weights
  has_weights <- !is.null(weights)
  
  # Priors
  priorLong           <- broadcast_prior(priorLong,           M)
  priorLong_aux       <- broadcast_prior(priorLong_aux,       M)
  priorLong_intercept <- broadcast_prior(priorLong_intercept, M)
  
  # Combine meta information
  meta_stuff <- nlist(is_jm,
                      id_var,
                      time_var,
                      epsilon,
                      auc_qnodes = 15L)
  
  #--------------------------
  # Longitudinal submodel(s)
  #--------------------------
  
  # info for separate longitudinal submodels
  y_mod <- xapply(formulaLong, dataLong, family, FUN = handle_y_mod)
  
  # construct single cnms list for all longitudinal submodels
  cnms <- get_common_cnms(y_mod, stub = stub)
  
  # construct single list with unique levels for each grouping factor
  flevels <- get_common_flevels(y_mod)
  
  # ensure id_var is a valid grouping factor in all submodels
  if (is_jm) {
    id_var  <- check_id_var (y_mod,   id_var)
    id_list <- check_id_list(y_mod,   id_var)
    weights <- check_weights(weights, id_var)
  }
  
  # observation weights
  y_weights <- lapply(y_mod, handle_weights, weights, id_var)

  #----------- Prior distributions -----------# 
  
  # valid prior distributions
  ok_dists <- nlist("normal", 
                    student_t = "t", 
                    "cauchy", 
                    "hs", 
                    "hs_plus", 
                    "laplace", 
                    "lasso") # disallow product normal
  ok_dists_int <- c(ok_dists[1:3])
  ok_dists_aux <- c(ok_dists[1:3], exponential = "exponential")
  ok_dists_cov <- c("decov", "lkj")
  
  y_vecs <- fetch(y_mod, "y", "y")     # used in autoscaling, response  vector
  x_mats <- fetch(y_mod, "x", "xtemp") # used in autoscaling, predictor matrix
  
  # note: *_user_prior_*_stuff objects are stored unchanged for constructing 
  # prior_summary, while *_prior_*_stuff objects are autoscaled
  
  # priors for longitudinal submodels
  y_links <- fetch(y_mod, "family", "link")
  y_nvars <- fetch(y_mod, "x", "K")
  y_user_prior_stuff <- y_prior_stuff <- 
    xapply(priorLong, 
           nvars = y_nvars, 
           link  = y_links,
           FUN   = handle_glm_prior, 
           args  = list(default_scale = 2.5, ok_dists = ok_dists))
  
  y_user_prior_intercept_stuff <- y_prior_intercept_stuff <- 
    xapply(priorLong_intercept,
           nvars = 1,
           link  = y_links, 
           FUN   = handle_glm_prior,
           args  = list(default_scale = 10, ok_dists = ok_dists_int))
  
  y_user_prior_aux_stuff <- y_prior_aux_stuff <- 
    xapply(priorLong_aux,
           nvars = 1,
           link  = NULL,
           FUN   = handle_glm_prior, 
           args  = list(default_scale = 5, ok_dists = ok_dists_aux))  
  
  b_user_prior_stuff <- b_prior_stuff <- 
    handle_cov_prior(prior_covariance, cnms = cnms, ok_dists = ok_dists_cov)
  
  # autoscaling of priors
  y_prior_stuff <- 
    xapply(y_prior_stuff, 
           response   = y_vecs, 
           predictors = x_mats, 
           family     = family, 
           FUN        = autoscale_prior)
  y_prior_intercept_stuff <- 
    xapply(y_prior_intercept_stuff, 
           response = y_vecs,
           family   = family, 
           FUN      = autoscale_prior)
  y_prior_aux_stuff <- 
    xapply(y_prior_aux_stuff, 
           response = y_vecs,
           family   = family, 
           FUN      = autoscale_prior)
  
  # autoscale priors for sd of random effects (for lkj prior only)
  if (b_prior_stuff$prior_dist_name == "lkj") { 
    b_prior_stuff <- 
      split_cov_prior(b_prior_stuff, cnms = cnms,
                      submodel_cnms = fetch(y_mod, "z", "group_cnms"))
    b_prior_stuff <- 
      xapply(names(cnms), FUN = function(nm) {
        z_mats <- fetch(y_mod, "z", "z", nm)
        xapply(b_prior_stuff[[nm]], 
               response   = y_vecs, 
               predictors = z_mats, 
               family     = family, 
               FUN        = autoscale_prior)})
  }
  
  #----------- Data for export to Stan -----------# 
  
  standata <- list(
    M           = ai(M), 
    has_weights = ai(!all(lapply(weights, is.null))),
    family      = fetch_array(y_mod, "family", "mvmer_family"),
    link        = fetch_array(y_mod, "family", "mvmer_link"),
    weights     = aa(numeric(0)), # not yet implemented
    prior_PD    = ai(prior_PD)
  )  
  
  # Dimensions
  standata$has_aux <- 
    fetch_array(y_mod, "has_aux", pad_length = 3)
  standata$resp_type <- 
    fetch_array(y_mod, "y", "resp_type", pad_length = 3)
  standata$intercept_type <- 
    fetch_array(y_mod, "intercept_type", "number", pad_length = 3)
  standata$yNobs <- 
    fetch_array(y_mod, "x", "N", pad_length = 3)
  standata$yNeta <- 
    fetch_array(y_mod, "x", "N", pad_length = 3) # same as Nobs for stan_mvmer
  standata$yK <- 
    fetch_array(y_mod, "x", "K", pad_length = 3)
  
  # Response vectors
  Y_integer <- fetch(y_mod, "y", "integer")
  standata$yInt1 <- if (M > 0) Y_integer[[1]] else aa(integer(0))  
  standata$yInt2 <- if (M > 1) Y_integer[[2]] else aa(integer(0))  
  standata$yInt3 <- if (M > 2) Y_integer[[3]] else aa(integer(0)) 
  
  Y_real <- fetch(y_mod, "y", "real")
  standata$yReal1 <- if (M > 0) Y_real[[1]] else aa(double(0)) 
  standata$yReal2 <- if (M > 1) Y_real[[2]] else aa(double(0)) 
  standata$yReal3 <- if (M > 2) Y_real[[3]] else aa(double(0)) 
  
  # Population level design matrices
  X <- fetch(y_mod, "x", "xtemp")
  standata$yX1 <- if (M > 0) X[[1]] else matrix(0,0,0)
  standata$yX2 <- if (M > 1) X[[2]] else matrix(0,0,0)
  standata$yX3 <- if (M > 2) X[[3]] else matrix(0,0,0)
  
  X_bar <- fetch(y_mod, "x", "x_bar")
  standata$yXbar1 <- if (M > 0) aa(X_bar[[1]]) else aa(double(0))
  standata$yXbar2 <- if (M > 1) aa(X_bar[[2]]) else aa(double(0))
  standata$yXbar3 <- if (M > 2) aa(X_bar[[3]]) else aa(double(0))
  
  # Data for group specific terms - group factor 1
  b1_varname <- names(cnms)[[1L]] # name of group factor 1
  b1_nvars <- fetch_(y_mod, "z", "nvars", b1_varname, 
                     null_to_zero = TRUE, pad_length = 3)
  b1_ngrps <- fetch_(y_mod, "z", "ngrps", b1_varname)
  if (!n_distinct(b1_ngrps) == 1L)
    stop("The number of groups for the grouping factor '", 
         b1_varname, "' should be the same in all submodels.")
  
  standata$bN1 <- b1_ngrps[[1L]] + 1L # add padding for _NEW_ group
  standata$bK1 <- sum(b1_nvars)
  standata$bK1_len <- aa(b1_nvars)
  standata$bK1_idx <- get_idx_array(b1_nvars)
  
  Z1 <- fetch(y_mod, "z", "z", b1_varname)
  Z1 <- lapply(Z1, transpose)
  Z1 <- lapply(Z1, convert_null, "matrix")
  standata$y1_Z1 <- if (M > 0) Z1[[1L]] else matrix(0,0,0)
  standata$y2_Z1 <- if (M > 1) Z1[[2L]] else matrix(0,0,0)
  standata$y3_Z1 <- if (M > 2) Z1[[3L]] else matrix(0,0,0)
  
  Z1_id <- fetch(y_mod, "z", "group_list", b1_varname)
  Z1_id <- lapply(Z1_id, groups)
  Z1_id <- lapply(Z1_id, convert_null, "arrayinteger")
  standata$y1_Z1_id <- if (M > 0) Z1_id[[1L]] else aa(integer(0))
  standata$y2_Z1_id <- if (M > 1) Z1_id[[2L]] else aa(integer(0))
  standata$y3_Z1_id <- if (M > 2) Z1_id[[3L]] else aa(integer(0))
  
  # Data for group specific terms - group factor 2
  if (length(cnms) > 1L) {
    # model has a second grouping factor
    b2_varname <- names(cnms)[[2L]] # name of group factor 2
    b2_nvars <- fetch_(y_mod, "z", "nvars", b2_varname, 
                       null_to_zero = TRUE, pad_length = 3)
    b2_ngrps <- fetch_(y_mod, "z", "ngrps", b2_varname)
    if (!n_distinct(b2_ngrps) == 1L)
      stop("The number of groups for the grouping factor '", 
           b2_varname, "' should be the same in all submodels.")
    standata$bN2 <- b2_ngrps[[1L]] + 1L # add padding for _NEW_ group
    standata$bK2 <- sum(b2_nvars)
    standata$bK2_len <- aa(b2_nvars)
    standata$bK2_idx <- get_idx_array(b2_nvars)
    
    Z2 <- fetch(y_mod, "z", "z", b2_varname)
    Z2 <- lapply(Z2, transpose)
    Z2 <- lapply(Z2, convert_null, "matrix")
    standata$y1_Z2 <- if (M > 0) Z2[[1L]] else matrix(0,0,0)
    standata$y2_Z2 <- if (M > 1) Z2[[2L]] else matrix(0,0,0)
    standata$y3_Z2 <- if (M > 2) Z2[[3L]] else matrix(0,0,0)
    
    Z2_id <- fetch(y_mod, "z", "group_list", b2_varname)
    Z2_id <- lapply(Z2_id, groups)
    Z2_id <- lapply(Z2_id, convert_null, "arrayinteger")
    standata$y1_Z2_id <- if (M > 0) Z2_id[[1L]] else aa(integer(0))
    standata$y2_Z2_id <- if (M > 1) Z2_id[[2L]] else aa(integer(0))
    standata$y3_Z2_id <- if (M > 2) Z2_id[[3L]] else aa(integer(0))
    
  } else {
    # no second grouping factor
    standata$bN2 <- 0L
    standata$bK2 <- 0L
    standata$bK2_len <- aa(rep(0,3L))
    standata$bK2_idx <- get_idx_array(rep(0,3L))
    standata$y1_Z2 <- matrix(0,0,0)
    standata$y2_Z2 <- matrix(0,0,0)
    standata$y3_Z2 <- matrix(0,0,0)
    standata$y1_Z2_id <- aa(integer(0))
    standata$y2_Z2_id <- aa(integer(0))
    standata$y3_Z2_id <- aa(integer(0))
  }
  
  # Priors
  standata$y_prior_dist_for_intercept <- 
    fetch_array(y_prior_intercept_stuff, "prior_dist")  
  standata$y_prior_mean_for_intercept <- 
    fetch_array(y_prior_intercept_stuff, "prior_mean")
  standata$y_prior_scale_for_intercept <- 
    fetch_array(y_prior_intercept_stuff, "prior_scale")
  standata$y_prior_df_for_intercept <- 
    fetch_array(y_prior_intercept_stuff, "prior_df")
  
  standata$y_prior_dist_for_aux <-
    fetch_array(y_prior_aux_stuff, "prior_dist")
  standata$y_prior_mean_for_aux <- 
    fetch_array(y_prior_aux_stuff, "prior_mean")
  standata$y_prior_scale_for_aux <- 
    fetch_array(y_prior_aux_stuff, "prior_scale")
  standata$y_prior_df_for_aux <- 
    fetch_array(y_prior_aux_stuff, "prior_df")
  
  standata$y_prior_dist <- 
    fetch_array(y_prior_stuff, "prior_dist", pad_length = 3)
  
  prior_mean <- fetch(y_prior_stuff, "prior_mean")
  standata$y_prior_mean1 <- if (M > 0) prior_mean[[1]] else aa(double(0))
  standata$y_prior_mean2 <- if (M > 1) prior_mean[[2]] else aa(double(0))
  standata$y_prior_mean3 <- if (M > 2) prior_mean[[3]] else aa(double(0))
  
  prior_scale <- fetch(y_prior_stuff, "prior_scale")
  standata$y_prior_scale1 <- if (M > 0) aa(prior_scale[[1]]) else aa(double(0))
  standata$y_prior_scale2 <- if (M > 1) aa(prior_scale[[2]]) else aa(double(0))
  standata$y_prior_scale3 <- if (M > 2) aa(prior_scale[[3]]) else aa(double(0))
  
  prior_df <- fetch(y_prior_stuff, "prior_df")
  standata$y_prior_df1 <- if (M > 0) prior_df[[1]] else aa(double(0))
  standata$y_prior_df2 <- if (M > 1) prior_df[[2]] else aa(double(0))
  standata$y_prior_df3 <- if (M > 2) prior_df[[3]] else aa(double(0))
  
  # hs priors only
  standata$y_global_prior_scale <- fetch_array(y_prior_stuff, "global_prior_scale") 
  standata$y_global_prior_df    <- fetch_array(y_prior_stuff, "global_prior_df")
  standata$y_slab_df            <- fetch_array(y_prior_stuff, "slab_df")
  standata$y_slab_scale         <- fetch_array(y_prior_stuff, "slab_scale")
  
  # Priors for group specific terms
  standata$t <- length(cnms)
  standata$p <- aa(sapply(cnms, length))
  standata$l <- aa(
    sapply(names(cnms), FUN = function(nm) {
      ngrps <- unique(fetch_(y_mod, "z", "ngrps", nm))
      ngrps + 1L # add padding for _NEW_ group
  }))
  standata$q <- sum(standata$p * standata$l)
  
  if (prior_covariance$dist == "decov") {
    
    # data for decov prior
    standata$prior_dist_for_cov     <- b_prior_stuff$prior_dist
    standata$b_prior_shape          <- b_prior_stuff$prior_shape
    standata$b_prior_scale          <- b_prior_stuff$prior_scale
    standata$b_prior_concentration  <- b_prior_stuff$prior_concentration
    standata$b_prior_regularization <- b_prior_stuff$prior_regularization
    standata$len_concentration      <- length(standata$b_prior_concentration)
    standata$len_regularization     <- length(standata$b_prior_regularization)
    standata$len_theta_L            <- sum(choose(standata$p, 2), standata$p)
    
    # pass empty lkj data
    standata$b1_prior_scale          <- aa(rep(0L, standata$bK1))
    standata$b2_prior_scale          <- aa(rep(0L, standata$bK2))
    standata$b1_prior_df             <- aa(rep(0L, standata$bK1))
    standata$b2_prior_df             <- aa(rep(0L, standata$bK2))
    standata$b1_prior_regularization <- 1.0
    standata$b2_prior_regularization <- 1.0   
    
  } else if (prior_covariance$dist == "lkj") {
    
    # data for lkj prior
    b1_prior_stuff          <- b_prior_stuff[[b1_varname]]
    b1_prior_dist           <- fetch_     (b1_prior_stuff, "prior_dist")
    b1_prior_scale          <- fetch_array(b1_prior_stuff, "prior_scale")
    b1_prior_df             <- fetch_array(b1_prior_stuff, "prior_df")
    b1_prior_regularization <- fetch_     (b1_prior_stuff, "prior_regularization")
    
    if (n_distinct(b1_prior_dist) > 1L)
      stop2("Bug found: covariance prior should be the same for all submodels.")
    if (n_distinct(b1_prior_regularization) > 1L)
      stop2("Bug found: prior_regularization should be the same for all submodels.")

    standata$prior_dist_for_cov <- unique(b1_prior_dist)
    standata$b1_prior_scale     <- b1_prior_scale
    standata$b1_prior_df        <- b1_prior_df
    standata$b1_prior_regularization <- if (length(b1_prior_regularization))
      unique(b1_prior_regularization) else 1.0
    
    if (standata$bK2 > 0) {
      # model has a second grouping factor
      b2_prior_stuff          <- b_prior_stuff[[b2_varname]]
      b2_prior_scale          <- fetch_array(b2_prior_stuff, "prior_scale")
      b2_prior_df             <- fetch_array(b2_prior_stuff, "prior_df")
      b2_prior_regularization <- fetch_     (b2_prior_stuff, "prior_regularization")
      standata$b2_prior_scale <- b2_prior_scale
      standata$b2_prior_df    <- b2_prior_df
      standata$b2_prior_regularization <- unique(b2_prior_regularization)
    } else {
      # model does not have a second grouping factor
      standata$b2_prior_scale          <- aa(double(0))
      standata$b2_prior_df             <- aa(double(0))
      standata$b2_prior_regularization <- 1.0
    }
    
    # pass empty decov data
    standata$len_theta_L            <- 0L
    standata$len_concentration      <- 0L
    standata$len_regularization     <- 0L
    standata$b_prior_shape          <- aa(rep(0L, standata$t))
    standata$b_prior_scale          <- aa(rep(0L, standata$t))
    standata$b_prior_concentration  <- aa(rep(0L, standata$len_concentration))
    standata$b_prior_regularization <- aa(rep(0L, standata$len_regularization))   
  }
  
  
  # Names for longitudinal submodel parameters
  y_nms_beta <- uapply(y_mod, get_beta_name)
  y_nms_int  <- uapply(y_mod, get_int_name)
  y_nms_aux  <- uapply(y_mod, get_aux_name)
  y_nms_ppd  <- uapply(y_mod, get_ppd_name)
  
  # Names for group specific coefficients ("b pars")
  b_nms <- get_ranef_nms(cnms, flevels)

  # Names for Sigma matrix
  y_nms_sigma <- get_Sigma_nms(cnms)
  
  #----------------
  # Event submodel
  #----------------

  if (is_jm) { # begin jm block

    # fit separate event submodel
    e_mod <- handle_e_mod(formula   = formulaEvent, 
                          data      = dataEvent, 
                          qnodes    = qnodes,
                          id_var    = id_var, 
                          y_id_list = id_list)
    
    basehaz <- e_mod$basehaz  
    
    # observation weights
    e_weights <- handle_weights(e_mod, weights, id_var)
  
    # check longitudinal observation times are not later than the event time
    lapply(dataLong, 
           FUN      = validate_observation_times,  
           exittime = e_mod$exittime, 
           id_var   = id_var, 
           time_var = time_var)
    
    #----------- Prior distributions -----------# 
    
    # valid prior distributions
    ok_e_aux_dists <- ok_dists[1:3]
  
    # note: *_user_prior_*_stuff objects are stored unchanged for constructing 
    # prior_summary, while *_prior_*_stuff objects are autoscaled
    
    # priors for event submodel
    e_user_prior_stuff <- e_prior_stuff <- 
      handle_glm_prior(priorEvent, 
                       nvars         = e_mod$K, 
                       default_scale = 2.5, 
                       link          = NULL, 
                       ok_dists      = ok_dists)  
    
    e_user_prior_intercept_stuff <- e_prior_intercept_stuff <- 
      handle_glm_prior(priorEvent_intercept, 
                       nvars         = 1, 
                       default_scale = 20,
                       link          = NULL, 
                       ok_dists      = ok_intercept_dists)  
    
    e_user_prior_aux_stuff <- e_prior_aux_stuff <-
      handle_glm_prior(priorEvent_aux, 
                       nvars         = basehaz$nvars,
                       default_scale = get_default_aux_scale(basehaz),
                       link          = NULL, 
                       ok_dists      = ok_e_aux_dists)

    # stop null priors if prior_PD is TRUE
    if (prior_PD) {
      if (is.null(prior))
        stop("'priorEvent' cannot be NULL if 'prior_PD' is TRUE")
      if (is.null(prior_intercept) && e_mod$has_intercept)
        stop("'priorEvent_intercept' cannot be NULL if 'prior_PD' is TRUE")
      if (is.null(prior_aux))
        stop("'priorEvent_aux' cannot be NULL if 'prior_PD' is TRUE")          
    }
    
    # autoscaling of priors
    e_prior_stuff           <- autoscale_prior(e_prior_stuff, predictors = e_mod$x)
    e_prior_intercept_stuff <- autoscale_prior(e_prior_intercept_stuff)
    e_prior_aux_stuff       <- autoscale_prior(e_prior_aux_stuff)
 
    #----------- Data for export to Stan -----------# 
    
    # dimensions
    standata$e_K              <- ai(e_mod$K)
    standata$nevent           <- ai(e_mod$nevent)
    standata$nrcens           <- ai(e_mod$nrcens)
    standata$nlcens           <- ai(e_mod$nlcens)
    standata$nicens           <- ai(e_mod$nicens)
    #standata$Npat             <- ai(e_mod$Npat)
    #standata$Nevents          <- ai(e_mod$Nevents)
    standata$qnodes           <- ai(e_mod$qnodes)
    standata$qrows            <- ai(e_mod$qrows)
    standata$qicens           <- ai(e_mod$qicens)
    standata$qdelayed         <- ai(e_mod$qdelayed)
    standata$e_has_intercept  <- ai(e_mod$has_intercept)
    
    # design matrices
    standata$epts             <- aa(e_mod$epts)
    standata$qpts             <- aa(e_mod$qpts)
    standata$qpts_upper       <- aa(e_mod$qpts_upper)
    standata$qpts_delayed     <- aa(e_mod$qpts_delayed)
    standata$qwts             <- aa(e_mod$qwts)
    standata$qwts_upper       <- aa(e_mod$qwts_upper)
    standata$qwts_delayed     <- aa(e_mod$qwts_delayed)
    standata$e_x_qpts         <- e_mod$x_qpts
    standata$e_x_qpts_upper   <- e_mod$x_qpts_upper
    standata$e_x_qpts_delayed <- e_mod$x_qpts_delayed
    standata$e_xbar           <- aa(e_mod$Xbar)
    standata$e_weights        <- aa(e_weights)
    standata$e_weights_rep    <- aa(rep(e_weights, times = qnodes))
    #standata$Npat_times_qnodes<- ai(e_mod$Npat * qnodes)
    
    # baseline hazard
    standata$basehaz_type       <- ai(basehaz$type)
    standata$basehaz_nvars      <- ai(basehaz$nvars)
    standata$basis_events       <- e_mod$basis_events
    standata$basis_qpts         <- e_mod$basis_qpts
    standata$basis_qpts_upper   <- e_mod$basis_qpts_upper
    standata$basis_qpts_delayed <- e_mod$basis_qpts_delayed
    standata$norm_const         <- e_mod$norm_const
    
    # priors
    standata$e_prior_dist              <- e_prior_stuff$prior_dist
    standata$e_prior_dist_for_intercept<- e_prior_intercept_stuff$prior_dist
    standata$e_prior_dist_for_aux      <- e_prior_aux_stuff$prior_dist
    
    # hyperparameters for event submodel priors
    standata$e_prior_mean               <- e_prior_stuff$prior_mean
    standata$e_prior_scale              <- e_prior_stuff$prior_scale
    standata$e_prior_df                 <- e_prior_stuff$prior_df
    standata$e_prior_mean_for_intercept <- c(e_prior_intercept_stuff$prior_mean)
    standata$e_prior_scale_for_intercept<- c(e_prior_intercept_stuff$prior_scale)
    standata$e_prior_df_for_intercept   <- c(e_prior_intercept_stuff$prior_df)
    standata$e_prior_scale_for_aux      <- e_prior_aux_stuff$prior_scale
    standata$e_prior_df_for_aux         <- e_prior_aux_stuff$prior_df
    standata$e_global_prior_scale       <- e_prior_stuff$global_prior_scale
    standata$e_global_prior_df          <- e_prior_stuff$global_prior_df
    standata$e_slab_df                  <- e_prior_stuff$slab_df
    standata$e_slab_scale               <- e_prior_stuff$slab_scale
    
    #-----------------------
    # Association structure
    #-----------------------
      
    # define the valid types of association structure
    # !! if order is changed here, then must also change standata$has_assoc !!
    ok_assoc <- c("null", 
                  "etavalue",
                  "etaslope", 
                  "etaauc", 
                  "muvalue", 
                  "muslope", 
                  "muauc", 
                  "shared_b", 
                  "shared_coef")
    
    ok_assoc_data         <- ok_assoc[c(2,3,5,6)] # ok to interact with covariates
    ok_assoc_interactions <- ok_assoc[c(2,5)]     # ok to interact across biomarkers
    ok_assoc_interval     <- ok_assoc[c(1:3,5,6)] # ok to use with interval censoring

    ok_assoc_list <- nlist(ok_assoc,
                           ok_assoc_data,
                           ok_assoc_interactions,
                           ok_assoc_interval)
    
    # check lag time is valid
    lag_assoc <- validate_lag_assoc(lag_assoc, M) 
    
    # return an array summarising the association structure
    assoc <- handle_assoc(assoc, lag_assoc, ok_assoc_list, meta_stuff)

    # for each submodel, identify grouping factors clustered within 'id_var' 
    # (i.e. lower level clustering)
    ok_assoc_grp <- c("sum",                      # valid inputs to 'grp_assoc' arg                  
                      "mean", 
                      "min", 
                      "max") 
    
    ok_assoc_with_grp <- c("etavalue",            # ok to use with non-NULL grp_assoc
                           "etavalue_data", 
                           "etaslope", 
                           "etaslope_data", 
                           "muvalue", 
                           "muvalue_data")  
    
    grp_basic <- xapply(FUN = get_basic_grp_info, y_mod = y_mod, id_var = id_var)
    grp_stuff <- xapply(FUN = get_extra_grp_info,
                        basic_info = grp_basic, 
                        flist = fetch(y_mod, "z", "group_list"),
                        args = nlist(id_var, grp_assoc, ok_assoc_grp))
    has_grp <- fetch_(grp_stuff, "has_grp")

    if (not.null(grp_assoc) && !any(has_grp))
      stop2("'grp_assoc' can only be specified when there is a grouping ",
            "factor clustered within patients.")  
    
    if (any(has_grp))
      validate_assoc_with_grp(grp_stuff, assoc, ok_assoc_with_grp)
    
    # design matrices for longitudinal submodel at the quadrature points
    auc_qnodes <- meta_stuff$auc_qnodes <- 15L
    a_mod <- xapply(FUN       = handle_assocmod, 
                    data      = dataLong, 
                    assoc     = apply(assoc, 2L, c), # converts array to list
                    y_mod     = y_mod,
                    grp_stuff = grp_stuff, 
                    args      = nlist(e_mod, meta_stuff))
        
    # number of association parameters
    a_K <- get_num_assoc_pars(assoc, a_mod)
    
    # use a stan_mvmer variational bayes model fit for:
    # - obtaining initial values for joint model parameters
    # - obtaining appropriate scaling for priors on association parameters
    dropargs  <- c("chains", "cores", "iter", "refresh", "test_grad", "control")
    init_dots <- list(...); for (i in dropargs) init_dots[[i]] <- NULL
    init_mod  <- stanmodels$mvmer
    init_data <- standata
    init_pars <- pars_to_monitor(standata, is_jm = FALSE)
    init_args <- nlist(object = init_mod, 
                       data   = init_data, 
                       pars   = init_pars, 
                       algorithm = "meanfield")
    init_args[names(init_dots)] <- init_dots
    utils::capture.output(init_fit <- do.call(rstan::vb, init_args))
    init_nms_all <- c(y_nms_int, 
                      y_nms_beta,
                      b_nms,
                      y_nms_aux, 
                      y_nms_sigma,
                      y_nms_ppd,
                      "log-posterior")
    init_fit  <- replace_stanfit_nms(init_fit, init_nms_all)
    init_mat  <- t(colMeans(am(init_fit))) # posterior means
    init_nms  <- collect_nms(colnames(init_mat), M, stub = "Long")
    init_beta <- lapply(1:M, function(m) init_mat[, init_nms$y[[m]]])
    init_b    <- lapply(1:M, function(m) {
      # can drop _NEW_ groups since they are not required for generating
      # the assoc_terms that are used in scaling the priors for 
      # the association parameters (ie. the Zt matrix returned by the 
      # function 'make_assoc_parts_for_stan' will not be padded).
      b <- init_mat[, init_nms$y_b[[m]]]
      b[!grepl("_NEW_", names(b), fixed = TRUE)]
    })
    
    if (is.character(init) && (init =="prefit"))
      init <- get_prefit_inits(init_fit, standata)
    
    #----------- Prior distributions -----------# 

    # Priors for association parameters
    e_user_prior_assoc_stuff <- e_prior_assoc_stuff <- 
      handle_glm_prior(priorEvent_assoc, 
                       nvars         = a_K, 
                       default_scale = 2.5,
                       link          = NULL, 
                       ok_dists      = ok_dists)  
    
    # Autoscaling of priors
    if (a_K) {
      e_prior_assoc_stuff <- autoscale_prior(e_prior_assoc_stuff, 
                                             family = family, 
                                             assoc  = assoc, 
                                             parts  = a_mod,
                                             beta   = init_beta, 
                                             b      = init_b)
    }
    
    #----------- Data for export to Stan -----------# 
 
    # dimensions   
    standata$assoc <- ai(a_K > 0L) # any association structure
    standata$a_K   <- ai(a_K)      # num association parameters
    
    # indicators for components required to build association terms
    standata$assoc_uses <- make_assoc_component_flags(assoc)
    
    # indicators for each possible type of association structure
    # !! must be careful with corresponding use of indexing in stan code !!
    # !! this is determined by the row ordering of the 'assoc' array     !!
    #   1  = ev
    #   2  = es
    #   3  = ea
    #   4  = mv
    #   5  = ms
    #   6  = ma
    #   7  = shared_b
    #   8  = shared_coef
    #   9  = ev_data
    #   10 = es_data
    #   11 = mv_data
    #   12 = ms_data
    #   13 = evev
    #   14 = evmv
    #   15 = mvev
    #   16 = mvmv
    standata$has_assoc <- make_assoc_type_flags(assoc)

    # data for calculating eta, slope, auc in GK quadrature
    standata <- standata_add_assoc_grp   (standata, a_mod, grp_stuff)
    standata <- standata_add_assoc_xz    (standata, a_mod, assoc, meta_stuff)
    standata <- standata_add_assoc_auc   (standata, a_mod, e_mod, meta_stuff)
    standata <- standata_add_assoc_extras(standata, a_mod, assoc, meta_stuff)
    
    # hyperparameters for assoc parameter priors
    standata$a_prior_dist         <- e_prior_assoc_stuff$prior_dist 
    standata$a_prior_mean         <- e_prior_assoc_stuff$prior_mean
    standata$a_prior_scale        <- aa(e_prior_assoc_stuff$prior_scale)
    standata$a_prior_df           <- e_prior_assoc_stuff$prior_df
    standata$a_global_prior_scale <- e_prior_assoc_stuff$global_prior_scale
    standata$a_global_prior_df    <- e_prior_assoc_stuff$global_prior_df
    standata$a_slab_df            <- e_prior_assoc_stuff$slab_df
    standata$a_slab_scale         <- e_prior_assoc_stuff$slab_scale
    
    # centering for association terms
    standata$a_xbar <- if (a_K) e_prior_assoc_stuff$a_xbar else numeric(0)    

  } # end jm block
  
  #---------------
  # Prior summary
  #---------------
  
  prior_info <- summarize_jm_prior(
    user_priorLong                      = y_user_prior_stuff,
    user_priorLong_intercept            = y_user_prior_intercept_stuff,
    user_priorLong_aux                  = y_user_prior_aux_stuff,
    user_prior_covariance               = prior_covariance,
    y_has_intercept                     = fetch_(y_mod, "x", "has_intercept"),
    y_has_predictors                    = fetch_(y_mod, "x", "K") > 0,
    adjusted_priorLong_scale            = fetch(y_prior_stuff, "prior_scale"),
    adjusted_priorLong_intercept_scale  = fetch(y_prior_intercept_stuff, "prior_scale"),
    adjusted_priorLong_aux_scale        = fetch(y_prior_aux_stuff, "prior_scale"),
    family                              = family, 
    if (is_jm) basehaz                  = basehaz,
    if (is_jm) user_priorEvent          = e_user_prior_stuff,
    if (is_jm) user_priorEvent_intercept= e_user_prior_intercept_stuff,
    if (is_jm) user_priorEvent_aux      = e_user_prior_aux_stuff,
    if (is_jm) user_priorEvent_assoc    = e_user_prior_assoc_stuff,
    if (is_jm) e_has_intercept          = standata$e_has_intercept,
    if (is_jm) e_has_predictors         = standata$e_K > 0,
    if (is_jm) has_assoc                = a_K > 0,
    if (is_jm) adjusted_priorEvent_scale           = e_prior_stuff$prior_scale,
    if (is_jm) adjusted_priorEvent_intercept_scale = e_prior_intercept_stuff$prior_scale,
    if (is_jm) adjusted_priorEvent_aux_scale       = e_prior_aux_stuff$prior_scale,
    if (is_jm) adjusted_priorEvent_assoc_scale     = e_prior_assoc_stuff$prior_scale)  
  
  #-----------
  # Fit model
  #-----------
  
  # obtain stan model code
  stanfit  <- if (is_jm) stanmodels$jm else stanmodels$mvmer
  
  # specify parameters for stan to monitor
  stanpars <- pars_to_monitor(standata, is_jm = is_jm)

  # report type of model to user
  txt1 <- if (M == 1) "uni"   else "multi"
  txt2 <- if (is_jm)  "joint" else "glmer"
  txt3 <- paste0("Fitting a ", txt1, "variate", txt2, "model.\n\n")
  txt4 <- "Please note the warmup may be much slower than later iterations!\n"
  
  # fit model using stan
  cat(txt3)
  if (algorithm == "sampling") { # mcmc
    cat(txt4)             
    args <- set_jm_sampling_args(
      object = stanfit,
      data   = standata,
      pars   = stanpars, 
      init   = init,
      cnms   = cnms,
      user_dots = list(...), 
      user_adapt_delta   = adapt_delta,
      user_max_treedepth = max_treedepth,
      show_messages = FALSE)
    stanfit <- do.call(rstan::sampling, args)
  } else { # meanfield or fullrank vb
    args <- nlist(
      object = stanfit,
      data   = standata,
      pars   = stanpars,
      algorithm
    )
    args[names(dots)] <- dots
    stanfit <- do.call(rstan::vb, args)
  }
  if (!isTRUE(check_stanfit(stanfit))) 
    return(standata)

  # Sigma values in stanmat
  if (prior_covariance$dist == "decov" && standata$len_theta_L)
    stanfit <- evaluate_Sigma(stanfit, cnms)
  
  if (is_jm) {
    e_nms_beta  <- get_beta_name (e_mod)
    e_nms_int   <- get_int_name  (e_mod, basehaz)
    e_nms_aux   <- get_aux_name  (e_mod, basehaz)
    e_nms_assoc <- get_assoc_name(a_mod, assoc)
  } else {
    e_nms_beta  <- NULL
    e_nms_int   <- NULL
    e_nms_aux   <- NULL
    e_nms_assoc <- NULL
  }
  
  # define new parameter names
  nms_all <- c(y_nms_int,
               y_nms_beta,
               e_nms_int,
               e_nms_beta,
               e_nms_assoc,
               b_nms,
               y_nms_aux,
               e_nms_aux,
               y_nms_sigma,
               y_nms_ppd, 
               "log-posterior")
  
  # substitute new parameter names into 'stanfit' object
  stanfit <- replace_stanfit_nms(stanfit, nms_all)
  
  # combine elements to add to returned structure
  if (!is_jm) e_mod <- a_mod <- assoc <- basehaz <- id_var <- grp_stuff <- NULL
  args <- nlist(.Data = stanfit,
                y_mod, 
                e_mod,
                a_mod,
                assoc,
                basehaz,
                prior_info, 
                id_var,
                cnms, 
                flevels,
                grp_stuff)
  
  do.call("structure", remove_null(args, recursive = FALSE))
}

