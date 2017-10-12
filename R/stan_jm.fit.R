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
stan_jm.fit <- function(formulaLong = NULL, dataLong = NULL, formulaEvent = NULL, 
                        dataEvent = NULL, time_var, id_var,  family = gaussian, 
                        assoc = "etavalue", lag_assoc = 0, grp_assoc, dataAssoc, 
                        epsilon = 1E-5, basehaz = c("weibull", "bs", "piecewise"), 
                        basehaz_ops, qnodes = 15, init = "prefit", weights, ...,					          
                        priorLong = normal(), priorLong_intercept = normal(), 
                        priorLong_aux = cauchy(0, 5), priorEvent = normal(), 
                        priorEvent_intercept = normal(), priorEvent_aux = cauchy(),
                        priorAssoc = normal(), prior_covariance = lkj(), prior_PD = FALSE, 
                        algorithm = c("sampling", "meanfield", "fullrank"), 
                        adapt_delta = NULL, QR = FALSE, sparse = FALSE) {
  
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
  
  if (missing(offset))      offset      <- NULL 
  if (missing(basehaz_ops)) basehaz_ops <- NULL
  if (missing(weights))     weights     <- NULL
  if (missing(id_var))      id_var      <- NULL
  if (missing(time_var))    time_var    <- NULL
  if (missing(grp_assoc))   grp_assoc   <- NULL
  if (missing(dataAssoc))   dataAssoc   <- NULL 
  
  if (!is.null(weights)) 
    stop("'weights' are not yet implemented.")
  if (!is.null(dataAssoc))
    stop("'dataAssoc' argument not yet implemented.")
  if (QR)               
    stop("'QR' decomposition is not yet implemented.")
  if (sparse)
    stop("'sparse' option is not yet implemented.")

  # Error if args not supplied together
  supplied_together(formulaLong, dataLong, error = TRUE)
  supplied_together(formulaEvent, dataEvent, error = TRUE)
  
  # Determine whether a joint longitudinal-survival model was specified
  is_jm <- supplied_together(formulaLong, formulaEvent)
  stub <- if (is_jm) "Long" else "y"

  if (is_jm && is.null(time_var))
    stop("'time_var' must be specified.")
  if (is_jm && is.null(id_var))
    stop("'id_var' must be specified.")
    
  # Formula & Data
  yF <- validate_arg(formulaLong, "formula"); M <- length(yF)
  yD <- validate_arg(dataLong, "data.frame", validate_length = M)  
  yD <- xapply(yF, yD, FUN = get_all_vars)
  if (is_jm) {
    yD <- check_vars_are_included(yD, id_var, time_var)
    eF <- formulaEvent
    eD <- as.data.frame(dataEvent)
    eD <- check_vars_are_included(eD, id_var)
    eD <- get_all_vars(eF, eD, dataEvent[[id_var]])
    names(eD) <- c(names(eD), id_var)
  }
  
  # Family
  ok_classes <- c("function", "family", "character")
  ok_families <- c("binomial", "gaussian", "Gamma", 
                   "inverse.gaussian", "poisson", "neg_binomial_2")
  family <- validate_arg(family, ok_classes, validate_length = M)
  family <- lapply(family, validate_famlink, ok_families)
  family <- lapply(family, append_mvmer_famlink)
  
  # Observation weights
  has_weights <- !is.null(weights)
  
  # Priors
  prior <- broadcast_prior(prior, M)
  prior_intercept <- broadcast_prior(prior_intercept, M)
  prior_aux <- broadcast_prior(prior_aux, M)
  
  #--------------------------------
  # Data for longitudinal submodel
  #--------------------------------
  
  # Fit separate longitudinal submodels
  y_mod <- xapply(yF, yD, family, FUN = handle_y_mod)
  
  # Construct single cnms list for all longitudinal submodels
  y_cnms  <- fetch(y_mod, "Z", "group_cnms")
  y_flist <- fetch(y_mod, "Z", "group_list")
  cnms <- get_common_cnms(y_cnms, stub = stub)
  cnms_nms <- names(cnms)
  if (length(cnms_nms) > 2L)
    stop("A maximum of 2 grouping factors are allowed.")
  
  # Ensure id_var is a valid grouping factor in all submodels
  if (is_jm) {
    id_var <- validate_grouping_factor(y_cnms, id_var)
    id_list <- check_id_list(id_var, y_flist)
  }
    
  #-------------------------
  # Data for event submodel
  #-------------------------

  if (is_jm) { # begin jm block
    
    # Fit separate event submodel
    e_mod <- handle_e_mod(formula = eF, data = eD, qnodes = qnodes, 
                          id_var = id_var, y_id_list = id_list)
    
    # Baseline hazard
    ok_basehaz <- nlist("weibull", "bs", "piecewise")
    basehaz <- handle_basehaz(basehaz, basehaz_ops, ok_basehaz = ok_basehaz, 
                              eventtime = e_mod$eventtime, status = e_mod_$status)
    e_mod$has_intercept <- (basehaz$type_name == "weibull")
    
    # Observation weights
    if (has_weights) {
      y_weights <- lapply(y_mod_stuff, handle_weights, weights, id_var)
      e_weights <- handle_weights(e_mod_stuff, weights, id_var)
    }
    
  } # end jm block
  
  #--------------------------------
  # Data for association structure
  #--------------------------------
    
  if (is_jm) { # begin jm block
    
    # Handle association structure
    # !! If order is changed here, then must also change standata$has_assoc !!
    ok_assoc <- c("null", "etavalue","etaslope", "etaauc", "muvalue", 
                  "muslope", "muauc", "shared_b", "shared_coef")
    ok_assoc_data <- ok_assoc[c(2:3,5:6)]
    ok_assoc_interactions <- ok_assoc[c(2,5)]
    
    lag_assoc <- validate_lag_assoc(lag_assoc, M)
    
    assoc <- mapply(assoc = assoc, y_mod = y_mod, lag = lag_assoc, 
                    FUN = validate_assoc, 
                    MoreArgs = list(ok_assoc = ok_assoc, ok_assoc_data = ok_assoc_data,
                                    ok_assoc_interactions = ok_assoc_interactions, 
                                    id_var = id_var, M = M))
    assoc <- check_order_of_assoc_interactions(assoc, ok_assoc_interactions)
    colnames(assoc) <- paste0("Long", 1:M)
    
    # For each submodel, identify any grouping factors that are
    # clustered within id_var (i.e. lower level clustering)
    ok_grp_assocs <- c("sum", "mean")
    clust_basic <- xapply(FUN = get_basic_clust, 
                          cnms  = y_cnms, flist = y_flist,
                          args = list(id_var = id_var))
    clust_stuff <- xapply(FUN = get_extra_clust, 
                          basic_info = clust_basic, flist = y_flist,
                          args = list(id_var = id_var, qnodes = qnodes, 
                                      grp_assoc = grp_assoc, 
                                      ok_grp_assocs = ok_grp_assocs))
    has_clust <- fetch_(clust_stuff, "has_clust")
    if (any(has_clust)) {
      clust_structure <- fetch(clust_stuff, "clust_ids")[has_clust]
      if (n_distinct(clust_structure) > 1L)
        stop("Any longitudinal submodels with a grouping factor clustered within ",
             "patients must use the same clustering structure; that is, the same ",
             "clustering variable and the same number of units clustered within a ",
             "given patient.")    
    } else if (!is.null(grp_assoc)) {
      stop("'grp_assoc' can only be specified when there is a grouping factor ",
           "clustered within patients.")  
    }    
    
    # Return design matrices for evaluating longitudinal submodel quantities
    # at the quadrature points
    auc_qnodes <- 15L
    a_mod_stuff <- mapply(handle_assocmod, 1:M, m_mc, dataLong, y_mod_stuff,
                          clust_stuff = clust_stuff, SIMPLIFY = FALSE, 
                          MoreArgs = list(id_list         = e_mod_stuff$cids, 
                                          times           = e_mod_stuff$cpts, 
                                          assoc           = assoc, 
                                          id_var          = id_var, 
                                          time_var        = time_var, 
                                          eps             = epsilon, 
                                          auc_qnodes      = auc_qnodes,
                                          dataAssoc       = dataAssoc,
                                          env             = calling_env))
    
    # Number of association parameters
    a_K <- get_num_assoc_pars(assoc, a_mod_stuff)
    
  } # end jm block
  
  #---------------------
  # Prior distributions
  #---------------------
  
  ok_dists <- nlist("normal", student_t = "t", "cauchy", "hs", "hs_plus", 
                    "laplace", "lasso")  # disallow product normal
  ok_intercept_dists <- ok_dists[1:3]
  ok_aux_dists <- c(ok_dists[1:3], exponential = "exponential")
  ok_covariance_dists <- c("decov", "lkj")
  
  # Note: *_user_prior_*_stuff objects are stored unchanged for constructing 
  # prior_summary, while *_prior_*_stuff objects are autoscaled

  # Priors for longitudinal submodel(s)
  y_user_prior_stuff <- y_prior_stuff <- 
    xapply(prior, nvars = fetch(y_mod, "X", "K"), link = fetch(y_mod, "family", "link"),
           FUN = handle_glm_prior, args = list(default_scale = 2.5, ok_dists = ok_dists))
  
  y_user_prior_intercept_stuff <- y_prior_intercept_stuff <- 
    xapply(prior_intercept, link = fetch(y_mod, "family", "link"), FUN = handle_glm_prior,
           args = list(nvars = 1, default_scale = 10, ok_dists = ok_intercept_dists))
  
  y_user_prior_aux_stuff <- y_prior_aux_stuff <- 
    xapply(prior_aux, FUN = handle_glm_prior, 
           args = list(nvars = 1, default_scale = 5, link = NULL, ok_dists = ok_aux_dists))  
  
  if (is_jm) { # begin jm block
    
    # Priors for event submodel
    e_user_prior_stuff <- e_prior_stuff <- 
      handle_glm_prior(priorEvent, nvars = e_mod$K, default_scale = 2.5, 
                       link = NULL, ok_dists = ok_dists)  
    
    e_user_prior_intercept_stuff <- e_prior_intercept_stuff <- 
      handle_glm_prior(priorEvent_intercept, nvars = 1, default_scale = 20,
                       link = NULL, ok_dists = ok_intercept_dists)  
    
    e_user_prior_aux_stuff <- e_prior_aux_stuff <-
      handle_glm_prior(priorEvent_aux, nvars = basehaz$df,
                       default_scale = if (basehaz$type_name == "weibull") 2 else 20,
                       link = NULL, ok_dists = ok_e_aux_dists)
    
    e_user_prior_assoc_stuff <- e_prior_assoc_stuff <- 
      handle_glm_prior(priorEvent_assoc, nvars = a_K, default_scale = 2.5,
                       link = NULL, ok_dists = ok_dists)  
    
  } # end jm block
  
  # Priors for covariance
  b_user_prior_stuff <- b_prior_stuff <- handle_cov_prior(
    prior_covariance, cnms = cnms, ok_dists = ok_covariance_dists)
     
  # Autoscaling of priors
  Y_vecs <- fetch(y_mod, "Y", "Y")
  X_mats <- fetch(y_mod, "X", "X")
  y_prior_stuff <- xapply(y_prior_stuff, response = Y_vecs, predictors = X_mats, 
                          family = family, FUN = autoscale_prior)
  y_prior_intercept_stuff <- xapply(y_prior_intercept_stuff, response = Y_vecs,
                                    family = family, FUN = autoscale_prior)
  y_prior_aux_stuff <- xapply(y_prior_aux_stuff, response = Y_vecs,
                              family = family, FUN = autoscale_prior)
  
  if (is_jm) { # begin jm block
    e_prior_stuff <- autoscale_prior(e_prior_stuff, predictors = e_mod$X$X)
    e_prior_intercept_stuff <- autoscale_prior(e_prior_intercept_stuff)
    e_prior_aux_stuff <- autoscale_prior(e_prior_aux_stuff)
    if (a_K)
      e_prior_assoc_stuff <- autoscale_prior(e_prior_assoc_stuff, a_mod_stuff, 
                                             assoc = assoc, family = family)
  } # end jm block
  
  if (b_prior_stuff$prior_dist_name == "lkj") { # autoscale priors for random effect sds
    b_prior_stuff <- split_cov_prior(b_prior_stuff, cnms = cnms, submodel_cnms = y_cnms)
    b_prior_stuff <- xapply(cnms_nms, 
                            FUN = function(nm) {
                              Z_mats <- fetch(y_mod, "Z", "Z", nm)
                              xapply(b_prior_stuff[[nm]], response = Y_vecs, predictors = Z_mats, 
                                     family = family, FUN = autoscale_prior)})
  }
  
  #-------------------------
  # Data for export to Stan
  #-------------------------
  
  standata <- list(  
    M            = as.integer(M),
    dense_X      = !sparse,
    special_case = as.integer(FALSE),
    has_offset   = as.integer(FALSE),
    has_weights  = as.integer(!all(lapply(weights, is.null))),
    family = fetch_array(y_mod, "family", "mvmer_family"),
    link   = fetch_array(y_mod, "family", "mvmer_link"),
    weights = as.array(numeric(0)), # not yet implemented
    prior_PD = as.integer(prior_PD)
  )
  
  # Dimensions
  standata$has_aux <- 
    fetch_array(y_mod, "has_aux", pad_length = 3)
  standata$resp_type <- 
    fetch_array(y_mod, "Y", "is_real", pad_length = 3)
  standata$intercept_type <- 
    fetch_array(y_mod, "intercept_type", "number", pad_length = 3)
  standata$yNobs <- 
    fetch_array(y_mod, "X", "N", pad_length = 3)
  standata$yNeta <- 
    fetch_array(y_mod, "X", "N", pad_length = 3) # same as Nobs for stan_mvmer
  standata$yK <- 
    fetch_array(y_mod, "X", "K", pad_length = 3)
  
  # Response vectors
  Y_integer <- fetch(y_mod, "Y", "integer")
  standata$yInt1 <- if (M > 0) Y_integer[[1]] else as.array(integer(0))  
  standata$yInt2 <- if (M > 1) Y_integer[[2]] else as.array(integer(0))  
  standata$yInt3 <- if (M > 2) Y_integer[[3]] else as.array(integer(0)) 
  
  Y_real <- fetch(y_mod, "Y", "real")
  standata$yReal1 <- if (M > 0) Y_real[[1]] else as.array(double(0)) 
  standata$yReal2 <- if (M > 1) Y_real[[2]] else as.array(double(0)) 
  standata$yReal3 <- if (M > 2) Y_real[[3]] else as.array(double(0)) 
  
  # Population level design matrices
  X <- fetch(y_mod, "X", "X")
  standata$yX1 <- if (M > 0) X[[1]] else matrix(0,0,0)
  standata$yX2 <- if (M > 1) X[[2]] else matrix(0,0,0)
  standata$yX3 <- if (M > 2) X[[3]] else matrix(0,0,0)
  
  X_bar <- fetch(y_mod, "X", "X_bar")
  standata$yXbar1 <- if (M > 0) as.array(X_bar[[1]]) else as.array(double(0))
  standata$yXbar2 <- if (M > 1) as.array(X_bar[[2]]) else as.array(double(0))
  standata$yXbar3 <- if (M > 2) as.array(X_bar[[3]]) else as.array(double(0))
  
  # Data for group specific terms - group factor 1
  b1_varname <- cnms_nms[[1L]] # name of group factor 1
  b1_nvars <- fetch_(y_mod, "Z", "nvars", b1_varname, 
                     null_to_zero = TRUE, pad_length = 3)
  b1_ngrps <- fetch_(y_mod, "Z", "ngrps", b1_varname)
  if (!n_distinct(b1_ngrps) == 1L)
    stop("The number of groups for the grouping factor '", 
         b1_varname, "' should be the same in all submodels.")
  
  standata$bN1 <- b1_ngrps[[1L]]
  standata$bK1 <- sum(b1_nvars)
  standata$bK1_len <- as.array(b1_nvars)
  standata$bK1_idx <- get_idx_array(b1_nvars)
  
  Z1 <- fetch(y_mod, "Z", "Z", b1_varname)
  #if (prior_covariance$dist == "lkj") {
  Z1 <- lapply(Z1, transpose)
  #}
  Z1 <- lapply(Z1, convert_null, "matrix")
  standata$y1_Z1 <- if (M > 0) Z1[[1L]] else matrix(0,0,0)
  standata$y2_Z1 <- if (M > 1) Z1[[2L]] else matrix(0,0,0)
  standata$y3_Z1 <- if (M > 2) Z1[[3L]] else matrix(0,0,0)
  
  Z1_id <- fetch(y_mod, "Z", "group_list", b1_varname)
  Z1_id <- lapply(Z1_id, groups)
  Z1_id <- lapply(Z1_id, convert_null, "arraydouble")
  standata$y1_Z1_id <- if (M > 0) Z1_id[[1L]] else as.array(double(0))
  standata$y2_Z1_id <- if (M > 1) Z1_id[[2L]] else as.array(double(0))
  standata$y3_Z1_id <- if (M > 2) Z1_id[[3L]] else as.array(double(0))
  
  # Data for group specific terms - group factor 2
  if (length(cnms) > 1L) {
    # model has a second grouping factor
    b2_varname <- cnms_nms[[2L]] # name of group factor 2
    b2_nvars <- fetch_(y_mod, "Z", "nvars", b2_varname, 
                       null_to_zero = TRUE, pad_length = 3)
    b2_ngrps <- fetch_(y_mod, "Z", "ngrps", b2_varname)
    if (!n_distinct(b2_ngrps) == 1L)
      stop("The number of groups for the grouping factor '", 
           b2_varname, "' should be the same in all submodels.")
    standata$bN2 <- b2_ngrps[[1L]]
    standata$bK2 <- sum(b2_nvars)
    standata$bK2_len <- as.array(b2_nvars)
    standata$bK2_idx <- get_idx_array(b2_nvars)
    
    Z2 <- fetch(y_mod, "Z", "Z", b2_varname)
    #if (prior_covariance$dist == "lkj") {
    Z2 <- lapply(Z2, transpose)
    #}
    Z2 <- lapply(Z2, convert_null, "matrix")
    standata$y1_Z2 <- if (M > 0) Z2[[1L]] else matrix(0,0,0)
    standata$y2_Z2 <- if (M > 1) Z2[[2L]] else matrix(0,0,0)
    standata$y3_Z2 <- if (M > 2) Z2[[3L]] else matrix(0,0,0)
    
    Z2_id <- fetch(y_mod, "Z", "group_list", b2_varname)
    Z2_id <- lapply(Z2_id, groups)
    Z2_id <- lapply(Z2_id, convert_null, "arraydouble")
    standata$y1_Z2_id <- if (M > 0) Z2[[1L]] else as.array(double(0))
    standata$y2_Z2_id <- if (M > 1) Z2[[2L]] else as.array(double(0))
    standata$y3_Z2_id <- if (M > 2) Z2[[3L]] else as.array(double(0))
    
  } else {
    # no second grouping factor
    standata$bN2 <- 0L
    standata$bK2 <- 0L
    standata$bK2_len <- as.array(rep(0,3L))
    standata$bK2_idx <- get_idx_array(rep(0,3L))
    standata$y1_Z2 <- matrix(0,0,0)
    standata$y2_Z2 <- matrix(0,0,0)
    standata$y3_Z2 <- matrix(0,0,0)
    standata$y1_Z2_id <- as.array(double(0))
    standata$y2_Z2_id <- as.array(double(0))
    standata$y3_Z2_id <- as.array(double(0))
  }
  
  # Priors for population level params
  standata$y_prior_dist <- 
    fetch_array(y_prior_stuff, "prior_dist", pad_length = 3)
  
  prior_mean <- fetch(y_prior_stuff, "prior_mean")
  standata$y_prior_mean1 <- if (M > 0) prior_mean[[1]] else as.array(double(0))
  standata$y_prior_mean2 <- if (M > 1) prior_mean[[2]] else as.array(double(0))
  standata$y_prior_mean3 <- if (M > 2) prior_mean[[3]] else as.array(double(0))
  
  prior_scale <- fetch(y_prior_stuff, "prior_scale")
  standata$y_prior_scale1 <- if (M > 0) as.array(prior_scale[[1]]) else as.array(double(0))
  standata$y_prior_scale2 <- if (M > 1) as.array(prior_scale[[2]]) else as.array(double(0))
  standata$y_prior_scale3 <- if (M > 2) as.array(prior_scale[[3]]) else as.array(double(0))
  
  prior_df <- fetch(y_prior_stuff, "prior_df")
  standata$y_prior_df1 <- if (M > 0) prior_df[[1]] else as.array(double(0))
  standata$y_prior_df2 <- if (M > 1) prior_df[[2]] else as.array(double(0))
  standata$y_prior_df3 <- if (M > 2) prior_df[[3]] else as.array(double(0))
  
  standata$y_global_prior_scale <- 
    fetch_array(y_prior_stuff, "global_prior_scale") # hs priors only
  standata$y_global_prior_df <- 
    fetch_array(y_prior_stuff, "global_prior_df") # hs priors only
  
  # Priors for intercepts 
  standata$y_prior_dist_for_intercept <- 
    fetch_array(y_prior_intercept_stuff, "prior_dist")  
  standata$y_prior_mean_for_intercept <- 
    fetch_array(y_prior_intercept_stuff, "prior_mean")
  standata$y_prior_scale_for_intercept <- 
    fetch_array(y_prior_intercept_stuff, "prior_scale")
  standata$y_prior_df_for_intercept <- 
    fetch_array(y_prior_intercept_stuff, "prior_df")
  
  # Priors for auxiliary params
  standata$y_prior_dist_for_aux <-
    fetch_array(y_prior_aux_stuff, "prior_dist")
  standata$y_prior_mean_for_aux <- 
    fetch_array(y_prior_aux_stuff, "prior_mean")
  standata$y_prior_scale_for_aux <- 
    fetch_array(y_prior_aux_stuff, "prior_scale")
  standata$y_prior_df_for_aux <- 
    fetch_array(y_prior_aux_stuff, "prior_df")
  
  # Priors for group specific terms
  standata$t <- length(cnms)
  standata$p <- as.array(sapply(cnms, length))
  standata$l <- as.array(sapply(cnms_nms, 
                                FUN = function(nm) unique(fetch_(y_mod, "Z", "ngrps", nm))))
  standata$q <- sum(standata$p * standata$l)
  if (prior_covariance$dist == "decov") {
    
    # data for decov prior
    standata$prior_dist_for_cov <- b_prior_stuff$prior_dist
    standata$b_prior_shape <- b_prior_stuff$prior_shape
    standata$b_prior_scale <- b_prior_stuff$prior_scale
    standata$b_prior_concentration <- b_prior_stuff$prior_concentration
    standata$b_prior_regularization <- b_prior_stuff$prior_regularization
    standata$len_concentration <- length(standata$b_prior_concentration)
    standata$len_regularization <- length(standata$b_prior_regularization)
    standata$len_theta_L <- sum(choose(standata$p, 2), standata$p)
    
    # pass empty lkj data
    standata$b1_prior_scale <- as.array(rep(0L, standata$bK1))
    standata$b2_prior_scale <- as.array(rep(0L, standata$bK2))
    standata$b1_prior_df <- as.array(rep(0L, standata$bK1))
    standata$b2_prior_df <- as.array(rep(0L, standata$bK2))
    standata$b1_prior_regularization <- 1.0
    standata$b2_prior_regularization <- 1.0   
    
  } else if (prior_covariance$dist == "lkj") {
    
    # data for lkj prior
    b1_prior_stuff <- b_prior_stuff[[b1_varname]]
    b1_prior_dist <- fetch_(b1_prior_stuff, "prior_dist")
    b1_prior_scale <- fetch_array(b1_prior_stuff, "prior_scale")
    b1_prior_df <- fetch_array(b1_prior_stuff, "prior_df")
    b1_prior_regularization <- fetch_(b1_prior_stuff, "prior_regularization")
    if (n_distinct(b1_prior_dist) > 1L)
      stop2("Bug found: covariance prior should be the same for all submodels.")
    if (n_distinct(b1_prior_regularization) > 1L) {
      stop2("Bug found: prior_regularization should be the same for all submodels.")
    }
    standata$prior_dist_for_cov <- unique(b1_prior_dist)
    standata$b1_prior_scale <- b1_prior_scale
    standata$b1_prior_df <- b1_prior_df
    standata$b1_prior_regularization <- if (length(b1_prior_regularization))
      unique(b1_prior_regularization) else 1.0
    
    if (standata$bK2 > 0) {
      # model has a second grouping factor
      b2_prior_stuff <- b_prior_stuff[[b2_varname]]
      b2_prior_scale <- fetch_array(b2_prior_stuff, "prior_scale")
      b2_prior_df    <- fetch_array(b2_prior_stuff, "prior_df")
      b2_prior_regularization <- fetch_(b2_prior_stuff, "prior_regularization")
      standata$b2_prior_scale <- b2_prior_scale
      standata$b2_prior_df    <- b2_prior_df
      standata$b2_prior_regularization <- unique(b2_prior_regularization)
    } else {
      # model does not have a second grouping factor
      standata$b2_prior_scale <- as.array(double(0))
      standata$b2_prior_df <- as.array(double(0))
      standata$b2_prior_regularization <- 1.0
    }
    
    # pass empty decov data
    standata$len_theta_L <- 0L
    standata$b_prior_shape <- as.array(rep(0L, standata$t))
    standata$b_prior_scale <- as.array(rep(0L, standata$t))
    standata$len_concentration <- 0L
    standata$len_regularization <- 0L
    standata$b_prior_concentration <- as.array(rep(0L, standata$len_concentration))
    standata$b_prior_regularization <- as.array(rep(0L, standata$len_regularization))   
  }
  
  #---------------
  # Prior summary
  #---------------
  
  prior_info <- summarize_jm_prior(
    user_priorLong = y_user_prior_stuff,
    user_priorLong_intercept = y_user_prior_intercept_stuff,
    user_priorLong_aux = y_user_prior_aux_stuff,
    user_priorEvent = e_user_prior_stuff,
    user_priorEvent_intercept = e_user_prior_intercept_stuff,
    user_priorEvent_aux = e_user_prior_aux_stuff,
    user_priorEvent_assoc = e_user_prior_assoc_stuff,
    user_prior_covariance = prior_covariance,
    y_has_intercept = fetch_(y_mod, "X", "has_intercept"),
    y_has_predictors = fetch_(y_mod, "X", "K") > 0,
    e_has_intercept = e_mod$has_intercept,
    e_has_predictors = e_mod$K > 0,
    has_assoc = a_K > 0,
    adjusted_priorLong_scale = fetch(y_prior_stuff, "prior_scale"),
    adjusted_priorLong_intercept_scale = fetch(y_prior_intercept_stuff, "prior_scale"),
    adjusted_priorLong_aux_scale = fetch(y_prior_aux_stuff, "prior_scale"),
    adjusted_priorEvent_scale = e_prior_stuff$prior_scale,
    adjusted_priorEvent_intercept_scale = e_prior_intercept_stuff$prior_scale,
    adjusted_priorEvent_aux_scale = e_prior_aux_stuff$prior_scale,
    adjusted_priorEvent_assoc_scale = e_prior_assoc_stuff$prior_scale,
    family = family, basehaz = basehaz
  )  
  
  #-----------
  # Fit model
  #-----------
  
  # call stan() to draw from posterior distribution
  stanfit <- if (is_jm) stanmodel$jm else stanmodels$mvmer
  pars <- c(if (M > 0 && standata$intercept_type[1]) "yAlpha1", 
            if (M > 1 && standata$intercept_type[2]) "yAlpha2", 
            if (M > 2 && standata$intercept_type[3]) "yAlpha3", 
            if (M > 0 && standata$yK[1]) "yBeta1",
            if (M > 1 && standata$yK[2]) "yBeta2",
            if (M > 2 && standata$yK[3]) "yBeta3",
            if (standata$e_has_intercept) "e_alpha",
            if (standata$e_K) "e_beta",
            if (standata$a_K) "a_beta",
            #if (standata$bK1 > 0) "bMat1",
            #if (standata$bK2 > 0) "bMat2",
            if (M > 0 && standata$has_aux[1]) "yAux1",
            if (M > 1 && standata$has_aux[2]) "yAux2",
            if (M > 2 && standata$has_aux[3]) "yAux3",
            if (length(standata$basehaz_X)) "e_aux",
            if (standata$prior_dist_for_cov == 2 && standata$bK1 > 0) "bCov1",
            if (standata$prior_dist_for_cov == 2 && standata$bK2 > 0) "bCov2",
            if (standata$prior_dist_for_cov == 1 && standata$len_theta_L) "theta_L",
            "mean_PPD")
  
  if (M == 1L) cat("Univariate", if (is_jm) "joint", "model specified\n") else 
  if (M  > 1L) cat("Multivariate", if (is_jm) "joint", "model specified\n")
  if (algorithm == "sampling") {
    cat("\nPlease note the warmup may be much slower than later iterations!\n")             
    sampling_args <- set_jm_sampling_args(
      object = stanfit,
      cnms = cnms,
      user_dots = list(...), 
      user_adapt_delta = adapt_delta,
      user_max_treedepth = max_treedepth,
      data = standata, 
      pars = pars, 
      init = init,
      show_messages = FALSE)
    stanfit <- do.call(sampling, sampling_args)
  } else {
    # meanfield or fullrank vb
    stanfit <- rstan::vb(stanfit, pars = pars, data = standata,
                         algorithm = algorithm, ...)    
  }
  check_stanfit(stanfit)
  
  # Names for pars
  y_intercept_nms <- uapply(1:M, function(m) {
    if (y_mod[[m]]$intercept_type$number > 0) 
      paste0(stub, m, "|(Intercept)") else NULL
  })
  y_beta_nms <- uapply(1:M, function(m) {
    if (!is.null(colnames(X[[m]]))) 
      paste0(stub, m, "|", colnames(X[[m]])) else NULL
  })
  y_aux_nms <- uapply(1:M, function(m) {
    famname_m <- family[[m]]$family
    if (is.gaussian(famname_m)) paste0(stub, m,"|sigma") else
      if (is.gamma(famname_m)) paste0(stub, m,"|shape") else
        if (is.ig(famname_m)) paste0(stub, m,"|lambda") else
          if (is.nb(famname_m)) paste0(stub, m,"|reciprocal_dispersion") else NULL
  })        
  
  if (is_jm) { # begin jm block
    
    e_intercept_nms <- if (e_mod$has_intercept) "Event|(Intercept)" else NULL
    e_beta_nms <- if (e_mod$K) paste0("Event|", colnames(e_mod$Xq)) else NULL  
    e_aux_nms <- if (basehaz$type == 1L) "Event|weibull-shape" else 
      paste0("Event|basehaz-coef", seq(basehaz$df))               
    e_assoc_nms <- character()  
    for (m in 1:M) {
      if (assoc["etavalue",         ][[m]]) e_assoc_nms <- c(e_assoc_nms, paste0("Assoc|Long", m,"|etavalue"))
      if (assoc["etavalue_data",    ][[m]]) e_assoc_nms <- c(e_assoc_nms, paste0("Assoc|Long", m,"|etavalue:", colnames(a_mod_stuff[[m]][["xq_data"]][["etavalue_data"]])))
      if (assoc["etavalue_etavalue",][[m]]) e_assoc_nms <- c(e_assoc_nms, paste0("Assoc|Long", m,"|etavalue:Long", assoc["which_interactions",][[m]][["etavalue_etavalue"]], "|etavalue"))
      if (assoc["etavalue_muvalue", ][[m]]) e_assoc_nms <- c(e_assoc_nms, paste0("Assoc|Long", m,"|etavalue:Long", assoc["which_interactions",][[m]][["etavalue_muvalue"]],  "|muvalue"))
      if (assoc["etaslope",         ][[m]]) e_assoc_nms <- c(e_assoc_nms, paste0("Assoc|Long", m,"|etaslope"))
      if (assoc["etaslope_data",    ][[m]]) e_assoc_nms <- c(e_assoc_nms, paste0("Assoc|Long", m,"|etaslope:", colnames(a_mod_stuff[[m]][["xq_data"]][["etaslope_data"]])))    
      if (assoc["etaauc",           ][[m]]) e_assoc_nms <- c(e_assoc_nms, paste0("Assoc|Long", m,"|etaauc"))
      if (assoc["muvalue",          ][[m]]) e_assoc_nms <- c(e_assoc_nms, paste0("Assoc|Long", m,"|muvalue"))
      if (assoc["muvalue_data",     ][[m]]) e_assoc_nms <- c(e_assoc_nms, paste0("Assoc|Long", m,"|muvalue:", colnames(a_mod_stuff[[m]][["xq_data"]][["muvalue_data"]])))    
      if (assoc["muvalue_etavalue", ][[m]]) e_assoc_nms <- c(e_assoc_nms, paste0("Assoc|Long", m,"|muvalue:Long", assoc["which_interactions",][[m]][["muvalue_etavalue"]], "|etavalue"))
      if (assoc["muvalue_muvalue",  ][[m]]) e_assoc_nms <- c(e_assoc_nms, paste0("Assoc|Long", m,"|muvalue:Long", assoc["which_interactions",][[m]][["muvalue_muvalue"]],  "|muvalue"))
      if (assoc["muslope",          ][[m]]) e_assoc_nms <- c(e_assoc_nms, paste0("Assoc|Long", m,"|muslope"))
      if (assoc["muslope_data",     ][[m]]) e_assoc_nms <- c(e_assoc_nms, paste0("Assoc|Long", m,"|muslope:", colnames(a_mod_stuff[[m]][["xq_data"]][["muslope_data"]])))    
      if (assoc["muauc",            ][[m]]) e_assoc_nms <- c(e_assoc_nms, paste0("Assoc|Long", m,"|muauc"))
    }
    if (sum(standata$size_which_b)) {
      temp_g_nms <- lapply(1:M, FUN = function(m) {
        all_nms <- paste0(paste0("Long", m, "|b["), y_mod[[m]]$Z$group_cnms[[id_var]], "]")
        all_nms[assoc["which_b_zindex",][[m]]]})
      e_assoc_nms <- c(e_assoc_nms, paste0("Assoc|", unlist(temp_g_nms)))
    }
    if (sum(standata$size_which_coef)) {
      temp_g_nms <- lapply(1:M, FUN = function(m) {
        all_nms <- paste0(paste0("Long", m, "|coef["), y_mod[[m]]$Z$group_cnms[[id_var]], "]")
        all_nms[assoc["which_coef_zindex",][[m]]]})
      e_assoc_nms <- c(e_assoc_nms, paste0("Assoc|", unlist(temp_g_nms)))
    }
    
  } # end jm block
  
  # Sigma values in stanmat, and Sigma names
  nc <- sapply(cnms, FUN = length)
  nms <- names(cnms) 
  Sigma_nms <- lapply(cnms, FUN = function(grp) {
    nm <- outer(grp, grp, FUN = paste, sep = ",")
    nm[lower.tri(nm, diag = TRUE)]
  })
  for (j in seq_along(Sigma_nms)) {
    Sigma_nms[[j]] <- paste0(nms[j], ":", Sigma_nms[[j]])
  }
  Sigma_nms <- unlist(Sigma_nms)
  
  if (prior_covariance$dist == "decov" && standata$len_theta_L) {
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
  } 
  
  new_names <- c(y_intercept_nms,
                 y_beta_nms,
                 e_intercept_nms,
                 e_beta_nms,
                 e_assoc_nms,                   
                 if (length(standata$q)) c(paste0("b[", b_nms, "]")),
                 y_aux_nms,
                 e_aux_nms,
                 paste0("Sigma[", Sigma_nms, "]"),
                 paste0(stub, 1:M, "|mean_PPD"), 
                 "log-posterior")
  stanfit@sim$fnames_oi <- new_names
  structure(stanfit, y_mod = y_mod, cnms = cnms)
}

