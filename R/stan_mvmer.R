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

#' Bayesian multivariate generalized linear models with correlated 
#' group-specific terms via Stan
#' 
#' Bayesian inference for multivariate GLMs with group-specific coefficients 
#' that are assumed to be correlated across the GLM submodels.
#' 
#' @export
#' @templateVar pkg stats
#' @templateVar pkgfun glm 
#' @templateVar rareargs na.action,contrasts
#' @template args-same-as-rarely
#' @template args-dots
#' @template args-prior_covariance
#' @template args-prior_PD
#' @template args-algorithm
#' @template args-adapt_delta
#' @template args-QR
#' @template args-sparse
#' 
#' @param formula A two-sided linear formula object describing both the 
#'   fixed-effects and random-effects parts of the longitudinal submodel  
#'   (see \code{\link[lme4]{glmer}} for details). For a multivariate GLM this 
#'   should be a list of such formula objects, with each element
#'   of the list providing the formula for one of the GLM submodels.
#' @param data A data frame containing the variables specified in
#'   \code{formula}. For a multivariate GLM, this can
#'   be either a single data frame which contains the data for all 
#'   GLM submodels, or it can be a list of data frames where each
#'   element of the list provides the data for one of the GLM submodels.
#' @param offset Not currently implemented. Same as \code{\link[stats]{glm}}.
#' @param family The family (and possibly also the link function) for the 
#'   GLM submodel(s). See \code{\link[lme4]{glmer}} for details. 
#'   If fitting a multivariate GLM, then this can optionally be a
#'   list of families, in which case each element of the list specifies the
#'   family for one of the GLM submodels. In other words, a different family
#'   can be specified for each GLM submodel. 
#' @param weights Same as in \code{\link[stats]{glm}},
#'   except that when fitting a multivariate GLM and a list of data frames 
#'   is provided in \code{data} then a corresponding list of weights 
#'   must be provided. If weights are 
#'   provided for one of the GLM submodels, then they must be provided for 
#'   all GLM submodels.
#' @param prior,prior_intercept,prior_aux Same as in \code{\link{stan_glmer}}
#'   except that for a multivariate GLM a list of priors can be provided for 
#'   any of \code{prior}, \code{prior_intercept} or \code{prior_aux} arguments. 
#'   That is, different priors can optionally be specified for each of the GLM  
#'   submodels. If a list is not provided, then the same prior distributions are 
#'   used for each GLM submodel. Note that the \code{"product_normal"} prior is
#'   not allowed for \code{stan_mvmer}.
#'   
#' @details The \code{stan_mvmer} function can be used to fit a multivariate
#'   generalized linear model (GLM) with group-specific terms. The model consists
#'   of distinct GLM submodels, each which contains group-specific terms; within
#'   a grouping factor (for example, patient ID) the grouping-specific terms are
#'   assumed to be correlated across the different GLM submodels. It is 
#'   possible to specify a different outcome type (for example a different
#'   family and/or link function) for each of the GLM submodels. \cr
#'   \cr
#'   Bayesian estimation of the model is performed via MCMC, in the same way as 
#'   for \code{\link{stan_glmer}}. Also, similar to \code{\link{stan_glmer}},
#'   an unstructured covariance matrix is used for the group-specific terms 
#'   within a given grouping factor, with priors on the terms of a decomposition
#'   of the covariance matrix.See \code{\link{priors}} for more information about 
#'   the priors distributions that are available for the covariance matrices, 
#'   the regression coefficients and the intercept and auxiliary parameters.
#'
#' @return A \link[=stanmvreg-objects]{stanmvreg} object is returned.
#' 
#' @seealso \code{\link{stan_glmer}}.
#'    
#' @examples
#' \donttest{
#' #####
#' # A multivariate GLM with two submodels. For the grouping factor 'id', the 
#' # group-specific intercept from the first submodel (logBili) is assumed to
#' # be correlated with the group-specific intercept and linear slope in the 
#' # second submodel (albumin)
#' f1 <- stan_mvmer(
#'         formula = list(
#'           logBili ~ year + (1 | id), 
#'           albumin ~ sex + year + (year | id)),
#'         data = pbcLong, 
#'         # this next line is only to keep the example small in size!
#'         chains = 1, cores = 1, seed = 12345, iter = 1000)
#' summary(f1) 
#' }
#' 
stan_mvmer <- function(formula, data, family = gaussian, weights, ...,				          
                       prior = normal(), prior_intercept = normal(), 
                       prior_aux = cauchy(0, 5),
                       prior_covariance = lkj(), prior_PD = FALSE, 
                       algorithm = c("sampling", "meanfield", "fullrank"), 
                       adapt_delta = NULL, QR = FALSE, sparse = FALSE) {
  
  #-----------------------------
  # Pre-processing of arguments
  #-----------------------------  
  
  algorithm <- match.arg(algorithm)
  if (missing(weights))
    weights <- NULL
  if (!is.null(weights)) 
    stop("Weights are not yet implemented for stan_mvmer")
  if (QR)               
    stop("QR decomposition not yet implemented for stan_jm")
  if (sparse)
    stop("'sparse' option is not yet implemented for stan_jm")
  
  # Formula
  yF <- validate_arg(formula, "formula"); M <- length(yF)
  
  # Data
  yD <- validate_arg(data, "data.frame", validate_length = M)  
  yD <- xapply(yF, yD, FUN = get_all_vars) # drop additional vars
  
  # Family
  ok_classes <- c("function", "family", "character")
  ok_families <- c("binomial", "gaussian", "Gamma", 
                   "inverse.gaussian", "poisson", "neg_binomial_2")
  family <- validate_arg(family, ok_classes, validate_length = M)
  family <- lapply(family, validate_famlink, ok_families)
  family <- lapply(family, append_mvmer_famlink)
  
  # Observation weights
  if (!is.null(weights)) {
    if (!is(weights, "list")) 
      weights <- rep(list(weights), M)
    weights <- lapply(weights, validate_weights)
  }
  
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
  cnms <- get_common_cnms(fetch(y_mod, "Z", "group_cnms"), stub = "y")
  cnms_nms <- names(cnms)
  if (length(cnms_nms) > 2L)
    stop("A maximum of 2 grouping factors are allowed.")
  
  #---------------------
  # Prior distributions
  #---------------------
  
  ok_dists <- nlist("normal", student_t = "t", "cauchy", "hs", "hs_plus", 
                    "laplace", "lasso")  # disallow product normal
  ok_intercept_dists <- ok_dists[1:3]
  ok_aux_dists <- c(ok_dists[1:3], exponential = "exponential")
  ok_covariance_dists <- c("decov", "lkj")
  
  # Note: y_user_prior_*_stuff objects are stored unchanged for constructing 
  # prior_summary, while y_prior_*_stuff objects are autoscaled
  y_user_prior_stuff <- y_prior_stuff <- xapply(
    prior, nvars = fetch(y_mod, "X", "K"), link = fetch(y_mod, "family", "link"),
    FUN = handle_glm_prior, args = list(default_scale = 2.5, ok_dists = ok_dists))
  
  y_user_prior_intercept_stuff <- y_prior_intercept_stuff <- xapply(
    prior_intercept, link = fetch(y_mod, "family", "link"), FUN = handle_glm_prior,
    args = list(nvars = 1, default_scale = 10, ok_dists = ok_intercept_dists))
  
  y_user_prior_aux_stuff <- y_prior_aux_stuff <- xapply(
    prior_aux, FUN = handle_glm_prior, 
    args = list(nvars = 1, default_scale = 5, link = NULL, ok_dists = ok_aux_dists))  
  
  b_user_prior_stuff <- b_prior_stuff <- handle_cov_prior(
      prior_covariance, cnms = cnms, ok_dists = ok_covariance_dists)
  b_prior_stuff <- split_cov_prior(
    b_prior_stuff, cnms = cnms, submodel_cnms = fetch(y_mod, "Z", "group_cnms"))
  
  # Autoscaling of priors
  response <- fetch(y_mod, "Y", "Y")
  Xpredictors <- fetch(y_mod, "X", "X")
  Zpredictors <- xapply(cnms_nms, FUN = function(nm) fetch(y_mod, "Z", "Z", nm))
  y_prior_stuff <- xapply(
    y_prior_stuff, response = response,
    predictors = Xpredictors, family = family,
    FUN = autoscale_prior)
  y_prior_intercept_stuff <- xapply(
    y_prior_intercept_stuff, response = response,
    family = family, FUN = autoscale_prior)
  y_prior_aux_stuff <- xapply(
    y_prior_aux_stuff, response = response,
    family = family, FUN = autoscale_prior)
  b_prior_stuff <- xapply(cnms_nms, FUN = function(nm) {
    xapply(b_prior_stuff[[nm]], predictors = Zpredictors[[nm]], 
           response = response, family = family, FUN = autoscale_prior)})
  
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
  
  Xbar <- fetch(y_mod, "X", "Xbar")
  standata$yXbar1 <- if (M > 0) as.array(Xbar[[1]]) else as.array(double(0))
  standata$yXbar2 <- if (M > 1) as.array(Xbar[[2]]) else as.array(double(0))
  standata$yXbar3 <- if (M > 2) as.array(Xbar[[3]]) else as.array(double(0))
  
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
  Z1 <- lapply(Z1, transpose)
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
    Z2 <- lapply(Z2, transpose)
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
  b1_prior_stuff <- b_prior_stuff[[b1_varname]]
  b_prior_dist <- unique(fetch_(b1_prior_stuff, "prior_dist"))
  if (!n_distinct(b_prior_dist) == 1L)
    stop2("Bug found: covariance prior should be the same for all submodels.")
  standata$prior_dist_for_covariance <- b_prior_dist
  
  if (b_prior_dist == 1L) {
    # decov prior
    
  } else if (b_prior_dist == 2L) {
    # lkj prior
    b1_prior_scale <- fetch_array(b1_prior_stuff, "prior_scale")
    b1_prior_df    <- fetch_array(b1_prior_stuff, "prior_df")
    b1_prior_regularization <- fetch_(b1_prior_stuff, "prior_regularization")
    standata$b1_prior_scale <- b1_prior_scale
    standata$b1_prior_df    <- b1_prior_df
    standata$b1_prior_regularization <- unique(b1_prior_regularization)
    if (!n_distinct(b1_prior_regularization) == 1L) {
      stop2("Bug found: prior_regularization should be the same in all ",
            "submodels with the b1 grouping factor.")
    }
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
  }
  
  #---------------
  # Prior summary
  #---------------
  

  #-----------
  # Fit model
  #-----------
  
  # call stan() to draw from posterior distribution
  stanfit <- stanmodels$mvmer
  pars <- c(if (M > 0 && standata$intercept_type[1]) "yAlpha1", 
            if (M > 1 && standata$intercept_type[2]) "yAlpha2", 
            if (M > 2 && standata$intercept_type[3]) "yAlpha3", 
            if (M > 0 && standata$yK[1]) "yBeta1",
            if (M > 1 && standata$yK[2]) "yBeta2",
            if (M > 2 && standata$yK[3]) "yBeta3",
            if (M > 0 && standata$has_aux[1]) "yAux1",
            if (M > 1 && standata$has_aux[2]) "yAux2",
            if (M > 2 && standata$has_aux[3]) "yAux3",
            if (standata$bK1 == 1) "bVec1",
            if (standata$bK1  > 1) "bMat1",
            if (standata$bK1  > 0) "bSd1",
            if (standata$bK1  > 1) "bCorr1",
            if (standata$bK2 == 1) "bVec2",
            if (standata$bK2  > 1) "bMat2",
            if (standata$bK2  > 0) "bSd2",
            if (standata$bK2  > 1) "bCorr2")
            #"mean_PPD")
  
  cat(paste0(if (M == 1L) "Uni" else "Multi", "variate model specified\n"))
  if (algorithm == "sampling") {
    cat("\nPlease note the warmup phase may be much slower than",
        "later iterations!\n")             
    sampling_args <- set_sampling_args(
      object = stanfit,
      prior = NULL,
      user_dots = list(...), 
      user_adapt_delta = adapt_delta,
      data = standata, 
      pars = pars, 
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
      paste0("y", m, "|(Intercept)") else NULL
  })
  y_beta_nms <- uapply(1:M, function(m) {
    if (!is.null(colnames(X[[m]]))) 
      paste0("y", m, "|", colnames(X[[m]])) else NULL
  })
  y_aux_nms <- uapply(1:M, function(m) {
    famname_m <- family[[m]]$family
    if (is.gaussian(famname_m)) paste0("y", m,"|sigma") else
      if (is.gamma(famname_m)) paste0("y", m,"|shape") else
        if (is.ig(famname_m)) paste0("y", m,"|lambda") else
          if (is.nb(famname_m)) paste0("y", m,"|reciprocal_dispersion") else NULL
  })        
  
  # Sigma values in stanmat, and Sigma names
  if (FALSE && standata$len_theta_L) {
    thetas <- extract(stanfit, pars = "theta_L", inc_warmup = TRUE, 
                      permuted = FALSE)
    nc <- sapply(cnms, FUN = length)
    nms <- names(cnms)
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
    Sigma_nms <- lapply(cnms, FUN = function(grp) {
      nm <- outer(grp, grp, FUN = paste, sep = ",")
      nm[lower.tri(nm, diag = TRUE)]
    })
    for (j in seq_along(Sigma_nms)) {
      Sigma_nms[[j]] <- paste0(nms[j], ":", Sigma_nms[[j]])
    }
    Sigma_nms <- unlist(Sigma_nms)
  }
  
  new_names <- c(y_intercept_nms,
                 y_beta_nms,
                 y_aux_nms,
                 #if (length(group)) c(paste0("b[", b_nms, "]")),
                 #if (standata$len_theta_L) paste0("Sigma[", Sigma_nms, "]"),
                 paste0("y", 1:M, "|mean_PPD"), 
                 "log-posterior")
  #stanfit@sim$fnames_oi <- new_names
  
  #n_grps <- standata$l - 1
  #names(n_grps) <- cnms_nms  # n_grps is num. of levels within each grouping factor
  #names(p) <- cnms_nms       # p is num. of variables within each grouping factor

  call <- match.call(expand.dots = TRUE)
  fit <- nlist(stanfit, formula, family, weights, M, cnms, ngrps,
               n_yobs = unlist(list_nms(fetch(y_mod, "X", "N"), M, stub = "y")), 
               data, algorithm, glmod = y_mod, 
               standata = NULL, terms = NULL, prior.info = NULL,
               stan_function = "stan_mvmer", call = match.call(expand.dots = TRUE))
  out <- stanmvreg(fit)
  return(out)
}


#------- internal

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
    decov_args <- prior
    prior_dist <- 1L
    prior_shape <- as.array(maybe_broadcast(decov_args$shape, t))
    prior_scale <- as.array(maybe_broadcast(decov_args$scale, t))
    prior_concentration <- 
      as.array(maybe_broadcast(decov_args$concentration, sum(p[p > 1])))
    prior_regularization <- 
      as.array(maybe_broadcast(decov_args$regularization, sum(p > 1)))
    prior_df <- NULL
  } else if (prior_dist_name == "lkj") {
    lkj_args <- prior
    prior_dist <- 2L
    prior_shape <- NULL
    prior_scale <- as.array(maybe_broadcast(lkj_args$scale, sum(p)))
    prior_concentration <- NULL
    prior_regularization <- 
      as.array(maybe_broadcast(lkj_args$regularization, sum(p > 1)))
    prior_df <- as.array(maybe_broadcast(lkj_args$df, sum(p)))
  } 
  
  nlist(prior_dist, prior_shape, prior_scale, prior_concentration, 
        prior_regularization, prior_df,
        prior_autoscale = isTRUE(prior$autoscale))
}  

# Construct a list with information on the glmer submodel
#
# @param formula The model formula for the glmer submodel
# @param data The data for the glmer submodel
# @param family The family object for the glmer submodel
# @return A named list with the following elements:
#   Y: named list with the reponse vector and related info
#   X: named list with the fe design matrix and related info
#   Z: named list with the re design matrices and related info
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
  
  # Binomial with >1 trials not allowed by stan_{mvmver,jm}
  is_binomial <- is.binomial(family$family)
  is_bernoulli <- is_binomial && NCOL(Y) == 1L && all(Y %in% 0:1)
  if (is_binomial && !is_bernoulli)
    STOP_binomial()
  
  # Various flags
  intercept_type <- check_intercept_type(X, family)
  has_aux <- check_for_aux(family)
  family <- append_mvmer_famlink(family, is_bernoulli)
  
  nlist(Y, X, Z, family, intercept_type, has_aux)
}

# Return the response vector
#
# @param formula The model formula
# @param model_frame The model frame
# @param family A family object
# @return A named list with: the outcome vector (Y), the outcome vector
#   if real or if integer, and a flag of whether the outcome type is real
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
# @return A named list with: the model matrix (possibly centred), the
#   predictor means, and a logical specifying whether the model formula 
#   included an intercept
make_X <- function(formula, model_frame = NULL, drop_intercept = TRUE, 
                   centre = TRUE) {
  fixed_form <- lme4::nobars(formula)
  fixed_terms <- terms(fixed_form)
  X <- model.matrix(fixed_terms, model_frame)
  has_intercept <- check_for_intercept(X, logical = TRUE)
  if (drop_intercept)
    X <- drop_intercept(X)
  if (centre) {
    if (!drop_intercept)
      stop2("Cannot centre 'x' without dropping the intercept.")
    Xbar <- colMeans(X)
    X <- sweep(X, 2, Xbar, FUN = "-")
  } else {
    Xbar <- rep(0, ncol(X))
  }
  # drop any column of x with < 2 unique values (empty interaction levels)
  sel <- (2 > apply(X, 2L, function(x) length(unique(x))))
  if (any(sel)) {
    warning("Dropped empty interaction levels: ",
            paste(colnames(X)[sel], collapse = ", "))
    X <- X[, !sel, drop = FALSE]
    Xbar <- Xbar[!sel]
  }    
  nlist(X, Xbar, has_intercept, N = NROW(X), K = NCOL(X))
}

# Return the design matrices for the group level terms
#
# @param formula The model formula
# @param model_frame The model frame
# @return A named list with: a list of design matrices for each of the 
#   grouping factors, a character vector with the name of each of the
#   grouping factors, and a list of list of vector giving the group IDs
#   corresponding to each row of the design matrices
make_Z <- function(formula, model_frame = NULL) {
  bars <- lme4::findbars(formula)
  if (length(bars) > 2L)
    stop2("A maximum of 2 grouping factors are allowed.")
  re_parts <- lapply(bars, split_at_bars)
  re_forms <- fetch(re_parts, "re_form")
  Z <- lapply(re_forms, model.matrix, model_frame)
  group_cnms <- lapply(Z, colnames)
  group_vars <- fetch(re_parts, "group_var")
  group_list <- lapply(group_vars, function(x) model_frame[[x]])
  nvars <- sapply(group_cnms, length)
  ngrps <- sapply(group_list, n_distinct)
  names(Z) <- names(group_cnms) <- names(group_list) <- 
    names(nvars) <- names(ngrps) <- group_vars
  nlist(Z, group_vars, group_cnms, group_list, nvars, ngrps)
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

# Split the random effects part of a model formula into
#   - the formula part (ie. the formula on the LHS of "|"), and 
#   - the name of the grouping factor (ie. the variable on the RHS of "|")
#
# @param x Random effects part of a model formula, as returned by lme4::findbars
# @return A named list with two elements: 
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

split_cov_prior <- function(prior_stuff, cnms, submodel_cnms) {
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
  new_prior_stuff
}

