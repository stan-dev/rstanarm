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
                       prior_covariance = decov(), prior_PD = FALSE, 
                       algorithm = c("sampling", "meanfield", "fullrank"), 
                       adapt_delta = NULL) {
  
  #-----------------------------
  # Pre-processing of arguments
  #-----------------------------  
  
  algorithm <- match.arg(algorithm)
  if (!missing(weights)) 
    stop("Weights are not yet implemented for stan_mvmer")
  
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
  y_mod <- xapply(yF, yD, family, FUN = handle_y_mod, args = list(id_var = id_var))
  
  # Construct single cnms list for all longitudinal submodels
  cnms <- get_common_cnms(fetch(y_mod_stuff, "cnms"), stub = "y")
  cnms_nms <- names(cnms)
  
  #---------------------
  # Prior distributions
  #---------------------
  
  ok_dists <- nlist("normal", student_t = "t", "cauchy", "hs", "hs_plus", 
                    "laplace", "lasso")  # disallow product normal
  ok_intercept_dists <- ok_dists[1:3]
  ok_aux_dists <- c(ok_dists[1:3], exponential = "exponential")
  
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
  
  # Minimum scaling of priors
  y_prior_stuff           <- Map(autoscale_prior, y_prior_stuff, y_mod_stuff, QR = QR, use_x = TRUE)
  y_prior_intercept_stuff <- Map(autoscale_prior, y_prior_intercept_stuff, y_mod_stuff, QR = QR, use_x = FALSE)
  y_prior_aux_stuff       <- Map(autoscale_prior, y_prior_aux_stuff, y_mod_stuff, QR = QR, use_x = FALSE)
  
  #-------------------------
  # Data for export to Stan
  #-------------------------
  
  standata <- list(  
    # dimensions
    dense_X      = !sparse,
    special_case = as.integer(FALSE),
    has_offset   = as.integer(FALSE),
    has_weights  = as.integer(!all(lapply(weights, is.null))),
    M            = as.integer(M),
    
    # data for longitudinal submodel(s)
    family = fetch_array(y_mod, "family", "mvmer_family"),
    link   = fetch_array(y_mod, "family", "mvmer_link"),
    intercept_type = fetch_array(y_mod, "intercept_type", "number"),
    has_aux = fetch_array(y_mod_stuff, "has_aux"),
    weights = as.array(numeric(0)), # not yet implemented

    # priors
    prior_dist               = fetch_array(y_prior_stuff, "prior_dist"), 
    prior_dist_for_intercept = fetch_array(y_prior_intercept_stuff, "prior_dist"),  
    prior_dist_for_aux       = fetch_array(y_prior_aux_stuff, "prior_dist"),  
    
    # hyperparameters for priors
    prior_mean                = fetch_array(y_prior_stuff,           "prior_mean"), 
    prior_scale               = fetch_array(y_prior_stuff,           "prior_scale"), 
    prior_df                  = fetch_array(y_prior_stuff,           "prior_df"), 
    prior_mean_for_intercept  = fetch_array(y_prior_intercept_stuff, "prior_mean"),
    prior_scale_for_intercept = fetch_array(y_prior_intercept_stuff, "prior_scale"), 
    prior_df_for_intercept    = fetch_array(y_prior_intercept_stuff, "prior_df"),  
    prior_mean_for_aux        = fetch_array(y_prior_aux_stuff,       "prior_mean"),
    prior_scale_for_aux       = fetch_array(y_prior_aux_stuff,       "prior_scale"),
    prior_df_for_aux          = fetch_array(y_prior_aux_stuff,       "prior_df"),
    global_prior_scale        = fetch_array(y_prior_stuff, "global_prior_scale"),
    global_prior_df           = fetch_array(y_prior_stuff, "global_prior_df"), 
    prior_PD = as.integer(prior_PD)
  )

  # Prior flag (same prior for all long submodel)
  standata$prior_special_case <- as.integer(
    (length(unique(standata$prior_dist)) == 1L) && all(standata$prior_dist %in% c(0,1,2)))
  
  # Dimensions
  standata$KM     <- fetch_array(y_mod_stuff, "K")
  standata$NM     <- fetch_array(y_mod_stuff, "N") 
  standata$NM_real<- fetch_array(y_mod_stuff, "real_N") 
  standata$NM_int <- fetch_array(y_mod_stuff, "int_N") 
  standata$K      <- as.integer(sum(standata$KM))
  standata$N      <- as.integer(sum(standata$NM))  
  standata$N_real <- as.integer(sum(fetch_(y_mod_stuff, "real_N"))) 
  standata$N_int  <- as.integer(sum(fetch_(y_mod_stuff, "int_N")))
  standata$N01    <- as.array(t(sapply(fetch(y_mod_stuff, "N01"), cbind))) 
  
  standata$yNobs <- fetch(y_mod, "X", "N")
  
  Y_real <- fetch(y_mod, "Y", "real")
  standata$yReal1 <- if (M > 0) Y_real[[1]] 
  standata$yReal2 <- if (M > 1) Y_real[[2]] 
  standata$yReal3 <- if (M > 2) Y_real[[3]] 

  Y_integer <- fetch(y_mod, "Y", "integer")
  standata$yInt1 <- if (M > 0) Y_integer[[1]] 
  standata$yInt2 <- if (M > 1) Y_integer[[2]] 
  standata$yInt3 <- if (M > 2) Y_integer[[3]] 
  
  Y_integer <- fetch(y_mod, "Y")
  standata$yReal1 <- if (M > 0) Y_real[[1]] 
  standata$yReal2 <- if (M > 1) Y_real[[2]] 
  standata$yReal3 <- if (M > 2) Y_real[[3]] 
  
  X <- fetch(y_mod, "X", "X")
  standata$yX1 <- if (M > 0) X[[1]] else matrix(0,0,0)
  standata$yX2 <- if (M > 1) X[[2]] else matrix(0,0,0)
  standata$yX3 <- if (M > 2) X[[3]] else matrix(0,0,0)
  
  Xbar <- fetch(y_mod, "X", "Xbar")
  standata$yXbar1 <- if (M > 0) Xbar[[1]] else double(0)
  standata$yXbar2 <- if (M > 1) Xbar[[2]] else double(0)
  standata$yXbar3 <- if (M > 2) Xbar[[3]] else double(0)
  
  
  # Combined response vector
  y_y <- fetch(y_mod_stuff, "y")
  y_is_real <- fetch_(y_mod_stuff, "is_real")
  standata$y_real <- as.array(as.numeric(unlist(y_y[y_is_real])))
  standata$y_int  <- as.array(as.integer(unlist(y_y[!y_is_real])))
  
  # Indexing for combined beta vector, response vector, design matrix, weights, etc
  standata$idx      <- get_idx_array(standata$NM)
  standata$idx_real <- get_idx_array(standata$NM_real)
  standata$idx_int  <- get_idx_array(standata$NM_int)
  standata$idx_K    <- get_idx_array(standata$KM)
  
  # Sum dimensions
  for (i in c("has_aux", paste0("has_intercept", c("", "_nob", "_lob", "_upb")))) 
    standata[[paste0("sum_", i)]] <- as.integer(sum(standata[[i]]))
  
  # Data for group-specific terms
  group <- lapply(y_mod_stuff, function(x) {
    pad_reTrms(Ztlist = x$Ztlist, 
               cnms   = x$cnms, 
               flist  = x$flist)})
  Z              <- fetch(group, "Z")
  y_cnms         <- fetch(group, "cnms")
  y_flist_padded <- fetch(group, "flist")
  t <- length(cnms_nms)   # num of unique grouping factors
  pmat <- matrix(0, t, M) # num of group-specific terms
  lmat <- matrix(0, t, M) # num of factor levels
  l <- c()
  for (i in 1:t) {
    for (j in 1:M) {
      pmat[i,j] <- length(y_cnms[[j]][[cnms_nms[i]]])
      lmat[i,j] <- nlevels(y_flist_padded[[j]][[cnms_nms[i]]])
    }
    l[i] <- max(lmat[i,])
    if (!all(lmat[i,] %in% c(0, l[i])))
      stop("The number of factor levels for each of the grouping factors ",
           "must be the same in each of the longitudinal submodels.")     
  }
  qmat <- l * pmat
  p  <- rowSums(pmat) # num group-specific terms for each grouping factor 
  q1 <- rowSums(qmat) # num group-specific coefs for each grouping factor
  q2 <- colSums(qmat) # num group-specific coefs for each submodel
  q  <- sum(qmat)     # total num group-specific coefs
  b_nms <- unlist(Map(make_b_nms, group, m = seq(M), stub = "y"))
  g_nms <- unlist(
    lapply(1:M, FUN = function(m) {
      lapply(1:length(group[[m]]$cnms), FUN = function(i) {
        paste(paste0("y", m), group[[m]]$cnms[[i]], names(group[[m]]$cnms)[i], sep = "|")
      })
    })
  )
  standata$t    <- as.integer(t)
  standata$pmat <- as.array(pmat)
  standata$p    <- as.array(p)
  standata$l    <- as.array(l)
  standata$qmat <- as.array(qmat)
  standata$q1   <- as.array(q1)
  standata$q2   <- as.array(q2)
  standata$q    <- as.integer(q)
  standata$len_theta_L <- sum(choose(p, 2), p)
  Zmerge <- Matrix::bdiag(Z)
  parts <- rstan::extract_sparse_parts(Zmerge)
  standata$num_non_zero <- as.integer(length(parts$w))
  standata$w <- parts$w
  standata$v <- parts$v
  standata$u <- as.array(parts$u)
  
  # Hyperparameters for decov prior
  if (prior_covariance$dist == "decov") {
    decov_args <- prior_covariance
    standata$shape <- as.array(maybe_broadcast(decov_args$shape, t))
    standata$scale <- as.array(maybe_broadcast(decov_args$scale, t))
    standata$len_concentration <- sum(p[p > 1])
    standata$concentration <- 
      as.array(maybe_broadcast(decov_args$concentration, sum(p[p > 1])))
    standata$len_regularization <- sum(p > 1)
    standata$regularization <- 
      as.array(maybe_broadcast(decov_args$regularization, sum(p > 1))) 
  }
  

  #---------------
  # Prior summary
  #---------------
  
  prior_info <- summarize_jm_prior(
    user_priorLong = y_user_prior_stuff,
    user_priorLong_intercept = y_user_prior_intercept_stuff,
    user_priorLong_aux = y_user_prior_aux_stuff,
    user_prior_covariance = prior_covariance,
    y_has_intercept = sapply(y_mod_stuff, `[[`, "has_intercept"),
    y_has_predictors = sapply(y_mod_stuff, function(x) x$K > 0),
    adjusted_priorLong_scale = fetch(y_prior_stuff, "prior_scale"),
    adjusted_priorLong_intercept_scale = fetch(y_prior_intercept_stuff, "prior_scale"),
    adjusted_priorLong_aux_scale = fetch(y_prior_aux_stuff, "prior_scale"),
    family = family, stub_for_names = "y"
  ) 
  names(prior_info) <- gsub("Long", "", names(prior_info), fixed = TRUE)
  
  #-----------
  # Fit model
  #-----------
  
  # call stan() to draw from posterior distribution
  stanfit <- stanmodels$mvmer
  pars <- c(if (standata$sum_has_intercept) "alpha", 
            if (standata$K) "beta",
            if (standata$q) "b",
            if (standata$sum_has_aux) "aux",
            if (standata$len_theta_L) "theta_L",
            "mean_PPD")
  
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
  y_nms <- unlist(lapply(1:M, function(m) 
    if (ncol(y_mod_stuff[[m]]$xtemp)) paste0("y", m, "|", colnames(y_mod_stuff[[m]]$xtemp))))
  y_int_nms <- unlist(lapply(1:M, function(m) 
    if (y_mod_stuff[[m]]$has_intercept) paste0("y", m, "|(Intercept)")))
  y_aux_nms <- character()  
  for (m in 1:M) {
    if (is.gaussian(y_mod_stuff[[m]]$famname))   y_aux_nms <- c(y_aux_nms, paste0("y", m,"|sigma"))
    else if (is.gamma(y_mod_stuff[[m]]$famname)) y_aux_nms <- c(y_aux_nms, paste0("y", m,"|shape"))
    else if (is.ig(y_mod_stuff[[m]]$famname))    y_aux_nms <- c(y_aux_nms, paste0("y", m,"|lambda"))
    else if (is.nb(y_mod_stuff[[m]]$famname))    y_aux_nms <- c(y_aux_nms, paste0("y", m,"|reciprocal_dispersion"))
  }  
  
  # Sigma values in stanmat, and Sigma names
  if (standata$len_theta_L) {
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
  
  new_names <- c(y_int_nms,
                 y_nms,
                 if (length(group)) c(paste0("b[", b_nms, "]")),
                 y_aux_nms,
                 if (standata$len_theta_L) paste0("Sigma[", Sigma_nms, "]"),
                 paste0("y", 1:M, "|mean_PPD"), 
                 "log-posterior")
  stanfit@sim$fnames_oi <- new_names
  
  n_grps <- standata$l - 1
  names(n_grps) <- cnms_nms  # n_grps is num. of levels within each grouping factor
  names(p) <- cnms_nms       # p is num. of variables within each grouping factor
  
  # Undo ordering of matrices if bernoulli
  y_mod_stuff <- lapply(y_mod_stuff, unorder_bernoulli)
  
  call <- match.call(expand.dots = TRUE)
  fit <- nlist(stanfit, family, formula, offset, weights, 
               M, cnms, n_yobs = 
                 unlist(list_nms(fetch(y_mod_stuff, "N"), M, stub = "y")), 
               n_grps, y_mod_stuff, y = 
                 list_nms(fetch(y_mod_stuff, "y"), M, stub = "y"),
               data, call, na.action, algorithm, glmod = fetch(y_mod_stuff, "mod"),
               standata = NULL, terms = NULL, prior.info = prior_info,
               stan_function = "stan_mvmer")
  out <- stanmvreg(fit)
  return(out)
}


#------- internal

# Validate the user specified family and link function and append the 
# family object with family and link information to be used by Stan
#
# @param family A family object
# @param supported_families A character vector of supported family names
# @return A family object appended with the numeric family and link used by Stan
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



handle_y_mod <- function(formula, data, family, id_var) {
  mf <- stats::model.frame(lme4::subbars(formula), data)
  if (!length(formula) == 3L)
    stop2("An outcome variable must be specified.")
  if (!id_var %in% colnames(data))
    stop2("'id_var' must appear in the data frame.")
  if (!id_var %in% colnames(mf))
    stop2("'id_var' must appear in the longitudinal submodel formula.")
  
  # Response vector, design matrices
  Y <- make_Y(formula, mf, family) 
  X <- make_X(formula, mf, drop_intercept = TRUE, centre = TRUE)
  Z <- make_Z(formula, mf) 
  
  # Binomial with >1 trials not allowed by stan_{mvmver,jm}
  is_binomial <- is.binomial(family$family)
  is_bernoulli <- is_binomial && (NCOL(Y) == 1L) && all(Y %in% 0:1)
  if (is_binomial && !is_bernoulli)
    STOP_binomial()
  
  # Various flags
  intercept_type <- check_intercept_type(X, family)
  has_aux <- check_for_aux(family)
  family <- append_mvmer_famlink(family, is_bernoulli)
  
  nlist(Y, X, Z, family, intercept_type, is_real, has_aux)
}

# Return the type of intercept required
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
  } else if (family == "gaussian") { # gaussian, !log
    type <- "no_bound"  
  } else { # gamma/inv-gaus/poisson/nb, !log 
    type <- "lower_bound"  
  }
  number <- switch(type, none = 0L, no_bound = 1L,
                   lower_bound = 2L, upper_bound = 3L)
  nlist(type, number) 
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
  if (subbars) lme4::subbars(fm) else fm
}


# Return the response vector
#
# @param formula The model formula
# @param data A data frame or model frame
# @param family A family object
# @return A named list with:
make_Y <- function(formula, data, family) {
  Y <- as.vector(model.response(mf))
  Y <- validate_glm_outcome_support(Y, family)
  is_real <- check_response_real(family)
  real <- if (is_real) Y else numeric(0) 
  int <- if (!is_real) Y else integer(0) 
  nlist(Y, real, int, is_real)
}

# Return the design matrix, possibly centred
#
# @param formula The model formula
# @param data A data frame or model frame
# @param drop_intercept Logical specifying whether to drop the intercept
#   from the returned model matrix
# @param centre Logical specifying whether to centre the predictors
# @return A named list with: the model matrix (possibly centred), the
#   predictor means, and a logical specifying whether the model formula 
#   included an intercept
make_X <- function(formula, data = NULL, drop_intercept = TRUE, 
                   centre = TRUE) {
  fixed_form <- lme4::nobars(formula)
  fixed_terms <- terms(fixed_form)
  X <- model.matrix(fixed_terms, data)
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
# @param A data frame or model frame
# @return A named list with: a list of design matrices for each of the 
#   grouping factors, a character vector with the name of each of the
#   grouping factors, and a list of list of vector giving the group IDs
#   corresponding to each row of the design matrices
make_Z <- function(formula, data = NULL) {
  re_parts <- lme4::findbars(formula)
  re_forms <- lapply(re_parts, reformulate_rhs, subbars = TRUE)
  if (length(re_forms) > 2L)
    stop2("A maximum of 2 grouping factors are allowed.")
  group_mat  <- lapply(re_forms, model.matrix, data)
  group_vars <- sapply(re_forms, get_group_factor)
  group_col  <- xapply(group_mat, group_vars,
    FUN = function(x, y) which(y %in% colnames(x)))
  group_list <- xapply(group_mat, group_col,
    FUN = function(x, y) x[, y])
  Z <- xapply(group_mat, group_col, 
    FUN = function(x, y) x[, -y, drop = FALSE])
  group_cnms <- lapply(Z, colnames)
  names(group_cnms) <- group_vars
  nlist(Z, group_list, group_vars, group_cnms)
}

# Return the name of the grouping factor in a random
# effects formula (ie. the variable on the RHS of "|")
#
# @param x The random effects formula (for one grouping factor)
# @return A character string, the name of the grouping factor
get_group_factor <- function(x) {
  terms <- strsplit(deparse(x), "\\s\\|\\s")
  if (!length(terms) == 2L)
    stop2("Could not parse the random effects formula.")
  terms[[2L]]
}




