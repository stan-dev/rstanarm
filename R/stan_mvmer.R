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
#' @importFrom lme4 lmerControl glmerControl
#' 
stan_mvmer <- function(formula, data, family = gaussian,
                       na.action = getOption("na.action", "na.omit"), weights, 
                       offset, contrasts, ...,				          
                       prior = normal(), prior_intercept = normal(), 
                       prior_aux = cauchy(0, 5),
                       prior_covariance = decov(), prior_PD = FALSE, 
                       algorithm = c("sampling", "meanfield", "fullrank"), 
                       adapt_delta = NULL, QR = FALSE, sparse = FALSE) {
  
  #-----------------------------
  # Pre-processing of arguments
  #-----------------------------  
  
  # Check for arguments not yet implemented
  if (!missing(weights)) 
    stop("Weights are not yet implemented for stan_mvmer")
  if (!missing(offset)) 
    stop("Offsets are not yet implemented for stan_mvmer")
  if (QR)               
    stop("QR decomposition not yet implemented for stan_mvmer")
  if (sparse)
    stop("'sparse' option is not yet implemented for stan_mvmer")
  if (missing(weights)) weights <- NULL
  if (missing(offset))  offset  <- NULL
  
  algorithm <- match.arg(algorithm)
  
  # Validate arguments
  formula <- validate_arg(formula, "formula")
  M <- length(formula)
  data <- validate_arg(data, "data.frame", null_ok = TRUE, 
                       validate_length = M, broadcast = TRUE)
  data <- lapply(data, as.data.frame)

  # Check family and link
  supported_families <- c("binomial", "gaussian", "Gamma", "inverse.gaussian",
                          "poisson", "neg_binomial_2")
  if (!is(family, "list")) {
    family <- rep(list(family), M) 
  } else if (!length(family) == M) {
    stop("family is a list of the incorrect length.")
  }
  family <- lapply(family, validate_family)
  fam <- lapply(family, function(x) 
    which(pmatch(supported_families, x$family, nomatch = 0L) == 1L))
  if (any(lapply(fam, length) == 0L)) 
    stop("'family' must be one of ", paste(supported_families, collapse = ", "))
  supported_links <- lapply(fam, function(x) supported_glm_links(supported_families[x]))
  link <- mapply(function(x, i) which(supported_links[[i]] == x$link),
                 family, seq_along(family), SIMPLIFY = TRUE)
  if (any(lapply(link, length) == 0L)) 
    stop("'link' must be one of ", paste(supported_links, collapse = ", "))
  
  # Matched call
  calling_env <- parent.frame()
  call <- match.call(expand.dots = TRUE)    
  mc   <- match.call(expand.dots = FALSE)
  mc$prior <- mc$prior_intercept <- mc$prior_covariance <-
    mc$prior_PD <- mc$algorithm <- mc$adapt_delta <-
    mc$... <- mc$QR <- mc$weights <- NULL 
  
  # Create call for longitudinal submodel  
  y_mc <- mc
  
  # Create call for each longitudinal submodel separately
  m_mc <- lapply(1:M, function(m, old_call, env) {
    new_call <- old_call
    fm     <- eval(old_call$formula, env)
    data   <- eval(old_call$data, env)
    family <- eval(old_call$family, env)
    new_call$formula <- if (is(fm, 'list')) fm[[m]] else fm
    new_call$data    <- if (is(data, 'list') && !inherits(data, 'data.frame')) data[[m]] else data
    new_call$family  <- if (is(family, 'list')) family[[m]] else family
    new_call
  }, old_call = y_mc, env = calling_env)
  
  # Is prior* already a list?
  prior           <- maybe_broadcast_priorarg(prior, M)
  prior_intercept <- maybe_broadcast_priorarg(prior_intercept, M)
  prior_aux       <- maybe_broadcast_priorarg(prior_aux, M)
  
  #--------------------------------
  # Data for longitudinal submodel
  #--------------------------------
  
  # Fit separate longitudinal submodels
  y_mod_stuff <- mapply(handle_glmod, m_mc, family, 
                        MoreArgs = list(supported_families = supported_families, 
                                        supported_links = supported_links, 
                                        sparse = sparse, env = calling_env), 
                        SIMPLIFY = FALSE)
  
  # Construct single cnms list for all longitudinal submodels
  cnms <- get_common_cnms(fetch(y_mod_stuff, "cnms"), stub = "y")
  cnms_nms <- names(cnms)
  
  # Prior weights
  if (!is.null(weights)) {
    if (!is(weights, "list")) weights <- rep(list(weights), M)
    weights <- lapply(weights, validate_weights)
  }
  
  #---------------------
  # Prior distributions
  #---------------------
  
  ok_dists <- nlist("normal", student_t = "t", "cauchy", "hs", "hs_plus", 
                    "laplace", "lasso")  # disallow product normal
  ok_intercept_dists <- ok_dists[1:3]
  ok_aux_dists <- c(ok_dists[1:3], exponential = "exponential")
  
  # Note: y_user_prior_*_stuff objects are stored unchanged for constructing 
  # prior_summary, while y_prior_*_stuff objects are autoscaled
  
  y_user_prior_stuff <- y_prior_stuff <- mapply(
    handle_glm_prior,
    prior,
    nvars = fetch(y_mod_stuff, "K"),
    link = fetch(family, "link"),
    MoreArgs = list(default_scale = 2.5, ok_dists = ok_dists), 
    SIMPLIFY = FALSE)
  
  y_user_prior_intercept_stuff <- y_prior_intercept_stuff <- mapply(
    handle_glm_prior,
    prior_intercept,
    link = fetch(family, "link"),
    MoreArgs = list(nvars = 1, default_scale = 10, ok_dists = ok_intercept_dists), 
    SIMPLIFY = FALSE)
  
  y_user_prior_aux_stuff <- y_prior_aux_stuff <- mapply(
    handle_glm_prior,
    prior_aux,
    MoreArgs = list(nvars = 1, default_scale = 5, link = NULL, ok_dists = ok_aux_dists), 
    SIMPLIFY = FALSE)  
  
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
    has_weights  = as.integer(!all(lapply(weights, is.null))),
    has_offset   = as.integer(!is.null(offset)),
    M            = as.integer(M),
    
    # data for longitudinal submodel(s)
    link = as.array(link),
    has_intercept     = fetch_array(y_mod_stuff, "has_intercept"),
    has_intercept_nob = fetch_array(y_mod_stuff, "has_intercept_unbound"),
    has_intercept_lob = fetch_array(y_mod_stuff, "has_intercept_lobound"),
    has_intercept_upb = fetch_array(y_mod_stuff, "has_intercept_upbound"),
    has_aux           = fetch_array(y_mod_stuff, "has_aux"),
    xbar              = fetch_array(y_mod_stuff, "xbar"),
    weights           = as.array(numeric(0)), # not yet implemented
    offset            = if (!is.null(offset))  stop("bug found. offset not yet implemented.") else double(0),
    
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
  
  # Not used
  standata$K_smooth   <- 0L
  standata$S          <- matrix(NA_real_, standata$N, 0L)
  standata$smooth_map <- integer(0)  
  
  # Design matrices
  X <- as.matrix(Matrix::bdiag(fetch(y_mod_stuff, "xtemp")))
  if (sparse) {
    parts <- extract_sparse_parts(X)
    standata$nnz_X <- length(parts$w)
    standata$w_X <- parts$w
    standata$v_X <- parts$v
    standata$u_X <- parts$u
    standata$X <- array(0, dim = c(0L, dim(X)))
  } else {
    standata$X <- array(X, dim = c(1L, dim(X)))
    standata$nnz_X <- 0L
    standata$w_X <- double(0)
    standata$v_X <- integer(0)
    standata$u_X <- integer(0)
  }  
  
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
  
  standata$family <- as.array(sapply(1:M, function(x) {
    return_fam <- switch(family[[x]]$family, 
                         gaussian = 1L, 
                         Gamma = 2L,
                         inverse.gaussian = 3L,
                         binomial = 5L,
                         poisson = 6L,
                         "neg_binomial_2" = 7L)
    if (y_mod_stuff[[x]]$is_bernoulli) return_fam <- 4L
    return_fam}))

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

