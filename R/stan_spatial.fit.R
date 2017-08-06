# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2013, 2014, 2015, 2016, 2017 Trustees of Columbia University
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

#' @export
#' 

stan_spatial.fit <- function(x, y, w,
                             trials = NULL,
                             family = gaussian(),
                             stan_function = c("stan_besag", "stan_bym"),
                             ...,
                             prior = normal(), prior_intercept = normal(),
                             prior_sigma = NULL, prior_tau = NULL, prior_nu = NULL,
                             prior_PD = FALSE,
                             algorithm = c("sampling", "optimizing", "meanfield", "fullrank"),
                             adapt_delta = NULL,
                             QR = FALSE) {
  # check that W is appropriate
  
  algorithm <- match.arg(algorithm)

  family <- validate_family(family)
  supported_families <- c("binomial", "gaussian", "poisson")
  fam <- which(pmatch(supported_families, family$family, nomatch = 0L) == 1L)
  if (!length(fam)) 
    stop("'family' must be one of ", paste(supported_families, collapse = ", "))
  
  supported_links <- supported_glm_links(supported_families[fam])
  link <- which(supported_links == family$link)
  if (!length(link)) 
    stop("'link' must be one of ", paste(supported_links, collapse = ", "))
  
  if (is.null(trials))
    trials <- rep(0,length(y))
  # Replace all this stuff with family = binomial(link = "logit") calls etc.
  if (family$family == "gaussian") {
    y_real <- y
    y_int <- rep(0, length(y))
    trials <- rep(0, length(y))
    family_num <- 1
  }
  else {
    y_real <- rep(0, length(y))
    y_int <- y
    if (family$family == "binomial") {
      family_num <- 3
      if (is.null(trials) | any(y > trials))
        stop("Outcome values must be less than or equal to the corresponding value in `trials`.")
        
    }
    else {  # poisson
      trials <- rep(0, length(y))
      family_num <- 2
    }
    if(!is.integer(y_int))
      stop("Outcome must an integer for count likelihoods.")
  }
  
  if (stan_function == "stan_besag")
    mod <- 1
  else if (stan_function == "stan_bym")
    mod <- 2
  
  sparse <- FALSE
  x_stuff <- center_x(x, sparse)
  for (i in names(x_stuff)) # xtemp, xbar, has_intercept
    assign(i, x_stuff[[i]])
  nvars <- ncol(xtemp)
  
  ok_dists <- nlist("normal", student_t = "t", "cauchy")
  ok_intercept_dists <- ok_dists
  ok_scale_dists <- nlist("normal", student_t = "t", "cauchy", "exponential")

  # Deal with prior
  prior_stuff <- handle_glm_prior(prior, nvars, family$link, default_scale = 2.5, 
                                  ok_dists = ok_dists)
  for (i in names(prior_stuff)) # prior_{dist, mean, scale, df, autoscale}
    assign(i, prior_stuff[[i]])
  
  # Deal with prior_intercept
  prior_intercept_stuff <- handle_glm_prior(prior_intercept, nvars = 1, 
                                            default_scale = 10, link = family$link,
                                            ok_dists = ok_intercept_dists)
  names(prior_intercept_stuff) <- paste0(names(prior_intercept_stuff), 
                                         "_for_intercept")
  for (i in names(prior_intercept_stuff)) # prior_{dist, mean, scale, df, autoscale}_for_intercept
    assign(i, prior_intercept_stuff[[i]])
  
  # Deal with prior_tau and prior_sigma
  if (stan_function == "stan_bym") {
    prior_sigma_stuff <- handle_glm_prior(prior_sigma, nvars = 1, family$link, default_scale = 1, 
                                          ok_dists = ok_scale_dists)
    names(prior_sigma_stuff) <- paste0(names(prior_sigma_stuff), 
                                           "_for_aux")
    for (i in names(prior_sigma_stuff)) # prior_{dist, mean, scale, df, autoscale}
      assign(i, prior_sigma_stuff[[i]])
    
    prior_tau_stuff <- list(alpha = prior_tau$alpha, beta = prior_tau$beta)
    prior_dist_for_tau <- 0
    prior_mean_for_tau <- 0
    prior_scale_for_tau <- 1
    prior_df_for_tau <- 1
  }
  else if (stan_function == "stan_besag") {
    prior_tau_stuff <- handle_glm_prior(prior_tau, nvars = 1, family$link, default_scale = 1, 
                                          ok_dists = ok_scale_dists)
    names(prior_tau_stuff) <- paste0(names(prior_tau_stuff), 
                                       "_for_tau")
    for (i in names(prior_tau_stuff)) # prior_{dist, mean, scale, df, autoscale}
      assign(i, prior_tau_stuff[[i]])
    prior_tau_stuff$alpha <- 0
    prior_tau_stuff$beta <- 0
    prior_dist_for_aux <- 0
    prior_mean_for_aux <- 0
    prior_scale_for_aux <- 1
    prior_df_for_aux <- 1
  }

  # Deal with prior_nu... needs more dist options
  if (family$family == "gaussian") {
    prior_nu_stuff <- handle_glm_prior(prior_nu, nvars = 1, family$link, default_scale = 1, 
                                    ok_dists = ok_dists)
    names(prior_nu_stuff) <- paste0(names(prior_nu_stuff), 
                                     "_for_nu")
    for (i in names(prior_nu_stuff)) # prior_{dist, mean, scale, df, autoscale}
      assign(i, prior_nu_stuff[[i]])
  }
  else {
    prior_dist_for_nu <- 0
    prior_mean_for_nu <- 0
    prior_scale_for_nu <- 1
    prior_df_for_nu <- 1
  }
  
  # QR decomposition for x
  if (QR) {
    if (nvars <= 1)
      stop("'QR' can only be specified when there are multiple predictors.")
    else {
      cn <- colnames(xtemp)
      decomposition <- qr(xtemp)
      sqrt_nm1 <- sqrt(nrow(xtemp) - 1L)
      Q <- qr.Q(decomposition)
      R_inv <- qr.solve(decomposition, Q) * sqrt_nm1
      xtemp <- Q * sqrt_nm1
      colnames(xtemp) <- cn
      xbar <- c(xbar %*% R_inv) 
    }
  }
  
  # pull out adjacency pairs from W
  adj_fun <- function(W) {
    W[upper.tri(W)] <- NA
    out <- which(W == 1, arr.ind = TRUE)
    return(out)
  }
  edges <- adj_fun(W)

  # need to use uncentered version
  standata <- nlist(N = nrow(xtemp),
                    K = ncol(xtemp),
                    edges = edges,
                    E_n = nrow(edges),
                    family = family_num,
                    link = link,
                    X = xtemp,
                    xbar = as.array(xbar),
                    y_real = y_real,
                    y_int = y_int,
                    trials = trials,
                    shape1_tau = c(prior_tau_stuff$alpha),
                    shape2_tau = c(prior_tau_stuff$beta),
                    prior_dist_for_intercept = prior_dist_for_intercept,
                    prior_dist = prior_dist,
                    prior_dist_tau = prior_dist_for_tau,
                    prior_dist_aux = prior_dist_for_aux,
                    prior_dist_nu = prior_dist_for_nu,
                    prior_mean_for_intercept = c(prior_mean_for_intercept),
                    prior_scale_for_intercept = c(prior_scale_for_intercept),
                    prior_df_for_intercept = c(prior_df_for_intercept),
                    prior_mean = as.array(prior_mean),
                    prior_scale = as.array(prior_scale),
                    prior_df = c(prior_df),
                    prior_mean_aux = c(prior_mean_for_aux),
                    prior_scale_aux = c(prior_scale_for_aux),
                    prior_df_aux = c(prior_df_for_aux),
                    prior_rate_aux = 0,  # FIX ME!!!
                    prior_mean_tau = c(prior_mean_for_tau),
                    prior_scale_tau = c(prior_scale_for_tau),
                    prior_df_tau = c(prior_df_for_tau),
                    prior_rate_tau = 0,  # FIX ME!!!
                    prior_mean_nu = c(prior_mean_for_nu),
                    prior_scale_nu = c(prior_scale_for_nu),
                    prior_df_nu = c(prior_df_for_nu),
                    has_intercept = has_intercept,
                    mod = mod)
  standata$X <- array(standata$X, dim = c(standata$N, standata$K))

  # create scaling_factor a la Dan Simpson
  create_scaling_factor <- function(dat, W) {
    #The ICAR precision matrix (note! This is singular)
    Q <-  diag(rowSums(W), dat$N) - W
    #Add a small jitter to the diagonal for numerical stability (optional but recommended)
    Q_pert <- Q + diag(dat$N) * max(diag(Q)) * sqrt(.Machine$double.eps)
    
    # Compute the diagonal elements of the covariance matrix subject to the 
    # constraint that the entries of the ICAR sum to zero.
    #See the function help for further details.
    Q_inv <- inla.qinv(Q_pert, constr=list(A = matrix(1,1,dat$N),e=0))
    Q_inv <- as.matrix(Q_inv)
    #Compute the geometric mean of the variances, which are on the diagonal of Q.inv
    scaling_factor <- exp(mean(log(diag(Q_inv))))
    return(scaling_factor)
  }
  standata$scaling_factor <- create_scaling_factor(standata, W)
  
  pars <- c(if (has_intercept) "alpha", "beta", "rho", if(mod == 2) c("tau"), if(family$family == "gaussian") "sigma",
            "mean_PPD", "psi")

  prior_info <- summarize_spatial_prior(
    user_prior = prior_stuff,
    user_prior_intercept = prior_intercept_stuff,
    user_prior_tau = prior_tau_stuff,
    user_prior_nu = prior_nu_stuff,
    user_prior_aux = prior_aux_stuff,
    has_intercept = has_intercept,
    has_predictors = nvars > 0,
    adjusted_prior_scale = prior_scale,
    adjusted_prior_intercept_scale = prior_scale_for_intercept,
    adjusted_prior_scale_tau = prior_scale_for_tau,
    adjusted_prior_scale_nu = prior_scale_for_nu,
    adjusted_prior_scale_aux = prior_scale_for_aux)

  stanfit <- stanmodels$spatial

  if (algorithm == "sampling") {
    sampling_args <- set_sampling_args(
      object = stanfit, 
      prior = prior, 
      user_dots = list(...), 
      user_adapt_delta = adapt_delta, 
      data = standata, 
      pars = pars, 
      show_messages = FALSE)
    stanfit <- do.call(sampling, sampling_args)
  }
  else {
    stop(paste("algorithm", algorithm, "is not supported."))
  }
  check_stanfit(stanfit)

  if (QR) {
    thetas <- extract(stanfit, pars = "beta", inc_warmup = TRUE, 
                      permuted = FALSE)
    betas <- apply(thetas, 1:2, FUN = function(theta) R_inv %*% theta)
    end <- tail(dim(betas), 1L)
    for (chain in 1:end) for (param in 1:nrow(betas)) {
      stanfit@sim$samples[[chain]][[has_intercept + param]] <-
        if (ncol(xtemp) > 1) betas[param, , chain] else betas[param, chain]
    } 
  }
  new_names <- c(if (has_intercept) "(Intercept)", 
                 colnames(xtemp), "rho", if(mod == 2) c("tau"),
                 if(family$family == "gaussian") "sigma", "mean_PPD", paste0("psi[", 1:standata$N, "]"), "log-posterior")
  stanfit@sim$fnames_oi <- new_names
  return(structure(stanfit, prior.info = prior_info))
}

# Summarize spatial prior

summarize_spatial_prior <- function(user_prior,
                                    user_prior_intercept,
                                    user_prior_tau,
                                    user_prior_nu,
                                    user_prior_aux,
                                    has_intercept,
                                    has_predictors,
                                    adjusted_prior_scale,
                                    adjusted_prior_intercept_scale,
                                    adjusted_prior_scale_tau,
                                    adjusted_prior_scale_nu,
                                    adjusted_prior_scale_aux) {
  rescaled_coef <-
    user_prior$prior_autoscale && has_predictors &&
    !is.na(user_prior$prior_dist_name) &&
    !all(user_prior$prior_scale == adjusted_prior_scale)
  rescaled_int <-
    user_prior_intercept$prior_autoscale_for_intercept && has_intercept &&
    !is.na(user_prior_intercept$prior_dist_name_for_intercept) &&
    (user_prior_intercept$prior_scale != adjusted_prior_intercept_scale)
  
  if (has_predictors && user_prior$prior_dist_name %in% "t") {
    if (all(user_prior$prior_df == 1)) {
      user_prior$prior_dist_name <- "cauchy"
    } else {
      user_prior$prior_dist_name <- "student_t"
    }
  }
  if (has_intercept &&
      user_prior_intercept$prior_dist_name_for_intercept %in% "t") {
    if (all(user_prior_intercept$prior_df_for_intercept == 1)) {
      user_prior_intercept$prior_dist_name_for_intercept <- "cauchy"
    } else {
      user_prior_intercept$prior_dist_name_for_intercept <- "student_t"
    }
  }
  prior_list <- list(
    prior = 
      if (!has_predictors) NULL else with(user_prior, list(
        dist = prior_dist_name,
        location = prior_mean,
        scale = prior_scale,
        adjusted_scale = if (rescaled_coef)
          adjusted_prior_scale else NULL,
        df = if (prior_dist_name %in% c("student_t", "hs", "hs_plus", 
                                        "lasso", "product_normal"))
          prior_df else NULL
      )),
    prior_intercept = 
      if (!has_intercept) NULL else with(user_prior_intercept, list(
        dist = prior_dist_name_for_intercept,
        location = prior_mean_for_intercept,
        scale = prior_scale_for_intercept,
        adjusted_scale = if (rescaled_int)
          adjusted_prior_intercept_scale else NULL,
        df = if (prior_dist_name_for_intercept %in% "student_t")
          prior_df_for_intercept else NULL
      ))#,
    #prior_aux = 
    #  if (!has_phi) NULL else with(user_prior_aux, list(
    #    dist = prior_dist_name_for_aux,
    #    location = if (!is.na(prior_dist_name_for_aux) && 
    #                   prior_dist_name_for_aux != "exponential")
    #      prior_mean_for_aux else NULL,
    #    scale = if (!is.na(prior_dist_name_for_aux) && 
    #                prior_dist_name_for_aux != "exponential")
    #      prior_scale_for_aux else NULL,
    #    df = if (!is.na(prior_dist_name_for_aux) && 
    #             prior_dist_name_for_aux %in% "student_t")
    #      prior_df_for_aux else NULL, 
    #    rate = if (!is.na(prior_dist_name_for_aux) && 
    #               prior_dist_name_for_aux %in% "exponential")
    #      1 / prior_scale_for_aux else NULL,
    #    aux_name = "phi"
    #  ))
  )
  return(prior_list)
}
