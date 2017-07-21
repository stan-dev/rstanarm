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
                             family = c("binomial", "poisson", "gaussian"),
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
  family <- match.arg(family)
  trials <- ifelse(is.null(trials), rep(0,length(y)), trials)
  # Replace all this stuff with family = binomial(link = "logit") calls etc.
  if (family == "gaussian") {
    y_real <- y
    y_int <- rep(0, length(y))
    link <- "identity"
    trials <- rep(0, length(y))
    family_num <- 1
  }
  else {
    y_real <- rep(0, length(y))
    y_int <- y
    if (family == "binomial") {
      link <- "logit"
      family_num <- 3
    }
    else {  # poisson
      link <- "log"
      trials <- rep(0, length(y))
      family_num <- 2
    }
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
  prior_stuff <- handle_glm_prior(prior, nvars, link, default_scale = 2.5, 
                                  ok_dists = ok_dists)
  for (i in names(prior_stuff)) # prior_{dist, mean, scale, df, autoscale}
    assign(i, prior_stuff[[i]])
  
  # Deal with prior_intercept
  prior_intercept_stuff <- handle_glm_prior(prior_intercept, nvars = 1, 
                                            default_scale = 10, link = link,
                                            ok_dists = ok_intercept_dists)
  names(prior_intercept_stuff) <- paste0(names(prior_intercept_stuff), 
                                         "_for_intercept")
  for (i in names(prior_intercept_stuff)) # prior_{dist, mean, scale, df, autoscale}_for_intercept
    assign(i, prior_intercept_stuff[[i]])
  
  # Deal with prior_tau and prior_sigma
  if (stan_function == "stan_bym") {
    prior_sigma_stuff <- handle_glm_prior(prior_sigma, nvars = 1, link, default_scale = 1, 
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
    prior_tau_stuff <- handle_glm_prior(prior_tau, nvars = 1, link, default_scale = 1, 
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
  if (family == "gaussian") {
    prior_nu_stuff <- handle_glm_prior(prior_nu, nvars = 1, link, default_scale = 1, 
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
                    link = 0,  # FIX ME!!!
                    X = if (has_intercept) x[,-1] else x ,  # use xtemp
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
  
  pars <- c(if (has_intercept) "alpha", "beta", "tau", if(mod == 2) c("sigma"), if(family == "gaussian") "nu",
            "mean_PPD", "psi")
  
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
                 colnames(xtemp), "tau", if(mod == 2) c("sigma"),
                 if(family == "gaussian") "nu", "mean_PPD", "log-posterior", paste0("psi[", 1:standata$N, "]"))
  stanfit@sim$fnames_oi <- new_names
  return(structure(stanfit))  # return(structure(stanfit, prior.info = prior_info))
}
