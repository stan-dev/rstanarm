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
                             stan_function,
                             ...,
                             prior = normal(), prior_intercept = normal(),
                             prior_sigma = normal(), prior_tau = beta(), prior_nu = NULL,
                             prior_PD = FALSE,
                             algorithm = c("sampling", "optimizing", "meanfield", "fullrank"),
                             adapt_delta = NULL,
                             QR = FALSE) {
  # check that W is appropriate
  
  algorithm <- match.arg(algorithm)
  family <- match.arg(family)
  
  if (family == "gaussian") {
    y_real <- y
    y_int <- rep(0, length(y))
    link <- "identity"
    trials <- rep(0, length(Y))
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
      trials <- rep(0, length(Y))
      family_num <- 2
    }
  }
  
  if (stan_function == "stan_icar")
    mod <- 1
  else if (stan_function == "stan_bym")
    mod <- 2
  
  sparse <- FALSE
  x_stuff <- center_x(x, sparse)
  for (i in names(x_stuff)) # xtemp, xbar, has_intercept
    assign(i, x_stuff[[i]])
  nvars <- ncol(xtemp)
  
  # temporarily suppress certain prior options
  # ok_dists <- nlist("normal", student_t = "t", "cauchy")
  ok_dists <- nlist("normal")
  ok_intercept_dists <- ok_dists

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
                                          ok_dists = "normal")
    for (i in names(prior_sigma_stuff)) # prior_{dist, mean, scale, df, autoscale}
      assign(i, prior_sigma_stuff[[i]])
    
    prior_tau_stuff <- list(alpha = prior_tau$alpha, beta = prior_tau$beta)
  }
  else if (stan_function == "stan_icar") {
    prior_tau_stuff <- list(alpha = 0, beta = 0)
    prior_sigma_stuff <- list(prior_mean = 0, prior_scale = 0)
  }
  
  # Deal with prior_nu
  if (family == "gaussian") {
    prior_nu_stuff <- handle_glm_prior(prior_nu, nvars = 1, link, default_scale = 1, 
                                    ok_dists = "normal")
    for (i in names(prior_nu_stuff)) # prior_{dist, mean, scale, df, autoscale}
      assign(i, prior_nu_stuff[[i]])
  }
  else {
    prior_nu_stuff <- list(prior_mean = 0, prior_scale = 0)
  }
  
  # QR decomposition for both x and z
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
                    X = if (has_intercept) x[,-1] else x ,  # use xtemp
                    y_real = y_real,
                    y_int = y_int,
                    trials = trials,
                    shape1_tau = c(prior_tau_stuff$alpha),
                    shape2_tau = c(prior_tau_stuff$beta),
                    loc_sigma = c(prior_sigma_stuff$prior_mean),
                    scale_sigma = c(prior_sigma_stuff$prior_scale),
                    loc_nu = c(prior_nu_stuff$prior_mean),
                    scale_nu = c(prior_nu_stuff$prior_scale),
                    loc_beta = as.array(prior_stuff$prior_mean),
                    scale_beta = as.array(prior_stuff$prior_scale),
                    loc_alpha = c(prior_intercept_stuff$prior_mean_for_intercept),
                    scale_alpha = c(prior_intercept_stuff$prior_scale_for_intercept),
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
  
  pars <- c(if (has_intercept) "alpha", "beta", if(mod == 2) c("sigma", "tau"), if(family == "gaussian") "nu",
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
                 colnames(xtemp), if(mod == 2) c("tau", "sigma"),
                 if(family == "gaussian") "nu", "mean_PPD", "log-posterior", paste0("psi[", 1:standata$N, "]"))
  stanfit@sim$fnames_oi <- new_names
  return(structure(stanfit))  # return(structure(stanfit, prior.info = prior_info))
}
