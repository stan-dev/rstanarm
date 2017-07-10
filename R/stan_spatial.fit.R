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
  
  if (stan_function == "stan_CARbym")
    mod <- 1
  else if(stan_function == "stan_CARleroux")
    mod <- 2
  
  sparse <- FALSE
  x_stuff <- center_x(x, sparse)
  for (i in names(x_stuff)) # xtemp, xbar, has_intercept
    assign(i, x_stuff[[i]])
  nvars <- ncol(xtemp)

  ok_dists <- nlist("normal", student_t = "t", "cauchy")
  ok_intercept_dists <- ok_dists

  prior_stuff <- handle_glm_prior(prior, nvars, link, default_scale = 2.5, 
                                  ok_dists = ok_dists)
  for (i in names(prior_stuff)) # prior_{dist, mean, scale, df, autoscale}
    assign(i, prior_stuff[[i]])
  
  prior_intercept_stuff <- handle_glm_prior(prior_intercept, nvars = 1, 
                                            default_scale = 10, link = link,
                                            ok_dists = ok_intercept_dists)
  names(prior_intercept_stuff) <- paste0(names(prior_intercept_stuff), 
                                         "_for_intercept")
  for (i in names(prior_intercept_stuff)) # prior_{dist, mean, scale, df, autoscale}_for_intercept
    assign(i, prior_intercept_stuff[[i]])
  
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
                    shape_tau = 1,
                    shape_sigma = 1,
                    shape_nu = 1,
                    scale_tau = 1,
                    scale_sigma = 1,
                    scale_nu = 1,
                    loc_beta = rep(0,ncol(xtemp)),
                    scale_beta = rep(1,ncol(xtemp)),
                    loc_alpha = 0,
                    scale_alpha = 1,
                    has_intercept = has_intercept,
                    mod = mod)
  
  pars <- c(if (has_intercept) "alpha", "beta", if(mod == 1) "sigma", "tau", if(family == "gaussian") "nu",
            "mean_PPD", if (mod < 2) "psi")
  
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
                 colnames(xtemp), "tau2", if(mod == 1) "sigma2",
                 if(family == "gaussian") "nu2", "mean_PPD", "log-posterior", paste0("psi[", 1:standata$N, "]"))
  stanfit@sim$fnames_oi <- new_names
  browser()
  return(structure(stanfit))  # return(structure(stanfit, prior.info = prior_info))
}
