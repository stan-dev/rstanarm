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

#' Workhorse function for CAR spatial models.
#' 
#' Both \code{stan_besag} and \code{stan_bym} call \code{stan_spatial.fit} to fit
#' the appropriate spatial model. See the documentation for either model for further
#' details on the arguments of \code{stan_spatial.fit}.
#' 
#' @export
#' 

stan_spatial.fit <- function(x, y, w,
                             trials = NULL,
                             family = gaussian(),
                             stan_function = c("stan_besag", "stan_bym"),
                             ...,
                             prior = normal(), prior_intercept = normal(),
                             prior_sigma = NULL, prior_rho = NULL, prior_tau = NULL,
                             prior_PD = FALSE,
                             algorithm = c("sampling", "meanfield", "fullrank"),
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
  
  # Deal with prior_rho and prior_sigma
  if (stan_function == "stan_bym") {
    has_tau <- 1
    prior_tau_stuff <- handle_glm_prior(prior_tau, nvars = 1, family$link, default_scale = 1, 
                                          ok_dists = ok_scale_dists)
    names(prior_tau_stuff) <- paste0(names(prior_tau_stuff), 
                                           "_for_tau")
    for (i in names(prior_tau_stuff)) # prior_{dist, mean, scale, df, autoscale}
      assign(i, prior_tau_stuff[[i]])
    # In BYM rho is constrained so appropriate prior is Beta dist.
    if (is.null(prior_rho)) {
      prior_dist_for_rho <- 0
      prior_rho$alpha <- 1
      prior_rho$beta <- 1
      prior_dist_name_for_rho <- NA
    }
    else {
      prior_dist_for_rho <- 1
      prior_dist_name_for_rho <- "beta"
    }

    prior_mean_for_rho <- 0
    prior_scale_for_rho <- 1
    prior_df_for_rho <- 1
    prior_rho_stuff <- list(prior_dist_name_for_rho = prior_dist_name_for_rho,
                            shape1 = prior_rho$alpha,
                            shape2 = prior_rho$beta,
                            prior_mean_for_rho = prior_mean_for_rho,
                            prior_scale_for_rho = prior_scale_for_rho,
                            prior_df_for_rho = prior_df_for_rho
                            )
  }
  else if (stan_function == "stan_besag") {
    has_tau <- 0
    prior_rho_stuff <- handle_glm_prior(prior_rho, nvars = 1, family$link, default_scale = 1, 
                                          ok_dists = ok_scale_dists)
    names(prior_rho_stuff) <- paste0(names(prior_rho_stuff), 
                                       "_for_rho")
    for (i in names(prior_rho_stuff)) # prior_{dist, mean, scale, df, autoscale}
      assign(i, prior_rho_stuff[[i]])

    prior_rho_stuff$shape1 <- 1
    prior_rho_stuff$shape2 <- 1
    prior_dist_for_tau <- 0
    prior_mean_for_tau <- 0
    prior_scale_for_tau <- 1
    prior_df_for_tau <- 1
  }

  # Deal with prior_sigma... needs more dist options
  if (family$family == "gaussian") {
    has_sigma <- 1
    prior_sigma_stuff <- handle_glm_prior(prior_sigma, nvars = 1, family$link, default_scale = 1, 
                                    ok_dists = ok_dists)
    names(prior_sigma_stuff) <- paste0(names(prior_sigma_stuff), 
                                     "_for_sigma")
    for (i in names(prior_sigma_stuff)) # prior_{dist, mean, scale, df, autoscale}
      assign(i, prior_sigma_stuff[[i]])
  }
  else {
    has_sigma <- 0
    prior_dist_for_sigma <- 0
    prior_mean_for_sigma <- 0
    prior_scale_for_sigma <- 1
    prior_df_for_sigma <- 1
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
                    shape1_rho = c(prior_rho_stuff$shape1),
                    shape2_rho = c(prior_rho_stuff$shape2),
                    prior_dist_for_intercept = prior_dist_for_intercept,
                    prior_dist = prior_dist,
                    prior_dist_rho = prior_dist_for_rho,
                    prior_dist_tau = prior_dist_for_tau,
                    prior_dist_sigma = prior_dist_for_sigma,
                    prior_mean_for_intercept = c(prior_mean_for_intercept),
                    prior_scale_for_intercept = c(prior_scale_for_intercept),
                    prior_df_for_intercept = c(prior_df_for_intercept),
                    prior_mean = as.array(prior_mean),
                    prior_scale = as.array(prior_scale),
                    prior_df = as.array(prior_df),
                    prior_mean_tau = c(prior_mean_for_tau),
                    prior_scale_tau = c(prior_scale_for_tau),
                    prior_df_tau = c(prior_df_for_tau),
                    prior_mean_rho = c(prior_mean_for_rho),
                    prior_scale_rho = c(prior_scale_for_rho),
                    prior_df_rho = c(prior_df_for_rho),
                    prior_mean_sigma = c(prior_mean_for_sigma),
                    prior_scale_sigma = c(prior_scale_for_sigma),
                    prior_df_sigma = c(prior_df_for_sigma),
                    has_intercept = has_intercept,
                    mod = mod)
  standata$X <- array(standata$X, dim = c(standata$N, standata$K))

  # create scaling_factor a la Dan Simpson
  create_scaling_factor <- function(dat) {
    edges <- dat$edges
    #Build the adjacency matrix
    adj.matrix <- Matrix::sparseMatrix(i=edges[,1],j=edges[,2],x=1,symmetric=TRUE)
    #The ICAR precision matrix (note! This is singular)
    Q <-  Matrix::Diagonal(dat$N, Matrix::rowSums(adj.matrix)) - adj.matrix
    
    #Add a small jitter to the diagonal for numerical stability (optional but recommended)
    Q_pert <- Q + Matrix::Diagonal(dat$N) * max(Matrix::diag(Q)) * sqrt(.Machine$double.eps)
    
    # Compute the diagonal elements of the covariance matrix subject to the 
    # constraint that the entries of the ICAR sum to zero.
    # See the function help for further details.
    Q_inv <- inla.qinv(Q_pert, constr=list(A = matrix(1,1,dat$N),e=0))
    
    # Compute the geometric mean of the variances, which are on the diagonal of Q.inv
    scaling_factor <- exp(mean(log(Matrix::diag(Q_inv))))
    return(scaling_factor)
  }
  
  if (stan_function == "stan_bym")
    standata$scaling_factor <- create_scaling_factor(standata)
  else
    standata$scaling_factor <- 0

  pars <- c(if (has_intercept) "alpha", "beta", "rho", if(mod == 2) c("tau"), if(family$family == "gaussian") "sigma",
            "mean_PPD", "psi")

  prior_info <- summarize_spatial_prior(
    user_prior = prior_stuff,
    user_prior_intercept = prior_intercept_stuff,
    user_prior_rho = prior_rho_stuff,
    user_prior_sigma = prior_sigma_stuff,
    user_prior_tau = prior_tau_stuff,
    has_intercept = has_intercept,
    has_predictors = nvars > 0,
    has_sigma = has_sigma,
    has_rho = 1,
    has_tau = has_tau,
    adjusted_prior_scale = prior_scale,
    adjusted_prior_intercept_scale = prior_scale_for_intercept,
    adjusted_prior_scale_rho = prior_scale_for_rho,
    adjusted_prior_scale_sigma = prior_scale_for_sigma,
    adjusted_prior_scale_tau = prior_scale_for_tau)

  stanfit <- stanmodels$spatial
  
  # n.b. optimizing is not supported
  if (algorithm == "optimizing") {
    out <-
      optimizing(stanfit,
                 data = standata,
                 draws = 1000,
                 constrained = TRUE,
                 hessian = TRUE,
                 ...)
    check_stanfit(out)
    out$par <- out$par[!grepl("(phi_raw|theta_raw)", names(out$par))]
    new_names <- names(out$par)
    mark <- grepl("^beta\\[[[:digit:]]+\\]$", new_names)
    if (QR && ncol(xtemp) > 1) {
      out$par[mark] <- R_inv %*% out$par[mark]
      out$theta_tilde[,mark] <- out$theta_tilde[, mark] %*% t(R_inv)
    }
    new_names[mark] <- colnames(xtemp)
    new_names[new_names == "alpha[1]"] <- "(Intercept)"
    names(out$par) <- new_names
    out$stanfit <- suppressMessages(sampling(stanfit, data = standata, chains = 0))
    return(structure(out, prior.info = prior_info))
  }
  else {
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
    else { # algorithm either "meanfield" or "fullrank"
      stanfit <- rstan::vb(stanfit, pars = pars, data = standata,
                           algorithm = algorithm, init = 0.001, ...)
      if (!QR) 
        recommend_QR_for_vb()
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
}

# Summarize spatial prior

summarize_spatial_prior <- function(user_prior,
                                    user_prior_intercept,
                                    user_prior_sigma,
                                    user_prior_rho,
                                    user_prior_tau,
                                    has_intercept,
                                    has_predictors,
                                    has_sigma,
                                    has_rho,
                                    has_tau,
                                    adjusted_prior_scale,
                                    adjusted_prior_intercept_scale,
                                    adjusted_prior_scale_sigma,
                                    adjusted_prior_scale_rho,
                                    adjusted_prior_scale_tau) {
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
  if (has_sigma &&
      user_prior_sigma$prior_dist_name_for_sigma %in% "t") {
    if (all(user_prior_sigma$prior_df_for_sigma == 1)) {
      user_prior_sigma$prior_dist_name_for_sigma <- "cauchy"
    } else {
      user_prior_sigma$prior_dist_name_for_sigma <- "student_t"
    }
  }

  if (has_rho &&
      user_prior_rho$prior_dist_name_for_rho %in% "t") {
    if (all(user_prior_rho$prior_df_for_rho == 1)) {
      user_prior_rho$prior_dist_name_for_rho <- "cauchy"
    } else if (has_rho && user_prior_rho$prior_dist_name_for_rho == "beta") {
      user_prior_rho$prior_dist_name_for_rho <- "beta"
    } else {
      user_prior_rho$prior_dist_name_for_rho <- "student_t"
    }
  }
  if (has_tau &&
      user_prior_tau$prior_dist_name_for_tau %in% "t") {
    if (all(user_prior_tau$prior_df_for_tau == 1)) {
      user_prior_tau$prior_dist_name_for_tau <- "cauchy"
    } else {
      user_prior_tau$prior_dist_name_for_tau <- "student_t"
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
      )),
    prior_sigma = 
      if (!has_sigma) NULL else with(user_prior_sigma, list(
        dist = prior_dist_name_for_sigma,
        location = prior_mean_for_sigma,
        scale = prior_scale_for_sigma,
        adjusted_scale = if (rescaled_int)
          adjusted_prior_sigma_scale else NULL,
        df = if (prior_dist_name_for_sigma %in% "student_t")
          prior_df_for_sigma else NULL
      )),
    prior_rho = 
      if (!has_rho) NULL else with(user_prior_rho, list(
        dist = prior_dist_name_for_rho,
        location = prior_mean_for_rho,
        scale = prior_scale_for_rho,
        shape1 = shape1,
        shape2 = shape2,
        adjusted_scale = if (rescaled_int)
          adjusted_prior_rho_scale else NULL,
        df = if (prior_dist_name_for_rho %in% "student_t")
          prior_df_for_rho else NULL
      )),
    prior_tau = 
      if (!has_tau) NULL else with(user_prior_tau, list(
        dist = prior_dist_name_for_tau,
        location = prior_mean_for_tau,
        scale = prior_scale_for_tau,
        adjusted_scale = if (rescaled_int)
          adjusted_prior_tau_scale else NULL,
        df = if (prior_dist_name_for_tau %in% "student_t")
          prior_df_for_tau else NULL
      ))
  )
  return(prior_list)
}
