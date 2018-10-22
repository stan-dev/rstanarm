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

#' Workhorse function for CAR models.
#' 
#' Both \code{stan_besag} and \code{stan_bym2} call \code{stan_spatial.fit} to
#' fit the appropriate spatial model. See the documentation for either modeling
#' function for further details on the arguments of \code{stan_spatial.fit}.
#' 
#' @export
#' 

stan_spatial.fit <- function(x, y, w,
                             stan_function = c("stan_besag", "stan_bym2"),
                             family = NULL,
                             trials = NULL,
                             order = c(1,2),
                             ...,
                             prior = normal(), prior_intercept = normal(),
                             prior_tau = normal(), prior_aux = NULL, prior_rho = NULL,
                             prior_PD = FALSE,
                             algorithm = c("sampling", "meanfield", "fullrank"),
                             adapt_delta = NULL,
                             QR = FALSE) {

  w[upper.tri(w)] <- 0
  
  # convert W to a sparse matrix if not already sparse.
  if(!is(w, "sparseMatrix"))
    w <- Matrix(w, sparse = TRUE)
  
  # pull out adjacency pairs from W
  edges <- summary(w)  # analagous to `which(w == 1, arr.ind = TRUE)` on dense matrix
  edges <- edges[,grep("^i$|^j$", colnames(edges))]
  
  algorithm <- match.arg(algorithm)

  family <- validate_family(family)
  supported_families <- c("binomial", "gaussian", "poisson", "Gamma", "neg_binomial_2")
  fam <- which(pmatch(supported_families, family$family, nomatch = 0L) == 1L)
  if (!length(fam)) 
    stop("'family' must be one of ", paste(supported_families, collapse = ", "))
  
  supported_links <- supported_glm_links(supported_families[fam])
  link <- which(supported_links == family$link)
  if (!length(link)) 
    stop("'link' must be one of ", paste(supported_links, collapse = ", "))
  
  family_num <- switch(family$family,
                       gaussian = 1,
                       poisson = 6,
                       neg_binomial_2 = 7,
                       binomial = 5,
                       Gamma = 2)
  
  # for when consistent-family-numbers gets merged
  # family_num <- switch(family$family,
  #                      gaussian = 1,
  #                      Gamma = 2,
  #                      inv_gaussian = 3,
  #                      beta = 4,
  #                      binomial = 5,
  #                      poisson = 6,
  #                      neg_binomial_2 = 7)
  
  if (family$family %in% c("gaussian", "Gamma")) {
    is_continuous <- TRUE
    y_real <- y
    y_int <- array(0, dim = c(0))
  }
  else {
    is_continuous <- FALSE
    y_real <- array(0, dim = c(0))
    y_int <- y
  }
  
  if (family$family %in% c("binomial", "poisson"))
    has_aux <- FALSE
  else
    has_aux <- TRUE
  
  if (family$family != "binomial")
    trials <- array(0, dim = c(0))
  
  if (family$family %in% c("binomial", "poisson", "neg_binomial_2")) {
    if(!is.integer(y_int))
      stop("Outcome must be an integer for count likelihoods.")
    if (family$family == "binomial" & (is.null(trials) | any(y > trials)))
      stop("Outcome values must be less than or equal to the corresponding value in `trials`.")
  }
  
  if (stan_function == "stan_besag")
    model_type <- 1
  else if (stan_function == "stan_bym")
    model_type <- 2
  else if (stan_function == "stan_bym2")
    model_type <- 3
  
  if (!(order %in% c(1,2)))
    stop("Argument 'order' must be 1 or 2.")
  
  x_stuff <- center_x(x, sparse = FALSE)
  for (i in names(x_stuff)) # xtemp, xbar, has_intercept
    assign(i, x_stuff[[i]])
  nvars <- ncol(xtemp)
  
  ok_dists <- nlist("normal", student_t = "t", "cauchy", "hs", "hs_plus", 
                    "laplace", "lasso", "product_normal")
  ok_intercept_dists <- ok_dists
  ok_scale_dists <- nlist("normal", student_t = "t", "cauchy", "exponential")
  
  # Deal with prior_intercept
  prior_intercept_stuff <- handle_glm_prior(prior_intercept, nvars = 1, 
                                            default_scale = 10, link = family$link,
                                            ok_dists = ok_intercept_dists)
  names(prior_intercept_stuff) <- paste0(names(prior_intercept_stuff), 
                                         "_for_intercept")
  for (i in names(prior_intercept_stuff)) # prior_{dist, mean, scale, df, autoscale}_for_intercept
    assign(i, prior_intercept_stuff[[i]])

  # Deal with prior
  prior_stuff <- handle_glm_prior(prior, nvars, family$link, default_scale = 2.5, 
                                  ok_dists = ok_dists)
  for (i in names(prior_stuff)) # prior_{dist, mean, scale, df, autoscale}
    assign(i, prior_stuff[[i]])
  
  # Deal with prior_tau
  prior_tau_stuff <- handle_glm_prior(prior_tau, nvars = 1, family$link, default_scale = 1, 
                                      ok_dists = ok_scale_dists)
  names(prior_tau_stuff) <- paste0(names(prior_tau_stuff), 
                                   "_for_tau")
  for (i in names(prior_tau_stuff)) # prior_{dist, mean, scale, df, autoscale}
    assign(i, prior_tau_stuff[[i]])
  
  # Deal with prior_rho
  if (stan_function == "stan_bym2") {
    has_rho <- 1
    if (is.null(prior_rho)) {
      prior_dist_for_rho <- 0
      prior_rho$alpha <- 0
      prior_rho$beta <- 0
      prior_dist_name_for_rho <- NA
    }
    else {
      prior_dist_for_rho <- 1
      prior_dist_name_for_rho <- "beta"
    }
    prior_rho_stuff <- handle_glm_prior(NULL, nvars = 1, family$link, default_scale = 1, 
                                        ok_dists = ok_scale_dists)
    names(prior_rho_stuff) <- paste0(names(prior_rho_stuff), 
                                     "_for_rho")
    for (i in names(prior_rho_stuff)) # prior_{dist, mean, scale, df, autoscale}
      assign(i, prior_rho_stuff[[i]])
    prior_rho_stuff$shape1 <- prior_rho$alpha
    prior_rho_stuff$shape2 <- prior_rho$beta
  }
  else if (stan_function == "stan_besag") {
    has_rho <- 0
    prior_dist_for_rho <- 0
    # prior_rho_stuff <- list(prior_dist_name_for_rho = NA)
    # prior_scale_for_rho <- 0
    # prior_rho_stuff$shape1 <- 0
    # prior_rho_stuff$shape2 <- 0
    prior_rho_stuff <- handle_glm_prior(NULL, nvars = 1, family$link, default_scale = 1, 
                                        ok_dists = ok_scale_dists)
    names(prior_rho_stuff) <- paste0(names(prior_rho_stuff), 
                                     "_for_rho")
    for (i in names(prior_rho_stuff)) # prior_{dist, mean, scale, df, autoscale}
      assign(i, prior_rho_stuff[[i]])
    prior_rho_stuff$shape1 <- 0
    prior_rho_stuff$shape2 <- 0
  }
  else if (stan_function == "stan_bym") {
    has_rho <- 1
    prior_rho_stuff <- handle_glm_prior(prior_rho, nvars = 1, family$link, default_scale = 1, 
                                        ok_dists = ok_scale_dists)
    names(prior_rho_stuff) <- paste0(names(prior_rho_stuff), 
                                     "_for_rho")
    for (i in names(prior_rho_stuff)) # prior_{dist, mean, scale, df, autoscale}
      assign(i, prior_rho_stuff[[i]])
    prior_rho_stuff$shape1 <- 0
    prior_rho_stuff$shape2 <- 0
  }
  
  # deal with auxiliary parameter
  if (has_aux) {
    prior_aux_stuff <- handle_glm_prior(prior_aux, nvars = 1, family$link, default_scale = 1, 
                                    ok_dists = ok_dists)
    names(prior_aux_stuff) <- paste0(names(prior_aux_stuff), 
                                     "_for_aux")
    for (i in names(prior_aux_stuff)) # prior_{dist, mean, scale, df, autoscale}
      assign(i, prior_aux_stuff[[i]])
  }
  else {
    prior_dist_for_aux <- 0
    prior_mean_for_aux <- 0
    prior_scale_for_aux <- 1
    prior_df_for_aux <- 1
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

  # need to use uncentered version
  standata <- nlist(N = nrow(xtemp),
                    K = ncol(xtemp),
                    edges = edges,
                    E_n = nrow(edges),
                    family = family_num,
                    link = link,
                    is_continuous = is_continuous,
                    has_aux = has_aux,
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
                    prior_dist_for_aux = prior_dist_for_aux,
                    prior_mean_for_intercept = c(prior_mean_for_intercept),
                    prior_scale_for_intercept = c(prior_scale_for_intercept),
                    prior_df_for_intercept = c(prior_df_for_intercept),
                    prior_mean = as.array(prior_mean),
                    prior_scale = as.array(prior_scale),
                    prior_df = as.array(prior_df),
                    prior_mean_rho = c(prior_mean_for_rho),
                    prior_scale_rho = c(prior_scale_for_rho),
                    prior_df_rho = c(prior_df_for_rho),
                    prior_mean_tau = c(prior_mean_for_tau),
                    prior_scale_tau = c(prior_scale_for_tau),
                    prior_df_tau = c(prior_df_for_tau),
                    prior_mean_for_aux = c(prior_mean_for_aux),
                    prior_scale_for_aux = c(prior_scale_for_aux),
                    prior_df_for_aux = c(prior_df_for_aux),
                    has_intercept = has_intercept,
                    model_type = model_type,
                    global_prior_df,
                    global_prior_df_for_intercept,
                    global_prior_scale,
                    global_prior_scale_for_intercept,
                    num_normals = if(prior_dist == 7) as.integer(prior_df) else integer(0))

  if (stan_function == "stan_bym2")
    standata$scaling_factor <- create_scaling_factor(standata)
  else
    standata$scaling_factor <- 0
  
  standata$order <- order
  if (order == 2) {
    Q <- Matrix::diag(Matrix::rowSums(w)) - w
    Q <- Q %*% Q
    sparse_stuff <- rstan::extract_sparse_parts(Q)
    standata$Q_n <- as.array(length(sparse_stuff$w), dim = 1)
    standata$w <- sparse_stuff$w
    standata$v <- sparse_stuff$v
    standata$u <- sparse_stuff$u
  }
  if (order == 1) {
    standata$Q_n <- array(0, dim = c(0))
    standata$w <- array(0, dim = c(0))
    standata$v <- array(0, dim = c(0))
    standata$u <- array(0, dim = c(0))
  }
  
  pars <- c(if (has_intercept) "alpha", "beta", if(model_type != 1) "rho", "tau", if(has_aux) "aux",
            "mean_PPD", "psi")
  
  switch_aux <- switch(family$family,
                       gaussian = "sigma",
                       poisson = NA,
                       neg_binomial_2 = "reciprocal_dispersion",
                       binomial = NA,
                       Gamma = "shape")

  prior_info <- summarize_spatial_prior(
    user_prior = prior_stuff,
    user_prior_intercept = prior_intercept_stuff,
    user_prior_rho = prior_rho_stuff,
    user_prior_aux = if (has_aux == 1) {prior_aux_stuff} else {NULL},
    user_prior_tau = prior_tau_stuff,
    has_intercept = has_intercept,
    has_predictors = nvars > 0,
    has_aux = has_aux,
    has_rho = has_rho,
    has_tau = 1,
    adjusted_prior_scale = prior_scale,
    adjusted_prior_intercept_scale = prior_scale_for_intercept,
    adjusted_prior_aux_scale = prior_scale_for_aux,
    adjusted_prior_scale_tau = prior_scale_for_tau,
    family = family)

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
                   colnames(xtemp),
                   if(model_type == 1) {"structured"},  # tau
                   if(model_type == 2) {c("structured", "unstructured")},  # rho, tau
                   if(model_type == 3) {c("mixing", "structured")},  # rho, tau
                   if(has_aux) switch_aux, "mean_PPD", paste0("psi[", 1:standata$N, "]"), "log-posterior")
    stanfit@sim$fnames_oi <- new_names
    return(structure(stanfit, prior.info = prior_info))
  }
}

# create scaling_factor a la Dan Simpson
create_scaling_factor <- function(dat) {
  edges <- dat$edges
  # Build the adjacency matrix
  adj.matrix <- Matrix::sparseMatrix(i=edges[,1],j=edges[,2],x=1,symmetric=TRUE)
  # The ICAR precision matrix (note! This is singular)
  Q <-  Matrix::Diagonal(dat$N, Matrix::rowSums(adj.matrix)) - adj.matrix
  # Add a small jitter to the diagonal for numerical stability (optional but recommended)
  Q_pert <- Q + Matrix::Diagonal(dat$N) * max(Matrix::diag(Q)) * sqrt(.Machine$double.eps)
  # Compute the diagonal elements of the covariance matrix subject to the 
  # constraint that the entries of the ICAR sum to zero.
  # See the function help for further details.
  # Q_inv <- INLA::inla.qinv(Q_pert, constr=list(A = matrix(1,1,dat$N),e=0))
  Q_inv <- qinv(Q_pert, A = matrix(1,1,dat$N))
  
  # Compute the geometric mean of the variances, which are on the diagonal of Q.inv
  scaling_factor <- exp(mean(log(Matrix::diag(Q_inv))))
  return(scaling_factor)
}

# qinv function (analagous to inla.qinv)

qinv <- function(Q, A = NULL) {
  # need to replace the line below with the sparse version, using recursions
  Sigma <- Matrix::solve(Q)
  if (is.null(A))
    return(Sigma)
  else {
    A <- matrix(1,1, nrow(Sigma))
    W <- Sigma %*% t(A)
    Sigma_const <- Sigma - W %*% solve(A %*% W) %*% t(W)
    return(Sigma_const)
  }
}

# Summarize spatial prior

summarize_spatial_prior <- function(user_prior,
                                    user_prior_intercept,
                                    user_prior_aux,
                                    user_prior_rho,
                                    user_prior_tau,
                                    has_intercept,
                                    has_predictors,
                                    has_aux,
                                    has_rho,
                                    has_tau,
                                    adjusted_prior_scale,
                                    adjusted_prior_intercept_scale,
                                    adjusted_prior_scale_rho,
                                    adjusted_prior_scale_tau,
                                    adjusted_prior_aux_scale,
                                    family) {
  rescaled_coef <-
    user_prior$prior_autoscale && has_predictors &&
    !is.na(user_prior$prior_dist_name) &&
    !all(user_prior$prior_scale == adjusted_prior_scale)
  rescaled_int <-
    user_prior_intercept$prior_autoscale_for_intercept && has_intercept &&
    !is.na(user_prior_intercept$prior_dist_name_for_intercept) &&
    (user_prior_intercept$prior_scale != adjusted_prior_intercept_scale)
  if (has_aux) {
    rescaled_aux <- user_prior_aux$prior_autoscale_for_aux &&
      !is.na(user_prior_aux$prior_dist_name_for_aux) &&
      (user_prior_aux$prior_scale_for_aux != adjusted_prior_aux_scale) 
  }
  
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
  if (has_aux &&
      user_prior_aux$prior_dist_name_for_aux %in% "t") {
    if (all(user_prior_aux$prior_df_for_aux == 1)) {
      user_prior_aux$prior_dist_name_for_aux <- "cauchy"
    } else {
      user_prior_aux$prior_dist_name_for_aux <- "student_t"
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
    prior_aux = 
      if (!has_aux) NULL else with(user_prior_aux, list(
        dist = prior_dist_name_for_aux,
        location = prior_mean_for_aux,
        scale = prior_scale_for_aux,
        adjusted_scale = if (rescaled_int)
          adjusted_prior_aux_scale else NULL,
        df = if (prior_dist_name_for_aux %in% "student_t")
          prior_df_for_aux else NULL
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
  aux_name <- .rename_aux(family)
  prior_list$prior_aux <- if (is.na(aux_name)) 
    NULL else with(user_prior_aux, list(
      dist = prior_dist_name_for_aux,
      location = if (!is.na(prior_dist_name_for_aux) && 
                     prior_dist_name_for_aux != "exponential")
        prior_mean_for_aux else NULL,
      scale = if (!is.na(prior_dist_name_for_aux) && 
                  prior_dist_name_for_aux != "exponential")
        prior_scale_for_aux else NULL,
      adjusted_scale = if (rescaled_aux)
        adjusted_prior_aux_scale else NULL,
      df = if (!is.na(prior_dist_name_for_aux) && 
               prior_dist_name_for_aux %in% "student_t")
        prior_df_for_aux else NULL, 
      rate = if (!is.na(prior_dist_name_for_aux) && 
                 prior_dist_name_for_aux %in% "exponential")
        1 / prior_scale_for_aux else NULL,
      aux_name = aux_name
    ))
  return(prior_list)
}
