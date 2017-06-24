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

#' @rdname stan_betareg
#' @export
#' @param z For \code{stan_betareg.fit}, a regressor matrix for \code{phi}.
#'   Defaults to an intercept only.
stan_betareg.fit <- function(x, y, z = NULL, 
                             weights = rep(1, NROW(x)), 
                             offset = rep(0, NROW(x)),
                             link = c("logit", "probit", "cloglog", 
                                      "cauchit", "log", "loglog"), 
                             link.phi = NULL, ...,
                             prior = normal(), 
                             prior_intercept = normal(),
                             prior_z = normal(), 
                             prior_intercept_z = normal(),
                             prior_phi = cauchy(0, 5),
                             prior_PD = FALSE, 
                             algorithm = c("sampling", "optimizing", 
                                           "meanfield", "fullrank"),
                             adapt_delta = NULL, 
                             QR = FALSE) {
  
  algorithm <- match.arg(algorithm)
  
  # determine whether the user has passed a matrix for the percision model (z)
  if (is.null(link.phi) && is.null(z)) {
    Z_true <- 0
    z <- model.matrix(y ~ 1)
  } else if (is.null(link.phi) && !(is.null(z))) {
    Z_true <- 1
    link.phi <- "log"
  } else {
    Z_true <- 1
  }

  # link for X variables
  link <- match.arg(link)
  supported_links <- c("logit", "probit", "cloglog", "cauchit", "log", "loglog")
  link_num <- which(supported_links == link)
  if (!length(link)) 
    stop("'link' must be one of ", paste(supported_links, collapse = ", "))
  
  # link for Z variables
  link.phi <- match.arg(link.phi, c(NULL, "log", "identity", "sqrt"))
  supported_phi_links <- c("log", "identity", "sqrt")
  link_num_phi <- which(supported_phi_links == link.phi)
  if (!length(link_num_phi)) 
    stop("'link' must be one of ", paste(supported_phi_links, collapse = ", "))
  if (Z_true == 0)
    link_num_phi <- 0

  # useless assignments to pass R CMD check
  has_intercept <- min_prior_scale <- 
    prior_df <- prior_df_for_intercept <- prior_df_for_intercept_z <- prior_df_z <-
    prior_dist <- prior_dist_for_intercept <- prior_dist_for_intercept_z <- prior_dist_z <-
    prior_mean <- prior_mean_for_intercept <- prior_mean_for_intercept_z <- prior_mean_z <-
    prior_scale <- prior_scale_for_intercept <- prior_scale_for_intercept_z <-
    prior_df_for_aux <- prior_dist_for_aux <- prior_mean_for_aux <- prior_scale_for_aux <-
    xbar <- xtemp <- prior_autoscale <- prior_autoscale_z <- global_prior_scale_z <- NULL

  sparse <- FALSE
  x_stuff <- center_x(x, sparse)
  for (i in names(x_stuff)) # xtemp, xbar, has_intercept
    assign(i, x_stuff[[i]])
  nvars <- ncol(xtemp)
  
  z_stuff <- center_x(z, sparse)
  ztemp <- z_stuff$xtemp
  zbar <- z_stuff$xbar
  has_intercept_z <- z_stuff$has_intercept
  nvars_z <- ncol(ztemp)
  if (Z_true == 0)
    has_intercept_z <- FALSE
  
  ok_dists <- nlist("normal", student_t = "t", "cauchy", "hs", "hs_plus",
                    "laplace", "lasso", "product_normal")
  ok_intercept_dists <- ok_dists[1:3]
  ok_aux_dists <- c(ok_dists[1:3], exponential = "exponential")
  
  # prior distributions (handle_glm_prior() from data_block.R)
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
  
  # prior distributions for parameters on z variables
  prior_stuff_z <- handle_glm_prior(prior_z, nvars_z, link = link.phi, 
                                    default_scale = 2.5, ok_dists = ok_dists)
  for (i in names(prior_stuff_z))
    assign(paste0(i,"_z"), prior_stuff_z[[i]])
  
  prior_intercept_stuff_z <- handle_glm_prior(prior_intercept_z, nvars = 1, 
                                              link = link.phi, default_scale = 10,
                                              ok_dists = ok_intercept_dists)
  names(prior_intercept_stuff_z) <- paste0(names(prior_intercept_stuff_z), 
                                           "_for_intercept")
  for (i in names(prior_intercept_stuff_z))
    assign(paste0(i, "_z"), prior_intercept_stuff_z[[i]])
  
  prior_aux <- prior_phi
  prior_aux_stuff <-
    handle_glm_prior(
      prior_aux,
      nvars = 1,
      default_scale = 5,
      link = NULL, # don't need to adjust scale based on logit vs probit
      ok_dists = ok_aux_dists
    )
  # prior_{dist, mean, scale, df, dist_name, autoscale}_for_aux
  names(prior_aux_stuff) <- paste0(names(prior_aux_stuff), "_for_aux")
  if (is.null(prior_aux)) {
    if (prior_PD)
      stop("'prior_aux' can't be NULL if 'prior_PD' is TRUE.")
    prior_aux_stuff$prior_scale_for_aux <- Inf
  }
  for (i in names(prior_aux_stuff)) 
    assign(i, prior_aux_stuff[[i]])
 
  if (nvars_z == 0) {
      prior_mean_z <- double()
      prior_scale_z <- double()
      prior_df_z <- integer()
  }

  # prior scaling (using sd of predictors)
  min_prior_scale <- 1e-12
  if (prior_dist > 0L && !QR && nvars != 0 && prior_autoscale) {
    prior_scale <- pmax(min_prior_scale, prior_scale / 
                          apply(xtemp, 2L, FUN = function(x) {
                            num.categories <- length(unique(x))
                            x.scale <- 1
                            if (num.categories == 2) {
                              x.scale <- diff(range(x))
                            } else if (num.categories > 2) {
                              x.scale <- sd(x)
                            }
                            return(x.scale)
                          }))
  }
  if (prior_dist_z > 0L && !QR && nvars_z != 0 && prior_autoscale_z) {
    prior_scale_z <- pmax(min_prior_scale, prior_scale_z / 
                            apply(ztemp, 2L, FUN = function(z) {
                              num.categories <- length(unique(z))
                              z.scale <- 1
                              if (num.categories == 2) {
                                z.scale <- diff(range(z))
                              } else if (num.categories > 2) {
                                z.scale <- sd(z)
                              }
                              return(z.scale)
                            }))
  }
  prior_scale <- as.array(pmin(.Machine$double.xmax, prior_scale))
  prior_scale_for_intercept <- 
    min(.Machine$double.xmax, prior_scale_for_intercept)
  if(nvars_z != 0) {
    prior_scale_z <- as.array(pmin(.Machine$double.xmax, prior_scale_z))
    prior_scale_for_intercept_z <- 
      min(.Machine$double.xmax, prior_scale_for_intercept_z)
  }

  # QR decomposition for both x and z
  if (QR) {
    if ((nvars <= 1 && nvars_z <= 1 && Z_true == 1) || 
        (nvars <= 1 && Z_true == 0))
      stop("'QR' can only be specified when there are multiple predictors.")
    if (nvars > 1) {
      cn <- colnames(xtemp)
      decomposition <- qr(xtemp)
      sqrt_nm1 <- sqrt(nrow(xtemp) - 1L)
      Q <- qr.Q(decomposition)
      R_inv <- qr.solve(decomposition, Q) * sqrt_nm1
      xtemp <- Q * sqrt_nm1
      colnames(xtemp) <- cn
      xbar <- c(xbar %*% R_inv) 
    }
    if (Z_true == 1 && nvars_z > 1) {
      cn_z <- colnames(ztemp)
      decomposition_z <- qr(ztemp)
      sqrt_nm1_z <- sqrt(nrow(ztemp) - 1L)
      Q_z <- qr.Q(decomposition_z)
      R_inv_z <- qr.solve(decomposition_z, Q_z) * sqrt_nm1_z
      ztemp <- Q_z * sqrt_nm1_z
      colnames(ztemp) <- cn_z
      zbar <- c(zbar %*% R_inv_z)
    }
  }
  
  # create entries in the data block of the .stan file
  standata <- nlist(
    N = nrow(xtemp), K = ncol(xtemp), 
    xbar = as.array(xbar), dense_X = !sparse,
    X = array(xtemp, dim = c(1L, dim(xtemp))),
    nnz_X = 0L, 
    w_X = double(), 
    v_X = integer(), 
    u_X = integer(), 
    y = y, lb_y = 0, ub_y = 1,
    prior_PD, has_intercept, family = 4L, link = link_num, 
    prior_dist, prior_mean, prior_scale = as.array(pmin(.Machine$double.xmax, prior_scale)), prior_df,
    prior_dist_for_intercept, prior_mean_for_intercept = c(prior_mean_for_intercept), 
    prior_scale_for_intercept = min(.Machine$double.xmax, prior_scale_for_intercept), 
    prior_df_for_intercept = c(prior_df_for_intercept),
    prior_dist_for_aux = prior_dist_for_aux,
    prior_scale_for_aux = prior_scale_for_aux %ORifINF% 0,
    prior_df_for_aux = c(prior_df_for_aux),
    prior_mean_for_aux = c(prior_mean_for_aux),
    prior_dist_for_smooth = 0L, prior_mean_for_smooth = array(NA_real_, dim = 0), 
    prior_scale_for_smooth = array(NA_real_, dim = 0),
    prior_df_for_smooth = array(NA_real_, dim = 0), K_smooth = 0L,
    S = matrix(NA_real_, nrow(xtemp), ncol = 0L), smooth_map = integer(),
    has_weights = length(weights) > 0, weights = weights,
    has_offset = length(offset) > 0, offset = offset,
    t = 0L, 
    p = integer(), 
    l = integer(), 
    q = 0L, 
    len_theta_L = 0L, shape = double(), scale = double(), 
    len_concentration = 0L, concentration = double(),
    len_regularization = 0L, regularization = double(),
    num_non_zero = 0L, 
    w = double(), 
    v = integer(), 
    u = integer(),
    special_case = 0L,
    z_dim = nvars_z,
    link_phi = link_num_phi,
    betareg_z = array(ztemp, dim = c(dim(ztemp))),
    has_intercept_z,
    zbar = array(zbar),
    prior_dist_z, prior_mean_z, prior_df_z,
    prior_scale_z = as.array(pmin(.Machine$double.xmax, prior_scale_z)),
    prior_dist_for_intercept_z, 
    prior_mean_for_intercept_z = c(prior_mean_for_intercept_z), 
    prior_df_for_intercept_z = c(prior_df_for_intercept_z),
    prior_scale_for_intercept_z = min(.Machine$double.xmax, prior_scale_for_intercept_z), 
    # for hs family priors
    global_prior_scale_z,
    # for product normal prior
    num_normals = if (prior_dist == 7) 
      as.array(as.integer(prior_df)) else integer(0),
    num_normals_z = if (prior_dist_z == 7) 
      as.array(as.integer(prior_df_z)) else integer(0)
    )

  # call stan() to draw from posterior distribution
  stanfit <- stanmodels$continuous
  if (Z_true == 1) {
    pars <- c(if (has_intercept) "alpha", 
              "beta", 
              "omega_int", 
              "omega", 
              "mean_PPD")
  } else {
    pars <- c(if (has_intercept) "alpha", 
              "beta", 
              "aux", 
              "mean_PPD")
  }
  
  prior_info <- summarize_betareg_prior(
    user_prior = prior_stuff,
    user_prior_intercept = prior_intercept_stuff,
    user_prior_z = prior_stuff_z,
    user_prior_intercept_z = prior_intercept_stuff_z,
    user_prior_aux = prior_aux_stuff,
    has_phi = !Z_true,
    has_intercept = has_intercept,
    has_intercept_z = has_intercept_z,
    has_predictors = nvars > 0,
    has_predictors_z = nvars_z > 0,
    adjusted_prior_scale = prior_scale,
    adjusted_prior_intercept_scale = prior_scale_for_intercept,
    adjusted_prior_scale_z = prior_scale_z,
    adjusted_prior_intercept_scale_z = prior_scale_for_intercept_z
  )
  

  if (algorithm == "optimizing") {
    out <-
      optimizing(stanfit,
                 data = standata,
                 draws = 1000,
                 constrained = TRUE,
                 ...)
    check_stanfit(out)
    out$par <- out$par[!grepl("eta_z", names(out$par))]
    out$theta_tilde <- out$theta_tilde[, !grepl("eta_z", colnames(out$theta_tilde))]
    new_names <- names(out$par)
    mark <- grepl("^beta\\[[[:digit:]]+\\]$", new_names)
    if (QR && ncol(xtemp) > 1) {
      out$par[mark] <- R_inv %*% out$par[mark]
      out$theta_tilde[,mark] <- out$theta_tilde[, mark] %*% t(R_inv)
    }
    new_names[mark] <- colnames(xtemp)
    new_names[new_names == "alpha[1]"] <- "(Intercept)"
    if (Z_true == 1) {
      new_names[new_names == "omega_int[1]"] <- "(phi)_(Intercept)"
      mark_z <- grepl("^omega\\[[[:digit:]]+\\]$", new_names)
      if (QR && ncol(ztemp) > 1) {
        out$par[mark_z] <- R_inv_z %*% out$par[mark_z]
        out$theta_tilde[,mark_z] <- out$theta_tilde[, mark_z] %*% t(R_inv_z)
      }
      new_names[mark_z] <- paste0("(phi)_", colnames(ztemp))
    } else {
      new_names[new_names == "aux"] <- "(phi)"
    }
    names(out$par) <- new_names
    colnames(out$theta_tilde) <- new_names
    out$stanfit <- suppressMessages(sampling(stanfit, data = standata, chains = 0))
    return(structure(out, prior.info = prior_info))
  } else {
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
    } else { # algorithm either "meanfield" or "fullrank"
      stanfit <- rstan::vb(stanfit, pars = pars, data = standata,
                           algorithm = algorithm, init = 0.001, ...)
      if (!QR) 
        recommend_QR_for_vb()
    }
    check_stanfit(stanfit)
    if (QR) {
      if (ncol(xtemp) > 1) {
        thetas <- extract(stanfit, pars = "beta", inc_warmup = TRUE, 
                          permuted = FALSE)
        betas <- apply(thetas, 1:2, FUN = function(theta) R_inv %*% theta)
        end <- tail(dim(betas), 1L)
        for (chain in 1:end) for (param in 1:nrow(betas)) {
          stanfit@sim$samples[[chain]][[has_intercept + param]] <-
            if (ncol(xtemp) > 1) betas[param, , chain] else betas[param, chain]
        } 
      }
      if (Z_true == 1 & ncol(ztemp) > 1) {
        thetas_z <- extract(stanfit, pars = "omega", 
                            inc_warmup = TRUE, permuted = FALSE)
        omegas <- apply(thetas_z, 1:2, FUN = function(theta) R_inv_z %*% theta)
        end_z <- tail(dim(omegas), 1L)
        for (chain_z in 1:end_z) for (param_z in 1:nrow(omegas)) {
          sel <- has_intercept + ncol(xtemp) + has_intercept_z + param_z
          stanfit@sim$samples[[chain_z]][[sel]] <-
            if (ncol(ztemp) > 1) omegas[param_z, , chain_z] else omegas[param_z, chain_z]
        }
      }
    }
    if (Z_true == 1) {
      new_names <- c(if (has_intercept) "(Intercept)", 
                     colnames(xtemp),
                     if (has_intercept_z) "(phi)_(Intercept)", 
                     paste0("(phi)_", colnames(ztemp)),
                     "mean_PPD", "log-posterior")
    } else {
      new_names <- c(if (has_intercept) "(Intercept)", 
                     colnames(xtemp), 
                     "(phi)",  
                     "mean_PPD", "log-posterior")
    }
    stanfit@sim$fnames_oi <- new_names
    return(structure(stanfit, prior.info = prior_info))
  }
}


# Create "prior.info" attribute needed for prior_summary()
#
# @param user_* The user's prior, prior_intercept, prior_covariance, and 
#   prior_options specifications. For prior and prior_intercept these should be
#   passed in after broadcasting the df/location/scale arguments if necessary.
# @param has_intercept T/F, does model have an intercept?
# @param has_predictors T/F, does model have predictors?
# @param adjusted_prior_* adjusted scales computed if prior_ops$scaled is TRUE
# @return A named list with components 'prior', 'prior_intercept', and possibly 
#   'prior_covariance', each of which itself is a list containing the needed
#   values for prior_summary.

summarize_betareg_prior <-
  function(user_prior,
           user_prior_intercept,
           user_prior_z,
           user_prior_intercept_z,
           user_prior_aux,
           has_phi,
           has_intercept, 
           has_intercept_z, 
           has_predictors,
           has_predictors_z,
           adjusted_prior_scale,
           adjusted_prior_intercept_scale,
           adjusted_prior_scale_z,
           adjusted_prior_intercept_scale_z) {
    rescaled_coef <-
      user_prior$prior_autoscale && has_predictors &&
      !is.na(user_prior$prior_dist_name) &&
      !all(user_prior$prior_scale == adjusted_prior_scale)
    rescaled_coef_z <-
      user_prior_z$prior_autoscale && has_predictors_z &&
      !is.na(user_prior_z$prior_dist_name) &&
      !all(user_prior_z$prior_scale == adjusted_prior_scale_z)
    rescaled_int <-
      user_prior_intercept$prior_autoscale_for_intercept && has_intercept &&
      !is.na(user_prior_intercept$prior_dist_name_for_intercept) &&
      (user_prior_intercept$prior_scale != adjusted_prior_intercept_scale)
    rescaled_int_z <-
      user_prior_intercept_z$prior_autoscale_for_intercept && has_intercept_z &&
      !is.na(user_prior_intercept_z$prior_dist_name_for_intercept) &&
      (user_prior_intercept_z$prior_scale != adjusted_prior_intercept_scale_z)
    
    
    if (has_predictors && user_prior$prior_dist_name %in% "t") {
      if (all(user_prior$prior_df == 1)) {
        user_prior$prior_dist_name <- "cauchy"
      } else {
        user_prior$prior_dist_name <- "student_t"
      }
    }
    if (has_predictors_z && user_prior_z$prior_dist_name %in% "t") {
      if (all(user_prior_z$prior_df == 1)) {
        user_prior_z$prior_dist_name <- "cauchy"
      } else {
        user_prior_z$prior_dist_name <- "student_t"
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
    if (has_intercept_z &&
        user_prior_intercept_z$prior_dist_name_for_intercept %in% "t") {
      if (all(user_prior_intercept_z$prior_df_for_intercept == 1)) {
        user_prior_intercept_z$prior_dist_name_for_intercept <- "cauchy"
      } else {
        user_prior_intercept_z$prior_dist_name_for_intercept <- "student_t"
      }
    }
    if (has_phi && user_prior_aux$prior_dist_name_for_aux %in% "t") {
      if (all(user_prior_aux$prior_df_for_aux == 1)) {
        user_prior_aux$prior_dist_name_for_aux <- "cauchy"
      } else {
        user_prior_aux$prior_dist_name_for_aux <- "student_t"
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
      prior_z = 
        if (!has_predictors_z) NULL else with(user_prior_z, list(
          dist = prior_dist_name,
          location = prior_mean,
          scale = prior_scale,
          adjusted_scale = if (rescaled_coef_z)
            adjusted_prior_scale_z else NULL,
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
      prior_intercept_z = 
        if (!has_intercept_z) NULL else with(user_prior_intercept_z, list(
          dist = prior_dist_name_for_intercept,
          location = prior_mean_for_intercept,
          scale = prior_scale_for_intercept,
          adjusted_scale = if (rescaled_int_z)
            adjusted_prior_intercept_scale_z else NULL,
          df = if (prior_dist_name_for_intercept %in% "student_t")
            prior_df_for_intercept else NULL
        )),
      prior_aux = 
      if (!has_phi) NULL else with(user_prior_aux, list(
          dist = prior_dist_name_for_aux,
          location = if (!is.na(prior_dist_name_for_aux) && 
                         prior_dist_name_for_aux != "exponential")
            prior_mean_for_aux else NULL,
          scale = if (!is.na(prior_dist_name_for_aux) && 
                      prior_dist_name_for_aux != "exponential")
            prior_scale_for_aux else NULL,
          df = if (!is.na(prior_dist_name_for_aux) && 
                   prior_dist_name_for_aux %in% "student_t")
            prior_df_for_aux else NULL, 
          rate = if (!is.na(prior_dist_name_for_aux) && 
                     prior_dist_name_for_aux %in% "exponential")
            1 / prior_scale_for_aux else NULL,
          aux_name = "phi"
        ))
    )
    return(prior_list)
  }
