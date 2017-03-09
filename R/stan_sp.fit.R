# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2013, 2014, 2015, 2016 Trustees of Columbia University
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
stan_sp.fit <- function(y, x, w, ..., sp_model,
                        prior_rho = beta(), prior_intercept = normal(), 
                        algorithm = c("sampling", "optimizing", "meanfield", "fullrank"),
                        adapt_delta = NULL, QR = FALSE, sparse = FALSE) {
  
  # check prior for spatial correlation parameter
  if(is.null(prior_rho))
    prior_rho = beta(1,1)  # kinda hacky
  else if(prior_rho$dist != "beta")
    stop("'prior_rho' can only be specified with 'beta()'")
  
  algorithm <- match.arg(algorithm)
  
  x_stuff <- center_x(x, sparse = FALSE)
  for (i in names(x_stuff)) # xtemp, xbar, has_intercept
    assign(i, x_stuff[[i]])
  nvars <- ncol(xtemp)

  ok_intercept_dists <- nlist("normal", student_t = "t", "cauchy")
  prior_intercept_stuff <- handle_glm_prior(
    prior_intercept,
    nvars = 1,
    default_scale = 10,
    link = gaussian()$link,
    ok_dists = ok_intercept_dists
  )
  names(prior_intercept_stuff) <- paste0(names(prior_intercept_stuff),"_for_intercept")
  for (i in names(prior_intercept_stuff))
    assign(i, prior_intercept_stuff[[i]])
  
  if(has_intercept == TRUE) {
    has_intercept <- 1
    pars <- c("alpha", "beta", "rho", "sigma", "mean_PPD")
  }
  else {
    has_intercept <- 0
    pars <- c("beta", "rho", "sigma", "mean_PPD")
  }
  
  if (QR) {
    if (ncol(xtemp) <= 1)
      stop("'QR' can only be specified when there are multiple predictors.")
    if (sparse)
      stop("'QR' and 'sparse' cannot both be TRUE.")
    cn <- colnames(xtemp)
    decomposition <- qr(xtemp)
    sqrt_nm1 <- sqrt(nrow(xtemp) - 1L)
    Q <- qr.Q(decomposition)
    R_inv <- qr.solve(decomposition, Q) * sqrt_nm1
    xtemp <- Q * sqrt_nm1
    colnames(xtemp) <- cn
    xbar <- c(xbar %*% R_inv)
  }
  
  # in stan you need to separate the intercept from the rest of the parameters
  # then if intercept == TRUE {intercept = intercept - xbar * beta}
  # if intercept == FALSE {eta = eta + xbar * beta} ??? (not sure)

  standata <- nlist(N = nrow(xtemp),
                    K = ncol(xtemp),
                    W = W,
                    X = xtemp,
                    y = y,
                    mod = if(sp_model == "lagsarlm"){1}else if(sp_model == "errorsarlm"){2},
                    shape1 = prior_rho$alpha,
                    shape2 = prior_rho$beta,
                    has_intercept = has_intercept,
                    xbar = xbar,
                    prior_dist_for_intercept = prior_dist_for_intercept,
                    prior_mean_for_intercept = c(prior_mean_for_intercept),
                    prior_scale_for_intercept = c(prior_scale_for_intercept),
                    prior_df_for_intercept = c(prior_df_for_intercept)
                    )
  
  stanfit <- stanmodels$spatial
  
  if(algorithm == "optimizing") {
    out <- optimizing(stanfit, data = standata, draws = 1000, constrained = TRUE, ...)
    check_stanfit(out)
    new_names <- names(out$par)
    mark <- grepl("^beta\\[[[:digit:]]+\\]$", new_names)
    if (QR) {
      out$par[mark] <- R_inv %*% out$par[mark]
      out$theta_tilde[,mark] <- out$theta_tilde[, mark] %*% t(R_inv)
    }
    new_names[mark] <- colnames(xtemp)
    new_names[new_names == "alpha[1]"] <- "(Intercept)"
    names(out$par) <- new_names
    colnames(out$theta_tilde) <- new_names
    out$stanfit <- suppressMessages(sampling(stanfit, data = standata, 
                                             chains = 0))
    return(structure(out))
  }
  else {
    if (algorithm == "sampling") {
      sampling_args <- set_sampling_args(
        object = stanfit, 
        prior = prior_rho, 
        user_dots = list(...), 
        user_adapt_delta = adapt_delta, 
        data = standata, 
        pars = pars, 
        show_messages = FALSE)
      stanfit <- do.call(sampling, sampling_args)
    }
    else {  # meanfield or fullrank vb
      stanfit <- rstan::vb(stanfit, pars = pars, data = standata,
                           algorithm = algorithm, ...)
      if (algorithm == "meanfield" && !QR) 
        msg_meanfieldQR()
    }
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
    check_stanfit(stanfit)
    if(has_intercept == 1)
      new_names <- c("(Intercept)", colnames(xtemp), "rho", "sigma", "mean_PPD", "log-posterior")
    else
      new_names <- c(colnames(xtemp), "rho", "sigma", "mean_PPD", "log-posterior")
    
    stanfit@sim$fnames_oi <- new_names
    return(structure(stanfit))  #  prior.info = prior_info
  }
}