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
                        prior_rho = beta(), prior_intercept, 
                        algorithm = c("sampling", "optimizing", "meanfield", "fullrank"),
                        adapt_delta = NULL, QR = FALSE, sparse = FALSE) {
  
  algorithm <- match.arg(algorithm)
  
  x_stuff <- center_x(x, sparse = FALSE)
  for (i in names(x_stuff)) # xtemp, xbar, has_intercept
    assign(i, x_stuff[[i]])
  nvars <- ncol(xtemp)
  
  if(has_intercept == TRUE) {
    has_intercept <- 1
    pars <- c("rho", "alpha", "beta", "mean_PPD")
  }
  else {
    has_intercept <- 0
    pars <- c("rho", "beta", "mean_PPD")
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
                    xbar = xbar
                    )
  
  stanfit <- stanmodels$spatial
  
  if(algorithm == "optimizing") {
    out <- optimizing(stanfit, data = standata, draws = 1000, constrained = TRUE, ...)
    check_stanfit(out)
  }
  else if(algorithm == "sampling") {
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
  else {  # algorithm == meanfield or fullrank
    stop("meanfield and fullrank algorithm not yet implemented.")
  }
  check_stanfit(stanfit)
  
  if(has_intercept == 1)
    new_names <- c("rho", "(Intercept)", colnames(xtemp), "mean_PPD", "log-posterior")
  else
    new_names <- c("rho", colnames(xtemp), "mean_PPD", "log-posterior")

  stanfit@sim$fnames_oi <- new_names
  return(structure(stanfit))  #  prior.info = prior_info
}