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
                        adapt_delta) {
  
  algorithm <- match.arg(algorithm)
  
  x_stuff <- center_x(x, sparse = FALSE)
  for (i in names(x_stuff)) # xtemp, xbar, has_intercept
    assign(i, x_stuff[[i]])
  nvars <- ncol(xtemp)

  standata <- nlist(N = nrow(X),
                    K = ncol(X),
                    W = W,
                    X = X,
                    y = y,
                    mod = if(sp_model == "lagsarlm"){1}else{2},
                    shape1 = prior_rho$alpha,
                    shape2 = prior_rho$beta
                    )
  
  stanfit <- stanmodels$spatial
  
  pars <- c("rho", "beta", "mean_PPD")
  
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
    stop("meanfield and fullrank algorithm not available.")
  }
  check_stanfit(stanfit)

  new_names <- c("rho", "(Intercept)", colnames(xtemp), "mean_PPD", "log-posterior")

  stanfit@sim$fnames_oi <- new_names
  return(structure(stanfit))  #  prior.info = prior_info
}