# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2015, 2016, 2017 Trustees of Columbia University
# Copyright (C) 2016, 2017, 2018 Sam Brilleman
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

# Function to create a stansurv object (fitted model object)
#
# @param object A list returned by a call to stan_surv
# @return A stansurv object
#
stansurv <- function(object) {
  
  alg     <- object$algorithm
  opt     <- alg == "optimizing"
  mcmc    <- alg == "sampling"
  stanfit <- object$stanfit
  basehaz <- object$basehaz
  K       <- NCOL(object$x)

  if (opt)
    stop2("Optimisation not implemented for 'stansurv' objects.")
  
  stan_summary <- make_stan_summary(stanfit)
  
  # number of parameters
  nvars  <- ncol(object$x) + has_intercept(basehaz) + basehaz$nvars
  
  # obtain medians
  coefs     <- stan_summary[seq(nvars), select_median(alg)] 
  coefs_nms <- rownames(stan_summary)[seq(nvars)]
  names(coefs) <- coefs_nms # ensure parameter names are retained
  
  # obtain standard errors and covariance matrix
  stanmat <- as.matrix(stanfit)[, seq(nvars), drop = FALSE]
  colnames(stanmat) <- coefs_nms   
  ses <- apply(stanmat, 2L, mad)
  covmat <- cov(stanmat)
  
  # for mcmc only
  if (mcmc) { 
    check_rhats(stan_summary[, "Rhat"])    # check rhats for all parameters
    runtime <- get_runtime(object$stanfit) # run time (in mins)
  }
  
  # return object of class 'stansurv'
  out <- nlist(
    coefficients  = coefs, 
    ses           = ses,
    covmat        = covmat,
    formula       = object$formula,
    has_tde       = object$has_tde,
    has_quadrature= object$has_quadrature,
    terms         = object$terms,
    data          = object$data,
    model_frame   = object$model_frame,
    x             = object$x,
    s_cpts        = object$s_cpts,
    entrytime     = object$t_beg, 
    eventtime     = object$t_end, 
    event         = object$event,      
    delayed       = object$delayed,    
    basehaz       = object$basehaz,
    nobs          = object$nobs,
    nevents       = object$nevents,
    nlcens        = object$nlcens,
    nrcens        = object$nrcens,
    nicens        = object$nicens,
    ncensor       = object$ncensor,
    ndelayed      = object$ndelayed,
    qnodes        = object$qnodes,
    prior.info    = object$prior_info,
    algorithm     = object$algorithm,
    stan_function = object$stan_function,
    call          = object$call,
    runtime       = if (mcmc) runtime else NULL,
    rstan_version    = utils::packageVersion("rstan"),
    rstanarm_version = utils::packageVersion("rstanarm"),
    stan_summary, 
    stanfit
  )
  out <- rm_null(out, recursive = FALSE)
  
  structure(out, class = c("stansurv", "stanreg"))
}


#---------- internal

# Return the model fitting time in minutes.
#
# @param stanfit An object of class 'stanfit'.
# @return A matrix of runtimes, stratified by warmup/sampling and chain/overall.
get_runtime <- function(stanfit) {
  tt <- rstan::get_elapsed_time(stanfit)
  tt <- round(tt / 60, digits = 1L)    # time per chain
  tt <- cbind(tt, total = rowSums(tt)) # time per chain & overall  
}
