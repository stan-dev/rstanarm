# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2015, 2016, 2017 Trustees of Columbia University
# Copyright (C) 2016, 2017 Sam Brilleman
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

# Function to create a stanmvreg object (fitted model object)
#
# @param object A list returned by a call to any of: stan_jm, stan_mvmer
# @return A stanmvreg object
#
stanmvreg <- function(object) {
  
  opt        <- object$algorithm == "optimizing"
  stanfit    <- object$stanfit
  M          <- object$M
  mvmer      <- is.mvmer(object)
  surv       <- is.surv(object)
  jm         <- is.jm(object)
  stub       <- if (jm) "Long" else "y"
  
  if (opt) {
    stop("Optimisation not implemented for stanmvreg objects.")
  } else {
    stan_summary <- make_stan_summary(stanfit)
    nms <- collect_nms(rownames(stan_summary), M, stub = get_stub(object))
    coefs <- list()
    ses <- list()
    
    # Coefs and SEs for longitudinal submodel(s)                    
    if (!is.null(object$y_mod_stuff)) {
      y_coefs <- lapply(1:M, function(m)
        stan_summary[c(nms$y[[m]], nms$y_b[[m]]), select_median(object$algorithm)])
      y_stanmat <- lapply(1:M, function(m) 
        as.matrix(stanfit)[, c(nms$y[[m]], nms$y_b[[m]]), drop = FALSE])
      y_ses <- lapply(y_stanmat, function(m) apply(m, 2L, mad))
      y_covmat <- lapply(y_stanmat, cov)
      for (m in 1:M) {
        rownames(y_covmat[[m]]) <- colnames(y_covmat[[m]]) <- 
          rownames(stan_summary)[c(nms$y[[m]], nms$y_b[[m]])]
      }
      # Remove padding
      coefs[1:M] <- list_nms(lapply(y_coefs, unpad_reTrms.default), M, stub = stub)
      ses[1:M]   <- list_nms(lapply(y_ses, unpad_reTrms.default), M, stub = stub)
    }

    # Coefs and SEs for event submodel    
    if (!is.null(object$e_mod_stuff)) {
      e_coefs <- stan_summary[c(nms$e, nms$a), select_median(object$algorithm)]        
      if (length(e_coefs) == 1L) 
        names(e_coefs) <- rownames(stan_summary)[c(nms$e, nms$a)[1L]]
      e_stanmat <- as.matrix(stanfit)[, c(nms$e, nms$a), drop = FALSE]
      e_ses <- apply(e_stanmat, 2L, mad)    
      e_covmat <- cov(e_stanmat)
      rownames(e_covmat) <- colnames(e_covmat) <- 
        rownames(stan_summary)[c(nms$e, nms$a)]
      coefs$Event <- e_coefs
      ses$Event <- e_ses
    }
    
    # Covariance matrix for fixed effects                    
    stanmat <- as.matrix(stanfit)[, c(nms$alpha, nms$beta), drop = FALSE]
    covmat <- cov(stanmat)
    
    if (object$algorithm == "sampling") { # for MCMC fits only
      # Check Rhats for all parameters
      check_rhats(stan_summary[, "Rhat"])    
      # Run time (mins)
      times <- round((rstan::get_elapsed_time(object$stanfit))/60, digits = 1)
      times <- cbind(times, total = rowSums(times))      
    } 
  }

  out <- nlist(
    coefficients      = coefs, 
    ses               = ses,
    covmat            = covmat,
    formula           = object$formula,
    prior.weights     = object$weights, 
    na.action         = object$na.action,
    prior.info        = object$prior.info,
    algorithm         = object$algorithm,
    call              = object$call,
    stan_function     = object$stan_function,
    runtime           = if (object$algorithm == "sampling") times else NULL,
    stan_summary, 
    stanfit
  )
  if (mvmer) {
    out$family <- object$family
    out$cnms <- object$cnms
    out$n_markers <- object$M
    out$n_yobs <- object$n_yobs
    out$n_grps <- object$n_grps
    out$y <- object$y
    out$data <- if (!jm) object$data else NULL
    out$glmod <- fetch(object$y_mod_stuff, "mod")
    out$glmod_stuff <- object$y_mod_stuff
  }
  if (surv) {
    out$n_subjects <- object$n_subjects
    out$n_events <- sum(object$d > 0)
    out$eventtime <- object$eventtime
    out$status <- object$d > 0
    out$basehaz <- object$basehaz
    out$data <- if (!jm) object$data else NULL
    out$coxmod <- object$e_mod_stuff$mod    
    out$coxmod_stuff <- object$e_mod_stuff    
  }
  if (jm) {
    out$time_var <- object$time_var
    out$id_var <- object$id_var
    out$quadnodes <- object$quadnodes
    out$assoc <- object$assoc
    out$epsilon <- object$epsilon    
    out$dataLong <- object$dataLong
    out$dataEvent <- object$dataEvent
    out$assocmod_stuff <- object$a_mod_stuff  
    out$fr <- object$fr
  }
  out <- Filter(function(x) !is.null(x), out)
  structure(out, class = c("stanmvreg", "stanreg", "lmerMod"))
}
