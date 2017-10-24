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
  is_mvmer   <- is.mvmer(object)
  is_surv    <- is.surv(object)
  is_jm      <- is.jm(object)
  stub       <- if (is_jm) "Long" else "y"
  
  if (opt) {
    stop("Optimisation not implemented for stanmvreg objects.")
  } else {
    stan_summary <- make_stan_summary(stanfit)
    nms <- collect_nms(rownames(stan_summary), M, stub = get_stub(object))
    coefs <- list()
    ses <- list()
    
    # Coefs and SEs for longitudinal submodel(s)                    
    if (is_mvmer) {
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
    if (is_surv) {
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
    formula       = list_nms(object$formula, M, stub),
    terms         = list_nms(object$terms, M, stub),
    coefficients  = coefs, 
    ses           = ses,
    covmat        = covmat,
    prior.weights = object$weights, 
    prior.info    = object$prior.info,
    algorithm     = object$algorithm,
    call          = object$call,
    stan_function = object$stan_function,
    runtime       = if (object$algorithm == "sampling") times else NULL,
    stan_summary, stanfit
  )
  if (is_mvmer) {
    out$cnms      <- object$cnms
    out$flevels   <- object$flevels
    out$n_markers <- object$M
    out$n_grps    <- object$n_grps
    out$n_yobs    <- list_nms(object$n_yobs, M, stub)
    out$family    <- list_nms(object$family, M, stub)
    out$glmod     <- list_nms(object$glmod, M, stub)
    out$data      <- if (!is_jm) list_nms(object$data, M, stub) else NULL
    classes <- c("stanmvreg", "stanreg", "lmerMod")
  }
  if (is_jm) {
    out$id_var    <- object$id_var
    out$time_var  <- object$time_var
    out$n_subjects<- object$n_subjects
    out$n_events  <- sum(object$survmod$status > 0)
    out$eventtime <- object$survmod$eventtime
    out$status    <- object$survmod$status > 0
    out$basehaz   <- object$basehaz
    out$survmod   <- object$survmod
    out$qnodes    <- object$qnodes
    out$epsilon   <- object$epsilon    
    out$assoc     <- object$assoc
    out$assocmod  <- list_nms(object$assocmod, M, stub) 
    out$dataLong  <- list_nms(object$dataLong, M, stub) 
    out$dataEvent <- object$dataEvent
    out$grp_stuff <- object$grp_stuff
    out$fr        <- object$fr
    classes <- c("stanjm", classes)
  }
  out <- rm_null(out, recursive = FALSE)
  structure(out, class = classes)
}
