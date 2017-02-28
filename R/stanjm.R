# Part of the rstanjm package
# Copyright (C) 2015, 2016 Trustees of Columbia University
# Copyright (C) 2016 Sam Brilleman
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

# Function to create a stanjm object (fitted model object)
#
# @param object A list returned by a call to stan_jm
# @return A stanjm object
#
stanjm <- function(object) {
  
  opt        <- object$algorithm == "optimizing"
  mer        <- rep(1L, object$M)
  stanfit    <- object$stanfit
  M          <- object$M

  if (opt) {
    stop("Optimisation not implemented for stan_jm")
  } else {
    stan_summary <- make_stan_summary(stanfit)
    nms <- collect_nms(rownames(stan_summary), M)
    
    # Coefs and SEs for longitudinal submodel(s)                    
    y_coefs <- lapply(1:M, function(m)
      stan_summary[c(nms$y[[m]], nms$y_b[[m]]), select_median(object$algorithm)])
    y_stanmat <- lapply(1:M, function(m) 
      as.matrix(stanfit)[, c(nms$y[[m]], nms$y_b[[m]]), drop = FALSE])
    y_ses <- lapply(y_stanmat, function(m) apply(m, 2L, mad))
    y_covmat <- lapply(y_stanmat, cov)
    for (m in 1:M) {
      rownames(y_covmat[[m]]) <- colnames(y_covmat[[m]]) <- rownames(stan_summary)[c(nms$y[[m]], nms$y_b[[m]])]
    }
 
    # Coefs and SEs for event submodel    
    e_coefs <- stan_summary[c(nms$e, nms$a), select_median(object$algorithm)]        
    if (length(e_coefs) == 1L) names(e_coefs) <- rownames(stan_summary)[c(nms$e, nms$a)[1L]]
    e_stanmat <- as.matrix(stanfit)[, c(nms$e, nms$a), drop = FALSE]
    e_ses <- apply(e_stanmat, 2L, mad)    
    e_covmat <- cov(e_stanmat)
    rownames(e_covmat) <- colnames(e_covmat) <- rownames(stan_summary)[c(nms$e, nms$a)]

    # Check Rhats for all parameters
    if (object$algorithm == "sampling") 
      check_rhats(stan_summary[, "Rhat"])
    
    # Covariance matrix for fixed effects                    
    stanmat <- as.matrix(stanfit)[, c(nms$alpha, nms$beta), drop = FALSE]
    covmat <- cov(stanmat)
  }
  
  # Linear predictor, fitted values
  y_eta <- lapply(1:M, function(m) linear_predictor.default(y_coefs[[m]], object$x[[m]], object$offset))
  y_mu  <- lapply(1:M, function(m) object$family[[m]]$linkinv(y_eta[[m]]))

  # Residuals
  y_tmp <- lapply(1:M, function(m) if (is.factor(object$y[[m]])) fac2bin(object$y[[m]]) else object$y[[m]])
  y_residuals <- lapply(1:M, function(m) y_tmp[[m]] - y_mu[[m]])

  # Observation labels
  y_nms      <- lapply(object$y, names)
  for (m in 1:M) {
    names(y_eta[[m]]) <- names(y_mu[[m]]) <- names(y_residuals[[m]]) <- y_nms[[m]]
  }

  # Remove padding
  y_coefs <- lapply(y_coefs, unpad_reTrms.default)
  y_ses   <- lapply(y_ses, unpad_reTrms.default)

  # Run time (mins)
  times <- round((rstan::get_elapsed_time(object$stanfit))/60, digits = 1)
  times <- cbind(times, total = rowSums(times))
  
  out <- nlist(
    coefficients      = list_nms(c(y_coefs, list(e_coefs)), M), 
    ses               = list_nms(c(y_ses, list(e_ses)), M),
    fitted.values     = list_nms(y_mu, M),
    linear.predictors = list_nms(y_eta, M),
    residuals         = list_nms(y_residuals, M), 
    covmat            = covmat,
    n_events          = sum(object$d > 0),
    n_markers         = object$M,
    n_subjects        = object$n_subjects,
    n_grps            = object$n_grps,
    n_yobs            = object$n_yobs,
    id_var            = object$id_var,
    time_var          = object$time_var,
    cnms              = object$cnms, 
    dataLong          = object$dataLong,
    dataEvent         = object$dataEvent,
    eventtime         = object$eventtime, 
    status            = object$d,   
    quadnodes         = object$quadnodes,
    y                 = object$y,
    fr                = list_nms(object$fr, M),
    #   offset = if (any(object$offset != 0)) object$offset else NULL,
    #   contrasts = object$contrasts, 
    prior.weights     = object$weights, 
    na.action         = object$na.action,
    call              = object$call, 
    formula           = list_nms(object$formula, M), 
    family            = list_nms(object$family, M), 
    assoc             = object$assoc,
    epsilon           = object$epsilon,
    basehaz           = object$basehaz,
    prior.info        = object$prior.info,
    algorithm         = object$algorithm,
    runtime           = times,
    glmod_stuff       = object$glmod_stuff,
    coxmod_stuff      = object$coxmod_stuff,    
    glmod             = object$glmod,
    coxmod            = object$coxmod,
    stan_summary, stanfit
  )
  
  structure(out, class = c("stanjm", "stanreg", "lmerMod"))
}
