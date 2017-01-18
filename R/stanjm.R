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
  family     <- object$family
  y          <- object$y
  x          <- object$x
  xq         <- object$xq
  xq_eps     <- object$xq_eps  
  e_x        <- object$e_x
  eventtime  <- object$eventtime
  d          <- object$d  
  dimensions <- object$dimensions
  y_cnms     <- object$y_cnms
  y_flist    <- object$y_flist
  y_nms      <- lapply(y, names)
  assoc      <- object$assoc
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
  y_eta <- lapply(1:M, function(m) linear_predictor.default(y_coefs[[m]], x[[m]], object$offset))
  y_mu  <- lapply(1:M, function(m) family[[m]]$linkinv(y_eta[[m]]))

  # Residuals
  y_tmp <- lapply(1:M, function(m) if (is.factor(y[[m]])) fac2bin(y[[m]]) else y[[m]])
  y_residuals <- lapply(1:M, function(m) y_tmp[[m]] - y_mu[[m]])
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
    coefficients = list_nms(c(y_coefs, list(e_coefs)), M), 
    ses = list_nms(c(y_ses, list(e_ses)), M),
    fitted.values = list_nms(y_mu, M),
    linear.predictors = list_nms(y_eta, M),
    residuals = list_nms(y_residuals, M), 
    covmat = covmat,
    n_markers = M,
    n_subjects = object$Npat,
    n_grps = object$n_grps,
    n_yobs = object$y_N,
    n_events = sum(d > 0),
    id_var = object$id_var,
    time_var = object$time_var,
    cnms = object$cnms, 
    #y_cnms = list_nms(y_cnms, M), 
    #y_flist = list_nms(y_flist, M),
    fr = list_nms(object$fr, M),
    x = list_nms(x, M), 
    xq = list_nms(xq, M), 
    xq_eps = if (sum(assoc$etaslope) || sum(assoc$muslope)) list_nms(xq_eps, M) else NULL,
    y = list_nms(y, M), 
    eventtime, 
    status = d,     
    quadpoints = object$quadpoints,     
    dataLong = object$dataLong,
    dataEvent = object$dataEvent,     
#    offset = if (any(object$offset != 0)) object$offset else NULL,
#    weights = object$weights, 
#    prior.weights = object$weights, 
#    contrasts = object$contrasts, 
    na.action = object$na.action,
    call = object$call, 
    formula = list_nms(object$formula, M), 
    family = list_nms(family, M), 
    base_haz = object$base_haz,
    assoc = assoc, # list_nms(assoc, M),
    quadnodes = object$quadnodes,
    algorithm = object$algorithm,
    prior.info = attr(stanfit, "prior.info"),
    times = times,
    stan_summary,  
    stanfit = if (opt) stanfit$stanfit else stanfit,
    glmod = object$glmod,
    coxmod = object$coxmod
  )
  if (opt) 
    out$asymptotic_sampling_dist <- stanmat
  
  structure(out, class = c("stanjm", "stanreg", "lmerMod"))
}


#-------- The following functions could be moved to misc.R

# Separates a names object into separate parts based on the longitudinal, 
# event, or association parameters.
# 
# @param x Character vector (often rownames(fit$stan_summary))
# @param M An integer specifying the number of longitudinal submodels.
# @param ... Arguments passed to grep
# @return A list with x separated out into those names corresponding
#   to parameters from the M longitudinal submodels, the event submodel
#   or association parameters.
collect_nms <- function(x, M, ...) {
  y <- lapply(1:M, function(m) grep(mod2rx(m), x, ...))
  y_extra <- lapply(1:M, function(m) 
    c(grep(paste0("^Long", m, "\\|sigma"), x, ...),
      grep(paste0("^Long", m, "\\|shape"), x, ...),
      grep(paste0("^Long", m, "\\|lambda"), x, ...),
      grep(paste0("^Long", m, "\\|overdispersion"), x, ...)))             
  y <- lapply(1:M, function(m) setdiff(y[[m]], y_extra[[m]]))
  e <- grep(mod2rx("^Event"), x, ...)     
  e_extra <- c(grep("^Event\\|weibull-shape|^Event\\|basehaz-coef", x, ...))         
  e <- setdiff(e, e_extra)
  a <- grep(mod2rx("^Assoc"), x, ...)
  b <- b_names2(x, ...)
  y_b <- lapply(1:M, function(m) b_names2(x, m, ...))
  alpha <- grep("^.{5}\\|\\(Intercept\\)", x, ...)      
  beta <- setdiff(c(unlist(y), e, a), alpha)  
  nlist(y, y_extra, y_b, e, e_extra, a, b, alpha, beta) 
}

# Grep for "b" parameters (ranef), can optionally be specified
# for a specific longitudinal submodel
#
# @param x Character vector (often rownames(fit$stan_summary))
# @param submodel Optional integer specifying which long submodel
# @param ... Passed to grep
b_names2 <- function(x, submodel = NULL, ...) {
  if (is.null(submodel)) {
    grep("^b\\[", x, ...)
  } else {
    grep(paste0("^b\\[Long", submodel, "\\|"), x, ...)
  }
}

# Converts "Long", "Event" or "Assoc" to the regular expression
# used at the start of variable names for the fitted joint model
#
# @param x The submodel for which the regular expression should be
#   obtained. Can be "Long", "Event", "Assoc", or an integer specifying
#   a specific longitudinal submodel.
mod2rx <- function(x) {
  if (x == "^Long") {
    c("^Long[1-9]\\|")
  } else if (x == "^Event") {
    c("^Event\\|")
  } else if (x == "^Assoc") {
    c("^Assoc\\|")
  } else if (x == "Long") {
    c("Long[1-9]\\|")
  } else if (x == "Event") {
    c("Event\\|")
  } else if (x == "Assoc") {
    c("Assoc\\|")
  } else {
    paste0("^Long", x, "\\|")
  }   
}