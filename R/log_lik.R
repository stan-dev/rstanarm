# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2015, 2016, 2017 Trustees of Columbia University
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

#' Pointwise log-likelihood matrix
#'
#' For models fit using MCMC only, the \code{log_lik} method returns the
#' \eqn{S} by \eqn{N} pointwise log-likelihood matrix, where \eqn{S} is the size
#' of the posterior sample and \eqn{N} is the number of data points, or in the
#' case of the \code{stanjm} method an \eqn{S} by \eqn{Npat} matrix where 
#' \eqn{Npat} is the number of individuals.
#'
#' @aliases log_lik
#' @export
#'
#' @templateVar stanregArg object
#' @template args-stanreg-object
#' @template args-dots-ignored
#' @param newdata An optional data frame of new data (e.g. holdout data) to use
#'   when evaluating the log-likelihood. See the description of \code{newdata}
#'   for \code{\link{posterior_predict}}.
#' @param offset A vector of offsets. Only required if \code{newdata} is
#'   specified and an \code{offset} was specified when fitting the model.
#'
#' @return For the \code{stanreg} method an \eqn{S} by \eqn{N} matrix, where 
#'   \eqn{S} is the size of the posterior sample and \eqn{N} is the number of 
#'   data points, or for the \code{stanjm} method an \eqn{S} by \eqn{Npat} 
#'   matrix where \eqn{Npat} is the number of individuals.
#'   
#'   
#' @examples 
#' \donttest{
#'  roaches$roach100 <- roaches$roach1 / 100
#'  fit <- stan_glm(
#'     y ~ roach100 + treatment + senior,
#'     offset = log(exposure2),
#'     data = roaches,
#'     family = poisson(link = "log"),
#'     prior = normal(0, 2.5),
#'     prior_intercept = normal(0, 10),
#'     iter = 500 # to speed up example
#'  )
#'  ll <- log_lik(fit)
#'  dim(ll)
#'  all.equal(ncol(ll), nobs(fit))
#'
#'  # using newdata argument
#'  nd <- roaches[1:2, ]
#'  nd$treatment[1:2] <- c(0, 1)
#'  ll2 <- log_lik(fit, newdata = nd, offset = c(0, 0))
#'  head(ll2)
#'  dim(ll2)
#'  all.equal(ncol(ll2), nrow(nd))
#' }
#'
log_lik.stanreg <- function(object, newdata = NULL, offset = NULL, ...) {
  if (!used.sampling(object))
    STOP_sampling_only("Pointwise log-likelihood matrix")
  newdata <- validate_newdata(newdata)
  calling_fun <- as.character(sys.call(-1))[1]
  args <- ll_args(object, newdata = newdata, offset = offset, 
                  reloo_or_kfold = calling_fun %in% c("kfold", "reloo"), 
                  ...)
  fun <- ll_fun(object)
  out <- vapply(
    seq_len(args$N),
    FUN = function(i) {
      as.vector(fun(
        i = i,
        data = args$data[i,, drop = FALSE],
        draws = args$draws
      ))
    },
    FUN.VALUE = numeric(length = args$S)
  )
  if (is.null(newdata)) colnames(out) <- rownames(model.frame(object))
  else colnames(out) <- rownames(newdata)
  return(out)
}

#' @rdname log_lik.stanreg
#' @export
#' @param newdataLong,newdataEvent Optional data frames containing new data 
#'   (e.g. holdout data) to use when evaluating the log-likelihood for a 
#'   model estimated using \code{\link{stan_jm}}. If the fitted model 
#'   was a multivariate joint model (i.e. more than one longitudinal outcome),
#'   then \code{newdataLong} is allowed to be a list of data frames. If supplying 
#'   new data, then \code{newdataEvent} should also include variables corresponding
#'   to the event time and event indicator. For more details, see the description 
#'   of \code{newdataLong} and \code{newdataEvent} for \code{\link{posterior_survfit}}.
#' 
log_lik.stanjm <- function(object, newdataLong = NULL, newdataEvent = NULL, ...) {
  if (!used.sampling(object))
    STOP_sampling_only("Pointwise log-likelihood matrix")
  validate_stanjm_object(object)
  M <- get_M(object)
  if (!identical(is.null(newdataLong), is.null(newdataEvent)))
    stop("Both newdataLong and newdataEvent must be supplied together.")
  if (!is.null(newdataLong)) {
    newdatas <- validate_newdatas(object, newdataLong, newdataEvent)
    newdataLong  <- newdatas[1:M]      
    newdataEvent <- newdatas[["Event"]]
  }
  pars <- extract_pars(object) # full array of draws
  data <- jm_data(object, newdataLong, newdataEvent) 
  calling_fun <- as.character(sys.call(-1))[1]
  reloo_or_kfold <- calling_fun %in% c("kfold", "reloo")
  val <- ll_jm(object, data, pars, reloo_or_kfold = reloo_or_kfold, ...)
  return(val)
}

# internal ----------------------------------------------------------------

# get log likelihood function for a particular model
# @param x stanreg object
# @return a function
ll_fun <- function(x, m = NULL) {
  validate_stanreg_object(x)
  f <- family(x, m = m)
  if (!is(f, "family") || is_scobit(x))
    return(.ll_polr_i)
  
  fun <- paste0(".ll_", f$family, "_i")
  get(fun, mode = "function")
}

# get arguments needed for ll_fun
# @param object stanreg object
# @param newdata same as posterior predict
# @param offset vector of offsets (only required if model has offset term and
#   newdata is specified)
# @param m Integer specifying which submodel for stanjm objects
# @param reloo_or_kfold logical. TRUE if ll_args is for reloo or kfold
# @param ... For models without group-specific terms (i.e., not stan_[g]lmer), 
#   if reloo_or_kfold is TRUE and 'newdata' is specified then ... is used to 
#   pass 'newx' and 'stanmat' from reloo or kfold (bypassing pp_data). This is a
#   workaround in case there are issues with newdata containing factors with
#   only a single level. Or for stanjm objects, then ... can be used to pass
#   'stanmat', which may be a matrix with a reduced number of draws (potentially
#   just a single MCMC draw).
# @return a named list with elements data, draws, S (posterior sample size) and
#   N = number of observations
ll_args <- function(object, ...) UseMethod("ll_args")
ll_args.stanreg <- function(object, newdata = NULL, offset = NULL, m = NULL, 
                            reloo_or_kfold = FALSE, ...) {
  validate_stanreg_object(object)
  f <- family(object, m = m)
  draws <- nlist(f)
  has_newdata <- !is.null(newdata)
  
  dots <- list(...)
  if (has_newdata && reloo_or_kfold && !is.mer(object)) {
    x <- dots$newx
    stanmat <- dots$stanmat
    y <- eval(formula(object)[[2L]], newdata)
  } else if (has_newdata) {
    ppdat <- pp_data(object, as.data.frame(newdata), offset = offset, m = m)
    tmp <- pp_eta(object, ppdat, m = m)
    eta <- tmp$eta
    stanmat <- tmp$stanmat
    x <- ppdat$x
    y <- eval(formula(object, m = m)[[2L]], newdata)
  } else {
    stanmat <- as.matrix.stanreg(object)
    x <- get_x(object, m = m)
    y <- get_y(object, m = m)
  }
  if (is.stanjm(object) && !is.null(dots$stanmat)) {
    stanmat <- dots$stanmat # potentially use a stanmat with a single draw
  }  
  
  if (is(f, "family") && !is_scobit(object)) {
    fname <- f$family
    if (!is.binomial(fname)) {
      data <- data.frame(y, x)
    } else {
      if (NCOL(y) == 2L) {
        trials <- rowSums(y)
        y <- y[, 1L]
      } else {
        trials <- 1
        if (is.factor(y)) 
          y <- fac2bin(y)
        stopifnot(all(y %in% c(0, 1)))
      }
      data <- data.frame(y, trials, x)
    }
    nms <- if (is.stanjm(object)) 
      collect_nms(colnames(stanmat), get_M(object)) else NULL  
    beta_sel <- if (is.null(nms)) seq_len(ncol(x)) else nms$y[[m]]
    draws$beta <- stanmat[, beta_sel, drop = FALSE]
    m_stub <- get_m_stub(m)
    if (is.gaussian(fname)) 
      draws$sigma <- stanmat[, paste0(m_stub, "sigma")]
    if (is.gamma(fname)) 
      draws$shape <- stanmat[, paste0(m_stub, "shape")]
    if (is.ig(fname)) 
      draws$lambda <- stanmat[, paste0(m_stub, "lambda")]
    if (is.nb(fname)) 
      draws$size <- stanmat[, paste0(m_stub, "reciprocal_dispersion")]
    if (is.beta(fname)) {
      draws$f_phi <- object$family_phi
      z_vars <- colnames(stanmat)[grepl("(phi)", colnames(stanmat))]
      if (length(z_vars) == 1 && z_vars == "(phi)") {
        draws$phi <- stanmat[, z_vars]
      } else {
        x_dat <- get_x(object)
        z_dat <- object$z
        colnames(x_dat) <- colnames(x_dat)
        colnames(z_dat) <- paste0("(phi)_", colnames(z_dat))
        data <- data.frame(y = get_y(object), cbind(x_dat, z_dat), check.names = FALSE)
        draws$phi <- stanmat[,z_vars]
      }
    }
  } else {
    stopifnot(is(object, "polr"))
    y <- as.integer(y)
    if (has_newdata) 
      x <- .validate_polr_x(object, x)
    data <- data.frame(y, x)
    draws$beta <- stanmat[, colnames(x), drop = FALSE]
    patt <- if (length(unique(y)) == 2L) "(Intercept)" else "|"
    zetas <- grep(patt, colnames(stanmat), fixed = TRUE, value = TRUE)
    draws$zeta <- stanmat[, zetas, drop = FALSE]
    draws$max_y <- max(y)
    if ("alpha" %in% colnames(stanmat)) {
      draws$alpha <- stanmat[, "alpha"]
      draws$f <- object$method
    }
  }
  
  data$offset <- if (has_newdata) offset else object$offset
  if (model_has_weights(object)) {
    if (is.stanjm(object)) 
      STOP_if_stanjm("posterior_survfit with weights")
    data$weights <- object$weights
  }
    
  if (is.mer(object)) {
    b_sel <- if (is.null(nms)) b_names(colnames(stanmat)) else nms$y_b[[m]]
    b <- stanmat[, b_sel, drop = FALSE]
    if (has_newdata) {
      Z_names <- ppdat$Z_names
      if (is.null(Z_names)) {
        b <- b[, !grepl("_NEW_", colnames(b), fixed = TRUE), drop = FALSE]
      } else {
        b <- pp_b_ord(b, Z_names)
      }
      if (is.null(ppdat$Zt)) z <- matrix(NA, nrow = nrow(x), ncol = 0)
      else z <- t(ppdat$Zt)
    } else {
      z <- get_z(object, m = m)
    }
    data <- cbind(data, as.matrix(z))
    draws$beta <- cbind(draws$beta, b)
  }
  
  nlist(data, draws, S = NROW(draws$beta), N = nrow(data))
}

# Alternative method for stanjm objects that allows data and pars to be
# passed directly, rather than constructed using pp_data within the ll_args
# method. This can be much faster when used in the MH algorithm within
# posterior_survfit, since it doesn't require repeated calls to pp_data.
#
# @param object A stanjm object
# @param data Output from jm_data
# @param pars Output from extract_pars
# @param m Integer specifying which submodel
# @param reloo_or_kfold logical. TRUE if ll_args is for reloo or kfold
ll_args.stanjm <- function(object, data, pars, m = 1, 
                           reloo_or_kfold = FALSE, ...) {
  validate_stanjm_object(object)
  if (model_has_weights(object))
    STOP_if_stanjm("posterior_survfit or log_lik with weights")
  f <- family(object, m = m)
  fname <- f$family
  draws <- nlist(f)
  stanmat <- pars$stanmat # potentially a stanmat with a single draw
  nms <- collect_nms(colnames(stanmat), get_M(object))
  y <- eval(formula(object, m = m)[[2L]], data$ndL[[m]])
  x <- data$yX[[m]]
  z <- t(data$yZt[[m]])
  if (!is.binomial(fname)) {
    dat <- data.frame(y, x)
  } else {
    if (NCOL(y) == 2L) {
      trials <- rowSums(y)
      y <- y[, 1L]
    } else {
      trials <- 1
      if (is.factor(y)) 
        y <- fac2bin(y)
      stopifnot(all(y %in% c(0, 1)))
    }
    dat <- data.frame(y, trials, x)
  }  
  dat <- cbind(dat, as.matrix(z))
  draws$beta <- stanmat[, nms$y[[m]], drop = FALSE]
  m_stub <- get_m_stub(m)
  if (is.gaussian(fname)) 
    draws$sigma <- stanmat[, paste0(m_stub, "sigma")]
  if (is.gamma(fname)) 
    draws$shape <- stanmat[, paste0(m_stub, "shape")]
  if (is.ig(fname)) 
    draws$lambda <- stanmat[, paste0(m_stub, "lambda")]
  if (is.nb(fname)) 
    draws$size <- stanmat[, paste0(m_stub, "reciprocal_dispersion")]
  b <- stanmat[, nms$y_b[[m]], drop = FALSE]
  b <- pp_b_ord(b, data$yZ_names[[m]])
  draws$beta <- cbind(draws$beta, b)
  nlist(data = dat, draws, S = NROW(draws$beta), N = nrow(dat))
}

# check intercept for polr models -----------------------------------------
# Check if a model fit with stan_polr has an intercept (i.e. if it's actually a 
# bernoulli model). If it doesn't have an intercept then the intercept column in
# x is dropped. This is only necessary if newdata is specified because otherwise
# the correct x is taken from the fitted model object.
.validate_polr_x <- function(object, x) {
  x0 <- get_x(object)
  has_intercept <- colnames(x0)[1L] == "(Intercept)" 
  if (!has_intercept && colnames(x)[1L] == "(Intercept)")
    x <- x[, -1L, drop = FALSE]
  x
}


# log-likelihood function helpers -----------------------------------------
.weighted <- function(val, w) {
  if (is.null(w)) {
    val
  } else {
    val * w
  } 
}

.xdata <- function(data) {
  sel <- c("y", "weights","offset", "trials")
  data[, -which(colnames(data) %in% sel)]
}
.mu <- function(data, draws) {
  eta <- as.vector(linear_predictor(draws$beta, .xdata(data), data$offset))
  draws$f$linkinv(eta)
}

# for stan_betareg only
.xdata_beta <- function(data) { 
  sel <- c("y", "weights","offset", "trials")
  data[, -c(which(colnames(data) %in% sel), grep("(phi)_", colnames(data), fixed = TRUE))]
}
.zdata_beta <- function(data) {
  sel <- c("y", "weights","offset", "trials")
  data[, grep("(phi)_", colnames(data), fixed = TRUE)]
}
.mu_beta <- function(data, draws) {
  eta <- as.vector(linear_predictor(draws$beta, .xdata_beta(data), data$offset))
  draws$f$linkinv(eta)
}
.phi_beta <- function(data, draws) {
  eta <- as.vector(linear_predictor(draws$phi, .zdata_beta(data), data$offset))
  draws$f_phi$linkinv(eta)
}

# log-likelihood functions ------------------------------------------------
.ll_gaussian_i <- function(i, data, draws) {
  val <- dnorm(data$y, mean = .mu(data,draws), sd = draws$sigma, log = TRUE)
  .weighted(val, data$weights)
}
.ll_binomial_i <- function(i, data, draws) {
  val <- dbinom(data$y, size = data$trials, prob = .mu(data,draws), log = TRUE)
  .weighted(val, data$weights)
}
.ll_poisson_i <- function(i, data, draws) {
  val <- dpois(data$y, lambda = .mu(data, draws), log = TRUE)
  .weighted(val, data$weights)
}
.ll_neg_binomial_2_i <- function(i, data, draws) {
  val <- dnbinom(data$y, size = draws$size, mu = .mu(data, draws), log = TRUE)
  .weighted(val, data$weights)
}
.ll_Gamma_i <- function(i, data, draws) {
  val <- dgamma(data$y, shape = draws$shape, 
                rate = draws$shape / .mu(data,draws), log = TRUE)
  .weighted(val, data$weights)
}
.ll_inverse.gaussian_i <- function(i, data, draws) {
  mu <- .mu(data, draws)
  val <- 0.5 * log(draws$lambda / (2 * pi)) - 
    1.5 * log(data$y) -
    0.5 * draws$lambda * (data$y - mu)^2 / 
    (data$y * mu^2)
  .weighted(val, data$weights)
}
.ll_polr_i <- function(i, data, draws) {
  eta <- linear_predictor(draws$beta, .xdata(data), data$offset)
  f <- draws$f
  J <- draws$max_y
  y_i <- data$y
  linkinv <- polr_linkinv(f)
  if (is.null(draws$alpha)) {
    if (y_i == 1) {
      val <- log(linkinv(draws$zeta[, 1] - eta))
    } else if (y_i == J) {
      val <- log1p(-linkinv(draws$zeta[, J-1] - eta))
    } else {
      val <- log(linkinv(draws$zeta[, y_i] - eta) - 
                   linkinv(draws$zeta[, y_i - 1L] - eta))
    }
  } else {
    if (y_i == 0) {
      val <- draws$alpha * log(linkinv(draws$zeta[, 1] - eta))
    } else if (y_i == 1) {
      val <- log1p(-linkinv(draws$zeta[, 1] - eta) ^ draws$alpha)
    } else {
      stop("Exponentiation only possible when there are exactly 2 outcomes.")
    }
  }
  .weighted(val, data$weights)
}
.ll_beta_i <- function(i, data, draws) {
  mu <- .mu_beta(data, draws)
  phi <- draws$phi
  if (length(grep("(phi)_", colnames(data), fixed = TRUE)) > 0) {
    phi <- .phi_beta(data, draws)
  }
  val <- dbeta(data$y, mu * phi, (1 - mu) * phi, log = TRUE)
  .weighted(val, data$weights)
}

# for stanjm objects only ------------------------------------------------

# Return log-likelihood for longitudinal submodel m
#
# @param object A stanjm object, or (when used in stan_jm function) a named list
#   with elements $basehaz, $family, $assoc
# @param data Output from jm_data
# @param pars Output from extract_pars
# @param m Integer specifying the longitudinal submodel
# @param user_b Parameters passed to ll_args and used to overide the selection
#   of b parameters which are normally taken from the stanmat draws
ll_long <- function(object, data, pars, m = 1, reloo_or_kfold = FALSE) {
  args <- ll_args(object, data, pars, m = m, reloo_or_kfold = reloo_or_kfold)
  fun  <- ll_fun(object, m = m)
  ll <- lapply(seq_len(args$N), function(j) as.vector(
    fun(i = j, data = args$data[j, , drop = FALSE], draws = args$draws)))
  ll <- do.call("cbind", ll)
  # return S * npat array by summing log-lik for y within each individual
  res <- t(apply(ll, 1L, function(row) tapply(row, data$flist[[m]], sum)))
  return(res) 
}

# Return survival probability or log-likelihood for event submodel
#
# @param object A stanjm object, or (when used in stan_jm function) a named list
#   with elements $basehaz, $family, $assoc
# @param data Output from jm_data
# @param pars Output from extract_pars
# @param one_draw A logical specifying whether the parameters provided in the 
#   pars argument are vectors for a single realisation of the parameter (e.g.
#   a single MCMC draw, or a posterior mean) (TRUE) or a stanmat array (FALSE)
# @param survprob A logical specifying whether to return the survival probability 
#   (TRUE) or the log likelihood for the event submodel (FALSE)
# @param An S by Npat matrix, or a length Npat vector, depending on the inputs
#   (where S is the size of the posterior sample and Npat is the number of 
#   individuals).
ll_event <- function(object, data, pars, one_draw = FALSE, survprob = FALSE) {
  basehaz <- object$basehaz
  family  <- object$family
  assoc   <- object$assoc          
  etimes  <- attr(data$assoc_parts, "etimes")
  estatus <- attr(data$assoc_parts, "estatus")
  qnodes  <- attr(data$assoc_parts, "qnodes")
  qtimes  <- attr(data$assoc_parts, "qtimes")
  qwts    <- attr(data$assoc_parts, "qwts") 
  times   <- c(etimes, qtimes)
  # To avoid an error in log(times) replace times equal to zero with a small 
  # non-zero value. Note that these times correspond to individuals where the,
  # event time (etimes) was zero, and therefore the cumhaz (at baseline) will 
  # be forced to zero for these individuals further down in the code anyhow.  
  times[times == 0] <- 1E-10 
  # Linear predictor for the event submodel
  e_eta <- linear_predictor(pars$ebeta, data$eXq) 
  if (one_draw) {
    aXq <- make_assoc_terms(parts = data$assoc_parts, assoc = assoc, 
                            family = family, beta = pars$beta, b = pars$b)
    e_eta <- e_eta + linear_predictor.default(pars$abeta, aXq)
  } else {
    aXq <- matrix(NA, NROW(data$eXq), NCOL(pars$abeta))
    for (s in 1:NROW(e_eta)) {
      abeta_s <- pars$abeta[s,]
      beta_s  <- lapply(pars$beta, function(x) x[s,])
      b_s     <- lapply(pars$b,    function(x) x[s,])
      aXq_s   <- make_assoc_terms(parts = data$assoc_parts, assoc = assoc, 
                                  family = family, beta = beta_s, b = b_s)
      e_eta[s,] <- e_eta[s,] + linear_predictor.default(abeta_s, aXq_s)
    }
  }
  # Baseline hazard
  if (basehaz$type_name == "weibull") { # pars$bhcoef == weibull shape
    log_basehaz <- as.vector(log(pars$bhcoef)) + 
      linear_predictor(pars$bhcoef - 1, log(times))
  } else if (basehaz$type_name == "bs") { # pars$bhcoef == spline coefs
    log_basehaz <- linear_predictor(pars$bhcoef, predict(basehaz$bs_basis, times))
  } else {
    stop("Not yet implemented for basehaz = ", basehaz$type_name)
  }  
  loghaz <- log_basehaz + e_eta # log haz at etimes (if not NULL) and qtimes  
  # Calculate survival prob or log_lik  
  if (one_draw) {
    qhaz <- tail(exp(loghaz), length(qtimes)) # haz at qtimes
    qwhaz <- qwts * qhaz
    splitting_vec <- rep(1:qnodes, each = data$Npat)
    cumhaz <- Reduce('+', split(qwhaz, splitting_vec))
    cumhaz[etimes == 0] <- 0 # force cumhaz to zero if last_time is zero
  } else {
    qhaz <- exp(loghaz[, tail(1:ncol(loghaz), length(qtimes)), drop = FALSE])
    qwhaz <- t(apply(qhaz, 1L, function(row) qwts * row))
    cumhaz <- Reduce('+', array2list(qwhaz, nsplits = qnodes))  
    cumhaz[, etimes == 0] <- 0 # force cumhaz to zero if last_time is zero
  }
  ll_survt <- -cumhaz
  if (survprob) { # return surv prob at time t (upper limit of integral)
    return(exp(ll_survt)) 
  } else { # return log_lik at event time
    if (is.null(etimes) || is.null(estatus))
      stop("'etimes' and 'estatus' cannot be NULL if 'survprob = FALSE'.")
    if (one_draw) { # return vector of length npat
      if (!length(loghaz) == (length(c(etimes, qtimes))))
        stop("Bug found: length of loghaz vector appears to be incorrect.")
      return(estatus * head(loghaz, length(etimes)) + ll_survt)
    } else { # return S * npat matrix
      if (!ncol(loghaz) == (length(c(etimes, qtimes))))
        stop("Bug found: number of cols in loghaz matrix appears to be incorrect.")
      eloghaz <- loghaz[, 1:length(etimes), drop = FALSE]
      ll_hazt <- t(apply(eloghaz, 1L, function(row) estatus * row))
      return(ll_hazt + ll_survt)
    }
  }
} 

# Return log likelihood for full joint model
#
# @param object A stanjm object, or (when used in stan_jm function) a named list
#   with elements $basehaz, $family, $assoc
# @param data Output from jm_data
# @param pars Output from extract_pars
# @param include_b A logical, if TRUE then the log likelihood for the random 
#   effects distribution is also included in the log likelihood calculation.
# @param sum A logical. If TRUE then the log likelihood is summed across all
#   individuals. If FALSE then the log likelihood is returned for each
#   individual (either as an S * Npat matrix, or a length Npat vector, depending
#   on the type of inputs to the pars argument).
# @param ... Arguments passed to ll_long. Can include 'reloo_or_kfold' which is
#   a logical specifying whether the function calling ll_jm was reloo or kfold.
# @return Either a matrix, a vector or a scalar, depending on the input types
#   and whether sum is set to TRUE.
ll_jm <- function(object, data, pars, include_b = FALSE, sum = FALSE, ...) {
  M <- get_M(object)
  # log-lik for longitudinal submodels
  ll_long <- lapply(1:M, function(m) ll_long(object, data, pars, m = m, ...))
  # log-lik for event submodel
  ll_event <- ll_event(object, data, pars)
  # log-lik for random effects model
  if (include_b) {
    if (length(object$cnms) > 1L)
      stop("Bug found: 'include_b' can only be TRUE when there is one grouping factor.")
    if ((data$Npat > 1) || (nrow(pars$stanmat) > 1L))
      stop("Bug found: 'include_b' can only be TRUE when 'data' is for one ",
           "individual, and stanmat is for a single draw.")
    id_var <- object$id_var
    nms    <- unlist(lapply(data$assoc_parts, function(x) x$mod_eta$Z_names))
    b      <- do.call("cbind", pars$b)
    b      <- as.vector(pp_b_ord(b, nms))
    mu     <- rep(0, length(b))
    Sigma  <- VarCorr(object, stanmat = pars$stanmat)[[id_var]]
    ll_b   <- mvtnorm::dmvnorm(b, mean = mu, sigma = Sigma, log = TRUE)
  } else ll_b <- NULL
  # check the dimensions of the various components
  if (is.matrix(ll_event)) { # S * Npat matrices
    mats <- unique(sapply(c(ll_long, list(ll_event)), is.matrix))
    dims <- unique(lapply(c(ll_long, list(ll_event)), dim)) 
    if ((length(dims) > 1L) || (length(mats) > 1L))
      stop("Bug found: elements of 'll_long' should be same class and ",
           "dimension as 'll_event'.")
    if (include_b && !identical(length(ll_b), ncol(ll_event)))
      stop("Bug found: length of 'll_b' should be equal to the number of ",
           "columns in 'll_event'.")
  } else { # length Npat vectors (ie, log-lik based on a single draw of pars)
    lens <- unique(sapply(c(ll_long, list(ll_event)), length))
    if (length(lens) > 1L)
      stop("Bug found: elements of 'll_long' should be same length as 'll_event'.")
    if (include_b && !identical(length(ll_b), length(ll_event)))
      stop("Bug found: length of 'll_b' should be equal to length of 'll_event'.")
  }  
  # sum the various components (long + event + random effects)
  val <- Reduce('+', c(ll_long, list(ll_event)))
  if (include_b && is.matrix(val)) {
    val <- t(apply(val, 1L, function(row) row + ll_b)) 
  } else if (include_b && is.vector(val)) {
    val <- val + ll_b
  } 
  # return log-lik for joint model
  if (!sum) 
    return(val)             # S * Npat matrix or length Npat vector
  else if (is.matrix(val)) 
    return(rowSums(val))    # length S vector
  else 
    return(sum(val))        # scalar 
}
