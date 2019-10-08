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
#' case of the \code{stanmvreg} method (when called on \code{\link{stan_jm}}
#' model objects) an \eqn{S} by \eqn{Npat} matrix where \eqn{Npat} is the number 
#' of individuals.
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
#' @return For the \code{stanreg} and \code{stanmvreg} methods an \eqn{S} by 
#'   \eqn{N} matrix, where \eqn{S} is the size of the posterior sample and 
#'   \eqn{N} is the number of data points. For the \code{stanjm} method 
#'   an \eqn{S} by \eqn{Npat} matrix where \eqn{Npat} is the number of individuals.
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
#'     iter = 500, # just to speed up example,
#'     refresh = 0
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
  newdata <- validate_newdata(object, newdata, m = NULL)
  calling_fun <- as.character(sys.call(-1))[1]
  dots <- list(...)
  if (is.stanmvreg(object)) {
    m <- dots[["m"]]
    if (is.null(m)) 
      STOP_arg_required_for_stanmvreg(m)
    if (!is.null(offset))
      stop2("'offset' cannot be specified for stanmvreg objects.")
  } else {
    m <- NULL
  }
  
  newdata <- validate_newdata(object, newdata = newdata, m = m)
  args <- ll_args.stanreg(object, newdata = newdata, offset = offset, 
                          reloo_or_kfold = calling_fun %in% c("kfold", "reloo"), 
                          ...)
  fun <- ll_fun(object, m = m)
  if (is_clogit(object)) {
    out <-
      vapply(
        seq_len(args$N),
        FUN.VALUE = numeric(length = args$S),
        FUN = function(i) {
          as.vector(fun(
            draws = args$draws,
            data_i = args$data[args$data$strata ==
                               levels(args$data$strata)[i], , drop = FALSE]
          ))
        }
      )
    return(out)
  } else {
    out <- vapply(
      seq_len(args$N),
      FUN = function(i) {
        as.vector(fun(
          data_i = args$data[i, , drop = FALSE],
          draws = args$draws
        ))
      },
      FUN.VALUE = numeric(length = args$S)
    )
  }
  if (is.null(newdata)) colnames(out) <- rownames(model.frame(object, m = m))
  else colnames(out) <- rownames(newdata)
  return(out)
}

#' @rdname log_lik.stanreg
#' @export
#' @templateVar mArg m
#' @template args-m
#'  
log_lik.stanmvreg <- function(object, m = 1, newdata = NULL, ...) {
  validate_stanmvreg_object(object)
  out <- log_lik.stanreg(object, newdata = newdata, m = m, ...)
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
#'   to the event time and event indicator as these are required for evaluating the
#'   log likelihood for the event submodel. For more details, see the description 
#'   of \code{newdataLong} and \code{newdataEvent} for \code{\link{posterior_survfit}}.
#' 
log_lik.stanjm <- function(object, newdataLong = NULL, newdataEvent = NULL, ...) {
  if (!used.sampling(object))
    STOP_sampling_only("Pointwise log-likelihood matrix")
  validate_stanjm_object(object)
  M <- get_M(object)
  if ("m" %in% names(list(...)))
    stop("'m' should not be specified for stan_jm objects since the ",
         "log-likelihood is calculated for the full joint model.")
  if (!identical(is.null(newdataLong), is.null(newdataEvent)))
    stop("Both newdataLong and newdataEvent must be supplied together.")
  if (!is.null(newdataLong)) {
    newdatas <- validate_newdatas(object, newdataLong, newdataEvent)
    newdataLong  <- newdatas[1:M]      
    newdataEvent <- newdatas[["Event"]]
  }
  pars <- extract_pars(object) # full array of draws
  data <- .pp_data_jm(object, newdataLong, newdataEvent) 
  calling_fun <- as.character(sys.call(-1))[1]
  reloo_or_kfold <- calling_fun %in% c("kfold", "reloo")
  val <- .ll_jm(object, data, pars, reloo_or_kfold = reloo_or_kfold, ...)
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
  else if (is_clogit(x)) 
    return(.ll_clogit_i)
  else if (is.nlmer(x)) 
    return(.ll_nlmer_i)
  
  fun <- paste0(".ll_", family(x, m = m)$family, "_i")
  get(fun, mode = "function")
}

# get arguments needed for ll_fun
# @param object stanreg object
# @param newdata same as posterior predict
# @param offset vector of offsets (only required if model has offset term and
#   newdata is specified)
# @param m Integer specifying which submodel for stanmvreg objects
# @param reloo_or_kfold logical. TRUE if ll_args is for reloo or kfold
# @param ... For models without group-specific terms (i.e., not stan_[g]lmer), 
#   if reloo_or_kfold is TRUE and 'newdata' is specified then ... is used to 
#   pass 'newx' and 'stanmat' from reloo or kfold (bypassing pp_data). This is a
#   workaround in case there are issues with newdata containing factors with
#   only a single level. Or for stanmvreg objects, then ... can be used to pass
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
  
  z_betareg <- NULL
  if (has_newdata && reloo_or_kfold && !is.mer(object)) {
    x <- dots$newx
    z_betareg <- dots$newz # NULL except for some stan_betareg models
    if (!is.null(z_betareg)) {
      z_betareg <- as.matrix(z_betareg)
    }
    stanmat <- dots$stanmat
    form <- as.formula(formula(object)) # in case formula is string
    y <- eval(form[[2L]], newdata)
  } else if (has_newdata) {
    ppdat <- pp_data(object, as.data.frame(newdata), offset = offset, m = m)
    pp_eta_dat <- pp_eta(object, ppdat, m = m)
    eta <- pp_eta_dat$eta
    stanmat <- pp_eta_dat$stanmat
    z_betareg <- ppdat$z_betareg
    x <- ppdat$x
    form <- as.formula(formula(object, m = m))
    y <- eval(form[[2L]], newdata)
  } else {
    stanmat <- as.matrix.stanreg(object)
    x <- get_x(object, m = m)
    y <- get_y(object, m = m)
  }
  if (is.stanmvreg(object) && !is.null(dots$stanmat)) {
    stanmat <- dots$stanmat # potentially use a stanmat with a single draw
  }  
  
  if (!is_polr(object)) { # not polr or scobit model
    fname <- f$family
    if (is.nlmer(object)) {
      draws <- list(mu = posterior_linpred(object, newdata = newdata),
                    sigma = stanmat[,"sigma"])
      data <- data.frame(y)
      data$offset <- if (has_newdata) offset else object$offset
      if (model_has_weights(object)) {
        data$weights <- object$weights
      }
      data$i_ <- seq_len(nrow(data))  # for nlmer need access to i inside .ll_nlmer_i
      return(nlist(data, draws, S = NROW(draws$mu), N = nrow(data)))
      
    } else if (!is.binomial(fname)) {
      data <- data.frame(y, x)
      if (!is.null(z_betareg)) {
        data <- cbind(data, z_betareg)
      }
    } else {
      if (NCOL(y) == 2L) {
        trials <- rowSums(y)
        y <- y[, 1L]
      } else if (is_clogit(object)) {
        if (has_newdata) strata <- eval(object$call$strata, newdata)
        else strata <- model.frame(object)[,"(weights)"]
        strata <- as.factor(strata)
        successes <- aggregate(y, by = list(strata), FUN = sum)$x
        formals(draws$f$linkinv)$g <- strata
        formals(draws$f$linkinv)$successes <- successes
        trials <- 1L
      } else {
        trials <- 1
        if (is.factor(y)) 
          y <- fac2bin(y)
        stopifnot(all(y %in% c(0, 1)))
      }
      data <- data.frame(y, trials, x)
    }
    nms <- if (is.stanmvreg(object)) 
      collect_nms(colnames(stanmat),
                  M = get_M(object), 
                  stub = get_stub(object)) else NULL  
    beta_sel <- if (is.null(nms)) seq_len(ncol(x)) else nms$y[[m]]
    draws$beta <- stanmat[, beta_sel, drop = FALSE]
    m_stub <- get_m_stub(m, stub = get_stub(object))
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
        if (has_newdata) {
          if (!is.null(z_betareg)) {
          colnames(data) <- c("y", colnames(get_x(object)), 
                              paste0("(phi)_", colnames(z_betareg)))
          }
        } else {
          x_dat <- get_x(object)
          z_dat <- as.matrix(object$z)
          colnames(x_dat) <- colnames(x_dat)
          colnames(z_dat) <- paste0("(phi)_", colnames(z_dat))
          data <- data.frame(y = get_y(object), cbind(x_dat, z_dat), check.names = FALSE)
        }
        draws$phi <- stanmat[,z_vars]
      }
    }
  } else {
    stopifnot(is_polr(object))
    y <- as.integer(y)
    if (has_newdata) {
      x <- .validate_polr_x(object, x)
    }
    data <- data.frame(y, x)
    draws$beta <- stanmat[, colnames(x), drop = FALSE]
    patt <- if (length(unique(y)) == 2L) "(Intercept)" else "|"
    zetas <- grep(patt, colnames(stanmat), fixed = TRUE, value = TRUE)
    draws$zeta <- stanmat[, zetas, drop = FALSE]
    draws$max_y <- max(y)
    if ("alpha" %in% colnames(stanmat)) { 
      stopifnot(is_scobit(object))
      # scobit
      draws$alpha <- stanmat[, "alpha"]
      draws$f <- object$method
    }
  }
  
  data$offset <- if (has_newdata) offset else object$offset
  if (model_has_weights(object)) {
    if (is.stanmvreg(object)) 
      STOP_if_stanmvreg("posterior_survfit with weights")
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
    data <- cbind(data, as.matrix(z)[1:NROW(x),, drop = FALSE])
    draws$beta <- cbind(draws$beta, b)
  }
  
  if (is_clogit(object)) {
    data$strata <- strata
    out <- nlist(data, draws, S = NROW(draws$beta), N = nlevels(strata))
  } else {
    out <- nlist(data, draws, S = NROW(draws$beta), N = nrow(data)) 
  }
  return(out)
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
  sel <- c("y", "weights","offset", "trials","strata")
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
.ll_gaussian_i <- function(data_i, draws) {
  val <- dnorm(data_i$y, mean = .mu(data_i, draws), sd = draws$sigma, log = TRUE)
  .weighted(val, data_i$weights)
}
.ll_binomial_i <- function(data_i, draws) {
  val <- dbinom(data_i$y, size = data_i$trials, prob = .mu(data_i, draws), log = TRUE)
  .weighted(val, data_i$weights)
}
.ll_clogit_i <- function(data_i, draws) {
  eta <- linear_predictor(draws$beta, .xdata(data_i), data_i$offset)
  denoms <- apply(eta, 1, log_clogit_denom, N_j = NCOL(eta), D_j = sum(data_i$y))
  rowSums(eta[,data_i$y == 1, drop = FALSE] - denoms)
}
.ll_poisson_i <- function(data_i, draws) {
  val <- dpois(data_i$y, lambda = .mu(data_i, draws), log = TRUE)
  .weighted(val, data_i$weights)
}
.ll_neg_binomial_2_i <- function(data_i, draws) {
  val <- dnbinom(data_i$y, size = draws$size, mu = .mu(data_i, draws), log = TRUE)
  .weighted(val, data_i$weights)
}
.ll_Gamma_i <- function(data_i, draws) {
  val <- dgamma(data_i$y, shape = draws$shape, 
                rate = draws$shape / .mu(data_i,draws), log = TRUE)
  .weighted(val, data_i$weights)
}
.ll_inverse.gaussian_i <- function(data_i, draws) {
  mu <- .mu(data_i, draws)
  val <- 0.5 * log(draws$lambda / (2 * pi)) - 
    1.5 * log(data_i$y) -
    0.5 * draws$lambda * (data_i$y - mu)^2 / 
    (data_i$y * mu^2)
  .weighted(val, data_i$weights)
}
.ll_polr_i <- function(data_i, draws) {
  eta <- linear_predictor(draws$beta, .xdata(data_i), data_i$offset)
  f <- draws$f
  J <- draws$max_y
  y_i <- data_i$y
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
  .weighted(val, data_i$weights)
}
.ll_beta_i <- function(data_i, draws) {
  mu <- .mu_beta(data_i, draws)
  phi <- draws$phi
  if (length(grep("(phi)_", colnames(data_i), fixed = TRUE)) > 0) {
    phi <- .phi_beta(data_i, draws)
  }
  val <- dbeta(data_i$y, mu * phi, (1 - mu) * phi, log = TRUE)
  .weighted(val, data_i$weights)
}
.ll_nlmer_i <- function(data_i, draws) {
  i_ <- data_i$i_
  val <- dnorm(data_i$y, mean = draws$mu[, i_], sd = draws$sigma, log = TRUE)
  .weighted(val, data_i$weights)
}

# log-likelihood functions for stanjm objects only ----------------------

# Alternative ll_args method for stanjm objects that allows data and pars to be
# passed directly, rather than constructed using pp_data within the ll_args
# method. This can be much faster when used in the MH algorithm within
# posterior_survfit, since it doesn't require repeated calls to pp_data.
#
# @param object A stanmvreg object
# @param data Output from .pp_data_jm
# @param pars Output from extract_pars
# @param m Integer specifying which submodel
# @param reloo_or_kfold logical. TRUE if ll_args is for reloo or kfold
ll_args.stanjm <- function(object, data, pars, m = 1, 
                           reloo_or_kfold = FALSE, ...) {
  validate_stanjm_object(object)
  if (model_has_weights(object))
    STOP_if_stanmvreg("posterior_survfit or log_lik with weights")
  f <- family(object, m = m)
  fname <- f$family
  draws <- nlist(f)
  stanmat <- pars$stanmat # potentially a stanmat with a single draw
  nms <- collect_nms(colnames(stanmat), get_M(object))
  if (is.jm(object)) {
    # for stan_jm models, log_lik is evaluated for the full
    # joint model, so data contains info on all submodels
    y <- data$y[[m]]
    x <- data$yX[[m]]
    z <- t(data$yZt[[m]])
    Z_names <- data$yZ_names[[m]]    
  } else { 
    # for stan_mvmer models, log_lik is only ever called for
    # one submodel at a time, so data is for one submodel
    y <- data$y
    x <- data$X
    z <- t(data$Zt)
    Z_names <- data$Z_names
  }
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
  b <- pp_b_ord(b, Z_names)
  draws$beta <- cbind(draws$beta, b)
  nlist(data = dat, draws, S = NROW(draws$beta), N = nrow(dat))
}

# Return log likelihood for full joint model
#
# @param object A stanmvreg object, or (when used in stan_jm function) a named list
#   with elements $basehaz, $family, $assoc
# @param data Output from .pp_data_jm
# @param pars Output from extract_pars
# @param include_long A logical, if TRUE then the log likelihood for the  
#   longitudinal submodels are included in the log likelihood calculation.
# @param include_b A logical, if TRUE then the log likelihood for the random 
#   effects distribution is also included in the log likelihood calculation.
# @param sum A logical. If TRUE then the log likelihood is summed across all
#   individuals. If FALSE then the log likelihood is returned for each
#   individual (either as an S * Npat matrix, or a length Npat vector, depending
#   on the type of inputs to the pars argument).
# @param ... Arguments passed to .ll_mvmer. Can include 'reloo_or_kfold' which is
#   a logical specifying whether the function calling ll_jm was reloo or kfold.
# @return Either a matrix, a vector or a scalar, depending on the input types
#   and whether sum is set to TRUE.
.ll_jm <- function(object, data, pars, include_long = TRUE, 
                   include_b = FALSE, sum = FALSE, ...) {
  
  M <- get_M(object)
  
  # Log likelihood for event submodel
  ll_event <- .ll_survival(object, data, pars)
  
  # Log likelihoods for longitudinal submodels
  if (include_long) {
    ll_long <- lapply(1:M, function(m) 
      .ll_long(object, data, pars, m = m, ...))
  }
  
  # Log likelihood for random effects submodel
  # NB this is only used in the Metropolis algorithm in 'posterior_survfit'
  #   when drawing random effects for new individuals. But it is not used
  #   in generating the pointwise log likelihood matrix under log_lik or loo.
  if (include_b) {
    if (length(object$cnms) > 2L)
      stop("Bug found: 'include_b' cannot be TRUE when there is more than ",
           "2 grouping factors.")
    if (length(object$cnms) == 2L && M > 1)
      stop("Bug found: 'include_b' cannot be TRUE when there is more than ",
           "one longitudinal submodel and more than one grouping factor.")
    if ((data$Npat > 1) || (nrow(pars$stanmat) > 1L))
      stop("Bug found: 'include_b' can only be TRUE when 'data' is for one ",
           "individual, and stanmat is for a single draw.")
    id_var <- object$id_var
    cnms   <- object$cnms
    Z_names <- fetch_(data$assoc_parts, "mod_eta", "Z_names")
    b <- do.call("cbind", pars$b)
    b <- as.vector(pp_b_ord(b, Z_names))
    Sigma_id <- VarCorr(object, stanmat = pars$stanmat)[[id_var]]
    if (length(cnms) > 1L) {
      b2_var <- grep(utils::glob2rx(id_var), names(cnms), 
                     value = TRUE, invert = TRUE)
      Sigma_b2 <- VarCorr(object, stanmat = pars$stanmat)[[b2_var]]
      Sigma_list <- rep(list(Sigma_b2), data$Ni)
      which_slot <- which(names(cnms) == b2_var)
      if (which_slot == 1L) {
        Sigma_bind <- c(Sigma_list, list(Sigma_id))
      } else {
        Sigma_bind <- c(list(Sigma_id), Sigma_list)
      }
      Sigma <- as.matrix(Matrix::bdiag(Sigma_bind))
    } else {
      Sigma <- Sigma_id
    }
    ll_b <- -0.5 * (c(determinant(Sigma, logarithm = TRUE)$modulus) + 
      (b %*% chol2inv(chol(Sigma)) %*% b)[1] + length(b) * log(2 * pi))
  } else {
    ll_b <- NULL
  }
  
  # Check the dimensions of the various components
  if (is.matrix(ll_event)) { # S * Npat matrices
    if (include_long) {
      mats <- unique(sapply(c(ll_long, list(ll_event)), is.matrix))
      dims <- unique(lapply(c(ll_long, list(ll_event)), dim)) 
      if ((length(dims) > 1L) || (length(mats) > 1L))
        stop("Bug found: elements of 'll_long' should be same class and ",
             "dimension as 'll_event'.")
    }
    if (include_b && !identical(length(ll_b), ncol(ll_event)))
      stop("Bug found: length of 'll_b' should be equal to the number of ",
           "columns in 'll_event'.")
  } else { # length Npat vectors (ie, log-lik based on a single draw of pars)
    if (include_long) {
      lens <- unique(sapply(c(ll_long, list(ll_event)), length))
      if (length(lens) > 1L)
        stop("Bug found: elements of 'll_long' should be same length as 'll_event'.")
    }
    if (include_b && !identical(length(ll_b), length(ll_event)))
      stop("Bug found: length of 'll_b' should be equal to length of 'll_event'.")
  }  
  
  # Sum the various components (long + event + random effects)
  if (include_long) {
    val <- Reduce('+', c(ll_long, list(ll_event)))
  } else {
    val <- ll_event
  }
  if (include_b && is.matrix(val)) {
    val <- sweep(val, 2L, ll_b, `+`) 
  } else if (include_b && is.vector(val)) {
    val <- val + ll_b
  }
  
  # Return log likelihood for joint model
  if (!sum) {
    return(val)             # S * Npat matrix or length Npat vector
  } else if (is.matrix(val)) {
    return(rowSums(val))    # length S vector
  } else {
    return(sum(val))        # scalar 
  }
}

# Return log-likelihood for longitudinal submodel m
#
# @param object A stanjm object.
# @param data Output from .pp_data_jm.
# @param pars Output from extract_pars.
# @param m Integer specifying the longitudinal submodel.
# @param reloo_or_kfold Logical specifying whether the call came from 
#   reloo or kfold.
# @return An S*Npat matrix.
.ll_long <- function(object, data, pars, m = 1, reloo_or_kfold = FALSE) {
  args <- ll_args.stanjm(object, data, pars, m = m, 
                         reloo_or_kfold = reloo_or_kfold)
  fun  <- ll_fun(object, m = m)
  ll <- lapply(seq_len(args$N), function(j) as.vector(
    fun(data_i = args$data[j, , drop = FALSE], draws = args$draws)))
  ll <- do.call("cbind", ll)
  # return S*Npat matrix by summing log-lik for y within each individual
  res <- apply(ll, 1L, function(row) tapply(row, data$flist[[m]], sum))
  res <- if (is.vector(res) & (args$S > 1L)) cbind(res) else t(res)
  return(res) 
}

# Return survival probability or log-likelihood for event submodel
#
# @param object A stanjm object.
# @param data Output from .pp_data_jm.
# @param pars Output from extract_pars.
# @param one_draw A logical specifying whether the parameters provided in the 
#   pars argument are vectors for a single realisation of the parameter (e.g.
#   a single MCMC draw, or a posterior mean) (TRUE) or a stanmat array (FALSE).
# @param survprob A logical specifying whether to return the survival probability 
#   (TRUE) or the log likelihood for the event submodel (FALSE).
# @param An S by Npat matrix, or a length Npat vector, depending on the inputs
#   (where S is the size of the posterior sample and Npat is the number of 
#   individuals).
.ll_survival <- function(object, data, pars, one_draw = FALSE, survprob = FALSE) {
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
  
  # Linear predictor for the survival submodel
  e_eta <- linear_predictor(pars$ebeta, data$eXq) 
  
  # Add on contribution from assoc structure
  if (length(pars$abeta)) {
    M <- get_M(object)
    # Temporary stop, until make_assoc_terms can handle it
    sel_stop <- grep("^shared", rownames(object$assoc))
    if (any(unlist(object$assoc[sel_stop,])))
      stop("'log_lik' cannot yet be used with shared_b or shared_coef ",
           "association structures.", call. = FALSE)
    pars$b <- lapply(1:M, function(m) {
      b_m <- pars$b[[m]]
      Z_names_m <- data$assoc_parts[[m]][["mod_eta"]][["Z_names"]]
      pp_b_ord(if (is.matrix(b_m)) b_m else t(b_m), Z_names_m)
    })
    if (one_draw) {
      aXq <- make_assoc_terms(parts = data$assoc_parts, assoc = assoc, 
                              family = family, beta = pars$beta, b = pars$b)
      e_eta <- e_eta + linear_predictor.default(pars$abeta, aXq)
    } else {
      aXq <- make_assoc_terms(parts = data$assoc_parts, assoc = assoc, 
                              family = family, beta = pars$beta, b = pars$b)
      for (k in 1:length(aXq)) {
        e_eta <- e_eta + sweep(aXq[[k]], 1L, pars$abeta[,k], `*`)
      }
    }    
  }
  
  # Log baseline hazard at etimes (if not NULL) and qtimes
  log_basehaz <- evaluate_log_basehaz(times = times, 
                                      basehaz = basehaz, 
                                      coefs = pars$bhcoef)
  
  # Log hazard at etimes (if not NULL) and qtimes
  log_haz <- log_basehaz + e_eta  
  
  # Extract log hazard at qtimes only
  if (is.vector(log_haz)) {
    q_log_haz <- tail(log_haz, length(qtimes))
  } else {
    sel_cols <- tail(1:ncol(log_haz), length(qtimes))
    q_log_haz <- log_haz[, sel_cols, drop = FALSE]
  }
  
  # Evaluate log survival
  log_surv <- evaluate_log_survival(log_haz = q_log_haz,
                                    qnodes = qnodes, qwts = qwts)
  
  # Force surv prob to 1 (ie. log surv prob to 0) if evaluating
  # at time t = 0; this avoids possible numerical errors
  log_surv[etimes == 0] <- 0
  
  # Possibly return surv prob at time t (upper limit of integral)
  if (survprob)
    return(exp(log_surv)) 
  
  # Otherwise return log likelihood at time t
  if (is.null(etimes) || is.null(estatus))
    stop("'etimes' and 'estatus' cannot be NULL if 'survprob = FALSE'.")
  times_length <- length(c(etimes, qtimes))
  if (one_draw) { # return vector of length npat
    if (!length(log_haz) == times_length)
      stop2("Bug found: length of log_haz vector is incorrect.")
    e_log_haz <- log_haz[1:length(etimes)]
    return(estatus * e_log_haz + log_surv)
  } else { # return S * npat matrix
    if (!ncol(log_haz) == times_length)
      stop2("Bug found: number of cols in log_haz matrix is incorrect.")
    e_log_haz <- log_haz[, 1:length(etimes), drop = FALSE]
    return(sweep(e_log_haz, 2L, estatus, `*`) + log_surv)
  }
} 

# Evaluate the log baseline hazard at the specified times
# given the vector or matrix of MCMC draws for the baseline
# hazard coeffients / parameters
#
# @param times A vector of times.
# @param basehaz A list with info about the baseline hazard.
# @param coefs A vector or matrix of parameter estimates (MCMC draws).
# @return A vector or matrix, depending on the input type of coefs.
evaluate_log_basehaz <- function(times, basehaz, coefs) {
  type <- basehaz$type_name
  if (type == "weibull") { 
    X  <- log(times) # log times
    B1 <- log(coefs) # log shape
    B2 <- coefs - 1  # shape - 1
    log_basehaz <- as.vector(B1) + linear_predictor(B2,X)
  } else if (type == "bs") { 
    X <- predict(basehaz$bs_basis, times) # b-spline basis
    B <- coefs                            # b-spline coefs
    log_basehaz <- linear_predictor(B,X)
  } else {
    stop2("Not yet implemented for basehaz = ", type)
  }
  log_basehaz
}

# Evaluate the log baseline hazard at the specified times
# given the vector or matrix of MCMC draws for the baseline
# hazard coeffients / parameters
#
# @param log_haz A vector containing the log hazard for each
#   individual, evaluated at each of the quadrature points. The
#   vector should be ordered such that the first N elements contain
#   the log_haz evaluated for each individual at quadrature point 1,
#   then the next N elements are the log_haz evaluated for each 
#   individual at quadrature point 2, and so on.
# @param qnodes Integer specifying the number of quadrature nodes
#   at which the log hazard was evaluated for each individual.
# @param qwts A vector of unstandardised GK quadrature weights.
# @return A vector or matrix of log survival probabilities.
evaluate_log_survival <- function(log_haz, qnodes, qwts) {
  UseMethod("evaluate_log_survival")
}

evaluate_log_survival.default <- function(log_haz, qnodes, qwts) {
  # convert log hazard to hazard
  haz <- exp(log_haz)
  # apply GK quadrature weights
  weighted_haz <- qwts * haz 
  # sum quadrature points for each individual to get cum_haz
  splitting_vec <- rep(1:qnodes, each = length(haz) / qnodes)
  cum_haz <- Reduce('+', split(weighted_haz, splitting_vec))
  # return: -cum_haz == log survival probability
  -cum_haz
}

evaluate_log_survival.matrix <- function(log_haz, qnodes, qwts) {
  # convert log hazard to hazard
  haz <- exp(log_haz)
  # apply GK quadrature weights
  weighted_haz <- sweep(haz, 2L, qwts, `*`)
  # sum quadrature points for each individual to get cum_haz
  cum_haz <- Reduce('+', array2list(weighted_haz, nsplits = qnodes))
  # return: -cum_haz == log survival probability
  -cum_haz
} 
