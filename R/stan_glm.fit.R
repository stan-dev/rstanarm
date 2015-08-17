# This file is part of rstanarm.
# Copyright 2013 Stan Development Team
# rstanarm is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# rstanarm is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with rstanarm.  If not, see <http://www.gnu.org/licenses/>.

#' @rdname stan_glm
#' @param prior.dist,prior.dist.for.intercept Character string indicating the 
#'   family of the prior distribution for the coefficients (see
#'   \code{\link{priors}}) or \code{NULL} to omit this prior
#' @param scaled Logical scalar indicating whether to rescale the predictors
#' @param prior.mean,prior.mean.for.intercept Numeric vector (possibly of
#'   length one) indicating the locations of the \code{prior.dist} and the
#'   \code{prior.dist.for.intercept} respectively
#' @param prior.scale,prior.scale.for.intercept Numeric vector (possibly of
#'   length one) indicating the scale of the \code{prior.dist} and the
#'   \code{prior.dist.for.intercept} respectively or \code{NULL}
#' @param prior.df,prior.df.for.intercept Numeric vector (possibly of length
#'   one) indicating the degrees of freedom of the \code{prior.dist} and
#'   the \code{prior.dist.for.intercept} in the Student t case
#' @param min.prior.scale Positive scalar indicating the smallest possible
#'   value to use for rescaling the predictors
#' @param prior.scale.for.dispersion Positive scalar indicating the scale
#'   parameter for the Cauchy prior on the dispersion parameter (if the
#'   model has one) or \code{NULL} to omit this prior
#' @param group A list, possibly of length zero, but otherwise having the
#'   structure of that produced by \code{\link[lme4]{mkReTrms}} to indicate
#'   the group-specific part of the model. In addition, this list must have
#'   elements for the \code{gamma_shape}, \code{scale}, \code{concentration}
#'   and \code{shape} components of a \code{\link{decov}} prior for the
#'   covariance matrices among the group-specific coefficients
#' @param prior.dist,prior.dist.for.intercept A character string, either 
#'   \code{"normal"} (the default), or \code{"t"} indicating the family of the 
#'   prior distribution for the coefficients, or \code{NULL} to omit this prior.
#' @param scaled A logical scalar indicating whether to rescale the predictors.
#' @param prior.mean,prior.mean.for.intercept A numeric vector (possibly of 
#'   length one) indicating the locations of the \code{prior.dist} and the 
#'   \code{prior.dist.for.intercept} respectively.
#' @param prior.scale,prior.scale.for.intercept A numeric vector (possibly of 
#'   length one) indicating the scale of the \code{prior.dist} and the 
#'   \code{prior.dist.for.intercept} respectively or \code{NULL}.
#' @param prior.df,prior.df.for.intercept A numeric vector (possibly of length 
#'   one) indicating the degrees of freedom of the \code{prior.dist} and the
#'   \code{prior.dist.for.intercept} in the Student t case.
#' @param min.prior.scale A positive scalar indicating the smallest possible 
#'   value to use for rescaling the predictors.
#' @param prior.scale.for.dispersion A positive scalar indicating the scale 
#'   parameter for the Cauchy prior on the dispersion parameter (if the model
#'   has one) or \code{NULL} to omit this prior.
#' @param group A list, possibly of length zero, but otherwise having the 
#'   structure of that produced by \code{\link[lme4]{mkReTrms}} to indicate the
#'   group-specific part of the model. In addition, this list must have elements
#'   for the \code{gamma_shape}, \code{scale}, \code{concentration} and
#'   \code{shape} components of a \code{\link{decov}} prior for the covariance
#'   matrices among the group-specific coefficients.
#' @export
#' 
stan_glm.fit <- function(x, y, weights = rep(1, NROW(x)), start = NULL, 
                         offset = rep(0, NROW(x)), family = gaussian(),
                         prior.dist = c("normal", "t", "horseshoe", 
                                        "horseshoe_plus"),
                         prior.dist.for.intercept = c("normal", "t"), 
                         scaled = TRUE, 
                         prior.mean = 0, prior.scale = NULL, prior.df = 1,
                         prior.mean.for.intercept = 0, 
                         prior.scale.for.intercept = NULL, 
                         prior.df.for.intercept = 1,
                         min.prior.scale = 1e-12, 
                         prior.scale.for.dispersion = 5, group = list(),
                         prior_PD = FALSE, algorithm = c("sampling", "optimizing"), 
                         ...) { # further arguments to sampling() or optimizing()
  
  if (is.character(family)) 
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) 
    family <- family()
  if(!is(family, "family"))
    stop("'family' must be a family")
  
  # these are from help(family)
  supported_families <- c("binomial", "gaussian", "Gamma",
                          "poisson", "Negative Binomial") # FIXME: add "inverse.gaussian"
  fam <- which(pmatch(supported_families, family$family, nomatch = 0L) == 1L)
  if(length(fam) == 0) stop(paste("'family' must be one of ", supported_families))
  
  # these are also from help(family)
  supported_links <- switch(supported_families[fam],
                            binomial = c("logit", "probit", "cauchit", "log", "cloglog"),
                            gaussian = c("identity", "log", "inverse"),
                            Gamma = c("inverse", "identity", "log"),
                            # inverse.gaussian = c("1/mu^2", "inverse", "identity", "log"),
                            "Negative Binomial" = , # intentional
                            poisson = c("log", "identity", "sqrt"),
                            stop("unsupported family"))
  link <- which(supported_links == family$link)
  if (length(link) == 0) stop(paste("'link' must be one of ", supported_links))
  
  if (family$family == "binomial") {
    if (NCOL(y) != 1L) {
      stopifnot(NCOL(y) == 2L)
      trials <- as.integer(y[, 1L] + y[, 2L])
      y <- as.integer(y[, 1L])
    } else {
      # convert factors to 0/1 using R's convention that first factor level is
      # treated as failure
      if (is.factor(y)) 
        y <- y != levels(y)[1L]
      y <- as.integer(y)
      if (!all(y %in% c(0L, 1L))) 
        stop("y values must be 0 or 1 for logistic regression")
    } 
  }
  
  x <- as.matrix(x)
  has_intercept <- colnames(x)[1L] == "(Intercept)"
  xtemp <- if (has_intercept) x[, -1L, drop=FALSE] else x
  xbar <- colMeans(xtemp)
  xtemp <- sweep(xtemp, 2, xbar, FUN = "-")
  nvars <- ncol(xtemp)
  
  # prior distributions
  if (!is.null(prior.dist)) {
    prior.dist <- match.arg(prior.dist)
    if (prior.dist == "normal" || prior.dist == "t") {
      prior.dist <- ifelse(prior.dist == "normal", 1L, 2L)
      prior.scale <- set_prior_scale(prior.scale, default = 2.5, link = family$link)
    }
    else prior.dist <- ifelse(prior.dist == "horseshoe", 3L, 4L)
    prior.df <- maybe_broadcast(prior.df, nvars)
    prior.df <- as.array(pmin(.Machine$double.xmax, prior.df))
    prior.mean <- maybe_broadcast(prior.mean, nvars)
    prior.mean <- as.array(prior.mean)
    prior.scale <- maybe_broadcast(prior.scale, nvars)
  }
  else {
    prior.dist <- 0L
    prior.mean <- as.array(rep(0, nvars))
    prior.scale <- prior.df <- as.array(rep(1, nvars))
  }
  if (!is.null(prior.dist.for.intercept)) {
    prior.dist.for.intercept <- match.arg(prior.dist.for.intercept)
    prior.dist.for.intercept <- ifelse(prior.dist.for.intercept == "normal", 1L, 2L)
    prior.scale.for.intercept <- set_prior_scale(prior.scale.for.intercept, 
                                                 default = 10, link = family$link)
    prior.df.for.intercept <- min(.Machine$double.xmax, prior.df.for.intercept)
  }
  else {
    prior.dist.for.intercept <- 0L
    prior.mean.for.intercept <- 0 
    prior.scale.for.intercept <- prior.df.for.intercept <- 1
  }
  
  is_gaussian <- family$family == "gaussian"
  if (scaled && prior.dist > 0L) {
    if (is_gaussian) {
      ss <- 2 * sd(y)
      prior.scale <- ss * prior.scale
      prior.scale.for.intercept <-  ss * prior.scale.for.intercept
    }
    prior.scale <- pmax(min.prior.scale, prior.scale /
                          apply(xtemp, 2, FUN = function(x) {
                            num.categories <- length(unique(x))
                            x.scale <- 1
                            if (num.categories == 2) x.scale <- diff(range(x))
                            else if (num.categories > 2) x.scale <- 2 * sd(x)
                          }))
  }
  prior.scale <- as.array(pmin(.Machine$double.xmax, prior.scale))
  priors.scale.for.intercept <- min(.Machine$double.xmax, prior.scale.for.intercept)
  
  is_bernoulli <- supported_families[fam] == "binomial" && all(y %in% 0:1)
  is_nb <- supported_families[fam] == "Negative Binomial"
  
  # create entries in the data block of the .stan file
  standata <- list(
    N = nrow(xtemp), K = ncol(xtemp), xbar = as.array(xbar), link = link,
    has_weights = as.integer(length(weights) > 0),
    has_offset = as.integer(length(offset) > 0),
    prior_dist = prior.dist, prior_mean = prior.mean, prior_scale = prior.scale, 
    prior_df = prior.df, prior_dist_for_intercept = prior.dist.for.intercept,
    prior_scale_for_intercept = prior.scale.for.intercept, 
    prior_mean_for_intercept = prior.mean.for.intercept,
    prior_df_for_intercept = prior.df.for.intercept,
    has_intercept = as.integer(has_intercept), prior_PD = as.integer(prior_PD))
  
  if (length(group) > 0) {
    if (!is_gaussian) stop("only Gaussian models are currently supported with group-specific terms")
    
#     Zind <- as.data.frame(summary(t(group$Zt)))
#     Zind$x <- as.integer(Zind$x)
#     Zind <- as.matrix(Zind)
    Z <- t(as.matrix(group$Zt))
    
    p <- sapply(group$cnms, FUN = length)
    l <- sapply(attributes(group$flist)$assign, function(i) nlevels(group$flist[,i]))
    t <- length(p)
    b_names <- unlist(lapply(1:t, FUN = function(i) rep(group$cnms[[i]], times = l[i])))
    b_names <- paste(b_names, colnames(Z))
    g_names <- unlist(lapply(1:t, FUN = function(i) {
      paste(group$cnms[[i]], names(group$cnms)[i], sep = "|")
    }))
    standata$t <- t
    standata$p <- as.array(p)
    standata$l <- as.array(l)
    standata$q <- ncol(Z)
    standata$Z <- Z
    standata$gamma_shape <- as.array(maybe_broadcast(group$decov$gamma_shape, t))
    standata$scale <- as.array(maybe_broadcast(group$decov$scale, t))
    standata$len_concentration <- sum(p[p > 1])
    standata$concentration <- as.array(maybe_broadcast(group$decov$concentration, 
                                                       sum(p[p > 1])))
    standata$len_shape <- sum(p > 1)
    standata$shape <- as.array(maybe_broadcast(group$decov$shape, sum(p > 1)))
  }
  else {
    standata$t <- 0L
    standata$p <- integer(0)
    standata$l <- integer(0)
    standata$q <- 0L
    standata$Z <- matrix(0, nrow = nrow(xtemp), ncol = 0L)
    standata$gamma_shape <- standata$scale <- standata$concentration <-
      standata$shape <- rep(0, 0)
    standata$len_concentration <- 0L
    standata$len_shape <- 0L
  }
  
  if (!is_bernoulli) {
    standata$X <- xtemp
    standata$y <- y
    standata$weights <- weights
    standata$offset <- offset
  }
  # call stan() to draw from posterior distribution
  if (supported_families[fam] == "gaussian") {
    standata$prior_scale_for_dispersion <- if (prior.scale.for.dispersion == Inf) 
      0 else prior.scale.for.dispersion
    stanfit <- get("stanfit_gaussian")
  }
  else if (supported_families[fam] == "binomial") {
    if (is_bernoulli) {
      y0 <- y == 0
      y1 <- y == 1
      standata$N <- c(sum(y0), sum(y1))
      standata$X0 <- xtemp[y0,]
      standata$X1 <- xtemp[y1,]
      if (length(weights) > 0) {
        standata$weights0 <- weights[y0]
        standata$weights1 <- weights[y1]
      }
      else {
        standata$weights0 <- double(0)
        standata$weights1 <- double(0)
      }
      if (length(offset) > 0) {
        standata$offset0 <- offset[y0]
        standata$offset1 <- offset[y1]
      }
      else {
        standata$offset0 <- double(0)
        standata$offset1 <- double(0)
      }
      stanfit <- get("stanfit_bernoulli")  
    }
    else {
      standata$trials <- trials
      stanfit <- get("stanfit_binomial")
    }
  }   
  else if (supported_families[fam] == "poisson") {
    standata$family <- 1L
    standata$prior_scale_for_dispersion <- if (prior.scale.for.dispersion == Inf) 
      0 else prior.scale.for.dispersion
    stanfit <- get("stanfit_count") 
  }
  else if (is_nb) {
    standata$family <- 2L
    standata$prior_scale_for_dispersion <- if (prior.scale.for.dispersion == Inf) 
      0 else prior.scale.for.dispersion
    stanfit <- get("stanfit_count") 
  }
  else stop(paste(family$family, "is not supported"))
  
  if (is.null(start)) start <- "random"
  else start <- as.list(start)
  
  pars <- c(if (has_intercept) "alpha", "beta", 
            if (length(group) > 0) c("b", "var_group"),
            if (is_gaussian) "sigma", if (is_nb) "theta",  "mean_PPD")
  algorithm <- match.arg(algorithm)
  if (algorithm == "sampling") {
    stanfit <- rstan::sampling(stanfit, pars = pars, data = standata, 
                               init = start, ...)
    new_names <- c(if (has_intercept) "(Intercept)", colnames(xtemp), 
                   if (length(group) > 0) c(paste0("b[", b_names, "]"),
                                            paste0("var[", g_names, "]")),
                   if (is_gaussian) "sigma", if (is_nb) "overdispersion", 
                   "mean_PPD", "log-posterior")
    stanfit@sim$fnames_oi <- new_names
    return(stanfit)
  }
  else if (algorithm == "optimizing") {
    out <- rstan::optimizing(stanfit, data = standata, init = start, hessian = TRUE)
    new_names <- c(if (has_intercept) "gamma", colnames(xtemp), 
                   if (is_gaussian) "sigma_unsaled", 
                   if (is_nb) "overdispersion",
                   if (length(group) > 0) c(paste0("u[", b_names, "]"),
                                            paste0("z_T[", 1:(sum(p[p > 2] - 1)), "]"),
                                            paste0("rho[", 1:standata$len_rho, "]"),
                                            paste0("zeta[", 1:standata$len_zeta, "]"),
                                            paste0("tau[", 1:length(group), "]")),
                   if (length(group) > 0) c(paste0("b[", b_names, "]"),
                                            paste0("var[", g_names, "]")),
                   if (is_gaussian) "sigma",
                   if (has_intercept) "(Intercept)", 
                   "mean_PPD", "log-likelihood")
    k <- ncol(out$hessian)
    names(out$par) <- new_names
    colnames(out$hessian) <- rownames(out$hessian) <- new_names[1:k]
    out$cov.scaled <- qr.solve(-out$hessian, diag(1, k , k))
    colnames(out$cov.scaled) <- rownames(out$cov.scaled) <- colnames(out$hessian)
    out$stanfit <- suppressMessages(rstan::sampling(stanfit, data = standata, chains = 0))
    return(out)
  }
}
