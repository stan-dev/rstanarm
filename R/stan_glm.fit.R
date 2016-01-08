# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2013, 2014, 2015 Trustees of Columbia University
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

#' @rdname stan_glm
#' @export
#' @param group A list, possibly of length zero (the default), but otherwise
#'   having the structure of that produced by \code{\link[lme4]{mkReTrms}} to
#'   indicate the group-specific part of the model. In addition, this list must
#'   have elements for the \code{regularization}, \code{concentration} 
#'   \code{shape}, and \code{scale} components of a \code{\link{decov}}
#'   prior for the covariance matrices among the group-specific coefficients.
#'   
stan_glm.fit <- function(x, y, weights = rep(1, NROW(x)), 
                         offset = rep(0, NROW(x)), family = gaussian(),
                         ...,
                         prior = normal(),
                         prior_intercept = normal(),
                         prior_ops = prior_options(),
                         group = list(),
                         prior_PD = FALSE, 
                         algorithm = c("sampling", "optimizing", "meanfield", "fullrank"), 
                         adapt_delta = NULL, QR = FALSE) {
  
  family <- validate_family(family)
  supported_families <- c("binomial", "gaussian", "Gamma", "inverse.gaussian",
                          "poisson", "neg_binomial_2")
  fam <- which(pmatch(supported_families, family$family, nomatch = 0L) == 1L)
  if (!length(fam)) 
    stop("'family' must be one of ", paste(supported_families, collapse = ", "))
  
  # these are also from help(family)
  supported_links <- switch(supported_families[fam],
                            binomial = c("logit", "probit", "cauchit", "log", "cloglog"),
                            gaussian = c("identity", "log", "inverse"),
                            Gamma = c("identity", "log", "inverse"),
                            inverse.gaussian = c("identity", "log", "inverse", "1/mu^2"),
                            "neg_binomial_2" = , # intentional
                            poisson = c("log", "identity", "sqrt"),
                            stop("unsupported family"))
  link <- which(supported_links == family$link)
  if (!length(link)) 
    stop("'link' must be one of ", paste(supported_links, collapse = ", "))
  
  if (is.binomial(family$family)) {
    if (NCOL(y) != 1L) {
      stopifnot(NCOL(y) == 2L)
      trials <- as.integer(y[, 1L] + y[, 2L])
      y <- as.integer(y[, 1L])
    } else if (all(weights == 1)) {
      # convert factors to 0/1 using R's convention that first factor level is
      # treated as failure
      if (is.factor(y)) 
        y <- y != levels(y)[1L]
      y <- as.integer(y)
      if (!all(y %in% c(0L, 1L))) 
        stop("y values must be 0 or 1 for bernoulli regression.")
    }
    else {
      if (!all(y >= 0 & y <= 1))
        stop("If weights are provided, then y values must be proportions ", 
             "between 0 and 1.")
      trials <- weights
    }
  }
  
  x <- as.matrix(x)
  has_intercept <- grepl("(Intercept", colnames(x)[1L], fixed = TRUE)
  xtemp <- if (has_intercept) x[, -1L, drop=FALSE] else x
  xbar <- colMeans(xtemp)
  xtemp <- sweep(xtemp, 2, xbar, FUN = "-")
  
  # drop any column of x with < 2 unique values (empty interaction levels)
  sel <- (2 > apply(xtemp, 2L, function(x) length(unique(x))))
  if (any(sel)) {
    warning("Dropped empty interaction levels: ",
            paste(colnames(xtemp)[sel], collapse = ", "))
    xtemp <- xtemp[, !sel, drop = FALSE]
    xbar <- xbar[!sel]
  }
  nvars <- ncol(xtemp)
  
  scaled <- prior_ops$scaled
  min_prior_scale <- prior_ops$min_prior_scale
  prior_scale_for_dispersion <- prior_ops$prior_scale_for_dispersion
  ok_dists <- nlist("normal", student_t = "t", "cauchy", "hs", "hs_plus")
  ok_intercept_dists <- ok_dists[1:3]
  
  # prior distributions
  if (!is.null(prior)) {
    if (!is.list(prior)) stop("'prior' should be a named list.")
    prior_dist <- prior$dist
    prior_scale <- prior$scale
    prior_mean <- prior$location
    prior_df <- prior$df
    prior_df[is.na(prior_df)] <- 1
    if (!prior_dist %in% unlist(ok_dists)) {
      stop("The prior distribution for the coefficients should be one of ",
           paste(names(ok_dists), collapse = ", "))
    }
    if (prior_dist == "normal" || prior_dist == "t") {
      prior_dist <- ifelse(prior_dist == "normal", 1L, 2L)
      prior_scale <- set_prior_scale(prior_scale, default = 2.5, 
                                     link = family$link)
    }
    else prior_dist <- ifelse(prior_dist == "hs", 3L, 4L)
    prior_df <- maybe_broadcast(prior_df, nvars)
    prior_df <- as.array(pmin(.Machine$double.xmax, prior_df))
    prior_mean <- maybe_broadcast(prior_mean, nvars)
    prior_mean <- as.array(prior_mean)
    prior_scale <- maybe_broadcast(prior_scale, nvars)
  }
  else {
    prior_dist <- 0L
    prior_mean <- as.array(rep(0, nvars))
    prior_scale <- prior_df <- as.array(rep(1, nvars))
  }
  if (!is.null(prior_intercept)) {
    if (!is.list(prior_intercept)) stop("'prior' should be a named list.")
    prior_dist_for_intercept <- prior_intercept$dist
    prior_scale_for_intercept <- prior_intercept$scale
    prior_mean_for_intercept <- prior_intercept$location
    prior_df_for_intercept <- prior_intercept$df 
    prior_df_for_intercept[is.na(prior_df_for_intercept)] <- 1
    
    if (!prior_dist_for_intercept %in% unlist(ok_intercept_dists)) {
      stop("The prior distribution for the intercept should be one of ",
           paste(names(ok_intercept_dists), collapse = ", "))
    }
    prior_dist_for_intercept <- 
      ifelse(prior_dist_for_intercept == "normal", 1L, 2L)
    prior_scale_for_intercept <- 
      set_prior_scale(prior_scale_for_intercept, default = 10, 
                      link = family$link)
    prior_df_for_intercept <- min(.Machine$double.xmax, prior_df_for_intercept)
  }
  else {
    prior_dist_for_intercept <- 0L
    prior_mean_for_intercept <- 0 
    prior_scale_for_intercept <- prior_df_for_intercept <- 1
  }
  
  is_gaussian <- is.gaussian(family$family)
  is_gamma <- is.gamma(family$family)
  is_ig <- is.ig(family$family)
  is_continuous <- is_gaussian || is_gamma || is_ig
  if (scaled && prior_dist > 0L) {
    if (is_gaussian) {
      ss <- 2 * sd(y)
      prior_scale <- ss * prior_scale
      prior_scale_for_intercept <-  ss * prior_scale_for_intercept
    }
    if (!QR) 
      prior_scale <- pmax(min_prior_scale, prior_scale / 
             apply(xtemp, 2L, FUN = function(x) {
               num.categories <- length(unique(x))
               x.scale <- 1
               if (num.categories == 2) x.scale <- diff(range(x))
               else if (num.categories > 2) x.scale <- 2 * sd(x)
               return(x.scale)
             }))
  }
  prior_scale <- as.array(pmin(.Machine$double.xmax, prior_scale))
  prior_scale_for_intercept <- min(.Machine$double.xmax, prior_scale_for_intercept)
  
  is_bernoulli <- is.binomial(supported_families[fam]) && all(y %in% 0:1)
  is_nb <- is.nb(supported_families[fam])

  if (QR) {
    if (ncol(xtemp) <= 1)
      stop("'QR' can only be specified when there are multiple predictors.")
    cn <- colnames(xtemp)
    decomposition <- qr(xtemp)
    sqrt_nm1 <- sqrt(nrow(xtemp) - 1L)
    Q <- qr.Q(decomposition)
    R_inv <- qr.solve(decomposition, Q) * sqrt_nm1
    xtemp <- Q * sqrt_nm1
    colnames(xtemp) <- cn
    xbar <- c(xbar %*% R_inv)
  }
  # create entries in the data block of the .stan file
  standata <- list(
    N = nrow(xtemp), K = ncol(xtemp), xbar = as.array(xbar), link = link,
    has_weights = as.integer(length(weights) > 0),
    has_offset = as.integer(length(offset) > 0),
    prior_dist = prior_dist, prior_mean = prior_mean, prior_scale = prior_scale, 
    prior_df = prior_df, prior_dist_for_intercept = prior_dist_for_intercept,
    prior_scale_for_intercept = prior_scale_for_intercept, 
    prior_mean_for_intercept = prior_mean_for_intercept,
    prior_df_for_intercept = prior_df_for_intercept,
    has_intercept = as.integer(has_intercept), prior_PD = as.integer(prior_PD))
  
  if (length(group)) {
    decov <- group$decov
    Z <- t(as.matrix(group$Zt))
    group <- pad_reTrms(Z = Z, cnms = group$cnms, flist = group$flist)
    Z <- group$Z
    p <- sapply(group$cnms, FUN = length)
    l <- sapply(attr(group$flist, "assign"), function(i) 
      nlevels(group$flist[[i]]))
    t <- length(p)
    group_nms <- names(group$cnms)
    b_names <- unlist(lapply(1:t, FUN = function(i) {
      nms_i <- paste(group$cnms[[i]], group_nms[i])
      if (length(nms_i) == 1) paste0(nms_i, ":", levels(group$flist[[i]]))
      else sapply(nms_i, paste0, ":", levels(group$flist[[i]]))
    }))
    g_names <- unlist(lapply(1:t, FUN = function(i) {
      paste(group$cnms[[i]], names(group$cnms)[i], sep = "|")
    }))
    standata$t <- t
    standata$p <- as.array(p)
    standata$l <- as.array(l)
    standata$q <- ncol(Z)
    standata$len_theta_L <- sum(choose(p,2), p)
    if (is_bernoulli) {
      parts0 <- extract_sparse_parts(Z[y == 0,, drop = FALSE])
      parts1 <- extract_sparse_parts(Z[y == 1,, drop = FALSE])
      standata$num_non_zero <- c(length(parts0$w), length(parts1$w))
      standata$w0 <- parts0$w
      standata$w1 <- parts1$w
      standata$v0 <- parts0$v
      standata$v1 <- parts1$v
      standata$u0 <- parts0$u
      standata$u1 <- parts1$u
    }
    else {
      parts <- extract_sparse_parts(Z)
      standata$num_non_zero <- length(parts$w)
      standata$w <- parts$w
      standata$v <- parts$v
      standata$u <- parts$u
    }
    standata$shape <- as.array(maybe_broadcast(decov$shape, t))
    standata$scale <- as.array(maybe_broadcast(decov$scale, t))
    standata$len_concentration <- sum(p[p > 1])
    standata$concentration <- as.array(maybe_broadcast(decov$concentration, 
                                                       sum(p[p > 1])))
    standata$len_regularization <- sum(p > 1)
    standata$regularization <- as.array(maybe_broadcast(
                                        decov$regularization, sum(p > 1)))
  }
  else {
    standata$t <- 0L
    standata$p <- integer(0)
    standata$l <- integer(0)
    standata$q <- 0L
    standata$len_theta_L <- 0L
    if (is_bernoulli) {
      standata$num_non_zero <- rep(0L, 2)
      standata$w0 <- standata$w1 <- double(0)
      standata$v0 <- standata$v1 <- integer(0)
      standata$u0 <- standata$u1 <- integer(0)
    }
    else {
      standata$num_non_zero <- 0L
      standata$w <- double(0)
      standata$v <- integer(0)
      standata$u <- integer(0)
    }
    standata$shape <- standata$scale <- standata$concentration <-
      standata$regularization <- rep(0, 0)
    standata$len_concentration <- 0L
    standata$len_regularization <- 0L
  }
  
  if (!is_bernoulli) {
    standata$X <- xtemp
    standata$y <- y
    standata$weights <- weights
    standata$offset <- offset
  }

  # call stan() to draw from posterior distribution
  if (is_continuous) {
    standata$prior_scale_for_dispersion <- prior_scale_for_dispersion %ORifINF% 0
    standata$family <- switch(family$family, 
                              gaussian = 1L, 
                              Gamma = 2L,
                              3L)
    stanfit <- stanmodels$continuous
  }
  else if (is.binomial(supported_families[fam])) {
    standata$prior_scale_for_dispersion <- 
      if (!length(group) || prior_scale_for_dispersion == Inf) 
        0 else prior_scale_for_dispersion
    standata$family <- 1L # not actually used
    if (is_bernoulli) {
      y0 <- y == 0
      y1 <- y == 1
      standata$N <- c(sum(y0), sum(y1))
      standata$X0 <- xtemp[y0,, drop = FALSE]
      standata$X1 <- xtemp[y1,, drop = FALSE]
      standata$Z0 <- standata$Z[y0,, drop = FALSE]
      standata$Z1 <- standata$Z[y1,, drop = FALSE]
      standata$Z <- NULL 
      if (length(weights)) {
        standata$weights0 <- weights[y0]
        standata$weights1 <- weights[y1]
      }
      else {
        standata$weights0 <- double(0)
        standata$weights1 <- double(0)
      }
      if (length(offset)) {
        standata$offset0 <- offset[y0]
        standata$offset1 <- offset[y1]
      }
      else {
        standata$offset0 <- double(0)
        standata$offset1 <- double(0)
      }
      stanfit <- stanmodels$bernoulli
    }
    else {
      standata$trials <- trials
      if (length(weights) & !all(weights == 1)) {
        standata$y <- round(y * trials)
        standata$weights <- double(0)
        standata$has_weights <- 0L
      }
      stanfit <- stanmodels$binomial
    }
  }   
  else if (is.poisson(supported_families[fam])) {
    standata$family <- 1L
    standata$prior_scale_for_dispersion <- prior_scale_for_dispersion %ORifINF% 0
    stanfit <- stanmodels$count 
  }
  else if (is_nb) {
    standata$family <- 2L
    standata$prior_scale_for_dispersion <- prior_scale_for_dispersion %ORifINF% 0
    stanfit <- stanmodels$count
  }
  else if (is_gamma) {
    # nothing
  }
  else stop(paste(family$family, "is not supported"))
  
  pars <- c(if (has_intercept) "alpha", "beta", 
            if (length(group)) "b",
            if (is_continuous | is_nb) "dispersion", "mean_PPD")
  algorithm <- match.arg(algorithm)
  if (algorithm == "optimizing") {
    out <- optimizing(stanfit, data = standata, 
                      draws = 1000, constrained = TRUE, ...)
    new_names <- names(out$par)
    mark <- grepl("^beta\\[[[:digit:]]+\\]$", new_names)
    if (QR) {
      out$par[mark] <- R_inv %*% out$par[mark]
      out$theta_tilde[,mark] <- out$theta_tilde[,mark] %*% t(R_inv)
    }
    new_names[mark] <- colnames(xtemp)
    new_names[new_names == "alpha[1]"] <- "(Intercept)"
    new_names[grepl("dispersion(\\[1\\])?$", new_names)] <- 
      if (is_gaussian) "sigma" else
        if (is_gamma) "shape" else
          if (is_ig) "lambda" else 
            if (is_nb) "overdispersion" else NA
    names(out$par) <- new_names
    colnames(out$theta_tilde) <- new_names
    out$stanfit <- suppressMessages(sampling(stanfit, data = standata, chains = 0))
    return(out)
  }
  else {
    if (algorithm == "sampling") {
      sampling_args <- set_sampling_args(
        object = stanfit, 
        prior = prior, 
        user_dots = list(...), 
        user_adapt_delta = adapt_delta, 
        data = standata, pars = pars, show_messages = FALSE)
      stanfit <- do.call(sampling, sampling_args)
    }
    else
      stanfit <- rstan::vb(stanfit, pars = pars, data = standata,
                           algorithm = algorithm, ...)
    if (QR) {
      thetas <- extract(stanfit, pars = "beta", inc_warmup = TRUE, permuted = FALSE)
      betas <- apply(thetas, 1:2, FUN = function(theta) R_inv %*% theta)
      end <- tail(dim(betas), 1)
      for (chain in 1:end) for (param in 1:nrow(betas)) {
        stanfit@sim$samples[[chain]][[has_intercept + param]] <-
          if (ncol(xtemp) > 1) betas[param,,chain] else betas[param,chain]
      }
    }
    new_names <- c(if (has_intercept) "(Intercept)", 
                   colnames(xtemp), 
                   if (length(group)) c(paste0("b[", b_names, "]")),
                   if (is_gaussian) "sigma", 
                   if (is_gamma) "shape", 
                   if (is_ig) "lambda",
                   if (is_nb) "overdispersion", 
                   "mean_PPD", 
                   "log-posterior")
    stanfit@sim$fnames_oi <- new_names
    return(stanfit)
  }
}


# Add extra level _NEW_ to each group
# 
# @param Z ranef indicator matrix
# @param cnms group$cnms
# @param flist group$flist
pad_reTrms <- function(Z, cnms, flist) {
  l <- sapply(attr(flist, "assign"), function(i) nlevels(flist[[i]]))
  p <- sapply(cnms, FUN = length)
  last <- cumsum(l * p)
  for (i in attr(flist, "assign")) {
    levels(flist[[i]]) <- c(levels(flist[[i]]), paste0("_NEW_", names(flist)[i]))
  }
  n <- nrow(Z)
  Z <- cbind(Z, matrix(0, nrow = n, ncol = p[length(p)], 
             dimnames = list(NULL, rep("_NEW_", p[length(p)]))))
  mark <- length(p) - 1L
  for (i in rev(head(last, -1))) {
    Z <- cbind(Z[,1:i, drop = FALSE],
               matrix(0, n, p[mark], dimnames = list(NULL, rep("_NEW_", p[mark]))),
               Z[,(i+1):ncol(Z), drop = FALSE])
    mark <- mark - 1L
  }
  return(nlist(Z, cnms, flist))
}

# Drop the extra reTrms from a matrix x
#
# @param x A matrix (e.g. the posterior sample or matrix of summary stats)
# @param columns Do the columns (TRUE) or rows (FALSE) correspond to the
#   variables?
unpad_reTrms <- function(x, columns = TRUE) {
  stopifnot(is.matrix(x))
  nms <- if (columns) colnames(x) else rownames(x)
  keep <- !grepl("_NEW_", nms, fixed = TRUE)
  if (columns) x[, keep, drop = FALSE] 
  else x[keep,, drop = FALSE]
}
