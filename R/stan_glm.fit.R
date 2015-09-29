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
#' @param group A list, possibly of length zero (the default), but otherwise
#'   having the structure of that produced by \code{\link[lme4]{mkReTrms}} to
#'   indicate the group-specific part of the model. In addition, this list must
#'   have elements for the \code{gamma_shape}, \code{scale},
#'   \code{concentration} and \code{shape} components of a \code{\link{decov}}
#'   prior for the covariance matrices among the group-specific coefficients.
#' @importFrom methods is   
#' @export
#' 
stan_glm.fit <- function(x, y, weights = rep(1, NROW(x)), 
                         offset = rep(0, NROW(x)), family = gaussian(),
                         ...,
                         prior = normal(),
                         prior_intercept = normal(),
                         prior_ops = prior_options(),
                         group = list(),
                         prior_PD = FALSE, 
                         algorithm = c("sampling", "optimizing")) {
  
  if (is.character(family)) 
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) 
    family <- family()
  if(!is(family, "family"))
    stop("'family' must be a family")
  
  # these are from help(family)
  supported_families <- c("binomial", "gaussian", "Gamma", "inverse.gaussian",
                          "poisson", "neg_binomial_2")
  fam <- which(pmatch(supported_families, family$family, nomatch = 0L) == 1L)
  if(length(fam) == 0) 
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
  if (length(link) == 0) 
    stop("'link' must be one of ", paste(supported_links, collapse = ", "))
  
  if (family$family == "binomial") {
    if (NCOL(y) != 1L) {
      stopifnot(NCOL(y) == 2L)
      trials <- as.integer(y[, 1L] + y[, 2L])
      y <- as.integer(y[, 1L])
    } else if (all(weights == 1)){
      # convert factors to 0/1 using R's convention that first factor level is
      # treated as failure
      if (is.factor(y)) 
        y <- y != levels(y)[1L]
      y <- as.integer(y)
      if (!all(y %in% c(0L, 1L))) 
        stop("y values must be 0 or 1 for bernoulli regression")
    }
    else {
      if (!all(y >= 0 & y <= 1))
        stop("y values must be between 0 and 1 for binomial regression")
      trials <- weights
    }
  }
  
  x <- as.matrix(x)
  has_intercept <- colnames(x)[1L] == "(Intercept)"
  xtemp <- if (has_intercept) x[, -1L, drop=FALSE] else x
  xbar <- colMeans(xtemp)
  xtemp <- sweep(xtemp, 2, xbar, FUN = "-")
  nvars <- ncol(xtemp)
  
  
  scaled <- prior_ops$scaled
  min_prior_scale <- prior_ops$min_prior_scale
  prior_scale_for_dispersion <- prior_ops$prior_scale_for_dispersion
  
  ok_dists <- nlist("normal", student_t = "t", "cauchy", 
                    "horseshoe", "horseshoe_plus")
  ok_intercept_dists <- ok_dists[1:3]
  
  # prior distributions
  if (!is.null(prior)) {
    prior_dist <- prior$dist
    prior_scale <- prior$scale
    prior_mean <- prior$location
    prior_df <- prior$df
    prior_df[is.na(prior_df)] <- 1
    if (!prior_dist %in% unlist(ok_dists)) {
      stop("The prior distribution for the coefficients should be one of ",
           paste(names(ok_dists), collapse = ", "), call. = FALSE)
    }
    if (prior_dist == "normal" || prior_dist == "t") {
      prior_dist <- ifelse(prior_dist == "normal", 1L, 2L)
      prior_scale <- set_prior_scale(prior_scale, default = 2.5, link = family$link)
    }
    else prior_dist <- ifelse(prior_dist == "horseshoe", 3L, 4L)
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
    prior_dist_for_intercept <- prior_intercept$dist
    prior_scale_for_intercept <- prior_intercept$scale
    prior_mean_for_intercept <- prior_intercept$location
    prior_df_for_intercept <- prior_intercept$df 
    prior_df_for_intercept[is.na(prior_df_for_intercept)] <- 1
    
    if (!prior_dist_for_intercept %in% unlist(ok_intercept_dists)) {
      stop("The prior distribution for the intercept should be one of ",
           paste(names(ok_intercept_dists), collapse = ", "), call. = FALSE)
    }
    prior_dist_for_intercept <- 
      ifelse(prior_dist_for_intercept == "normal", 1L, 2L)
    prior_scale_for_intercept <- 
      set_prior_scale(prior_scale_for_intercept, default = 10, link = family$link)
    prior_df_for_intercept <- min(.Machine$double.xmax, prior_df_for_intercept)
  }
  else {
    prior_dist_for_intercept <- 0L
    prior_mean_for_intercept <- 0 
    prior_scale_for_intercept <- prior_df_for_intercept <- 1
  }
  
  is_gaussian <- family$family == "gaussian"
  is_gamma <- family$family == "Gamma"
  is_ig <- family$family == "inverse.gaussian"
  is_continuous <- is_gaussian || is_gamma || is_ig
  if (scaled && prior_dist > 0L) {
    if (is_gaussian) {
      ss <- 2 * sd(y)
      prior_scale <- ss * prior_scale
      prior_scale_for_intercept <-  ss * prior_scale_for_intercept
    }
    prior_scale <- pmax(min_prior_scale, prior_scale /
                          apply(xtemp, 2, FUN = function(x) {
                            num.categories <- length(unique(x))
                            x.scale <- 1
                            if (num.categories == 2) x.scale <- diff(range(x))
                            else if (num.categories > 2) x.scale <- 2 * sd(x)
                          }))
  }
  prior_scale <- as.array(pmin(.Machine$double.xmax, prior_scale))
  prior_scale_for_intercept <- min(.Machine$double.xmax, prior_scale_for_intercept)
  
  is_bernoulli <- supported_families[fam] == "binomial" && all(y %in% 0:1)
  is_nb <- supported_families[fam] == "neg_binomial_2"
  
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
    Z <- t(as.matrix(group$Zt))
    p <- sapply(group$cnms, FUN = length)
    l <- sapply(attributes(group$flist)$assign, function(i) nlevels(group$flist[,i]))
    t <- length(p)
    group_nms <- names(group$cnms)
    b_names <- unlist(lapply(1:t, FUN = function(i) {
      nms_i <- paste(group$cnms[[i]], group_nms[i])
      rep(nms_i, times = l[i])
    }))
    b_names <- paste(b_names, colnames(Z), sep = ":")
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
  if (is_continuous) {
    standata$prior_scale_for_dispersion <- if (prior_scale_for_dispersion == Inf) 
      0 else prior_scale_for_dispersion
    standata$family <- switch(family$family, 
                              gaussian = 1L, 
                              Gamma = 2L,
                              3L)
    stanfit <- get("stanfit_continuous")
  }
  else if (supported_families[fam] == "binomial") {
    standata$prior_scale_for_dispersion <- 
      if (!length(group) || prior_scale_for_dispersion == Inf) 
        0 else prior_scale_for_dispersion
  
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
      stanfit <- get("stanfit_bernoulli")  
    }
    else {
      standata$trials <- trials
      if (length(weights) & !all(weights == 1)) {
        standata$y <- round(y * trials)
        standata$weights <- double(0)
        standata$has_weights <- 0L
      }
      stanfit <- get("stanfit_binomial")
    }
  }   
  else if (supported_families[fam] == "poisson") {
    standata$family <- 1L
    standata$prior_scale_for_dispersion <- if (prior_scale_for_dispersion == Inf) 
      0 else prior_scale_for_dispersion
    stanfit <- get("stanfit_count") 
  }
  else if (is_nb) {
    standata$family <- 2L
    standata$prior_scale_for_dispersion <- if (prior_scale_for_dispersion == Inf) 
      0 else prior_scale_for_dispersion
    stanfit <- get("stanfit_count") 
  }
  else if (is_gamma) {
    # nothing
  }
  else stop(paste(family$family, "is not supported"))
  
  pars <- c(if (has_intercept) "alpha", "beta", 
            if (length(group)) "b",
            if (is_continuous) "dispersion", if (is_nb) "theta",  "mean_PPD")
  algorithm <- match.arg(algorithm)
  if (algorithm == "optimizing") {
    out <- rstan::optimizing(stanfit, data = standata, hessian = TRUE, ...)
    k <- ncol(out$hessian)
    rownames(out$hessian) <- colnames(out$hessian) <- head(names(out$par), k)
    new_names <- names(out$par)
    new_names[grepl("^beta\\[[[:digit:]]+\\]$", new_names)] <- colnames(xtemp)
    new_names[new_names == "alpha[1]"] <- "(Intercept)"
    new_names[new_names == "dispersion"] <- if (is_gaussian) "sigma" else
                                            if (is_gamma) "scale" else
                                            if (is_ig) "lambda" else NA
    if (is_nb) new_names[new_names == "theta[1]"] <- "overdispersion"
    if (length(group)) {
      new_names[grepl("^b\\[[[:digit:]]+\\]$", new_names)] <- paste0("b[", b_names, "]")
      # new_names[grepl("^var_group\\[[[:digit:]]+\\]$", new_names)] <- paste0("var[", g_names, "]")
      # rename theta_L ?
    }
    names(out$par) <- new_names
    out$cov.scaled <- qr.solve(-out$hessian, diag(1, k, k))
    colnames(out$cov.scaled) <- rownames(out$cov.scaled) <- colnames(out$hessian)
    out$stanfit <- suppressMessages(rstan::sampling(stanfit, data = standata, chains = 0))
    return(out)
  }
  else {
    if ("control" %in% names(list(...))) {
      stanfit <- rstan::sampling(stanfit, data = standata, pars = pars, 
                                 show_messages = FALSE, ...)
    }
    else stanfit <- rstan::sampling(stanfit, data = standata, pars = pars,
                                    control = stan_control, show_messages = FALSE, ...)
    
    # else
    #   stanfit <- rstan::vb(stanfit, pars = pars, data = standata, 
    #                        algorithm = algorithm, ...)
    new_names <- c(if (has_intercept) "(Intercept)", colnames(xtemp), 
                   if (length(group)) c(paste0("b[", b_names, "]")),
                                        # paste0("var[", g_names, "]")),
                   if (is_gaussian) "sigma", if (is_gamma) "shape", 
                   if (is_ig) "lambda",
                   if (is_nb) "overdispersion", 
                   "mean_PPD", "log-posterior")
    stanfit@sim$fnames_oi <- new_names
    return(stanfit)
  }
}
