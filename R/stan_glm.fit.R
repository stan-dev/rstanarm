# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2013, 2014, 2015, 2016 Trustees of Columbia University
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
#' @param prior_ops Deprecated. See \link{rstanarm-deprecated} for details.
#' @param group A list, possibly of length zero (the default), but otherwise
#'   having the structure of that produced by \code{\link[lme4]{mkReTrms}} to
#'   indicate the group-specific part of the model. In addition, this list must
#'   have elements for the \code{regularization}, \code{concentration} 
#'   \code{shape}, and \code{scale} components of a \code{\link{decov}}
#'   prior for the covariance matrices among the group-specific coefficients.
#' @importFrom lme4 mkVarCorr
stan_glm.fit <- function(x, y, weights = rep(1, NROW(x)), 
                         offset = rep(0, NROW(x)), family = gaussian(),
                         ...,
                         prior = normal(),
                         prior_intercept = normal(),
                         prior_dispersion = cauchy(0, 5),
                         prior_ops = NULL,
                         group = list(),
                         prior_PD = FALSE, 
                         algorithm = c("sampling", "optimizing", 
                                       "meanfield", "fullrank"), 
                         adapt_delta = NULL, QR = FALSE, sparse = FALSE) {
  
  # prior_ops deprecated but make sure it still works until 
  # removed in future release
  if (!is.null(prior_ops)) {
    tmp <- .support_deprecated_prior_options(prior, prior_intercept, 
                                             prior_dispersion, prior_ops)
    prior <- tmp[["prior"]]
    prior_intercept <- tmp[["prior_intercept"]]
    prior_dispersion <- tmp[["prior_dispersion"]]
    prior_ops <- NULL
  }
  
  algorithm <- match.arg(algorithm)
  family <- validate_family(family)
  supported_families <- c("binomial", "gaussian", "Gamma", "inverse.gaussian",
                          "poisson", "neg_binomial_2")
  fam <- which(pmatch(supported_families, family$family, nomatch = 0L) == 1L)
  if (!length(fam)) 
    stop("'family' must be one of ", paste(supported_families, collapse = ", "))
  
  supported_links <- switch(
    supported_families[fam],
    binomial = c("logit", "probit", "cauchit", "log", "cloglog"),
    gaussian = c("identity", "log", "inverse"),
    Gamma = c("identity", "log", "inverse"),
    inverse.gaussian = c("identity", "log", "inverse", "1/mu^2"),
    "neg_binomial_2" = , # intentional
    poisson = c("log", "identity", "sqrt"),
    stop("unsupported family")
  )
  link <- which(supported_links == family$link)
  if (!length(link)) 
    stop("'link' must be one of ", paste(supported_links, collapse = ", "))
  
  if (binom_y_prop(y, family, weights))
    stop("To specify 'y' as proportion of successes and 'weights' as ",
         "number of trials please use stan_glm rather than calling ",
         "stan_glm.fit directly.")
  if (is.binomial(family$family)) {
    if (NCOL(y) == 1L) {
      if (is.numeric(y) || is.logical(y)) 
        y <- as.integer(y)
      if (is.factor(y)) 
        y <- fac2bin(y)
      if (!all(y %in% c(0L, 1L))) 
        stop("y values must be 0 or 1 for bernoulli regression.")
    } else {
      if (!isTRUE(NCOL(y) == 2L))
        stop("y should either be a vector or a matrix 1 or 2 columns.")
      trials <- as.integer(y[, 1L] + y[, 2L])
      y <- as.integer(y[, 1L])
    }
  }

  # useless assignments to pass R CMD check
  has_intercept <- 
    prior_df <- prior_df_for_intercept <- prior_df_for_dispersion <-
    prior_dist <- prior_dist_for_intercept <- prior_dist_for_dispersion <- 
    prior_mean <- prior_mean_for_intercept <- prior_mean_for_dispersion <- 
    prior_scale <- prior_scale_for_intercept <- prior_scale_for_dispersion <- 
    prior_autoscale <- prior_autoscale_for_intercept <-
    global_prior_scale <- global_prior_df <- NULL
  
  x_stuff <- center_x(x, sparse)
  for (i in names(x_stuff)) # xtemp, xbar, has_intercept
    assign(i, x_stuff[[i]])
  nvars <- ncol(xtemp)

  ok_dists <- nlist("normal", student_t = "t", "cauchy", "hs", "hs_plus")
  ok_intercept_dists <- ok_dists[1:3]
  ok_dispersion_dists <- c(ok_dists[1:3], exponential = "exponential")
  
  # prior distributions
  prior_stuff <- handle_glm_prior(
    prior,
    nvars,
    link = family$link,
    default_scale = 2.5,
    ok_dists = ok_dists
  )
  # prior_{dist, mean, scale, df, dist_name, autoscale}, global_prior_df, global_prior_scale
  for (i in names(prior_stuff))
    assign(i, prior_stuff[[i]])
  
  prior_intercept_stuff <- handle_glm_prior(
    prior_intercept,
    nvars = 1,
    default_scale = 10,
    link = family$link,
    ok_dists = ok_intercept_dists
  )
  # prior_{dist, mean, scale, df, dist_name, autoscale}_for_intercept
  names(prior_intercept_stuff) <- paste0(names(prior_intercept_stuff), "_for_intercept")
  for (i in names(prior_intercept_stuff))
    assign(i, prior_intercept_stuff[[i]])
  
  prior_dispersion_stuff <-
    handle_glm_prior(
      prior_dispersion,
      nvars = 1,
      default_scale = 5,
      link = NULL, # don't need to adjust scale based on logit vs probit
      ok_dists = ok_dispersion_dists
    )
  # prior_{dist, mean, scale, df, dist_name, autoscale}_for_dispersion
  names(prior_dispersion_stuff) <- paste0(names(prior_dispersion_stuff), "_for_dispersion")
  for (i in names(prior_dispersion_stuff)) 
    assign(i, prior_dispersion_stuff[[i]])
  
  famname <- supported_families[fam]
  is_bernoulli <- is.binomial(famname) && all(y %in% 0:1)
  is_nb <- is.nb(famname)
  is_gaussian <- is.gaussian(famname)
  is_gamma <- is.gamma(famname)
  is_ig <- is.ig(famname)
  is_continuous <- is_gaussian || is_gamma || is_ig
  
  # require intercept for certain family and link combinations
  if (!has_intercept) {
    linkname <- supported_links[link]
    needs_intercept <- !is_gaussian && linkname == "identity" ||
      is_gamma && linkname == "inverse" ||
      is.binomial(famname) && linkname == "log"
    if (needs_intercept)
      stop("To use this combination of family and link ", 
           "the model must have an intercept.")
  }

  if (prior_dist > 0L) {
    if (is_gaussian) {
      ss <- 2 * sd(y)
      if (prior_autoscale) 
        prior_scale <- ss * prior_scale
      if (prior_autoscale_for_intercept && prior_dist_for_intercept > 0L) 
        prior_scale_for_intercept <-  ss * prior_scale_for_intercept
    }
    if (!QR && prior_autoscale) {
      min_prior_scale <- 1e-12 # used to be set in prior_options()
      prior_scale <- pmax(min_prior_scale, prior_scale / 
             apply(xtemp, 2L, FUN = function(x) {
               num.categories <- length(unique(x))
               x.scale <- 1
               if (num.categories == 2) x.scale <- diff(range(x))
               else if (num.categories > 2) x.scale <- 2 * sd(x)
               return(x.scale)
             }))
    }
  }
  prior_scale <- 
    as.array(pmin(.Machine$double.xmax, prior_scale))
  prior_scale_for_intercept <- 
    min(.Machine$double.xmax, prior_scale_for_intercept)
  
  if (QR) {
    if (ncol(xtemp) <= 1)
      stop("'QR' can only be specified when there are multiple predictors.")
    if (sparse)
      stop("'QR' and 'sparse' cannot both be TRUE")
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
  standata <- nlist(
    N = nrow(xtemp),
    K = ncol(xtemp),
    xbar = as.array(xbar),
    dense_X = !sparse,
    link,
    has_weights = length(weights) > 0,
    has_offset = length(offset) > 0,
    has_intercept,
    prior_PD,
    prior_dist,
    prior_mean,
    prior_scale,
    prior_df,
    prior_dist_for_intercept,
    prior_scale_for_intercept = c(prior_scale_for_intercept),
    prior_mean_for_intercept = c(prior_mean_for_intercept),
    prior_df_for_intercept = c(prior_df_for_intercept), 
    global_prior_df, global_prior_scale, # for hs prior
    has_intercept, prior_PD
    # mean,df,scale for dispersion added below depending on family
  )

  # make a copy of user specification before modifying 'group' (used for keeping
  # track of priors)
  user_covariance <- if (!length(group)) NULL else group[["decov"]]
  
  if (length(group)) {
    check_reTrms(group)
    decov <- group$decov
    Z <- t(group$Zt)
    group <- pad_reTrms(Ztlist = group$Ztlist, cnms = group$cnms, flist = group$flist)
    Z <- group$Z
    p <- sapply(group$cnms, FUN = length)
    l <- sapply(attr(group$flist, "assign"), function(i) 
      nlevels(group$flist[[i]]))
    t <- length(l)
    b_nms <- make_b_nms(group)
    g_nms <- unlist(lapply(1:t, FUN = function(i) {
      paste(group$cnms[[i]], names(group$cnms)[i], sep = "|")
    }))
    standata$t <- t
    standata$p <- as.array(p)
    standata$l <- as.array(l)
    standata$q <- ncol(Z)
    standata$len_theta_L <- sum(choose(p, 2), p)
    if (is_bernoulli) {
      parts0 <- extract_sparse_parts(Z[y == 0, , drop = FALSE])
      parts1 <- extract_sparse_parts(Z[y == 1, , drop = FALSE])
      standata$num_non_zero <- c(length(parts0$w), length(parts1$w))
      standata$w0 <- parts0$w
      standata$w1 <- parts1$w
      standata$v0 <- parts0$v
      standata$v1 <- parts1$v
      standata$u0 <- parts0$u
      standata$u1 <- parts1$u
    } else {
      parts <- extract_sparse_parts(Z)
      standata$num_non_zero <- length(parts$w)
      standata$w <- parts$w
      standata$v <- parts$v
      standata$u <- parts$u
    }
    standata$shape <- as.array(maybe_broadcast(decov$shape, t))
    standata$scale <- as.array(maybe_broadcast(decov$scale, t))
    standata$len_concentration <- sum(p[p > 1])
    standata$concentration <- 
      as.array(maybe_broadcast(decov$concentration, sum(p[p > 1])))
    standata$len_regularization <- sum(p > 1)
    standata$regularization <- 
      as.array(maybe_broadcast(decov$regularization, sum(p > 1)))
    
  } else { # !length(group)
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
    } else {
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
    if (sparse) {
      parts <- extract_sparse_parts(xtemp)
      standata$nnz_X <- length(parts$w)
      standata$w_X <- parts$w
      standata$v_X <- parts$v
      standata$u_X <- parts$u
      standata$X <- array(0, dim = c(0L, dim(xtemp)))
    } else {
      standata$X <- array(xtemp, dim = c(1L, dim(xtemp)))
      standata$nnz_X <- 0L
      standata$w_X <- double(0)
      standata$v_X <- integer(0)
      standata$u_X <- integer(0)
    }
    standata$y <- y
    standata$weights <- weights
    standata$offset <- offset
  }

  # call stan() to draw from posterior distribution
  if (is_continuous) {
    standata$prior_scale_for_dispersion <- prior_scale_for_dispersion %ORifINF% 0
    standata$prior_df_for_dispersion <- c(prior_df_for_dispersion)
    standata$prior_mean_for_dispersion <- c(prior_mean_for_dispersion)
    standata$family <- switch(family$family, 
                              gaussian = 1L, 
                              Gamma = 2L,
                              3L)
    stanfit <- stanmodels$continuous
  } else if (is.binomial(famname)) {
    standata$prior_scale_for_dispersion <- 
      if (!length(group) || prior_scale_for_dispersion == Inf) 
        0 else prior_scale_for_dispersion
    standata$prior_mean_for_dispersion <- 0
    standata$prior_df_for_dispersion <- 0
    standata$family <- 1L # not actually used
    if (is_bernoulli) {
      y0 <- y == 0
      y1 <- y == 1
      standata$N <- c(sum(y0), sum(y1))
      if (sparse) {
        standata$X0 <- array(0, dim = c(0L, sum(y0), ncol(xtemp)))
        standata$X1 <- array(0, dim = c(0L, sum(y1), ncol(xtemp)))
        parts0 <- extract_sparse_parts(xtemp[y0, , drop = FALSE])
        standata$nnz_X0 <- length(parts0$w)
        standata$w_X0 = parts0$w
        standata$v_X0 = parts0$v
        standata$u_X0 = parts0$u
        parts1 <- extract_sparse_parts(xtemp[y1, , drop = FALSE])
        standata$nnz_X1 <- length(parts1$w)
        standata$w_X1 = parts1$w
        standata$v_X1 = parts1$v
        standata$u_X1 = parts1$u
      } else {
        standata$X0 <- array(xtemp[y0, , drop = FALSE], dim = c(1, sum(y0), ncol(xtemp)))
        standata$X1 <- array(xtemp[y1, , drop = FALSE], dim = c(1, sum(y1), ncol(xtemp)))
        standata$nnz_X0 = 0L 
        standata$w_X0 = double(0)
        standata$v_X0 = integer(0)
        standata$u_X0 = integer(0)
        standata$nnz_X1 = 0L 
        standata$w_X1 = double(0)
        standata$v_X1 = integer(0)
        standata$u_X1 = integer(0)
      }
      if (length(weights)) { 
        # nocov start
        # this code is unused because weights are interpreted as number of 
        # trials for binomial glms
        standata$weights0 <- weights[y0]
        standata$weights1 <- weights[y1]
        # nocov end
      } else {
        standata$weights0 <- double(0)
        standata$weights1 <- double(0)
      }
      if (length(offset)) {
        # nocov start
        standata$offset0 <- offset[y0]
        standata$offset1 <- offset[y1]
        # nocov end
      } else {
        standata$offset0 <- double(0)
        standata$offset1 <- double(0)
      }
      stanfit <- stanmodels$bernoulli
    } else {
      standata$trials <- trials
      stanfit <- stanmodels$binomial
    }
  } else if (is.poisson(famname)) {
    standata$family <- 1L
    standata$prior_scale_for_dispersion <- prior_scale_for_dispersion %ORifINF% 0
    standata$prior_mean_for_dispersion <- 0
    standata$prior_df_for_dispersion <- 0
    stanfit <- stanmodels$count 
  } else if (is_nb) {
    standata$family <- 2L
    standata$prior_scale_for_dispersion <- prior_scale_for_dispersion %ORifINF% 0
    standata$prior_df_for_dispersion <- c(prior_df_for_dispersion)
    standata$prior_mean_for_dispersion <- c(prior_mean_for_dispersion)
    stanfit <- stanmodels$count
  } else if (is_gamma) {
    # nothing
  } else { # nocov start
    # family already checked above
    stop(paste(famname, "is not supported."))
  } # nocov end
  
  prior_info <- summarize_glm_prior(
    user_prior = prior_stuff,
    user_prior_intercept = prior_intercept_stuff,
    user_prior_dispersion = prior_dispersion_stuff,
    user_prior_covariance = user_covariance,
    has_intercept = has_intercept,
    has_predictors = nvars > 0,
    adjusted_prior_scale = prior_scale,
    adjusted_prior_intercept_scale = prior_scale_for_intercept,
    family = family
  )
  
  pars <- c(if (has_intercept) "alpha", 
            "beta", 
            if (length(group)) "b",
            if (is_continuous | is_nb) "dispersion",
            if (standata$len_theta_L) "theta_L",
            "mean_PPD")
  if (algorithm == "optimizing") {
    out <- optimizing(stanfit, data = standata, 
                      draws = 1000, constrained = TRUE, ...)
    new_names <- names(out$par)
    mark <- grepl("^beta\\[[[:digit:]]+\\]$", new_names)
    if (QR) {
      out$par[mark] <- R_inv %*% out$par[mark]
      out$theta_tilde[,mark] <- out$theta_tilde[, mark] %*% t(R_inv)
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
    out$stanfit <- suppressMessages(sampling(stanfit, data = standata, 
                                             chains = 0))
    return(structure(out, prior.info = prior_info))
    
  } else {
    if (algorithm == "sampling") {
      sampling_args <- set_sampling_args(
        object = stanfit, 
        prior = prior, 
        user_dots = list(...), 
        user_adapt_delta = adapt_delta, 
        data = standata, 
        pars = pars, 
        show_messages = FALSE)
      stanfit <- do.call(sampling, sampling_args)
    } else {
      # meanfield or fullrank vb
      stanfit <- rstan::vb(stanfit, pars = pars, data = standata,
                           algorithm = algorithm, ...)
      if (algorithm == "meanfield" && !QR) 
        msg_meanfieldQR()
    }
    if (QR) {
      thetas <- extract(stanfit, pars = "beta", inc_warmup = TRUE, 
                        permuted = FALSE)
      betas <- apply(thetas, 1:2, FUN = function(theta) R_inv %*% theta)
      end <- tail(dim(betas), 1L)
      for (chain in 1:end) for (param in 1:nrow(betas)) {
        stanfit@sim$samples[[chain]][[has_intercept + param]] <-
          if (ncol(xtemp) > 1) betas[param, , chain] else betas[param, chain]
      }
    }
    if (standata$len_theta_L) {
      thetas <- extract(stanfit, pars = "theta_L", inc_warmup = TRUE, 
                        permuted = FALSE)
      cnms <- group$cnms
      nc <- sapply(cnms, FUN = length)
      nms <- names(cnms)
      Sigma <- apply(thetas, 1:2, FUN = function(theta) {
        Sigma <- mkVarCorr(sc = 1, cnms, nc, theta, nms)
        unlist(sapply(Sigma, simplify = FALSE, 
                      FUN = function(x) x[lower.tri(x, TRUE)]))
      })
      l <- length(dim(Sigma))
      end <- tail(dim(Sigma), 1L)
      shift <- grep("^theta_L", names(stanfit@sim$samples[[1]]))[1] - 1L
      if (l == 3) for (chain in 1:end) for (param in 1:nrow(Sigma)) {
        stanfit@sim$samples[[chain]][[shift + param]] <- Sigma[param, , chain] 
      }
      else for (chain in 1:end) {
        stanfit@sim$samples[[chain]][[shift + 1]] <- Sigma[, chain]
      }
      Sigma_nms <- lapply(cnms, FUN = function(grp) {
        nm <- outer(grp, grp, FUN = paste, sep = ",")
        nm[lower.tri(nm, diag = TRUE)]
      })
      for (j in seq_along(Sigma_nms)) {
        Sigma_nms[[j]] <- paste0(nms[j], ":", Sigma_nms[[j]])
      }
      Sigma_nms <- unlist(Sigma_nms)
    }
    new_names <- c(if (has_intercept) "(Intercept)", 
                   colnames(xtemp), 
                   if (length(group)) c(paste0("b[", b_nms, "]")),
                   if (is_gaussian) "sigma", 
                   if (is_gamma) "shape", 
                   if (is_ig) "lambda",
                   if (is_nb) "overdispersion", 
                   if (standata$len_theta_L) paste0("Sigma[", Sigma_nms, "]"),
                   "mean_PPD", 
                   "log-posterior")
    stanfit@sim$fnames_oi <- new_names
    return(structure(stanfit, prior.info = prior_info))
  }
}


# Add extra level _NEW_ to each group
# 
# @param Ztlist ranef indicator matrices
# @param cnms group$cnms
# @param flist group$flist
#' @importFrom Matrix rBind
pad_reTrms <- function(Ztlist, cnms, flist) {
  stopifnot(is.list(Ztlist))
  l <- sapply(attr(flist, "assign"), function(i) nlevels(flist[[i]]))
  p <- sapply(cnms, FUN = length)
  n <- ncol(Ztlist[[1]])
  for (i in attr(flist, "assign")) {
    if (grepl("^Xr", names(p)[i])) next
    levels(flist[[i]]) <- c(gsub(" ", "_", levels(flist[[i]])), 
                            paste0("_NEW_", names(flist)[i]))
  }
  for (i in 1:length(p)) {
    if (grepl("^Xr", names(p)[i])) next
    Ztlist[[i]] <- if (getRversion() < "3.2.0") {
      rBind( Ztlist[[i]], Matrix(0, nrow = p[i], ncol = n, sparse = TRUE))
    } else {
      rbind2(Ztlist[[i]], Matrix(0, nrow = p[i], ncol = n, sparse = TRUE))
    }
  }
  Z <- t(do.call(rbind, args = Ztlist))
  return(nlist(Z, cnms, flist))
}

# Drop the extra reTrms from a matrix x
#
# @param x A matrix or array (e.g. the posterior sample or matrix of summary
#   stats)
# @param columns Do the columns (TRUE) or rows (FALSE) correspond to the 
#   variables?
unpad_reTrms <- function(x, ...) UseMethod("unpad_reTrms")
unpad_reTrms.default <- function(x, ...) {
  if (is.matrix(x) || is.array(x))
    return(unpad_reTrms.array(x, ...))
  keep <- !grepl("_NEW_", names(x), fixed = TRUE)
  x[keep]
}

unpad_reTrms.array <- function(x, columns = TRUE, ...) {
  ndim <- length(dim(x))
  if (ndim > 3)
    stop("'x' should be a matrix or 3-D array")
  
  nms <- if (columns) 
    last_dimnames(x) else rownames(x)
  keep <- !grepl("_NEW_", nms, fixed = TRUE)
  if (length(dim(x)) == 2) {
    x_keep <- if (columns) 
      x[, keep, drop = FALSE] else x[keep, , drop = FALSE]
  } else {
    x_keep <- if (columns) 
      x[, , keep, drop = FALSE] else x[keep, , , drop = FALSE]
  }
  return(x_keep)
}

make_b_nms <- function(group) {
  group_nms <- names(group$cnms)
  b_nms <- character()
  for (i in seq_along(group$cnms)) {
    nm <- group_nms[i]
    nms_i <- paste(group$cnms[[i]], nm)
    if (length(nms_i) == 1) {
      b_nms <- c(b_nms, paste0(nms_i, ":", levels(group$flist[[nm]])))
    } else {
      b_nms <- c(b_nms, c(t(sapply(nms_i, paste0, ":", levels(group$flist[[nm]])))))
    }
  }
  return(b_nms)  
}


# Create "prior.info" attribute needed for prior_summary()
#
# @param user_* The user's prior, prior_intercept, prior_covariance, and 
#   prior_dispersion specifications. For prior and prior_intercept these should be
#   passed in after broadcasting the df/location/scale arguments if necessary.
# @param has_intercept T/F, does model have an intercept?
# @param has_predictors T/F, does model have predictors?
# @param adjusted_prior_* adjusted scales computed if using autoscaled priors
# @param family Family object.
# @return A named list with components 'prior', 'prior_intercept', and possibly 
#   'prior_covariance' and 'prior_dispersion' each of which itself is a list
#   containing the needed values for prior_summary.
summarize_glm_prior <-
  function(user_prior,
           user_prior_intercept,
           user_prior_dispersion,
           user_prior_covariance,
           has_intercept, 
           has_predictors,
           adjusted_prior_scale,
           adjusted_prior_intercept_scale, 
           family) {
    rescaled_coef <-
      user_prior$prior_autoscale && 
      has_predictors &&
      !is.na(user_prior$prior_dist_name) &&
      !all(user_prior$prior_scale == adjusted_prior_scale)
    rescaled_int <-
      user_prior_intercept$prior_autoscale_for_intercept &&
      has_intercept &&
      !is.na(user_prior_intercept$prior_dist_name_for_intercept) &&
      (user_prior_intercept$prior_scale_for_intercept != adjusted_prior_intercept_scale)
    
    if (has_predictors && user_prior$prior_dist_name %in% "t") {
      if (all(user_prior$prior_df == 1)) {
        user_prior$prior_dist_name <- "cauchy"
      } else {
        user_prior$prior_dist_name <- "student_t"
      }
    }
    if (has_intercept &&
        user_prior_intercept$prior_dist_name_for_intercept %in% "t") {
      if (all(user_prior_intercept$prior_df_for_intercept == 1)) {
        user_prior_intercept$prior_dist_name_for_intercept <- "cauchy"
      } else {
        user_prior_intercept$prior_dist_name_for_intercept <- "student_t"
      }
    }
    if (user_prior_dispersion$prior_dist_name_for_dispersion %in% "t") {
      if (all(user_prior_dispersion$prior_df_for_dispersion == 1)) {
        user_prior_dispersion$prior_dist_name_for_dispersion <- "cauchy"
      } else {
        user_prior_dispersion$prior_dist_name_for_dispersion <- "student_t"
      }
    }
    prior_list <- list(
      prior = 
        if (!has_predictors) NULL else with(user_prior, list(
          dist = prior_dist_name,
          location = prior_mean,
          scale = prior_scale,
          adjusted_scale = if (rescaled_coef)
            adjusted_prior_scale else NULL,
          df = if (prior_dist_name %in% c("student_t", "hs", "hs_plus"))
            prior_df else NULL
        )),
      prior_intercept = 
        if (!has_intercept) NULL else with(user_prior_intercept, list(
          dist = prior_dist_name_for_intercept,
          location = prior_mean_for_intercept,
          scale = prior_scale_for_intercept,
          adjusted_scale = if (rescaled_int)
            adjusted_prior_intercept_scale else NULL,
          df = if (prior_dist_name_for_intercept %in% "student_t")
            prior_df_for_intercept else NULL
        ))
    )
    if (length(user_prior_covariance))
      prior_list$prior_covariance <- user_prior_covariance
    
    dispersion_name <- .rename_dispersion(family)
    prior_list$prior_dispersion <- if (is.na(dispersion_name)) 
      NULL else with(user_prior_dispersion, list(
        dist = prior_dist_name_for_dispersion,
        location = if (prior_dist_name_for_dispersion != "exponential")
          prior_mean_for_dispersion else NULL,
        scale = if (prior_dist_name_for_dispersion != "exponential")
          prior_scale_for_dispersion else NULL,
        df = if (prior_dist_name_for_dispersion %in% "student_t")
          prior_df_for_dispersion else NULL, 
        rate = if (prior_dist_name_for_dispersion %in% "exponential")
          1 / prior_scale_for_dispersion else NULL,
        dispersion_name = dispersion_name
      ))
      
    return(prior_list)
  }

.rename_dispersion <- function(family) {
  fam <- family$family
  if (is.gaussian(fam)) "sigma" else
    if (is.gamma(fam)) "shape" else
      if (is.ig(fam)) "lambda" else 
        if (is.nb(fam)) "overdispersion" else NA
}
