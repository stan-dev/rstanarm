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

stanglm <-
  function(formula, family = gaussian(), data, weights, subset, 
           na.action, start = NULL, etastart, mustart, offset, control = list(), 
           model = TRUE, method = "glm.fit", x = FALSE, y = TRUE, contrasts = NULL, 
           # above arguments from glm(), below arguments from arm:::bayesglm()
           prior.mean = 0, prior.scale = NULL, prior.df = 1,
           prior.mean.for.intercept = 0, prior.scale.for.intercept = NULL,
           prior.df.for.intercept = 1, min.prior.scale = 1e-12, scaled = TRUE,
           prior.scale.for.dispersion = 5, ...) { # further arguments to stan()

  # Parse like glm()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action", 
               "etastart", "mustart", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  Y <- model.response(mf, "any")
  if (length(dim(Y)) == 1L) {
    nm <- rownames(Y)
    dim(Y) <- NULL
    if (!is.null(nm)) 
      names(Y) <- nm
  }
  if (!is.empty.model(mt)) X <- model.matrix(mt, mf, contrasts) 
  else X <- matrix(, NROW(Y), 0L)
  weights <- as.vector(model.weights(mf))
  if (!is.null(weights) && !is.numeric(weights)) stop("'weights' must be a numeric vector")
  if (!is.null(weights) && any(weights < 0))     stop("negative weights not allowed")
  if ( is.null(weights)) weights <- rep(1.0, NROW(Y))
  offset <- as.vector(model.offset(mf))
  if (!is.null(offset)) {
    if (length(offset) != NROW(Y)) 
      stop(gettextf("number of offsets is %d should equal %d (number of observations)", 
                    length(offset), NROW(Y)), domain = NA)
  }
  else offset <- rep(0, nrow(X))
  
  stanglm.fit(X, Y, weights, start, etastart, mustart, offset, family, list(), TRUE,
              prior.mean, prior.scale, prior.df,
              prior.mean.for.intercept, prior.scale.for.intercept,
              prior.df.for.intercept, min.prior.scale, 
              scaled, prior.scale.for.dispersion, ...)
}

stanglm.fit <-
  function(x, y, weights = rep(1, NROW(x)), start = NULL, etastart = NULL, 
           mustart = NULL, offset = rep(0, NROW(x)), family = gaussian(), 
           control = list(), intercept = TRUE,
           # above arguments from glm(), below arguments from arm:::bayesglm()
           prior.mean = 0, prior.scale = NULL, prior.df = 1,
           prior.mean.for.intercept = 0, prior.scale.for.intercept = NULL,
           prior.df.for.intercept = 1, min.prior.scale = 1e-12, scaled = TRUE,
           prior.scale.for.dispersion = 5, ...) { # further arguments to stan()
  
  if(!is(family, "family")) stop("'family' must be a family")
  x <- as.matrix(x)
  
  # these are from help(family)
  supported_families <- c("binomial", "gaussian", "Gamma", 
                          "inverse.gaussian","poisson")
  fam <- which(supported_families == family$family)
  if(length(fam) == 0) {
    stop("'family' must be one of ", supported_families)
  }
  
  # these are also from help(family)
  supported_links <- switch(supported_families[fam],
                            binomial = c("logit", "probit", "cauchit", "log", "cloglog"),
                            gaussian = c("identity", "log", "inverse"),
                            Gamma = c("inverse", "identity", "log"),
                            inverse.gaussian = c("1/mu^2", "inverse", "identity", "log"),
                            poisson = c("log", "identity", "sqrt"),
                            stop("unsupported family")
  )
  link <- which(supported_links == family$link) 
  if(length(link) == 0) {
    stop("'link' must be one of ", supported_links)
  }
  
  # process hyperprior values like arm:::bayesglm() does
  if(is.null(prior.scale)) {
    prior.scale <- 2.5
    if(family$link == "probit") prior.scale <- prior.scale * dnorm(0) / dlogis(0)
  }
  if(is.null(prior.scale.for.intercept)) {
    prior.scale.for.intercept <- 10
    if(family$link == "probit") {
      prior.scale.for.intercept <- prior.scale.for.intercept * dnorm(0) / dlogis(0)
    }
  }
  nvars <- ncol(x) - 1
  if(length(prior.df) == 1)    prior.df <- rep(prior.df, nvars)
  prior.df <- as.array(pmin(.Machine$double.xmax, prior.df))
  prior.df.for.intercept <- min(.Machine$double.xmax, prior.df.for.intercept)
  if(length(prior.mean)  == 1) prior.mean  <- rep(prior.mean,  nvars)
  prior.mean <- as.array(prior.mean)
  if(length(prior.scale) == 1) prior.scale <- rep(prior.scale, nvars)
  if(scaled) {
    if(family$family == "gaussian") {
      prior.scale <- prior.scale * 2 * sd(y)
      prior.scale.for.intercept <- prior.scale.for.intercept * 2 * sd(y)
    }
    prior.scale <- pmax(min.prior.scale, prior.scale / 
                          apply(x[,-1,drop=FALSE], 2, FUN = function(x) {
                            num.categories <- length(unique(x))
                            x.scale <- 1
                            if(num.categories == 2)     x.scale <- diff(range(x))
                            else if(num.categories > 2) x.scale <- 2 * sd(x)
                          }))
  }
  prior.scale <- as.array(pmin(.Machine$double.xmax, prior.scale))
  priors.scale.for.intercept <- min(.Machine$double.xmax, prior.scale.for.intercept)
  
  # create entries in the data {} block of the .stan file
  standata <- list(N = nrow(x),
                   K = ncol(x),
                   X = x,
                   y = y,
                   family = fam,
                   link = link,
                   has_weights = as.integer(!all(weights == 1)),
                   weights = weights,
                   has_offset = as.integer(!all(offset == 0)),
                   offset = offset,
                   prior_scale = prior.scale,
                   prior_scale_for_intercept = prior.scale.for.intercept,
                   prior_mean = prior.mean,
                   prior_mean_for_intercept = prior.mean.for.intercept,
                   prior_df = prior.df,
                   prior_df_for_intercept = prior.df.for.intercept)
  
  if(family$family == "gaussian") {
    standata$prior_scale_for_dispersion <- prior.scale.for.dispersion
  }
  
  # call stan() to draw from posterior distribution
  if(supported_families[fam] == "gaussian") stanfit <- stanfit_gaussian
  else if(supported_families[fam] %in% c("binomial", "poisson")) stanfit <- stanfit_discrete
  else stop("model not supported yet") # FIXME           
  
  if(is.null(start)) start <- "random"
  else start <- as.list(start)
  stanfit <- rstan:::stan(fit = stanfit, data = standata, init = start, ...)
  betas <- grepl("beta[", dimnames(stanfit)$parameters, fixed = TRUE)
  stanfit@sim$fnames_oi[betas] <- colnames(x)
  return(stanfit)
}
