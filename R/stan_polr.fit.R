# This file is part of rstanarm.
# Copyright 2015 Stan Development Team
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

#' @rdname stan_polr
#' @export
#' @param x A design matrix.
#' @param y A response variable, which must be a (preferably ordered) factor.
#' @param wt A numeric vector (possibly \code{NULL}) of observation weights.
#' @param offset A numeric vector (possibly \code{NULL}) of offsets.
#' 
#' @importFrom utils head tail
stan_polr.fit <- function (x, y, wt = NULL, offset = NULL, 
                           method = c("logistic", "probit", "loglog", 
                                      "cloglog", "cauchit"), ...,
                           prior = R2(stop("'location' must be specified")), 
                           prior_counts = dirichlet(1), prior_PD = FALSE, 
                           algorithm = c("sampling", "optimizing")) {
  algorithm <- match.arg(algorithm)
  method <- match.arg(method)
  link <- which(c("logistic", "probit", "loglog", "cloglog", "cauchit") == method)
  if (!is.factor(y)) stop("'y' must be a factor")
  y_lev <- levels(y)
  J <- length(y_lev)
  y <- as.integer(y)
  if (colnames(x)[1] == "(Intercept)") x <- x[,-1,drop=FALSE]
  xbar <- as.array(colMeans(x))
  X <- sweep(x, 2, xbar, FUN = "-")
  s_X <- as.array(apply(X, 2, sd))
  has_weights <- length(wt) > 0 && !all(wt == 1)
  if (!has_weights) weights <- double(0)
  has_offset <- length(offset) > 0 && !all(offset == 0)
  if (!has_offset) offset <- double(0)

  if (length(prior) > 0) {
    shape <- make_eta(prior$location, prior$what, K = ncol(x))
    prior_dist <- 1L
  }
  else {
    shape <- 0
    prior_dist <- 0L
  }
  
  if (length(prior_counts) == 0) prior_counts <- rep(1,J)
  else prior_counts <- maybe_broadcast(prior_counts$concentration, J)
  
  N <- nrow(X)
  K <- ncol(X)
  
  standata <- nlist(J, N, K, X, xbar, s_X, y, prior_PD, link, 
                    has_weights, weights, has_offset, offset,
                    prior_dist, shape, prior_counts)
  if (prior_dist == 1) {
    L <- t(chol(cor(x)))
    dim(L) <- c(1L, dim(L))
    halfK <- K / 2
    R2 <- as.array(halfK / (halfK + standata$shape))
    z_beta <- NULL
  }
  else {
    L <- array(0, dim = c(0L, K, K))
    R2 <- array(0, dim = 0L)
    z_beta <- rep(0, K)
  }
  pi <- table(y) / N               
  start <- function(chain_id) {
    list(pi = pi, L = L, R2 = R2, z_beta = z_beta)
  }

  stanfit <- stanmodels$polr
  if (algorithm == "optimizing") {
    standata$do_residuals <- 0L
    out <- optimizing(stanfit, data = standata, hessian = TRUE, init = start)
    new_names <- c(paste0("pi[", y_lev, "]"), paste0("z_beta[", colnames(x), "]"),
                   "Delta_y", colnames(x), paste0("cutpoints[", y_lev[-1], "]"),
                   paste(head(y_lev, -1), tail(y_lev, -1), sep = "|"),
                   paste("mean_PPD", y_lev, sep = ":"))
    names(out$par) <- new_names
    K <- ncol(out$hessian)
    out$cov.scaled <- qr.solve(-out$hessian, diag(1, K , K))
    colnames(out$cov.scaled) <- rownames(out$cov.scaled)
    out$stanfit <- suppressMessages(sampling(stanfit, data = standata, chains = 0))
    return(out)
  }
  else {
    if (J > 2) pars <- c("beta", "zeta", "mean_PPD")
    else       pars <- c("zeta", "beta", "mean_PPD")
    standata$do_residuals <- J > 2
    if ("control" %in% names(list(...))) {
      stanfit <- sampling(stanfit, data = standata, pars = pars, 
                                 init = start, show_messages = FALSE, ...)
    }
    else stanfit <- sampling(stanfit, data = standata, pars = pars, init = start,
                                    control = stan_control, show_messages = FALSE, ...)
    
    # else 
    #   stanfit <- vb(stanfit, pars = pars, data = standata, algorithm = algorithm, ...)
    if (J > 2)
      new_names <- c(colnames(x), paste(head(y_lev, -1), tail(y_lev, -1), sep = "|"),
                     paste("mean_PPD", y_lev, sep = ":"), "log-posterior")
    else new_names <- c("(Intercept)", colnames(x), "mean_PPD", "log-posterior")
    stanfit@sim$fnames_oi <- new_names
    return(stanfit)
  }
}
