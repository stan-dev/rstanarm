# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2017 Trustees of Columbia University
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

#' Conditional logistic (clogit) models via Stan
#'
#' clogit
#'
#' @export
#' @templateVar pkg survival
#' @templateVar pkgfun clogit
#' @templateVar sameargs model,offset
#' @templateVar rareargs na.action,contrasts
#' @templateVar fun stan_clogit
#' @templateVar fitfun stan_glm.fit
#' @template return-stanreg-object
#' @template see-also
#' @template args-priors
#' @template args-prior_PD
#' @template args-algorithm
#' @template args-adapt_delta
#' @template args-QR
#' @template args-sparse
#' 
#' @param formula,data,subset,na.action Same as for \code{\link[lme4]{glmer}}, 
#'   except that any intercept included in the formula will be dropped. \emph{We
#'   strongly advise against omitting the \code{data} argument}. Unless 
#'   \code{data} is specified (and is a data frame) many post-estimation 
#'   functions (including \code{update}, \code{loo}, \code{kfold}) are not 
#'   guaranteed to work properly.
#' @param strata A factor indicating the groups in the data where the number of
#'   successes (possibly one) is fixed by the research design. It may be useful
#'   to use \code{\link{interaction}} or \code{\link[survival]{strata}} to create
#'   this factor
#' 
#' @details The \code{stan_clogit} function is mostly similar in syntax to 
#'   \code{\link[survival]{clogit}} but rather than performing maximum likelihood 
#'   estimation of generalized linear models, full Bayesian estimation is 
#'   performed (if \code{algorithm} is \code{"sampling"}) via MCMC. The Bayesian
#'   model adds priors (independent by default) on the coefficients of the GLM.
#'   
#' @seealso The Bernoulli vignette
#' 
#' @examples
#' stan_clogit(case ~ spontaneous + induced, strata = stratum,
#'             data = infert[order(infert$stratum, !infert$case),], QR = TRUE)
#'
stan_clogit <- function(formula, data, subset, na.action = NULL, ..., 
                        strata, prior = normal(), 
                        prior_covariance = decov(), prior_PD = FALSE, 
                        algorithm = c("sampling", "optimizing", 
                                      "meanfield", "fullrank"),
                        adapt_delta = NULL, QR = FALSE, sparse = FALSE) {
  
  algorithm <- match.arg(algorithm)
  data <- validate_data(data, if_missing = environment(formula))
  call <- match.call(expand.dots = TRUE)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action", "strata"), 
             table = names(mf), nomatch = 0L)
  mf <- mf[c(1L, m)]
  names(mf)[length(mf)] <- "weights"
  has_bars <- length(findbars(formula)) > 0
  if (has_bars) {
    if (is.null(prior_covariance))
      stop("'prior_covariance' can't be NULL.", call. = FALSE)
    mf[[1L]] <- quote(lme4::glFormula)
    mf$control <- make_glmerControl()
    glmod <- eval(mf, parent.frame())
    X <- glmod$X
    mf <- glmod$fr
    Y <- mf[, as.character(glmod$formula[2L])]
    group <- glmod$reTrms
    group$strata <- glmod$strata <- as.factor(mf[,"(weights)"])
    group$decov <- prior_covariance
  }
  else {
    validate_glm_formula(formula)
    mf[[1L]] <- as.name("model.frame")
    mf$drop.unused.levels <- TRUE
    mf <- eval(mf, parent.frame())
    group <- list(strata = as.factor(mf[,"(weights)"]))
    mt <- attr(mf, "terms")
    X <- model.matrix(mt, mf, contrasts)
    Y <- array1D_check(model.response(mf, type = "any"))
  }
  ord <- order(group$strata, !Y)
  if (!identical(ord, 1:NROW(Y)))
    stop("data must be sorted by 'stratum' and successes before failures within 'stratum'")
  offset <- model.offset(mf) %ORifNULL% double(0)
  weights <- double(0)
  mf <- check_constant_vars(mf)
  mt <- attr(mf, "terms")
  if (is.empty.model(mt))
    stop("No intercept or predictors specified.", call. = FALSE)
  xint <- match("(Intercept)", colnames(X), nomatch = 0L)
  if (xint > 0L) {
    X <- X[, -xint, drop = FALSE]
    attr(mt, "intercept") <- 0L
  }
  f <- binomial(link = "logit")
  stanfit <- stan_glm.fit(x = X, y = Y, weights = weights, 
                          offset = offset, family = f,
                          prior = prior, 
                          prior_PD = prior_PD, 
                          algorithm = algorithm, adapt_delta = adapt_delta, 
                          group = group, QR = QR, sparse = sparse, ...)
  f$link <- "clogit"
  f$linkinv <- function(eta, g = group$strata, 
                        successes = aggregate(Y, by = list(g), FUN = sum)$x) {
    denoms <- unlist(lapply(1:length(successes), FUN = function(j) {
      mark <- g == levels(g)[j]
      log_clogit_denom(sum(mark), successes[j], eta[mark])
    }))
    exp(eta - denoms[as.integer(g)])
  }
  f$linkfun <- log
  f$mu.eta <- function(eta) stop("'mu.eta' should not have been called")
  fit <- nlist(stanfit, algorithm, family = f, formula, data, offset, weights,
               x = X, y = Y, model = mf,  terms = mt, call, 
               na.action = attr(mf, "na.action"), 
               contrasts = attr(X, "contrasts"), 
               modeling_function = "stan_clogit", if(has_bars) glmod)
  out <- stanreg(fit)
  out$xlevels <- .getXlevels(mt, mf)
  class(out) <- c(class(out), if(has_bars) "lmerMod", "clogit")
  return(out)
}

log_clogit_denom <- function(N_j, D_j, eta_j) {
  if (D_j == 1 && N_j == NROW(eta_j)) return(log_sum_exp(eta_j));
  if (D_j == 0) return(0)
  if (N_j == D_j) {
    if (D_j == 1) return(eta_j[N_j])
    return(sum(eta_j[(N_j - 1):(N_j + 1)]))
  }
  else {
    N_jm1 <- N_j - 1
    return( log_sum_exp2(log_clogit_denom(N_jm1, D_j, eta_j),
            log_clogit_denom(N_jm1, D_j - 1, eta_j) + eta_j[N_j]) )
  }
}
