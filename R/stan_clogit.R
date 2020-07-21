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

#' Conditional logistic (clogit) regression models via Stan
#'
#' \if{html}{\figure{stanlogo.png}{options: width="25px" alt="http://mc-stan.org/about/logo/"}}
#' A model for case-control studies with optional prior distributions for the
#' coefficients, intercept, and auxiliary parameters.
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
#' @template args-dots
#' 
#' @param formula,data,subset,na.action Same as for \code{\link[lme4]{glmer}},
#'   except that any global intercept included in the formula will be dropped.
#'   \emph{We strongly advise against omitting the \code{data} argument}. Unless
#'   \code{data} is specified (and is a data frame) many post-estimation
#'   functions (including \code{update}, \code{loo}, \code{kfold}) are not
#'   guaranteed to work properly.
#' @param strata A factor indicating the groups in the data where the number of 
#'   successes (possibly one) is fixed by the research design. It may be useful 
#'   to use \code{\link{interaction}} or \code{\link[survival]{strata}} to
#'   create this factor. However, the \code{strata} argument must not rely on
#'   any object besides the \code{data} \code{\link{data.frame}}.
#' @param prior_covariance Cannot be \code{NULL} when lme4-style group-specific
#'   terms are included in the \code{formula}. See \code{\link{decov}} for
#'   more information about the default arguments. Ignored when there are no
#'   group-specific terms.
#' 
#' @details The \code{stan_clogit} function is mostly similar in syntax to 
#'   \code{\link[survival]{clogit}} but rather than performing maximum
#'   likelihood estimation of generalized linear models, full Bayesian
#'   estimation is performed (if \code{algorithm} is \code{"sampling"}) via
#'   MCMC. The Bayesian model adds priors (independent by default) on the
#'   coefficients of the GLM.
#'   
#'   The \code{data.frame} passed to the \code{data} argument must be sorted by 
#'   the variable passed to the \code{strata} argument.
#'   
#'   The \code{formula} may have group-specific terms like in
#'   \code{\link{stan_glmer}} but should not allow the intercept to vary by the
#'   stratifying variable, since there is no information in the data with which
#'   to estimate such deviations in the intercept.
#'   
#' @seealso The vignette for Bernoulli and binomial models, which has more
#'   details on using \code{stan_clogit}.
#'   \url{http://mc-stan.org/rstanarm/articles/}
#' 
#' @examples
#' dat <- infert[order(infert$stratum), ] # order by strata
#' post <- stan_clogit(case ~ spontaneous + induced + (1 | education), 
#'                     strata = stratum,
#'                     data = dat,
#'                     subset = parity <= 2,
#'                     QR = TRUE,
#'                     chains = 2, iter = 500) # for speed only
#'
#' nd <- dat[dat$parity > 2, c("case", "spontaneous", "induced", "education", "stratum")]
#' # next line would fail without case and stratum variables                                 
#' pr <- posterior_epred(post, newdata = nd) # get predicted probabilities
#' 
#' # not a random variable b/c probabilities add to 1 within strata
#' all.equal(rep(sum(nd$case), nrow(pr)), rowSums(pr)) 
#'             
#' @importFrom lme4 findbars
stan_clogit <- function(formula, data, subset, na.action = NULL, ..., 
                        strata, prior = normal(autoscale=TRUE), 
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
  mf$data <- data
  err <- try(eval(mf$weights, data, enclos = NULL), silent = TRUE)
  if (inherits(err, "try-error")) {
    stop("the 'stratum' argument must be evaluatable solely within 'data'")
  }
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
  
  ord <- order(group$strata)
  if (any(diff(ord) <= 0)) {
    stop("Data must be sorted by 'strata' (in increasing order).")
  }
  offset <- model.offset(mf) %ORifNULL% double(0)
  weights <- double(0)
  mf <- check_constant_vars(mf)
  mt <- attr(mf, "terms")
  if (is.empty.model(mt))
    stop("Predictors specified.", call. = FALSE)
  xint <- match("(Intercept)", colnames(X), nomatch = 0L)
  if (xint > 0L) {
    X <- X[, -xint, drop = FALSE]
    # I cannot remember why I was calling drop.terms() to get rid of the intercept
    # mt <- drop.terms(mt, dropx = xint)
    attr(mt, "intercept") <- 0L
  }
  f <- binomial(link = "logit")
  stanfit <- stan_glm.fit(x = X, y = Y, weights = weights, 
                          offset = offset, family = f,
                          prior = prior, 
                          prior_PD = prior_PD, 
                          algorithm = algorithm, adapt_delta = adapt_delta, 
                          group = group, QR = QR, sparse = sparse, ...)
  if (algorithm != "optimizing" && !is(stanfit, "stanfit")) return(stanfit)  
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
               stan_function = "stan_clogit", 
               glmod = if(has_bars) glmod)
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
