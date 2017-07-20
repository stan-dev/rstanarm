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

#' Multinomial logit regression models via Stan
#'
#' Multinomial logit modeling with optional prior distributions for the 
#' coefficients and intercepts.
#'
#' @export
#' @templateVar armRef (Ch. 3-6)
#' @templateVar pkg mlogit
#' @templateVar pkgfun mlogit
#' @templateVar sameargs model,offset,weights 
#' @templateVar rareargs na.action
#' @templateVar fun stan_mlogit
#' @templateVar fitfun stan_glm.fit
#' @template return-stanreg-object
#' @template return-stanfit-object
#' @template see-also
#' @template args-formula-data-subset
#' @template args-same-as
#' @template args-same-as-rarely
#' @template args-x-y
#' @template args-dots
#' @template args-priors
#' @template args-prior_PD
#' @template args-algorithm
#' @template args-adapt_delta
#' @template args-QR
#' 
#' 
#' @details Add details
#'   
#' @seealso The vignette for \code{stan_mlogit}.
#' #' 
#' @examples 
#' ### Add 5 second example
stan_mlogit <-
  function(formula, # figure out which arguments of mlogit::mlogit we want to support
           data,
           subset,
           weights,
           na.action,
           alt.subset = NULL,
           reflevel = NULL,
           ...,
           prior = normal(),
           prior_intercept = normal(),
           prior_PD = FALSE,
           algorithm = c("sampling", "optimizing", "meanfield", "fullrank"),
           adapt_delta = NULL,
           QR = FALSE,
           sparse = FALSE) {
    
    if (!requireNamespace("mlogit", quietly = TRUE))
      stop("Please install the mlogit package before using 'stan_mlogit'.")
    
    data <- validate_data(data, if_missing = environment(formula))
    mc <- match.call(expand.dots = FALSE)

    # NULLify any Stan specific arguments in mc
    mc$prior <- mc$prior_intercept <- mc$prior_PD <- mc$algorithm <-
      mc$adapt_delta <- mc$QR <- mc$sparse <- NULL
    
    mc[[1L]] <- quote(mlogit::mlogit)
    mc$estimate <- FALSE
    mf <- eval(mc, parent.frame())
    mt <- terms(mf)
    mf <- check_constant_vars(mf)
    Y <- array1D_check(model.response(mf, type = "any"))
    X <- model.matrix(mlogit::mFormula(formula), data = mf)
    weights <- validate_weights(as.vector(model.weights(mf)))
    offset <- double() # no offset allowed in mlogit::mlogit
    
    algorithm <- match.arg(algorithm)

    stop("modify stan_glm.fit.R to handle a multinomial logit model")
    # you can probably use exec/binomial.stan for this purpose with some additions
    stanfit <- 
      stan_glm.fit(x = X, y = Y, 
                   weights = weights, offset = offset,
                   ...,
                   prior = prior,
                   prior_intercept = prior_intercept, 
                   algorithm = algorithm, adapt_delta = adapt_delta, 
                   QR = QR, sparse = sparse)
    fit <- 
      nlist(stanfit, algorithm, data, offset, weights,
            x = X, y = Y,
            formula, model = mf, terms = mt, call = match.call(),
            na.action = attr(mf, "na.action"), contrasts = attr(X, "contrasts"), 
            stan_function = "stan_mlogit")
    out <- stanreg(fit)
    out$xlevels <- lapply(mf[,-1], FUN = function(x) {
      xlev <- if (is.factor(x) || is.character(x)) levels(x) else NULL
      xlev[!vapply(xlev, is.null, NA)]
    })
    structure(out, class = c("stanreg", "mlogit"))
  }
