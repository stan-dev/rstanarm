# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2015, 2016 Trustees of Columbia University
#  Copyright (C) 1995-2015 The R Core Team
#  Copyright (C) 1998 B. D. Ripley
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

#' @rdname stan_lm
#' @export
#' @param projections For \code{stan_aov}, a logical scalar (defaulting to
#'   \code{FALSE}) indicating whether \code{\link[stats]{proj}} should be called
#'   on the fit.
#' @examples
#' \dontrun{
#' stan_aov(yield ~ block + N*P*K, data = npk, contrasts = "contr.poly",
#'          prior = R2(0.5), seed = 12345) 
#' }
#'             
stan_aov <- function(formula, data = NULL, projections = FALSE,
                     contrasts = NULL, ...,
                     prior = R2(stop("'location' must be specified")), 
                     prior_PD = FALSE, 
                     algorithm = c("sampling", "meanfield", "fullrank"), 
                     adapt_delta = NULL) {
    # parse like aov() does
    Terms <- if(missing(data)) 
      terms(formula, "Error") else terms(formula, "Error", data = data)
    indError <- attr(Terms, "specials")$Error
    ## NB: this is only used for n > 1, so singular form makes no sense
    ## in English.  But some languages have multiple plurals.
    if(length(indError) > 1L)
        stop(sprintf(ngettext(length(indError),
                              "there are %d Error terms: only 1 is allowed",
                              "there are %d Error terms: only 1 is allowed"),
                     length(indError)), domain = NA)
    lmcall <- Call <- match.call()
    ## need rstanarm:: for non-standard evaluation
    lmcall[[1L]] <- quote(rstanarm::stan_lm)
    lmcall$singular.ok <- FALSE
    if (projections) 
      qr <- lmcall$qr <- TRUE
    lmcall$projections <- NULL
    if (is.null(indError)) {
        ## no Error term
        fit <- eval(lmcall, parent.frame())
        fit$terms <- Terms
        fit$qr <- qr(model.matrix(Terms, data = fit$data))
        R <- qr.R(fit$qr)
        beta <- extract(fit$stanfit, pars = "beta", permuted = FALSE)
        pnames <- dimnames(beta)$parameters
        rownames(R) <- colnames(R)
        R <- R[pnames, pnames, drop = FALSE]
        effects <- apply(beta, 1:2, FUN = function(x) R %*% x)
        effects <- aperm(effects, c(2,3,1))
        fit$effects <- effects
        class(fit) <- c("stanreg", "aov", "lm")
        if (projections) 
          fit$projections <- proj(fit)
        fit$call <- Call
        return(fit)
    } else { # nocov start
        stop("Error terms not supported yet")
        if(pmatch("weights", names(match.call()), 0L))
            stop("weights are not supported in a multistratum aov() fit")
        ##  Helmert contrasts can be helpful: do we want to force them?
        ##  this version does for the Error model.
        opcons <- options("contrasts")
        options(contrasts = c("contr.helmert", "contr.poly"))
        on.exit(options(opcons))
        allTerms <- Terms
        errorterm <-  attr(Terms, "variables")[[1 + indError]]
        eTerm <- deparse(errorterm[[2L]], width.cutoff = 500L, backtick = TRUE)
        intercept <- attr(Terms, "intercept")
        ecall <- lmcall
        ecall$formula <-
            as.formula(paste(deparse(formula[[2L]], width.cutoff = 500L,
                                     backtick = TRUE), "~", eTerm,
                             if(!intercept) "- 1"),
                       env = environment(formula))

        ecall$method <- "qr"
        ecall$qr <- TRUE
        ecall$contrasts <- NULL
        er.fit <- eval(ecall, parent.frame())
        options(opcons)
        nmstrata <- attr(terms(er.fit), "term.labels")
        ## remove backticks from simple labels for strata (only)
        nmstrata <- sub("^`(.*)`$", "\\1", nmstrata)
        nmstrata <- c("(Intercept)", nmstrata)
        qr.e <- er.fit$qr
        rank.e <- er.fit$rank
        if(rank.e < NROW(er.fit$coefficients))
            warning("Error() model is singular")
        qty <- er.fit$residuals
        maov <- is.matrix(qty)
        asgn.e <- er.fit$assign[qr.e$pivot[1L:rank.e]]
        ## we want this to label the rows of qtx, not cols of x.
        maxasgn <- length(nmstrata) - 1L
        nobs <- NROW(qty)
	len <- if(nobs > rank.e) {
	    asgn.e[(rank.e+1):nobs] <- maxasgn + 1L
	    nmstrata <- c(nmstrata, "Within")
	    maxasgn + 2L
	} else maxasgn + 1L
        result <- setNames(vector("list", len), nmstrata)
        lmcall$formula <- form <-
            update(formula, paste(". ~ .-", deparse(errorterm, width.cutoff = 500L, backtick = TRUE)))
        Terms <- terms(form)
        lmcall$method <- "model.frame"
        mf <- eval(lmcall, parent.frame())
        xlev <- .getXlevels(Terms, mf)
        resp <- model.response(mf)
        qtx <- model.matrix(Terms, mf, contrasts)
        cons <- attr(qtx, "contrasts")
        dnx <- colnames(qtx)
        asgn.t <- attr(qtx, "assign")
        if(length(wts <- model.weights(mf))) {
            wts <- sqrt(wts)
            resp <- resp * wts
            qtx <- qtx * wts
        }
        qty <- as.matrix(qr.qty(qr.e, resp))
        if((nc <- ncol(qty)) > 1L) {
            dny <- colnames(resp)
            if(is.null(dny)) dny <- paste0("Y", 1L:nc)
            dimnames(qty) <- list(seq(nrow(qty)), dny)
        } else dimnames(qty) <- list(seq(nrow(qty)), NULL)
        qtx <- qr.qty(qr.e, qtx)
        dimnames(qtx) <- list(seq(nrow(qtx)) , dnx)
        for(i in seq_along(nmstrata)) {
            select <- asgn.e == (i-1L)
            ni <- sum(select)
            if(!ni) next
            ## helpful to drop constant columns.
            xi <- qtx[select, , drop = FALSE]
            cols <- colSums(xi^2) > 1e-5
            if(any(cols)) {
                xi <- xi[, cols, drop = FALSE]
                attr(xi, "assign") <- asgn.t[cols]
                fiti <- lm.fit(xi, qty[select,,drop=FALSE])
                fiti$terms <- Terms
            } else {
                y <- qty[select,,drop=FALSE]
                fiti <- list(coefficients = numeric(), residuals = y,
                             fitted.values = 0 * y, weights = wts, rank = 0L,
                             df.residual = NROW(y))
            }
            if(projections) fiti$projections <- proj(fiti)
            class(fiti) <- c(if(maov) "maov", "aov", oldClass(er.fit))
            result[[i]] <- fiti
        }
        ## drop empty strata
        result <- result[!sapply(result, is.null)]
        class(result) <- c("aovlist", "listof")
        if(qr) attr(result, "error.qr") <- qr.e
        attr(result, "call") <- Call
        if(length(wts)) attr(result, "weights") <- wts
        attr(result, "terms") <- allTerms
        attr(result, "contrasts") <- cons
        attr(result, "xlevels") <- xlev
        result
    } # nocov end
}
