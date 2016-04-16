# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2015, 2016 Trustees of Columbia University
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

#' Bayesian generalized linear additive models with group-specific terms via
#' Stan
#' 
#' Bayesian inference for GAMMs with flexible priors.
#' 
#' @export
#' @templateVar fun stan_gamm4
#' @templateVar pkg gamm4
#' @templateVar pkgfun gamm4
#' @template return-stanreg-object
#' @template see-also
#' @template args-priors
#' @template args-prior_PD
#' @template args-algorithm
#' @template args-adapt_delta
#' @template args-QR
#' 
#' @param formula,random,family,data,knots,drop.unused.levels Same as for 
#'   \code{\link[gamm4]{gamm4}}.
#' @param subset,weights,na.action Same as \code{\link[stats]{glm}}, 
#'   but rarely specified.
#' @param ... Further arguments passed to \code{\link[rstan]{sampling}} (e.g. 
#'   \code{iter}, \code{chains}, \code{cores}, etc.) or to
#'   \code{\link[rstan]{vb}} (if \code{algorithm} is \code{"meanfield"} or
#'   \code{"fullrank"}).
#' @param prior_covariance Cannot be \code{NULL}; see \code{\link{decov}} for
#'   more information about the default arguments.
#'
#' @details The \code{stan_gamm4} function is similar in syntax to 
#'   \code{\link[gamm4]{gamm4}}, which accepts a syntax that is similar to (but 
#'   not quite as extensive as) that for \code{\link[mgcv]{gamm}} and converts 
#'   it internally into the syntax accepted by \code{\link[lme4]{glmer}}. But 
#'   rather than performing (restricted) maximum likelihood estimation, the 
#'   \code{stan_gamm4} function utilizes MCMC to perform Bayesian estimation. 
#'   The Bayesian model adds independent priors on the common regression 
#'   coefficients (in the same way as \code{\link{stan_glm}}) and priors on the 
#'   terms of a decomposition of the covariance matrices of the group-specific 
#'   parameters, including the smooths. Estimating these models via MCMC avoids
#'   the optimization issues that often crop up with GAMMs and provides better
#'   estimates for the uncertainty in the parameter estimates. 
#'   
#'   See \code{\link[gamm4]{gamm4}} for more information about the model
#'   specicification and \code{\link{priors}} for more information about the
#'   priors.
#' @examples
#' # see example(gamm4, package = "gamm4") but prefix gamm4() calls with stan_
#' 
stan_gamm4 <- function(formula, random = NULL, family = gaussian(), data = list(), 
                       weights = NULL, subset = NULL, na.action, knots = NULL, 
                       drop.unused.levels = TRUE, ..., 
                       prior = normal(), prior_intercept = normal(),
                       prior_ops = prior_options(),
                       prior_covariance = decov(), prior_PD = FALSE, 
                       algorithm = c("sampling", "meanfield", "fullrank"), 
                       adapt_delta = NULL, QR = FALSE) {

  mc <- match.call()
  glmod <- gamm4_to_glmer(formula, random, family, data, weights, subset, 
                          na.action, knots, drop.unused.levels)

  X <- glmod$X
  y <- glmod$fr[, as.character(glmod$formula[2L])]
  if (is.matrix(y) && ncol(y) == 1L) y <- as.vector(y)
  offset <- model.offset(glmod$fr) %ORifNULL% double(0)
  weights <- validate_weights(glmod$fr$weights)
  
  if (is.null(prior))
    prior <- list()
  if (is.null(prior_intercept)) 
    prior_intercept <- list()
  if (!length(prior_ops)) 
    prior_ops <- list(scaled = FALSE, prior_scale_for_dispersion = Inf)
  

  group <- glmod$reTrms
  group$decov <- prior_covariance
  algorithm <- match.arg(algorithm)
  
  stanfit <- stan_glm.fit(x = X, y = y, weights = weights,
                          offset = offset, family = family,
                          prior = prior, prior_intercept = prior_intercept,
                          prior_ops = prior_ops, prior_PD = prior_PD, 
                          algorithm = algorithm, adapt_delta = adapt_delta,
                          group = group, QR = QR, ...)
  
  Z <- pad_reTrms(Z = t(group$Zt), cnms = group$cnms, 
                  flist = group$flist)$Z
  colnames(Z) <- b_names(names(stanfit), value = TRUE)
  fit <- nlist(stanfit, family, formula, offset, weights, 
               x = if (getRversion() < "3.2.0") cBind(X, Z) else cbind2(X, Z), 
               y = y, data, call = mc, terms = NULL, model = NULL, 
               prior.info = get_prior_info(call, formals()),
               algorithm, glmod = glmod)
  out <- stanreg(fit)
  class(out) <- c(class(out), "lmerMod")
  return(out)
  # TODO: maybe convert back to gam parameterization?
}

#' @importFrom lme4 findbars mkReTrms
#' @importFrom mgcv interpret.gam gamm new.name
gamm4_to_glmer <- function(formula, random = NULL, family = gaussian(), data = list(), 
                           weights = NULL, subset = NULL, na.action, knots = NULL, 
                           drop.unused.levels = TRUE) {

  if (!is.null(random)) {
    if (!inherits(random,"formula")) stop("`random' must be a lme4-style formula")
    random.vars <- all.vars(random)
    random_list <- findbars(random)
    names(random_list) <- names(mkReTrms(random_list, data, drop.unused.levels)$cnms)
    for (i in seq_along(random_list)) { 
      random_list[[i]] <- as.formula(paste("~", random_list[i]))
    }
  } 
  else random_list <- random.vars <- NULL

  # create model frame.....
  gp <- mgcv::interpret.gam(formula) # interpret the formula 
  # call mgcv::gamm() with zero iterations to set everything (?) up
  # calling gamm4::gamm4() is too slow because we cannot do zero iterations
  settings <- list(maxIter = 0, msMaxIter = 0, niterEM = 0, msMaxEval = 0, 
                   gradHess = FALSE, apVar = FALSE, opt = "optim", 
                   optimMethod = "BFGS", natural = FALSE)
  call <- match.call(expand.dots = TRUE)
  mc <- match.call(expand.dots = FALSE)
  mc[[1]] <- quote(mgcv::gamm)
  mc$control <- settings
  mc$random <- random_list
  if (is.null(weights)) mc$weights <- NULL
  if (is.null(subset))  mc$subset  <- NULL
  mc$niterPQL <- 0L
  mc$verbosePQL <- FALSE
  GAM <- eval(mc, parent.frame())
  G <- GAM$gam
  mf <- G$model
  X <- mf$X
  r.name <- grep("^Xr", colnames(mf), value = TRUE)
  G$random <- sapply(r.name, simplify = FALSE, FUN = function(r) mf[[r]])
  
  Terms <- attr(mf,"terms")
  
  ## summarize the *raw* input variables
  ## note can't use get_all_vars here -- buggy with matrices
  vars <- all.vars(gp$fake.formula[-2]) ## drop response here
  inp <- parse(text = paste("list(", paste(vars, collapse = ","),")"))
  dl <- eval(inp, data, parent.frame())
  names(dl) <- vars ## list of all variables needed

  ## lmer offset handling work around...
  ## variables not in mf raw -- can cause lmer problem
  mvars <- vars[!(vars %in% names(mf))]
  if (length(mvars) > 0) ## append raw versions to mf
    for (i in 1:length(mvars)) mf[[mvars[i]]] <- dl[[mvars[i]]]
  
  rm(dl) ## save space 
  
  n.sr <- length(G$random) # number of random smooths (i.e. s(...,fx=FALSE,...) terms)
  
  if (is.null(random) && n.sr == 0) 
    stop("gamm4-style models must have at least 1 smooth with unknown smoothing parameter or at least one other random effect")
  
  offset.name <- attr(mf,"names")[attr(attr(mf,"terms"),"offset")]
  
  yname <- new.name("y",names(mf))
  # eval(parse(text=paste("mf$",yname,"<-mf$y",sep="")))
  mf[[yname]] <- mf$y
  Xname <- new.name("X",names(mf))
  # eval(parse(text=paste("mf$",Xname,"<-mf$X",sep="")))
  mf[[Xname]] <- mf$X
  lme4.formula <- paste(yname,"~", Xname, "-1")
  if (length(offset.name))
    lme4.formula <- paste(lme4.formula,"+",offset.name)
  
  ## Basic trick is to call (g)lFormula to set up model, with simple i.i.d. dummy random effects for the 
  ## penalized component of each smooth. This results in columns of Z being produced for these dummy's,
  ## which can be over-written with the right thing. NOTE: that lambdat could also be modified, I think!!
  
  ## Add the random effect dummy variables for the smooth
  if (n.sr) for (i in 1:n.sr) { # adding the constructed variables to the model frame avoiding name duplication
    mf[[r.name[i]]] <- factor(rep(1:ncol(G$random[[i]]), length=nrow(G$random[[i]])))
    lme4.formula <- paste(lme4.formula,"+ (1|",r.name[i],")")
  }
  
  if (!is.null(random)) ## append the regular random effects
    lme4.formula <- paste(lme4.formula, "+" , substring(deparse(random), first = 2))
  
  lme4.formula <- as.formula(lme4.formula)
  
  b <- glFormula(lme4.formula, data = mf, family = family, weights = G$w,
                 control = glmerControl(check.nlev.gtreq.5 = "ignore",
                                        check.nlev.gtr.1 = "stop",
                                        check.nobs.vs.rankZ = "ignore",
                                        check.nobs.vs.nlev = "ignore",
                                        check.nobs.vs.nRE = "ignore"))
  
  if (n.sr) { ## Fabian Scheipl's trick of overwriting dummy slots revised for new structure
    tn <- names(b$reTrms$cnms) ## names associated with columns of Z (same order as Gp)
    ind <- 1:length(tn)
    sn <- names(G$random) ## names of smooth random components
    s_labels <- sapply(G$smooth, FUN = function(s) s$label)
    for (i in 1:n.sr) { ## loop through random effect smooths
      k <- ind[sn[i] == tn] ## which term should contain G$random[[i]] 
      ii <- (b$reTrms$Gp[k]+1):b$reTrms$Gp[k+1]
      b$reTrms$Zt[ii,] <- as(t(G$random[[i]]),"dgCMatrix")
      b$reTrms$cnms[[k]] <- s_labels[[i]]
    }
  }
  return(b)
}
