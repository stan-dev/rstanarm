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
#' @template args-sparse
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
#'   \code{\link[gamm4]{gamm4}} in the \pkg{gamm4} package, which accepts a syntax 
#'   that is similar to (but not quite as extensive as) that for 
#'   \code{\link[mgcv]{gamm}} in the \pkg{mgcv} package and converts 
#'   it internally into the syntax accepted by \code{\link[lme4]{glmer}} in the
#'   \pkg{lme4} package. But rather than performing (restricted) maximum likelihood 
#'   estimation, the \code{stan_gamm4} function utilizes MCMC to perform Bayesian 
#'   estimation. The Bayesian model adds independent priors on the common regression 
#'   coefficients (in the same way as \code{\link{stan_glm}}) and priors on the 
#'   terms of a decomposition of the covariance matrices of the group-specific 
#'   parameters, including the smooths. Estimating these models via MCMC avoids
#'   the optimization issues that often crop up with GAMMs and provides better
#'   estimates for the uncertainty in the parameter estimates. 
#'   
#'   See \code{\link[gamm4]{gamm4}} for more information about the model
#'   specicification and \code{\link{priors}} for more information about the
#'   priors.
#' @references 
#' Crainiceanu, C., Ruppert D., and Wand, M. (2005). Bayesian Analysis for 
#' Penalized Spline Regression Using WinBUGS. 
#' \emph{Journal of Statistical Software}. \strong{14}(14), 1--22.
#' \url{https://www.jstatsoft.org/article/view/v014i14}
#' @examples
#' # from example(gamm4, package = "gamm4"), prefixing gamm4() calls with stan_
#' \donttest{
#' dat <- mgcv::gamSim(1,n=400,scale=2) ## simulate 4 term additive truth
#' ## Now add 20 level random effect `fac'...
#' dat$fac <- fac <- as.factor(sample(1:20,400,replace=TRUE))
#' dat$y <- dat$y + model.matrix(~fac-1)%*%rnorm(20)*.5
#'
#' br <- stan_gamm4(y~s(x0)+x1+s(x2),data=dat,random=~(1|fac), QR = TRUE)
#' br
#'
#' ##########################
#' ## Poisson example GAMM...
#' ##########################
#' ## simulate data...
#' x <- runif(100)
#' fac <- sample(1:20,100,replace=TRUE)
#' eta <- x^2*3 + fac/20; fac <- as.factor(fac)
#' y <- rpois(100,exp(eta))
#'
#' ## fit model and examine it...
#' bp <- stan_gamm4(y~s(x), family=poisson, random=~(1|fac),
#'                  data = data.frame(x, y, fac))
#' bp
#'
#' #################################################################
#' ## Add a factor to the linear predictor, to be modelled as random
#' ## and make response Poisson.
#' #################################################################
#' g <- as.factor(sample(1:20,400,replace=TRUE))
#' dat$f <- dat$f + model.matrix(~ g-1)%*%rnorm(20)*2
#' dat$y <- rpois(400,exp(dat$f/7))
#' dat$g <- g
#'
#' b2 <- stan_gamm4(y~s(x0)+s(x1)+s(x2)+s(x3), family=poisson,
#'                  data=dat, random=~(1|g), QR = TRUE)
#' b2
#'
#'
#' ##################################
#' # Multivariate varying coefficient
#' ## Start by simulating data...
#'
#' f0 <- function(x, z, sx = 0.3, sz = 0.4) {
#'   (pi^sx * sz) * (1.2 * exp(-(x - 0.2)^2/sx^2 - 
#'   (z - 0.3)^2/sz^2) + 0.8 * exp(-(x - 0.7)^2/sx^2 - (z - 0.8)^2/sz^2))
#' }
#' f1 <- function(x2) 2 * sin(pi * x2)
#' f2 <- function(x2) exp(2 * x2) - 3.75887
#' f3 <- function (x2) 
#'   0.2 * x2^11 * (10 * (1 - x2))^6 + 10 * (10 * x2)^3 * (1 - x2)^10
#'
#' n <- 1000
#'
#' ## first set up a continuous-within-group effect...
#'
#' g <- factor(sample(1:50,n,replace=TRUE)) ## grouping factor
#' x <- runif(n)                       ## continuous covariate
#' X <- model.matrix(~g-1)
#' mu <- X%*%rnorm(50)*.5 + (x*X)%*%rnorm(50)
#'
#' ## now add nested factors...
#' a <- factor(rep(1:20,rep(50,20)))
#' b <- factor(rep(rep(1:25,rep(2,25)),rep(20,50)))
#' Xa <- model.matrix(~a-1)
#' Xb <- model.matrix(~a/b-a-1)
#' mu <- mu + Xa%*%rnorm(20) + Xb%*%rnorm(500)*.5
#'
#' ## finally simulate the smooth terms
#' v <- runif(n); w <- runif(n); z <- runif(n)
#' r <- runif(n)
#' mu <- mu + f0(v,w)*z*10 + f3(r) 
#'
#' y <- mu + rnorm(n)*2 ## response data
#'
#'
#' br <- stan_gamm4(y ~ s(v,w,by=z) + s(r,k=20,bs="cr"), random = ~ (1|a/b),
#'                  data = data.frame(y, v, w, z, r, a, b), QR = TRUE)
#' br
#'
#' ## now fit the full model
#'
#' br <- stan_gamm4(y ~ s(v,w,by=z) + s(r,k=20,bs="cr"),
#'                  random = ~ (x+0|g) + (1|g) + (1|a/b),
#'                  data = data.frame(y, v, w, z, r, a, b), QR = TRUE)
#'
#' br
#'
#' ## try a Poisson example, based on the same linear predictor...
#'
#' lp <- mu/5
#' y <- rpois(exp(lp),exp(lp)) ## simulated response
#'
#'
#' br <- stan_gamm4(y ~ s(v,w,by=z) + s(r,k=20,bs="cr"),
#'                  family=poisson, random = ~ (1|a/b),
#'                  data = data.frame(y, v, w, z, r, a, b), QR = TRUE)
#'
#' ## and now fit full version
#'
#' br <- stan_gamm4(y ~ s(v,w,by=z) + s(r,k=20,bs="cr"),
#'                  family=poisson,random = ~ (x|g) + (1|a/b),
#'                  data = data.frame(y, v, w, z, r, a, b), QR = TRUE)
#' br
#' 
#'
#' ####################################
#' # Different smooths of x2 depending 
#' # on factor `fac'...
#' ####################################
#' dat <- mgcv::gamSim(4)
#'
#' br <- stan_gamm4(y ~ fac+s(x2,by=fac)+s(x0), data=dat, QR = TRUE)
#' summary(br)
#'
#' ######################################
#' ## A "signal" regression example, in
#' ## which a univariate response depends
#' ## on functional predictors.
#' ######################################
#'
#' ## simulate data first....
#'
#' rf <- function(x=seq(0,1,length=100)) {
#'  ## generates random functions...
#'  m <- ceiling(runif(1)*5) ## number of components
#'  f <- x*0;
#'  mu <- runif(m,min(x),max(x));sig <- (runif(m)+.5)*(max(x)-min(x))/10
#'  for (i in 1:m) f <- f+ dnorm(x,mu[i],sig[i])
#'  f
#' }
#'
#' x <- seq(0,1,length=100) ## evaluation points
#'
#'
#' ## simulate 200 functions and store in rows of L...
#' L <- matrix(NA,200,100) 
#' for (i in 1:200) L[i,] <- rf()  ## simulate the functional predictors
#'
#' f2 <- function(x) { ## the coefficient function
#'  (0.2*x^11*(10*(1-x))^6+10*(10*x)^3*(1-x)^10)/10 
#' }
#'
#' f <- f2(x) ## the true coefficient function
#'
#' y <- L%*%f + rnorm(200)*20 ## simulated response data
#'
#' ## Now fit the model E(y) = L%*%f(x) where f is a smooth function.
#' ## The summation convention is used to evaluate smooth at each value
#' ## in matrix X to get matrix F, say. Then rowSum(L*F) gives E(y).
#'
#' ## create matrix of eval points for each function. Note that
#' ## `smoothCon' is smart and will recognize the duplication...
#' X <- matrix(x,200,100,byrow=TRUE) 
#'
#' br <- stan_gamm4(y~s(X,by=L,k=20), QR = TRUE)
#' br
#' }
stan_gamm4 <- function(formula, random = NULL, family = gaussian(), data = list(), 
                       weights = NULL, subset = NULL, na.action, knots = NULL, 
                       drop.unused.levels = TRUE, ..., 
                       prior = normal(), prior_intercept = normal(),
                       prior_ops = prior_options(),
                       prior_covariance = decov(), prior_PD = FALSE, 
                       algorithm = c("sampling", "meanfield", "fullrank"), 
                       adapt_delta = NULL, QR = FALSE, sparse = FALSE) {

  mc <- match.call(expand.dots = FALSE)
  family <- validate_family(family)
  glmod <- suppressWarnings(gamm4_to_glmer(formula, random, family, data, weights, 
                                           subset, na.action, knots, drop.unused.levels))

  colnames(glmod$X) <- gsub("^X\\.0", "", colnames(glmod$X))
  colnames(glmod$X) <- gsub("Fx1$",   "", colnames(glmod$X))
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
  
  Z <- pad_reTrms(Ztlist = group$Ztlist, cnms = group$cnms, 
                  flist = group$flist)$Z
  colnames(Z) <- b_names(names(stanfit), value = TRUE)
  fit <- nlist(stanfit, family, formula, offset, weights, 
               x = if (getRversion() < "3.2.0") cBind(X, Z) else cbind2(X, Z), 
               y = y, data, call = mc, terms = NULL, model = NULL, 
               prior.info = get_prior_info(call, formals()),
               algorithm, glmod = glmod)
  out <- stanreg(fit)
  class(out) <- c(class(out), "gamm4", "lmerMod")
  return(out)
  # TODO: maybe convert back to gam parameterization?
}

