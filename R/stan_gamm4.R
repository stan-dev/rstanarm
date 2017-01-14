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
#' @template args-prior_aux
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
#'   estimation. The Bayesian model adds priors on the common regression 
#'   coefficients (in the same way as \code{\link{stan_glm}}) and priors on the 
#'   terms of a decomposition of the covariance matrices of the group-specific 
#'   parameters, including the smooths. Estimating these models via MCMC avoids
#'   the optimization issues that often crop up with GAMMs and provides better
#'   estimates for the uncertainty in the parameter estimates. 
#'   
#'   See \code{\link[gamm4]{gamm4}} for more information about the model
#'   specicification and \code{\link{priors}} for more information about the
#'   priors. If \code{random = NULL}, the output is a subset of that produced by
#'   \code{\link[mgcv]{gam}} in the sense that there are several estimated components
#'   for each smooth term. However, the parameterization used to estimate the model
#'   is different and corresponds to the parameterization in 
#'   \code{\link[gamm4]{gamm4}} where is smooth term is decomposed into a linear
#'   and a non-linear part. If \code{prior} is not \code{NULL}, then the number 
#'   of parameters to place priors on is equal to the number of linear terms in
#'   the \code{formula}. The prior on the non-linear part of each smooth term is
#'   handled by the \code{\link{decov}} function. If \code{random} is not \code{NULL},
#'   then there are additional group-specific terms whose priors are also handled
#'   by the \code{\link{decov}} function and whose posterior medians can be extracted
#'   by calling \code{\link[lme4]{ranef}}.
#'   
#'   The \code{plot_nonlinear} function creates a ggplot object with one facet for
#'   each smooth function specified in the call to \code{stan_gamm4} in the case
#'   where all smooths are univariate. A subset of the smooth functions can be 
#'   specified using the \code{smooths} argument, which is necessary to plot a
#'   bivariate smooth or to exclude the bivariate smooth and plot the univariate
#'   ones. In the bivariate case, a plot is produced using 
#'   \code{\link[ggplot2]{geom_contour}}. In the univariate case, the resulting
#'   plot is conceptually similar to \code{\link[mgcv]{plot.gam}} except the 
#'   outer lines here demark the edges of posterior uncertainty intervals 
#'   (credible intervals) rather than confidence intervals and the inner line
#'   is the posterior median of the function rather than the function implied
#'   by a point estimate. To change the colors used in the plot see 
#'   \code{\link[bayesplot]{color_scheme_set}}.
#'   
#' @references 
#' Crainiceanu, C., Ruppert D., and Wand, M. (2005). Bayesian analysis for 
#' penalized spline regression using WinBUGS. \emph{Journal of Statistical
#' Software}. \strong{14}(14), 1--22. 
#' \url{https://www.jstatsoft.org/article/view/v014i14}
#' 
#' @examples
#' # from example(gamm4, package = "gamm4"), prefixing gamm4() call with stan_
#' \donttest{
#' dat <- mgcv::gamSim(1, n = 400, scale = 2) ## simulate 4 term additive truth
#' ## Now add 20 level random effect `fac'...
#' dat$fac <- fac <- as.factor(sample(1:20, 400, replace = TRUE))
#' dat$y <- dat$y + model.matrix(~ fac - 1) %*% rnorm(20) * .5
#'
#' br <- stan_gamm4(y ~ s(x0) + x1 + s(x2), data = dat, random = ~ (1 | fac), 
#'                  QR = TRUE, chains = 1, iter = 200) # for example speed
#' print(br)
#' plot_nonlinear(br)
#' plot_nonlinear(br, smooths = "s(x0)", alpha = 2/3)
#' }
#' 
#' @importFrom lme4 getME
stan_gamm4 <- function(formula, random = NULL, family = gaussian(), data = list(), 
                       weights = NULL, subset = NULL, na.action, knots = NULL, 
                       drop.unused.levels = TRUE, ..., 
                       prior = normal(), prior_intercept = normal(),
                       prior_aux = cauchy(0, 5),
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
  if (is.null(prior_aux)) 
    prior_aux <- list()

  group <- glmod$reTrms
  group$decov <- prior_covariance
  algorithm <- match.arg(algorithm)
  
  stanfit <- stan_glm.fit(x = X, y = y, weights = weights,
                          offset = offset, family = family,
                          prior = prior, prior_intercept = prior_intercept,
                          prior_aux = prior_aux, prior_PD = prior_PD, 
                          algorithm = algorithm, adapt_delta = adapt_delta,
                          group = group, QR = QR, ...)

  if (nrow(glmod$raw_X) < ncol(glmod$raw_X)) {
    warning("cannot reparameterize output, returning a stanfit object with reduced functionaliy")
    return(stanfit)
  }
  # Convert back to gam parameterization
  genuine <- lme4::findbars(random)
  nrows <- sapply(group$Ztlist, FUN = nrow)
  Z <- pad_reTrms(Ztlist = group$Ztlist, cnms = group$cnms, 
                  flist = group$flist)$Z
  colnames(Z) <- b_names(names(stanfit), value = TRUE)
  if (length(genuine)) {
    mark <- grep(" Xr", colnames(Z), fixed = TRUE)
    Z_nl <- Z[,mark, drop = FALSE]
  }
  else Z_nl <- Z
  XZ <- if (getRversion() < "3.2.0") cBind(X, Z_nl) else cbind2(X, Z_nl)
  arr <- as.array(stanfit)
  XtX <- crossprod(glmod$raw_X)
  for (chain in 1:ncol(stanfit)) {
    end <- tail(grep("b[", names(stanfit@sim$samples[[chain]]), 
                     fixed = TRUE, value = FALSE), 1)
    if (length(genuine)) {
      for (start in rev(grep(" Xr", colnames(Z), fixed = TRUE, invert = TRUE))) {
        stanfit@sim$samples[[chain]][end] <- stanfit@sim$samples[[chain]][start]
        names(stanfit@sim$samples[[chain]])[end] <-
        names(stanfit@sim$samples[[chain]])[start]  
        end <- end - 1L
      }
    }
    eta <- linear_predictor(arr[, chain, colnames(XZ)], x = XZ)
    beta <- qr.solve(XtX, t(glmod$raw_X) %*% t(eta))
    for (j in 1:nrow(beta)) {
      stanfit@sim$samples[[chain]][[j]] <- beta[j,]
      names(stanfit@sim$samples[[chain]])[j] <- rownames(beta)[j]
    }
  }
  X <- glmod$raw_X
  stanfit@par_dims$beta <- ncol(X) - stanfit@par_dims$alpha
  if (length(genuine)) Z <- Z[,-mark,drop = FALSE]
  else Z <- Matrix::Matrix(nrow = nrow(Z), ncol = 0, sparse = TRUE)
  stanfit@par_dims$b <- ncol(Z)
  XZ <- if (getRversion() < "3.2.0") cBind(X, Z) else cbind2(X, Z)
  stanfit@sim$fnames_oi[1:ncol(XZ)] <- colnames(XZ)
  fit <- nlist(stanfit, family, formula, offset, weights, 
               x = XZ, y = y, data, terms = NULL, model = NULL, 
               call = match.call(expand.dots = TRUE),
               algorithm, glmod = glmod)
  out <- stanreg(fit)
  class(out) <- c(class(out), "gamm4", "lmerMod")
  return(out)
}

#' @rdname stan_gamm4
#' @export
#' @param x An object produced by \code{stan_gamm4}.
#' @param smooths An optional character vector specifying a subset of the smooth
#'   functions specified in the call to \code{stan_gamm4}. The default is
#'   include all smooth terms.
#' @param prob For univarite smooths, a scalar between 0 and 1 governing the
#'   width of the uncertainty interval.
#' @param facet_args An optional named list of arguments passed to 
#'   \code{\link[ggplot2]{facet_wrap}} (other than the \code{facets} argument).
#' @param alpha,size For univariate smooths, passed to 
#'   \code{\link[ggplot2]{geom_ribbon}}. For bivariate smooths, \code{size/2} is
#'   passed to \code{\link[ggplot2]{geom_contour}}.
#'   
#' @return \code{plot_nonlinear} returns a ggplot object.
#' 
#' @importFrom ggplot2 aes_ aes_string facet_wrap ggplot geom_contour geom_line geom_ribbon labs scale_color_gradient2
#' 
plot_nonlinear <- function(x, smooths, ..., 
                           prob = 0.9, facet_args = list(), 
                           alpha = 1, size = 0.75) {
  validate_stanreg_object(x)
  if (!is(x, "gamm4"))
    stop("Plot only available for models fit using the stan_gamm4 function.")
  
  scheme <- bayesplot::color_scheme_get()
  
  XZ <- x$x
  XZ <- XZ[,!grepl("_NEW_", colnames(XZ), fixed = TRUE)]
  labels <- sapply(x$glmod$smooths, "[[", "label")
  xnames <- sapply(x$glmod$smooths, "[[", "vn")
  names(xnames) <- labels
  if (!missing(smooths)) {
    found <- smooths %in% labels
    if (all(!found)) {
      stop("All specified terms are invalid. Valid terms are: ", 
              paste(grep(",", labels, fixed = TRUE, value = TRUE, invert = TRUE), 
                    collapse = ", "))
    } else if (any(!found)) {
      warning("The following specified terms were not found and ignored: ", 
              paste(smooths[!found], collapse = ", "))
    }
    labels <- smooths[found]
    if (!is.matrix(xnames)) xnames <- xnames[found]
  }
  
  B <- as.matrix(x)[, colnames(XZ), drop = FALSE]
  original <- x$glmod$model
  
  bivariate <- any(grepl(",", labels, fixed = TRUE))
  if (bivariate) {
    if (length(labels) > 1)
      stop("Multivariate functions can only be plotted one at a time; specify 'smooths'.")
    if (length(xnames) > 2)
      stop("Only univariate and bivariate functions can be plotted currently.")
    xrange <- range(original[, xnames[1]])
    yrange <- range(original[, xnames[2]])
    xz <- expand.grid(seq(from = xrange[1], to = xrange[2], length.out = 100),
                      seq(from = yrange[1], to = yrange[2], length.out = 100))
    colnames(xz) <- xnames[1:2]
    plot_data <- data.frame(x = xz[, 1], y = xz[, 2])
    for (i in colnames(original)) {
      if (i %in% colnames(xz)) next
      xz[[i]] <- 1
    }
    XZ <- mgcv::predict.gam(mgcv::gam(formula(x), data = x$data), 
                            newdata = xz, type = "lpmatrix")
    incl <- grepl(labels, colnames(B), fixed = TRUE)
    b <- B[, incl, drop = FALSE]
    xz <- XZ[, grepl(labels, colnames(XZ), fixed = TRUE), drop = FALSE]
    plot_data$z <- apply(linear_predictor.matrix(b, xz), 2, FUN = median)
    return(
      ggplot(plot_data, aes_(x = ~x, y = ~y, z = ~z)) + 
             geom_contour(aes_string(color = "..level.."), size = size/2) + 
             labs(x = xnames[1], y = xnames[2]) + 
             scale_color_gradient2(low = scheme[[1]],
                                   mid = scheme[[3]], 
                                   high = scheme[[6]]) +
             bayesplot::theme_default()
    )
  }
  
  df_list <- lapply(labels, FUN = function(term) {
    incl <- grepl(term, colnames(B), fixed = TRUE)
    xz <- XZ[, incl, drop = FALSE]
    b <- B[, incl, drop = FALSE]
    x <- original[, xnames[term]]
    xz <- xz[order(x), , drop=FALSE]
    f <- linear_predictor.matrix(b, xz)
    data.frame(
      predictor = sort(x),
      lower  = apply(f, 2, quantile, probs = (1 - prob) / 2),
      upper  = apply(f, 2, quantile, probs = prob + (1 - prob) / 2),
      middle = apply(f, 2, median),
      term = term
    )
  })
  plot_data <- do.call(rbind, df_list)
  
  facet_args[["facets"]] <- ~ term
  if (is.null(facet_args[["scales"]]))
    facet_args[["scales"]] <- "free"
  if (is.null(facet_args[["strip.position"]]))
    facet_args[["strip.position"]] <- "left"
  
  ggplot(plot_data, aes_(x = ~ predictor)) + 
    geom_ribbon(aes_(ymin = ~ lower, ymax = ~ upper), 
                fill = scheme[[1]], color = scheme[[2]],
                alpha = alpha, size = size) + 
    geom_line(aes_(y = ~ middle), color = scheme[[5]], 
              size = 0.75 * size, lineend = "round") + 
    labs(y = NULL) + 
    do.call(facet_wrap, facet_args) + 
    bayesplot::theme_default()
}
