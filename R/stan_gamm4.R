# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2016 Simon N. Wood
# Copyright (C) 2015, 2016, 2017 Trustees of Columbia University
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

#' Bayesian generalized linear additive models with optional group-specific
#' terms via Stan
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
#' @template args-prior_smooth
#' @template args-prior_PD
#' @template args-algorithm
#' @template args-adapt_delta
#' @template args-QR
#' @template args-sparse
#' 
#' @param formula,random,family,data,knots,drop.unused.levels Same as for 
#'   \code{\link[gamm4]{gamm4}}. \emph{We strongly advise against
#'   omitting the \code{data} argument}. Unless \code{data} is specified (and is
#'   a data frame) many post-estimation functions (including \code{update},
#'   \code{loo}, \code{kfold}) are not guaranteed to work properly.
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
#'   \code{\link[gamm4]{gamm4}} in the \pkg{gamm4} package. But rather than performing 
#'   (restricted) maximum likelihood estimation with the \pkg{lme4} package,
#'   the \code{stan_gamm4} function utilizes MCMC to perform Bayesian 
#'   estimation. The Bayesian model adds priors on the common regression 
#'   coefficients (in the same way as \code{\link{stan_glm}}), priors on the 
#'   standard deviations of the smooth terms, and a prior on the decomposition
#'   of the covariance matrices of any group-specific parameters (as in 
#'   \code{\link{stan_glmer}}). Estimating these models via MCMC avoids
#'   the optimization issues that often crop up with GAMMs and provides better
#'   estimates for the uncertainty in the parameter estimates. 
#'   
#'   See \code{\link[gamm4]{gamm4}} for more information about the model
#'   specicification and \code{\link{priors}} for more information about the
#'   priors on the main coefficients. The \code{formula} should include at least
#'   one smooth term, which can be specified in any way that is supported by the
#'   \code{\link[mgcv]{jagam}} function in the \pkg{mgcv} package. The 
#'   \code{prior_smooth} argument should be used to specify a prior on the unknown
#'   standard deviations that govern how smooth the smooth function is. The
#'   \code{prior_covariance} argument can be used to specify the prior on the
#'   components of the covariance matrix for any (optional) group-specific terms.
#'   The \code{\link[gamm4]{gamm4}} function in the \pkg{gamm4} package uses
#'   group-specific terms to implement the departure from linearity in the smooth
#'   terms, but that is not the case for \code{stan_gamm4} where the group-specific
#'   terms are exactly the same as in \code{\link{stan_glmer}}.
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
#'                  chains = 1, iter = 200) # for example speed
#' print(br)
#' plot_nonlinear(br)
#' plot_nonlinear(br, smooths = "s(x0)", alpha = 2/3)
#' }
#' 
stan_gamm4 <- function(formula, random = NULL, family = gaussian(), data, 
                       weights = NULL, subset = NULL, na.action, knots = NULL, 
                       drop.unused.levels = TRUE, ..., 
                       prior = normal(), prior_intercept = normal(),
                       prior_smooth = exponential(autoscale = FALSE), 
                       prior_aux = cauchy(0, 5),
                       prior_covariance = decov(), prior_PD = FALSE, 
                       algorithm = c("sampling", "meanfield", "fullrank"), 
                       adapt_delta = NULL, QR = FALSE, sparse = FALSE) {

  data <- validate_data(data, if_missing = list())
  family <- validate_family(family)
  
  if (!is.null(random)) {
    fake.formula <- as.character(mgcv::interpret.gam(formula)$fake.formula)
    form <- paste(fake.formula[2], fake.formula[1], fake.formula[3],
                  "+", random[2], collapse = " ")
    glmod <- lme4::glFormula(as.formula(form), data, family = gaussian,
                             subset, weights, na.action,
                             control = make_glmerControl())
    data <- glmod$fr
    weights <- validate_weights(glmod$fr$weights)
  }
  else {
    weights <- validate_weights(weights)
    glmod <- NULL
  }
  
  jd <- mgcv::jagam(formula = formula, family = gaussian(), data = data,
                    file = tempfile(fileext = ".jags"), weights = NULL,
                    na.action = na.action, offset = NULL, knots = knots,
                    drop.unused.levels = drop.unused.levels, diagonalize = TRUE)

  y <- jd$jags.data$y
  # there is no offset allowed by gamm4::gamm4
  offset <- validate_offset(as.vector(model.offset(jd$pregam$model)), y = y)
  X <- jd$jags.data$X
  mark <- which(colnames(X) != "")
  colnames(X) <- colnames(jd$pregam$X) <- jd$pregam$term.names
  S <- lapply(jd$pregam$smooth, FUN = function(s) {
    ranks <- s$rank
    start <- s$first.para
    out <- list()
    for (r in seq_along(ranks)) {
      end <- start + ranks[r] - 1L
      out[[r]] <- X[,start:end, drop = FALSE]
      start <- end + 1L
    }
    return(out)
  })
  if (any(sapply(S, length) > 1)) S <- unlist(S, recursive = FALSE)
  names(S) <- names(jd$pregam$sp)
  X <- X[,mark, drop = FALSE]
  X <- c(list(X), S)
  
  if (is.null(prior)) prior <- list()
  if (is.null(prior_intercept)) prior_intercept <- list()
  if (is.null(prior_aux)) prior_aux <- list()
  if (is.null(prior_smooth)) prior_smooth <- list()
  
  if (is.null(random)) {
    group <- list()
    prior_covariance <- list()
  }
  else {
    group <- glmod$reTrms
    group$decov <- prior_covariance
  }
  algorithm <- match.arg(algorithm)
  
  stanfit <- stan_glm.fit(x = X, y = y, weights = weights,
                          offset = offset, family = family,
                          prior = prior, prior_intercept = prior_intercept,
                          prior_aux = prior_aux, prior_smooth = prior_smooth,
                          prior_PD = prior_PD, algorithm = algorithm, 
                          adapt_delta = adapt_delta, group = group, QR = QR, ...)
  if (family$family == "Beta regression") family$family <- "beta"
  X <- do.call(cbind, args = X)
  if (is.null(random)) Z <- Matrix::Matrix(nrow = length(y), ncol = 0, sparse = TRUE)
  else {
    Z <- pad_reTrms(Ztlist = group$Ztlist, cnms = group$cnms, 
                    flist = group$flist)$Z
    colnames(Z) <- b_names(names(stanfit), value = TRUE)
  }
  if (getRversion() < "3.2.0") XZ <- cBind(X, Z) 
  else XZ <- cbind2(X, Z)
  
  # make jam object with point estimates, see ?mgcv::sim2jam
  mat <- as.matrix(stanfit)
  mark <- 1:ncol(X)
  jd$pregam$Vp <- cov(mat[,mark, drop = FALSE])
  jd$pregam$coefficients <- colMeans(mat[,mark, drop = FALSE])
  jd$pregam$sig2 <- if ("sigma" %in% colnames(mat)) mean(mat[,"sigma"]) else 1
  eta <- X %*% t(mat[,mark,drop = FALSE])
  mu <- rowMeans(family$linkinv(eta))
  eta <- rowMeans(eta)
  w <- as.numeric(jd$pregam$w * family$mu.eta(eta) ^ 2 / family$variance(mu))
  XWX <- t(X) %*% (w * X)
  jd$pregam$edf <- rowSums(jd$pregam$Vp * t(XWX)) / jd$pregam$sig2
  class(jd$pregam) <- c("jam", "gam")
  fit <- nlist(stanfit, family, formula, offset, weights, 
               x = XZ, y = y, data, terms = jd$pregam$terms, 
               model = if (is.null(random)) jd$pregam$model else glmod$fr,
               call = match.call(expand.dots = TRUE),
               algorithm, glmod = glmod,
               stan_function = "stan_gamm4")
  out <- stanreg(fit)
  out$jam <- jd$pregam
  class(out) <- c(class(out), "gamm4", if (!is.null(glmod)) "lmerMod")
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
  on.exit("try plot(x$jam) instead")
  scheme <- bayesplot::color_scheme_get()
  
  XZ <- x$x
  XZ <- XZ[,!grepl("_NEW_", colnames(XZ), fixed = TRUE)]
  labels <- sapply(x$jam$smooth, "[[", "label")
  xnames <- sapply(x$jam$smooth, "[[", "vn")
  names(x$jam$smooth) <- labels
  names(xnames) <- labels
  fs <- sapply(x$jam$smooth, FUN = "inherits", what = "fs.interaction")
  
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
    fs <- fs[found]
    if (!is.matrix(xnames)) xnames <- xnames[found]
  }
  else smooths <- 1:length(labels)
  
  B <- as.matrix(x)[, colnames(XZ), drop = FALSE]
  original <- x$jam$model
  
  bivariate <- any(grepl(",", labels, fixed = TRUE))
  if (bivariate && !any(fs)) {
    if (length(labels) > 1) {
      on.exit(NULL)
      stop("Multivariate functions can only be plotted one at a time; specify 'smooths'.")
    }
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
    XZ <- predict(x$jam, newdata = xz, type = "lpmatrix")
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
  
  df_list <- lapply(x$jam$smooth[smooths], FUN = function(s) {
    incl <- s$first.para:s$last.para
    b <- B[, incl, drop = FALSE]
    if (inherits(s, "fs.interaction")) { # see mgcv:::plot.fs.interaction
      xx <- original[,s$base$term]
      fac <- original[,s$fterm]
      out <- by(data.frame(fac, xx), list(fac), FUN = function(df) {
        df <- df[order(df[,2]),]
        names(df) <- c(s$fterm, s$base$term)
        xz <- mgcv::PredictMat(s, df)
        f <- linear_predictor.matrix(b, xz)
        data.frame(
          predictor = df[,2],
          lower  = apply(f, 2, quantile, probs = (1 - prob) / 2),
          upper  = apply(f, 2, quantile, probs = prob + (1 - prob) / 2),
          middle = apply(f, 2, median),
          term = paste(s$label, df[,1], sep = ".")
        )
      })
      do.call(rbind, args = out)
    }
    else {
      xz <- XZ[, incl, drop = FALSE]
      x <- original[, s$term]
      ord <- order(x)
      x <- x[ord]
      xz <- xz[ord, , drop=FALSE]
      if (!is.null(s$by.level)) {
        fac <- original[,s$by][ord]
        mark <- fac == s$by.level
        x <- x[mark]
        xz <- xz[mark, , drop = FALSE]
      }
      f <- linear_predictor.matrix(b, xz)
      data.frame(
        predictor = x,
        lower  = apply(f, 2, quantile, probs = (1 - prob) / 2),
        upper  = apply(f, 2, quantile, probs = prob + (1 - prob) / 2),
        middle = apply(f, 2, median),
        term = s$label
      )
    }
  })
  plot_data <- do.call(rbind, df_list)
  
  facet_args[["facets"]] <- ~ term
  if (is.null(facet_args[["scales"]]))
    facet_args[["scales"]] <- "free"
  if (is.null(facet_args[["strip.position"]]))
    facet_args[["strip.position"]] <- "left"

  on.exit(NULL)  
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
