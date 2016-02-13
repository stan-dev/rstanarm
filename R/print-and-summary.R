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

.printfr <- function(x, digits, ...) {
  print(format(round(x, digits), nsmall = digits), quote = FALSE, ...)
}
.median_and_madsd <- function(x) {
  cbind(Median = apply(x, 2, median), MAD_SD = apply(x, 2, mad))
}

#' Print method for stanreg objects
#' 
#' The \code{print} method for stanreg objects displays a compact summary of the
#' fitted model. See the Details section below for a description of the printed 
#' output. For additional summary statistics and diagnostics use the 
#' \code{\link[=summary.stanreg]{summary}} method.
#' 
#' @export
#' @method print stanreg
#' @templateVar stanregArg x
#' @template args-stanreg-object
#' @param digits Number of digits to use for formatting numbers.
#' @param ... Ignored.
#' @return Returns \code{x}, invisibly.
#' @details 
#' \subsection{Point estimates}{
#' Regardless of the estimation algorithm, point estimates are medians computed 
#' from simulations. For models fit using MCMC (\code{"sampling"}) the posterior
#' sample is used. For optimization (\code{"optimizing"}), the simulations are
#' generated from the asymptotic Gaussian sampling distribution of the
#' parameters. For the \code{"meanfield"} and \code{"fullrank"} variational
#' approximations, draws from the variational approximation to the posterior are
#' used. In all cases, the point estimates reported are the same as the values
#' returned by \code{\link[=coef.stanreg]{coef}}.
#' }
#' \subsection{Uncertainty estimates}{
#' The standard deviations reported (labeled MAD_SD in the print output) are 
#' computed from the same set of draws described above and are proportional to 
#' the median absolute deviation (\code{\link[stats]{mad}}) from the median. 
#' Compared to the raw posterior standard deviation, the MAD_SD will be more 
#' robust for long-tailed distributions. These are the same as the values 
#' returned by \code{\link[=se.stanreg]{se}}.
#' }
#' \subsection{Additional output}{
#' For models fit using MCMC or a variational approximation, the median and 
#' MAD_SD are also reported for mean_PPD, the sample average (\eqn{X = 
#' \bar{X}}{X = xbar}) posterior predictive distribution of the outcome.
#' 
#' For GLMs with group-specific terms (see \code{\link{stan_glmer}}) the printed 
#' output also shows point estimates of the standard deviations of the group 
#' effects (and correlations if there are both intercept and slopes that vary by
#' group).
#' 
#' For analysis of variance models (see \code{\link{stan_aov}}) models, an
#' ANOVA-like table is also displayed.
#' }
#' 
#' @seealso \code{\link{summary.stanreg}}, \code{\link{stanreg-methods}}
#' 
print.stanreg <- function(x, digits = 1, ...) {
  print(x$call)
  cat("\nEstimates:\n")
  
  mer <- is.mer(x)
  ord <- is(x, "polr") && !("(Intercept)" %in% rownames(x$stan_summary))
  if (!used.optimizing(x)) {
    mat <- as.matrix(x$stanfit) # don't used as.matrix.stanreg method b/c want access to mean_PPD
    nms <- setdiff(rownames(x$stan_summary), "log-posterior")
    if (mer) 
      nms <- setdiff(nms, grep("^b\\[", nms, value = TRUE))
    if (ord) {
      cut_nms <- grep("|", nms, fixed = TRUE, value = TRUE)
      nms <- setdiff(nms, cut_nms)
      cut_mat <- mat[, cut_nms, drop = FALSE]
      cut_estimates <- .median_and_madsd(cut_mat)
    }
    ppd_nms <- grep("^mean_PPD", nms, value = TRUE)
    nms <- setdiff(nms, ppd_nms)
    coef_mat <- mat[, nms, drop = FALSE]
    ppd_mat <- mat[, ppd_nms, drop = FALSE]
    estimates <- .median_and_madsd(coef_mat)
    ppd_estimates <- .median_and_madsd(ppd_mat)
    
    .printfr(estimates, digits, ...)
    if (ord) {
      cat("\nCutpoints:\n")
      .printfr(cut_estimates, digits, ...)
    }
    if (mer) {
      cat("\nError terms:\n")
      print(VarCorr(x), digits = digits + 1, ...)
      cat("Num. levels:", 
          paste(names(ngrps(x)), unname(ngrps(x)), collapse = ", "), "\n")
    }
    cat("\nSample avg. posterior predictive \ndistribution of y (X = xbar):\n")
    .printfr(ppd_estimates, digits, ...)
    
  } else { 
    # used optimization
    nms <- names(x$coefficients)
    famname <- family(x)$family
    if (is.gaussian(famname)) {
      nms <- c(nms, "sigma")
    } else if (is.gamma(famname)) {
      nms <- c(nms, "shape")
    } else if (is.ig(famname)) {
      nms <- c(nms, "lambda")
    } else if (is.nb(famname)) {
      nms <- c(nms, "overdispersion")
    }
    nms <- c(nms, grep("^mean_PPD", rownames(x$stan_summary), value = TRUE))
    estimates <- x$stan_summary[nms,1:2]
    .printfr(estimates, digits, ...)
  }
  
  if (is(x, "aov")) {
    labels <- attributes(x$terms)$term.labels
    patterns <- gsub(":", ".*:", labels)
    dnms <- dimnames(extract(x$stanfit, pars = "beta", 
                             permuted = FALSE))$parameters
    groups <- sapply(patterns, simplify = FALSE, FUN = grep, x = dnms)
    names(groups) <- gsub(".*", "", names(groups), fixed = TRUE)
    groups <- groups[sapply(groups, length) > 0]
    effects_dim <- dim(x$effects)
    effects <- x$effects^2
    effects <- sapply(groups, FUN = function(i) {
      apply(effects[, , i, drop = FALSE], 1:2, mean)
    })
    dim(effects) <- c(effects_dim[-3], ncol(effects))
    dim(effects) <- c(nrow(effects) * ncol(effects), dim(effects)[3])
    colnames(effects) <- paste("Mean Sq", names(groups))
    cat("\nANOVA-like table:\n")
    anova_table <- .median_and_madsd(effects)
    .printfr(anova_table, digits, ...)
  }
  
  invisible(x)
}


#' Summary method for stanreg objects
#' 
#' Summaries of parameter estimates and MCMC convergence diagnostics 
#' (Monte Carlo error, effective sample size, Rhat).
#' 
#' @export
#' @method summary stanreg
#' 
#' @templateVar stanregArg object
#' @template args-stanreg-object
#' @template args-regex-pars
#' 
#' @param ... Currently ignored.
#' @param pars An optional character vector specifying a subset of parameters to
#'   display. Parameters can be specified by name or several shortcuts can be 
#'   used. Using \code{pars="beta"} will restrict the displayed parameters to 
#'   only the regression coefficients (without the intercept). \code{"alpha"} 
#'   can also be used as a shortcut for \code{"(Intercept)"}. If the model has 
#'   varying intercepts and/or slopes they can be selected using \code{pars = 
#'   "varying"}. See Examples.
#' @param probs For models fit using MCMC or one of the variational algorithms, 
#'   an optional numeric vector of probabilities passed to 
#'   \code{\link[stats]{quantile}}.
#' @param digits Number of digits to use for formatting numbers when printing. 
#'   When calling \code{summary}, the value of digits is stored as the 
#'   \code{"print.digits"} attribute of the returned object.
#'   
#' @return The \code{summary} method returns an object of class 
#'   \code{"summary.stanreg"}, which is a matrix of summary statistics and 
#'   diagnostics, with attributes storing information for use by the
#'   \code{print} method. The \code{print} method for \code{summary.stanreg}
#'   objects is called for its side effect and does not return anything. The 
#'   \code{as.data.frame} method for \code{summary.stanreg} objects converts the
#'   matrix to a data.frame, preserving row and column names but dropping the 
#'   \code{print}-related attributes.
#' 
#' @seealso \code{\link{print.stanreg}}, \code{\link{stanreg-methods}}
#' 
#' @examples
#' summary(example_model, probs = c(0.1, 0.9))
#' 
#' # These produce the same output for this example, 
#' # but the second method can be used for any model
#' summary(example_model, pars = c("(Intercept)", "size", 
#'                                 paste0("period", 2:4)))
#' summary(example_model, pars = c("alpha", "beta"))
#' 
#' # Only show parameters varying by group
#' summary(example_model, pars = "varying") 
#' 
#' @importMethodsFrom rstan summary
summary.stanreg <- function(object, pars = NULL, regex_pars = NULL, 
                            probs = NULL, ..., digits = 1) {
  mer <- is.mer(object)
  pars <- collect_pars(object, pars, regex_pars)
  if (!used.optimizing(object)) {
    args <- list(object = object$stanfit)
    if (!is.null(probs)) 
      args$probs <- probs
    out <- do.call("summary", args)$summary
    
    if (is.null(pars) && used.variational(object))
      out <- out[!rownames(out) %in% "log-posterior", , drop = FALSE]
    if (!is.null(pars)) {
      pars2 <- NA
      if ("alpha" %in% pars) 
        pars2 <- c(pars2, "(Intercept)")
      if ("beta" %in% pars) {
        beta_nms <- if (mer) 
          names(fixef(object)) else names(object$coefficients)
        pars2 <- c(pars2, setdiff(beta_nms, "(Intercept)"))
      }
      if ("b" %in% pars) {
        if (mer) {
          pars2 <- c(pars2, b_names(rownames(object$stan_summary), value = TRUE))
        } else {
          warning("No group-specific parameters. 'varying' ignored.", 
                  call. = FALSE) 
        }
      }
      pars2 <- c(pars2, setdiff(pars, c("alpha", "beta", "varying")))
      pars <- pars2[!is.na(pars2)]
      out <- out[rownames(out) %in% pars, , drop = FALSE]
    }
    out <- out[!grepl(":_NEW_", rownames(out), fixed = TRUE), , drop = FALSE]
    stats <- colnames(out)
    if ("n_eff" %in% stats)
      out[, "n_eff"] <- round(out[, "n_eff"])
    if ("se_mean" %in% stats) # So people don't confuse se_mean and sd
      colnames(out)[stats %in% "se_mean"] <- "mcse"
    
  } else { # used optimization
    if (!is.null(probs)) 
      warning("'probs' ignored if for models fit using optimization.",
              call. = FALSE)
    if (is.null(pars)) {
      famname <- family(object)$family
      mark <- names(object$coefficients)
      if (is.gaussian(famname)) 
        mark <- c(mark, "sigma")
      if (is.nb(famname)) 
        mark <- c(mark, "overdispersion") 
    } else {
      mark <- NA
      if ("alpha" %in% pars) 
        mark <- c(mark, "(Intercept)")
      if ("beta" %in% pars) 
        mark <- c(mark, setdiff(names(object$coefficients), "(Intercept)"))
      mark <- c(mark, setdiff(pars, c("alpha", "beta")))
      mark <- mark[!is.na(mark)]
    }
    out <- object$stan_summary[mark, , drop=FALSE]
  }
  structure(out, 
            call = object$call, 
            algorithm = object$algorithm,
            posterior_sample_size = posterior_sample_size(object),
            nobs = nobs(object),
            ngrps = if (mer) ngrps(object) else NULL,
            print.digits = digits, 
            class = "summary.stanreg")
}

#' @rdname summary.stanreg
#' @export
#' @method print summary.stanreg
#'
#' @param x An object of class \code{"summary.stanreg"}.
print.summary.stanreg <- function(x, digits = max(1, attr(x, "print.digits")), 
                                  ...) {
  atts <- attributes(x)
  print(atts$call)
  cat("\nAlgorithm:", atts$algorithm)
  if (!is.null(atts$posterior_sample_size) && atts$algorithm == "sampling")
    cat("\nPosterior sample size:", atts$posterior_sample_size)
  cat("\nObservations:", atts$nobs)
  if (!is.null(atts$ngrps))
    cat("\nGroups:", paste(names(atts$ngrps), unname(atts$ngrps), 
                           collapse = ", "))
  
  cat("\n\nEstimates:\n")
  sel <- which(colnames(x) %in% c("mcse", "n_eff", "Rhat"))
  if (!length(sel)) {
    .printfr(x, digits)
  } else {
    xtemp <- x[, -sel, drop = FALSE]
    colnames(xtemp) <- paste(" ", colnames(xtemp))
    .printfr(xtemp, digits)
    cat("\nDiagnostics:\n")
    mcse_rhat <- format(round(x[, c("mcse", "Rhat"), drop = FALSE], digits), 
                        nsmall = digits)
    n_eff <- format(x[, "n_eff", drop = FALSE], drop0trailing = TRUE)
    print(cbind(mcse_rhat, n_eff), quote = FALSE)
    cat("\nFor each parameter, mcse is Monte Carlo standard error, ", 
        "n_eff is a crude measure of effective sample size, ", 
        "and Rhat is the potential scale reduction factor on split chains", 
        " (at convergence Rhat=1).\n", sep = '')
  }
  invisible(x)
}

#' @rdname summary.stanreg
#' @method as.data.frame summary.stanreg
#' @export
as.data.frame.summary.stanreg <- function(x, ...) {
  as.data.frame(unclass(x), ...)
}
