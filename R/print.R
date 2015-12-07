.printfr <- function(x, digits, ...) {
  print(format(round(x, digits), nsmall = digits), quote = FALSE, ...)
}
.median_and_madsd <- function(x) {
  cbind(Median = apply(x, 2, median), MAD_SD = apply(x, 2, mad))
}

#' Print method for stanreg objects
#' 
#' The \code{print} method for stanreg objects displays a compact summary of the
#' fitted model. For a more detailed summary of the fitted model use the
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
#' The estimates reported by \code{print} depend on the modeling function used
#' and the estimation algorithm. For all models, regression coefficients are 
#' summarized by:
#' \describe{
#' \item{Point estimates}{
#'  If \code{algorithm='sampling'}, posterior medians are used as point
#'  estimates. If \code{algorithm='optimizing'}, the point estimates are also
#'  medians, but they are computed from 1000 draws from the asymptotic sampling
#'  distribution of the parameters. In all cases, these are the same as the
#'  estimates returned by calling \code{\link[=coef.stanreg]{coef}}.
#' }
#' \item{Uncertainty estimates}{
#'  The 'standard errors' reported (which we denote MAD_SD) are proportional to
#'  the median absolute deviation (\code{\link[stats]{mad}}) from the posterior
#'  median (if \code{algorithm='sampling'}) or the median of the asymptotic
#'  sampling distribution (if \code{algorithm='optimizing'}). These are the same
#'  as the standard errors returned by calling \code{\link[=se.stanreg]{se}}.}
#' }
#' 
#' For models fit using MCMC (\code{algorithm='sampling'}) the median and MAD_SD
#' are also reported for \code{mean_PPD}, the sample average (\eqn{X =
#' \bar{X}}{X = xbar}) \code{\link[=posterior_predict]{posterior predictive
#' distribution}} of the outcome.
#' 
#' For \code{\link[=stan_glmer]{GLMs with group-specific terms}} the printed 
#' output also shows point estimates of the standard deviations of the group 
#' effects (and correlations if there are both intercept and slopes that vary by
#' group).
#' 
#' For \code{\link[=stan_aov]{analysis of variance}} models, an ANOVA-like table
#' is also displayed.
#' 
#' @seealso \code{\link{summary.stanreg}}, \code{\link{stanreg-methods}}
#' 
print.stanreg <- function(x, digits = 1, ...) {
  print(x$call)
  cat("\nEstimates:\n")
  
  mer <- is(x, "lmerMod")
  ord <- is(x, "polr") && !("(Intercept)" %in% rownames(x$stan_summary))
  if (!used.optimizing(x)) {
    # don't used as.matrix.stanreg method b/c want access to mean_PPD
    mat <- as.matrix(x$stanfit)
    nms <- setdiff(rownames(x$stan_summary), "log-posterior")
    if (mer) nms <- setdiff(nms, grep("^b\\[", nms, value = TRUE))
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
  }
  else { # used optimizing
    nms <- names(x$coefficients)
    if (ord) 
      nms <- c(nms, grep("|", rownames(x$stan_summary), 
                         fixed = TRUE, value = TRUE))
    else {
      famname <- x$family$family
      if (is.gaussian(famname))
        nms <- c(nms, "sigma")
      else if (is.gamma(famname))
        nms <- c(nms, "shape")
      else if (is.ig(famname))
        nms <- c(nms, "lambda")
      else if (is.nb(famname))
        nms <- c(nms, "overdispersion")
    }
    nms <- c(nms, grep("^mean_PPD", rownames(x$stan_summary), value = TRUE))
    estimates <- x$stan_summary[nms,1:2]
    .printfr(estimates, digits, ...)
  }
  
  if (is(x, "aov")) {
    labels <- attributes(x$terms)$term.labels
    patterns <- gsub(":", ".*:", labels)
    groups <- sapply(patterns, simplify = FALSE, FUN = grep, 
                     x = dimnames(extract(x$stanfit, pars = "beta", permuted = FALSE))$parameters)
    names(groups) <- gsub(".*", "", names(groups), fixed = TRUE)
    groups <- groups[sapply(groups, length) > 0]
    effects_dim <- dim(x$effects)
    effects <- x$effects^2
    effects <- sapply(groups, FUN = function(i) apply(effects[,,i,drop = FALSE], 1:2, mean))
    dim(effects) <- c(effects_dim[-3], ncol(effects))
    dim(effects) <- c(nrow(effects) * ncol(effects), dim(effects)[3])
    colnames(effects) <- paste("Mean Sq", names(groups))
    cat("\nANOVA-like table:\n")
    anova_table <- .median_and_madsd(effects)
    .printfr(anova_table, digits, ...)
  }
  return(invisible(x))
}
