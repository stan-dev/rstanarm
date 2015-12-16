.posterior_sample_size <- function(x) {
  stopifnot(is.stanreg(x))
  if (used.sampling(x)) return(sum(x$stanfit@sim$n_save - x$stanfit@sim$warmup2))
  else if (used.variational(x)) return(x$stanfit@sim$n_save) 
  else return(NULL)
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
#' 
#' @param ... Ignored.
#' @param pars Optional character vector specifying a subset of parameters to 
#'   display. Parameters can be specified by name or several shortcuts can be 
#'   used. Using \code{pars="beta"} will restrict the displayed parameters to 
#'   just the regression coefficients only. \code{"alpha"} can also be used as a
#'   shortcut for \code{"(Intercept)"}. If the model has varying intercepts 
#'   and/or slopes they can be selected using \code{pars = "varying"}. See 
#'   Examples.
#' @param probs For models fit using MCMC, an optional numeric vector of
#'   probabilities specifying which \code{\link[stats]{quantile}}s to display.
#' @param digits Number of digits to use for formatting numbers.
#' 
#' @seealso \code{\link{print.stanreg}}, \code{\link{stanreg-methods}}
#' 
#' @examples
#' summary(example_model, probs = c(0.1, 0.9))
#' 
#' # Only show varying intercept and slope terms
#' summary(example_model, pars = "varying") 
#' 
#' # These produce the same output for this example, 
#' # but the second method can be used for any model
#' summary(example_model, pars = c("(Intercept)", "size", 
#'                                 paste0("period", 2:4)))
#' summary(example_model, pars = c("alpha", "beta"))
#' 
#' @importMethodsFrom rstan summary
summary.stanreg <- function(object, ..., pars, probs, digits = 1) {
  if (object$stanfit@mode == 0) {
    if (!missing(pars)) pars[pars == "varying"] <- "b"
    out <- summary(object$stanfit, pars, probs)$summary
    if ("n_eff" %in% colnames(out))
      out[, "n_eff"] <- round(out[, "n_eff"])
    if ("se_mean" %in% colnames(out)) # So people don't confuse se_mean and sd
      colnames(out)[which(colnames(out) == "se_mean")] <- "mcse" 
  }
  else {
    if (missing(pars)) {
      mark <- names(object$coefficients)
      if (is.gaussian(object$family$family)) mark <- c(mark, "sigma")
      else if (is.nb(object$family$family)) mark <- c(mark, "overdispersion") 
    } else {
      mark <- NA
      if ("alpha" %in% pars) mark <- c(mark, "(Intercept)")
      if ("beta" %in% pars) mark <- c(mark, setdiff(names(object$coefficients), "(Intercept)"))
      mark <- c(mark, setdiff(pars, c("alpha", "beta")))
      mark <- mark[!is.na(mark)]
    }
    out <- object$stan_summary[mark,, drop=FALSE]
  }
  structure(out, 
            call = object$call, 
            algorithm = object$algorithm,
            posterior_sample_size = .posterior_sample_size(object),
            nobs = nobs(object),
            ngrps = if (!is.null(object$glmod)) ngrps(object) else NULL,
            print.digits = digits, 
            class = "summary.stanreg")
}


#' @method print summary.stanreg
#' @export
print.summary.stanreg <- function(x, ...) {
  dots <- list(...)
  atts <- attributes(x)
  digits <- dots$digits %ORifNULL% atts$print.digits
  print(atts$call)
  cat("\nAlgorithm:", atts$algorithm)
  if (!is.null(atts$posterior_sample_size) && atts$algorithm == "sampling")
    cat("\nPosterior sample size:", atts$posterior_sample_size)
  cat("\nObservations:", atts$nobs)
  if (!is.null(atts$ngrps))
    cat("\nGroups:", paste(names(atts$ngrps), unname(atts$ngrps), collapse = ", "))
  
  cat("\n\nEstimates:\n")
  sel <- which(colnames(x) %in% c("mcse", "n_eff", "Rhat"))
  if (!length(sel)) {
    print(format(round(x, digits), nsmall = digits), quote = FALSE, ...)
  } else {
    xtemp <- x[, -sel, drop = FALSE]
    colnames(xtemp) <- paste(" ", colnames(xtemp))
    print(format(round(xtemp, digits), nsmall = digits), quote = FALSE, ...)
    cat("\nDiagnostics:\n")
    mcse_rhat <- format(round(x[, c("mcse", "Rhat"), drop=FALSE], digits), nsmall = digits)
    n_eff <- format(x[, "n_eff", drop=FALSE], drop0trailing = TRUE)
    print(cbind(mcse_rhat, n_eff), quote = FALSE)
    cat("\nFor each parameter, mcse is Monte Carlo standard error, ", 
        "n_eff is a crude measure of effective sample size, ", 
        "and Rhat is the potential scale reduction factor on split chains", 
        " (at convergence Rhat=1).\n", sep = '')
  }
  invisible(x)
}
