#' Summary method for stanreg objects
#' 
#' Summaries of parameter estimates and MCMC convergence diagnostics 
#' (Monte Carlo error, effective sample size, Rhat).
#' 
#' @method summary stanreg
#' @export
#' 
#' @param object A fitted model object returned by one of the \pkg{rstanarm} 
#'   modeling functions. This will be a list with class 'stanreg' as well as at 
#'   least one of 'lm', 'glm', 'polr', 'lmerMod', or 'aov'.
#' @param ... Ignored.
#' @param pars Optional character vector specifying a subset of parameters to 
#'   display. If \code{pars} is missing all parameters are used. Parameters can 
#'   be specified by name or several shortcuts can be used. Using 
#'   \code{pars="beta"} will restrict the displayed parameters to just the 
#'   regression coefficients (excluding the intercept) only. \code{"alpha"} can
#'   also be used as a shortcut for \code{"(Intercept)"}. If the model has 
#'   varying intercepts and/or slopes, these can be selected using \code{pars = 
#'   "varying"}. See Examples.
#' @param probs For models fit using MCMC, an optional numeric vector of
#'   probabilities specifying which \code{\link[stats]{quantile}}s to display.
#' @param digits Number of digits to use for formatting numbers.
#' 
#' @examples
#' \dontrun{
#' summary(stan_glm(mpg ~ wt, data = mtcars), probs = c(0.1, 0.9))
#' 
#' fit <- stan_glmer(mpg ~ wt + (1|cyl) + (1+wt|gear), data = mtcars)
#' summary(fit)
#' 
#' # Only show varying intercept and slope terms
#' summary(fit, pars = "varying") 
#' 
#' # These produce the same output for this example, 
#' # but the second method can be used for any model
#' summary(fit, pars = c("(Intercept)", "wt")) 
#' summary(fit, pars = c("alpha", "beta"))
#' }   
#' 
summary.stanreg <- function(object, ..., pars, probs, digits = 1) {
  # use RStan's summary just as placeholder. we should replace this with our own
  # summary 
  
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
      mark <- na.omit(mark)
    }
    out <- object$stan_summary[mark,, drop=FALSE]
  }
  SS <- if (object$algorithm != "sampling") 
    NULL else sum(object$stanfit@sim$n_save - object$stanfit@sim$warmup2)
  groups <- if (!is.null(object$glmod)) ngrps(object) else NULL
  
  structure(out, 
            call = object$call, 
            algorithm = object$algorithm,
            nobs = NROW(get_y(object)),
            posterior_sample_size = SS, 
            ngrps = groups,
            print.digits = digits, 
            class = "summary.stanreg")
}


#' @method print summary.stanreg
#' @export
print.summary.stanreg <- function(x, ...) {
  dots <- list(...)
  atts <- attributes(x)
  digits <- dots$digits %ORifNULL% atts$print.digits
  cat("\n")
  print(atts$call)
  cat("\nAlgorithm:", atts$algorithm)
  if (!is.null(atts$posterior_sample_size))
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
