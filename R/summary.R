#' @method summary stanreg
#' @export
summary.stanreg <- function(object, ..., digits = 1) {
  # use RStan's summary just as placeholder. we should replace this with our own
  # summary 
  if (object$stanfit@mode == 0) {
    out <- summary(object$stanfit, ...)$summary
    if ("n_eff" %in% colnames(out))
      out[, "n_eff"] <- round(out[, "n_eff"])
    if ("se_mean" %in% colnames(out)) # So people don't confuse se_mean and sd
      colnames(out)[which(colnames(out) == "se_mean")] <- "mcse" 
  }
  else {
    mark <- names(object$coefficients)
    if (is.gaussian(object$family$family)) mark <- c(mark, "sigma")
    else if (is.nb(object$family$family)) mark <- c(mark, "overdispersion")
    out <- object$stan_summary[mark,,drop=FALSE]
  }
  SS <- if (object$algorithm != "sampling") 
    NULL else sum(object$stanfit@sim$n_save - object$stanfit@sim$warmup2)
  
  structure(out, call = object$call, algorithm = object$algorithm,
            posterior_sample_size = SS, print.digits = digits, 
            class = "summary.stanreg")
}


#' @method print summary.stanreg
#' @export
print.summary.stanreg <- function(x, ...) {
  dots <- list(...)
  digits <- dots$digits %ORifNULL% attr(x, "print.digits")
  cat("\n")
  print(attr(x, "call"))
  cat("\nAlgorithm:", attr(x, "algorithm"))
  if (!is.null(attr(x, "posterior_sample_size")))
    cat("\nPosterior sample size:", attr(x, "posterior_sample_size"))
  
  cat("\n\nEstimates:\n")
  sel <- which(colnames(x) %in% c("n_eff", "Rhat"))
  if (!length(sel)) {
    print(format(round(x, digits), nsmall = digits), quote = FALSE, ...)
  } else {
    xtemp <- x[, -sel, drop = FALSE]
    colnames(xtemp) <- paste(" ", colnames(xtemp))
    print(format(round(xtemp, digits), nsmall = digits), quote = FALSE, ...)
    cat("\n")
    Rhat <- format(round(x[, "Rhat", drop=FALSE], digits), nsmall = digits)
    n_eff <- format(x[, "n_eff", drop=FALSE], drop0trailing = TRUE)
    print(cbind(n_eff, Rhat), quote = FALSE)
    cat("\nFor each parameter, n_eff is a crude measure of effective sample size,\n", 
        "and Rhat is the potential scale reduction factor on split chains (at \n",
        "convergence, Rhat=1).\n", sep = '')
  }
  invisible(x)
}
