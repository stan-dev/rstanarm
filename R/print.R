#' @method print stanreg
#' @export
print.stanreg <- function(x, digits = 1, ...) {
  mer <- is(x, "lmerMod")
  if (x$algorithm != "optimizing") {
    nms <- setdiff(rownames(x$stan_summary), "log-posterior")
    if (mer) nms <- setdiff(nms, grep("^b\\[", nms, value = TRUE))
    mat <- as.matrix(x$stanfit)[,nms,drop=FALSE]
    estimates <- cbind(Median = apply(mat, 2, median), 
                       MAD_SD = apply(mat, 2, mad))
  }
  else {
    nms <- names(x$coefficients)
    if (is(x, "polr")) 
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
    # if (mer) nms <- setdiff(nms, grep("^b\\[", nms, value = TRUE))
    estimates <- x$stan_summary[nms,1:2]
  }
  
  print(x$call)
  cat("\nAlgorithm:", x$algorithm)
  cat("\n\n")
  print(format(round(estimates, digits), nsmall = digits), quote = FALSE, ...)
  
  if (mer) {
    cat("\nGroups:\n")
    print(ngrps(x))
    cat("\nError terms:\n")
    print(VarCorr(x), digits = digits + 1, ...)
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
    cat("\nANOVA-like table\n")
    anova_table <- cbind(Median = apply(effects, 2, median),
                         MAD_SD = apply(effects, 2, mad))
    print(format(round(anova_table, digits), nsmall = digits), quote = FALSE, ...)
  }
  return(invisible(NULL))
}
