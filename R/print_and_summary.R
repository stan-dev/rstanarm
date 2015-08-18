print.stanreg <- function(x, ...) {
  # use RStan's print just as placeholder. we should replace this with our own
  # print method
  if (x$stanfit@mode == 0) print(x$stanfit, pars = "lp__", include = FALSE, ...)
  else if (is.null(x$family)) {
    mark <- c(names(x$coefficients), 
              grep("|", rownames(x$stan_summary), fixed = TRUE, value = TRUE))
    print(x$stan_summary[mark,,drop = FALSE], ...)
  }
  else {
    mark <- names(x$coefficients)
    if (x$family$family == "gaussian") mark <- c(mark, "sigma")
    else if (x$family$family == "Negative Binomial") mark <- c(mark, "overdispersion")
    print(x$stan_summary[mark,,drop=FALSE], ...)
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
    dimnames(effects) <- list(iterations = NULL, chains = paste0("chain", 1:ncol(x$stanfit), sep = ":"),
                              parameters = paste("Mean Sq", names(groups)))
    cat("ANOVA-like table\n")
    rstan::monitor(effects, warmup = 0, print = TRUE, digits_summary = 2)
  }
}

summary.stanreg <- function(object, ..., digits = 2) {
  # use RStan's summary just as placeholder. we should replace this with our own
  # summary 
  if (object$stanfit@mode == 0) {
    out <- round(summary(object$stanfit, ...)$summary, digits = digits)
    if ("n_eff" %in% colnames(out))
      out[, "n_eff"] <- round(out[, "n_eff"])
    if ("se_mean" %in% colnames(out)) # So people don't confuse se_mean and sd
      colnames(out)[which(colnames(out) == "se_mean")] <- "mcse" 
    out
  }
  else {
    mark <- names(object$coefficients)
    if (object$family$family == "gaussian") mark <- c(mark, "sigma")
    else if (object$family$family == "Negative Binomial") mark <- c(mark, "overdispersion")
    object$stan_summary[mark,,drop=FALSE]
  }
}
