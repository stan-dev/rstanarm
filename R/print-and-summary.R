# Part of the rstanarm package for estimating model parameters
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
#' 
#' For joint longitudinal and time-to-event (see \code{\link{stan_jm}}) models
#' the estimates are presented separately for each of the distinct submodels.  
#' }
#' 
#' @seealso \code{\link{summary.stanreg}}, \code{\link{stanreg-methods}}
#' 
print.stanreg <- function(x, digits = 1, ...) {
  cat(x$modeling_function)
  cat("\n family: ", family_plus_link(x))
  cat("\n formula:", formula_string(formula(x)))
  
  cat("\n------\n")
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
    if (mer)
      estimates <- estimates[!grepl("^Sigma\\[", rownames(estimates)),, drop=FALSE]
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
      nms <- c(nms, "reciprocal_dispersion")
    } else if (is.beta(famname)) {}
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

  cat("\n------\n")
  cat("For info on the priors used see help('prior_summary.stanreg').")
  
  invisible(x)
}

#' @rdname print.stanreg
#' @export
#' @method print stanjm
print.stanjm <- function(x, digits = 3, ...) {
  print(x$call) 
  
  M <- x$n_markers
  link    <- sapply(1:M, function(m) x$family[[m]]$link)
  
  mat <- as.matrix(x$stanfit)
  nms <- collect_nms(rownames(x$stan_summary), M, value = TRUE)
  
  # Estimates table for longitudinal submodel(s)
  for (m in 1:M) {
    terms_m <- terms(x)[[m]]
    sel <- attr(terms_m, "response")
    yvar <- rownames(attr(terms_m, "factors"))[sel]
    cat(paste0("\nLongitudinal submodel", if (M > 1) paste0(" ", m), 
               ": ", yvar,"\n"))
    coef_mat <- mat[, c(nms$y[[m]], nms$y_extra[[m]]), drop = FALSE]
    
    # Calculate median and MAD
    estimates <- .median_and_madsd(coef_mat)
    
    # Add column with eform
    if (link[m] %in% c("log", "logit")) 
      estimates <- cbind(estimates, 
                         "exp(Median)" = c(exp(estimates[nms$y[[m]], "Median"]), 
                                           rep(NA, length(nms$y_extra[[m]]))))
    
    # Print estimates
    rownames(estimates) <- 
      gsub(paste0("^Long", m, "\\|"), "", rownames(estimates))     
    .printfr(estimates, digits, ...)
  }
  
  # Estimates table for event submodel
  cat("\nEvent submodel:\n")   
  coef_mat <- mat[, c(nms$e, nms$a, nms$e_extra), drop = FALSE]
  
  # Calculate median and MAD
  estimates <- .median_and_madsd(coef_mat)
  
  # Add column with eform
  estimates <- cbind(estimates, 
                     "exp(Median)" = c(exp(estimates[c(nms$e, nms$a), "Median"]), 
                                       rep(NA, length(nms$e_extra))))
  
  rownames(estimates) <- gsub("^Event\\|", "", rownames(estimates))  
  rownames(estimates) <- gsub("^Assoc\\|", "", rownames(estimates))   
  .printfr(estimates, digits, ...)
  
  # Estimates table for group-level random effects
  cat("\nGroup-level random effects:\n") 
  print(VarCorr(x), digits = digits + 1, ...)
  cat("Num. levels:", paste(names(ngrps(x)), unname(ngrps(x)), 
                            collapse = ", "), "\n")  
  
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
#'   "varying"}.
#'   
#'   In addition, for \code{stanjm} objects there are some additional shortcuts 
#'   available. Using \code{pars = "long"} will display the 
#'   parameter estimates for the longitudinal submodels only (excluding group-specific
#'   pparameters, but including auxiliary parameters).
#'   Using \code{pars = "event"} will display the 
#'   parameter estimates for the event submodel only, including any association
#'   parameters. 
#'   Using \code{pars = "assoc"} will display only the 
#'   association parameters. 
#'   Using \code{pars = "fixef"} will display all fixed effects, but not
#'   the random effects or the auxiliary parameters. 
#'    \code{pars} and \code{regex_pars} are set to \code{NULL} then all 
#'   fixed effect regression coefficients are selected, as well as any 
#'   auxiliary parameters and the log posterior.   
#'   
#'   If \code{pars} is \code{NULL} all parameters are selected for a \code{stanreg}
#'   object, while for a \code{stanjm} object all 
#'   fixed effect regression coefficients are selected as well as any 
#'   auxiliary parameters and the log posterior. See 
#'   \strong{Examples}.
#' @param probs For models fit using MCMC or one of the variational algorithms, 
#'   an optional numeric vector of probabilities passed to 
#'   \code{\link[stats]{quantile}}.
#' @param digits Number of digits to use for formatting numbers when printing. 
#'   When calling \code{summary}, the value of digits is stored as the 
#'   \code{"print.digits"} attribute of the returned object.
#'   
#' @return The \code{summary} method returns an object of class 
#'   \code{"summary.stanreg"} (or \code{"summary.stanjm"}, inheriting 
#'   \code{"summary.stanreg"}), which is a matrix of 
#'   summary statistics and 
#'   diagnostics, with attributes storing information for use by the
#'   \code{print} method. The \code{print} method for \code{summary.stanreg} or
#'   \code{summary.stanjm} objects is called for its side effect and just returns 
#'   its input. The \code{as.data.frame} method for \code{summary.stanreg} 
#'   objects converts the matrix to a data.frame, preserving row and column 
#'   names but dropping the \code{print}-related attributes.
#' 
#' @seealso \code{\link{prior_summary}} to extract or print a summary of the 
#'   priors used for a particular model.
#' 
#' @examples
#' if (!exists("example_model")) example(example_model) 
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
#' as.data.frame(summary(example_model, pars = "varying"))
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
      pars <- allow_special_parnames(object, pars)
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
        mark <- c(mark, "reciprocal_dispersion")
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
            modeling_function = object$modeling_function,
            family = family_plus_link(object),
            formula = formula(object),
            posterior_sample_size = posterior_sample_size(object),
            nobs = nobs(object),
            ngrps = if (mer) ngrps(object) else NULL,
            print.digits = digits, 
            priors = object$prior.info,
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
  cat("\nModel Info:\n")
  cat("\n function: ", atts$modeling_function)
  cat("\n family:   ", atts$family)
  cat("\n formula:  ", formula_string(atts$formula))
  cat("\n algorithm:", atts$algorithm)
  cat("\n priors:   ", "see help('prior_summary')")
  if (!is.null(atts$posterior_sample_size) && atts$algorithm == "sampling")
    cat("\n sample:   ", atts$posterior_sample_size, "(posterior sample size)")
  cat("\n num obs:  ", atts$nobs)
  if (!is.null(atts$ngrps))
    cat("\n groups:   ", paste0(names(atts$ngrps), " (", unname(atts$ngrps), ")",
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

#' @rdname summary.stanreg
#' @export
#' @method summary stanjm
summary.stanjm <- function(object, pars = NULL, regex_pars = NULL, 
                           probs = NULL, ..., digits = 3) {
  pars <- collect_pars(object, pars, regex_pars)
  M <- object$n_markers
  
  # Outcome variable for each submodel
  y_vars <- sapply(1:M, function(m, object) {
    terms_m <- terms(object)[[m]]
    sel <- attr(terms_m, "response")
    ret <- rownames(attr(terms_m, "factors"))[sel]
  }, object = object)
  
  # Family and link for each submodel
  fam <- sapply(object$family, function(x) 
    paste0(x$family, " (", x$link, ")")) 
  
  # Construct summary table  
  args <- list(object = object$stanfit)
  if (!is.null(probs)) 
    args$probs <- probs
  out <- do.call("summary", args)$summary
  
  nms <- collect_nms(rownames(object$stan_summary), M, value = TRUE)
  if (!is.null(pars)) {
    pars2 <- NA     
    if ("alpha" %in% pars) pars2 <- c(pars2, nms$alpha)
    if ("beta" %in% pars) pars2 <- c(pars2, nms$beta)
    if ("long" %in% pars) pars2 <- c(pars2, unlist(nms$y), unlist(nms$y_extra))
    if ("event" %in% pars) pars2 <- c(pars2, nms$e, nms$a, nms$e_extra)
    if ("assoc" %in% pars) pars2 <- c(pars2, nms$a)      
    if ("fixef" %in% pars) pars2 <- c(pars2, unlist(nms$y), nms$e, nms$a)
    if ("b" %in% pars) pars2 <- c(pars2, nms$b)
    pars2 <- c(pars2, setdiff(pars, 
                              c("alpha", "beta", "varying", "b",
                                "long", "event", "assoc", "fixef")))
    pars <- pars2[!is.na(pars2)]
  } else {
    pars <- c(unlist(nms$y), unlist(nms$y_extra),
              nms$e, nms$a, nms$e_extra, "log-posterior")
  }
  
  out <- out[rownames(out) %in% pars, , drop = FALSE]
  out <- out[!grepl(":_NEW_", rownames(out), fixed = TRUE), , drop = FALSE]
  stats <- colnames(out)
  if ("n_eff" %in% stats)
    out[, "n_eff"] <- round(out[, "n_eff"])
  if ("se_mean" %in% stats) # So people don't confuse se_mean and sd
    colnames(out)[stats %in% "se_mean"] <- "mcse"
  
  # Reorder rows of output table
  nms_tmp <- rownames(out)  
  nms_tmp_y <- lapply(1:M, function(m) 
    grep(paste0("^Long", m, "\\|"), nms_tmp, value = TRUE))
  nms_tmp_e <- grep("^Event\\|", nms_tmp, value = TRUE)
  nms_tmp_a <- grep("^Assoc\\|", nms_tmp, value = TRUE)
  nms_tmp_b <- b_names(nms_tmp, value = TRUE)
  nms_tmp_lp <- grep("^log-posterior$", nms_tmp, value = TRUE)
  out <- out[c(unlist(nms_tmp_y), nms_tmp_e, nms_tmp_a, nms_tmp_b, nms_tmp_lp), , drop = FALSE]
  
  # Run times
  times <- round((rstan::get_elapsed_time(object$stanfit))/60, digits = 1)
  times <- cbind(times, total = rowSums(times))
  
  # Output object
  structure(out, 
            call = object$call, 
            algorithm = object$algorithm,
            n_markers = object$n_markers,
            n_subjects = object$n_subjects,
            n_grps = object$n_grps,
            n_events = object$n_events,
            n_yobs = object$n_yobs,
            y_vars = y_vars,
            id_var = object$id_var,
            time_var = object$time_var,
            family = fam,
            basehaz = object$basehaz$type_name,
            posterior_sample_size = posterior_sample_size(object),
            times = times,
            print.digits = digits, 
            class = c("summary.stanjm", "summary.stanreg"))
}

#' @rdname summary.stanreg
#' @export
#' @method print summary.stanjm
print.summary.stanjm <- function(x, digits = max(1, attr(x, "print.digits")), 
                                 ...) {
  atts <- attributes(x)
  M <- atts$n_markers
  
  print(atts$call)
  
  cat(paste0("\n", if (M == 1) "Uni" else "Multi", 
             "variate joint model, consisting of:")) 
  for (m in 1:M) {
    cat(paste0("\n  Outcome variable", 
               if (M > 1) paste0(" (Long ", m, ")"), 
               ": ", atts$y_vars[m]))
  }
  for (m in 1:M) {
    cat(paste0("\n  Family", 
               if (M > 1) paste0(" (Long ", m, ")"), 
               ": ", atts$family[m]))
  }
  cat(paste0("\n  Baseline hazard: ", atts$basehaz))  
  cat(paste0("\n  Clustering variables: ", paste(names(atts$n_grps), sep = ",")))
  if (!is.null(atts$n_subjects))
    cat(paste0("\n  Num. subjects (", atts$id_var, "): ", atts$n_subjects))
  cat(paste0("\n  Num. events: ", atts$n_events, " (", 
             round(100 * atts$n_events/atts$n_subjects, 1), "%)"))
  cat("\n  Num. long observations: ")
  cat(paste0(atts$n_yobs, if (M > 1) paste0(" (Long ", 1:M, ")"), collapse = ", "))
  cat("\n  Posterior sample size:", atts$posterior_sample_size, "MCMC iterations")
  
  cat("\n\nTime taken for sampling (mins):\n")
  print(atts$times)
  
  cat("\nEstimates:\n")
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


# internal ----------------------------------------------------------------
.printfr <- function(x, digits, ...) {
  print(format(round(x, digits), nsmall = digits), quote = FALSE, ...)
}
.median_and_madsd <- function(x) {
  cbind(Median = apply(x, 2, median), MAD_SD = apply(x, 2, mad))
}

# Allow "alpha", "beta", "varying" as shortcuts 
#
# @param object stanreg object
# @param pars result of calling collect_pars(object, pars, regex_pars)
allow_special_parnames <- function(object, pars) {
  pars[pars == "varying"] <- "b"
  pars2 <- NA
  if ("alpha" %in% pars)
    pars2 <- c(pars2, "(Intercept)")
  if ("beta" %in% pars) {
    beta_nms <- if (is.mer(object))
      names(fixef(object)) else names(object$coefficients)
    pars2 <- c(pars2, setdiff(beta_nms, "(Intercept)"))
  }
  if ("b" %in% pars) {
    if (is.mer(object)) {
      pars2 <-
        c(pars2, b_names(rownames(object$stan_summary), value = TRUE))
    } else {
      warning("No group-specific parameters. 'varying' ignored.",
              call. = FALSE)
    }
  }
  pars2 <- c(pars2, setdiff(pars, c("alpha", "beta", "varying")))
  pars2[!is.na(pars2)]
}

# Family name with link in parenthesis 
# @param x stanreg object
family_plus_link <- function(x) {
  fam <- family(x)
  if (is.character(fam)) {
    stopifnot(identical(fam, x$method))
    fam <- paste0("ordered [", fam, "]")
  } else if (inherits(x, "betareg")) {
    fam <- paste0("beta [",
                  x$family$link,
                  ", link.phi=",
                  x$family_phi$link,
                  "]")
  } else {
    fam <- paste0(fam$family, " [", fam$link, "]")
  }
  return(fam)
}

# @param formula formula object
formula_string <- function(formula, break_and_indent = TRUE) {
  coll <- if (break_and_indent) "--MARK--" else " "
  char <- gsub("\\s+", " ", paste(deparse(formula), collapse = coll))
  if (!break_and_indent)
    return(char)
  gsub("--MARK--", "\n\t  ", char, fixed = TRUE)
}

