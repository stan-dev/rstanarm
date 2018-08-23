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
#' fitted model. See the \strong{Details} section below for descriptions of the
#' different components of the printed output. For additional summary statistics
#' and diagnostics use the \code{\link[=summary.stanreg]{summary}} method.
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
#' \subsection{Uncertainty estimates (MAD_SD)}{
#' The standard deviations reported (labeled \code{MAD_SD} in the print output)
#' are computed from the same set of draws described above and are proportional
#' to the median absolute deviation (\code{\link[stats]{mad}}) from the median.
#' Compared to the raw posterior standard deviation, the MAD_SD will be
#' more robust for long-tailed distributions. These are the same as the values
#' returned by \code{\link[=se.stanreg]{se}}.
#' }
#' \subsection{Additional output}{
#' \itemize{
#' \item The median and MAD_SD are also reported for \code{mean_PPD}, the sample
#' average posterior predictive distribution of the outcome. This is useful as a
#' quick diagnostic. A useful heuristic is to check if \code{mean_PPD} is
#' plausible when compared to \code{mean(y)}. If it is plausible then this does
#' \emph{not} mean that the model is good in general (only that it can reproduce
#' the sample mean), however if \code{mean_PPD} is implausible then it is a sign
#' that something is wrong (severe model misspecification, problems with the
#' data, computational issues, etc.).
#' 
#' \item For GLMs with group-specific terms (see \code{\link{stan_glmer}}) the printed 
#' output also shows point estimates of the standard deviations of the group 
#' effects (and correlations if there are both intercept and slopes that vary by
#' group).
#' 
#' \item For analysis of variance models (see \code{\link{stan_aov}}) models, an
#' ANOVA-like table is also displayed.
#' 
#' \item For joint longitudinal and time-to-event (see \code{\link{stan_jm}}) models
#' the estimates are presented separately for each of the distinct submodels.  
#' }
#' }
#' 
#' @seealso \code{\link{summary.stanreg}}, \code{\link{stanreg-methods}}
#' 
print.stanreg <- function(x, digits = 1, ...) {
  cat(x$stan_function)
  surv <- is.surv(x)
  if (surv) {
    cat("\n baseline hazard:", basehaz_string(x$basehaz)) 
    cat("\n formula:        ", formula_string(formula(x)))
    cat("\n observations:   ", x$nobs) 
    cat("\n events:         ", x$nevents, percent_string(x$nevents, x$nobs))
    cat("\n censored:       ", x$ncensor, percent_string(x$ncensor, x$nobs))
    cat("\n delayed entry:  ", yes_no_string(x$ndelayed))
  } else {
    cat("\n family:      ", family_plus_link(x))
    cat("\n formula:     ", formula_string(formula(x)))
    cat("\n observations:", nobs(x))
    if (isTRUE(x$stan_function %in% c("stan_glm", "stan_glm.nb", "stan_lm", "stan_aov")))
      cat("\n predictors:  ", length(coef(x)))
  }
  
  cat("\n------\n")

  mer <- is.mer(x)
  gamm <- isTRUE(x$stan_function == "stan_gamm4")
  ord <- is_polr(x) && !("(Intercept)" %in% rownames(x$stan_summary))

  aux_nms <- .aux_name(x)
  
  if (!used.optimizing(x)) {
    
    if (isTRUE(x$stan_function %in% c("stan_lm", "stan_aov"))) {
      aux_nms <- c("R2", "log-fit_ratio", aux_nms)
    }
    mat <- as.matrix(x$stanfit) # don't used as.matrix.stanreg method b/c want access to mean_PPD
    nms <- setdiff(rownames(x$stan_summary), c("log-posterior", aux_nms))
    
    if (gamm) {
      smooth_sd_nms <- grep("^smooth_sd\\[", nms, value = TRUE)
      nms <- setdiff(nms, smooth_sd_nms)
      smooth_sd_mat <- mat[, smooth_sd_nms, drop = FALSE]
      smooth_sd_estimates <- .median_and_madsd(smooth_sd_mat)
    }
    if (mer) {
      nms <- setdiff(nms, grep("^b\\[", nms, value = TRUE))
    }
    if (ord) {
      cut_nms <- grep("|", nms, fixed = TRUE, value = TRUE)
      nms <- setdiff(nms, cut_nms)
      cut_mat <- mat[, cut_nms, drop = FALSE]
      cut_estimates <- .median_and_madsd(cut_mat)
    }
    
    ppd_nms <- grep("^mean_PPD", nms, value = TRUE)
    nms <- setdiff(nms, ppd_nms)
    coef_mat <- mat[, nms, drop = FALSE]
    estimates <- .median_and_madsd(coef_mat)
    
    if (mer) {
      estimates <- estimates[!grepl("^Sigma\\[", rownames(estimates)),, drop=FALSE]
    }
    .printfr(estimates, digits, ...)
    
    if (length(aux_nms)) {
      aux_estimates <- .median_and_madsd(mat[, aux_nms, drop=FALSE])
      cat("\nAuxiliary parameter(s):\n")
      .printfr(aux_estimates, digits, ...)
    }
    if (ord) {
      cat("\nCutpoints:\n")
      .printfr(cut_estimates, digits, ...)
    }
    if (gamm) {
      cat("\nSmoothing terms:\n")
      .printfr(smooth_sd_estimates, digits, ...)
    }
    if (mer) {
      cat("\nError terms:\n")
      print(VarCorr(x), digits = digits + 1, ...)
      cat("Num. levels:", 
          paste(names(ngrps(x)), unname(ngrps(x)), collapse = ", "), "\n")
    }
    if (is(x, "aov")) {
      print_anova_table(x, digits, ...)
    }
    if (!no_mean_PPD(x) && !is_clogit(x)) {
      ppd_mat <- mat[, ppd_nms, drop = FALSE]
      ppd_estimates <- .median_and_madsd(ppd_mat)
      
      cat("\nSample avg. posterior predictive distribution of y:\n")
      .printfr(ppd_estimates, digits, ...)
    }
    
  } else { 
    # used optimization
    nms <- names(x$coefficients)
    ppd_nms <- grep("^mean_PPD", rownames(x$stan_summary), value = TRUE)
    
    estimates <- x$stan_summary[nms, 1:2, drop=FALSE]
    .printfr(estimates, digits, ...)
    
    if (length(aux_nms)) {
      cat("\nAuxiliary parameter(s):\n")
      .printfr(x$stan_summary[aux_nms, 1:2, drop=FALSE], digits, ...)
    }
    if (length(ppd_nms) && !no_mean_PPD(x)) {
      cat("\nSample avg. posterior predictive distribution of y:\n")
      .printfr(x$stan_summary[ppd_nms, 1:2, drop=FALSE], digits, ...)
    }
  }
  
  cat("\n------\n")
  cat("* For help interpreting the printed output see ?print.stanreg\n")
  cat("* For info on the priors used see ?prior_summary.stanreg\n")
  
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
#'   In addition, for \code{stanmvreg} objects there are some additional shortcuts 
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
#'   object, while for a \code{stanmvreg} object all 
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
#'   \code{"summary.stanreg"} (or \code{"summary.stanmvreg"}, inheriting 
#'   \code{"summary.stanreg"}), which is a matrix of 
#'   summary statistics and 
#'   diagnostics, with attributes storing information for use by the
#'   \code{print} method. The \code{print} method for \code{summary.stanreg} or
#'   \code{summary.stanmvreg} objects is called for its side effect and just returns 
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
  surv <- is.surv(object)
  mer  <- is.mer(object)
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
  
  is_glm <- 
    isTRUE(object$stan_function %in% c("stan_glm", "stan_glm.nb", "stan_lm"))
  
  structure(
    out,
    call          = object$call,
    algorithm     = object$algorithm,
    stan_function = object$stan_function,
    family        = family_plus_link(object),
    formula       = formula(object),
    basehaz       = if (surv) basehaz_string(object$basehaz) else NULL,
    posterior_sample_size = posterior_sample_size(object),
    nobs          = nobs(object),
    npreds        = if (is_glm) length(coef(object)) else NULL,
    ngrps         = if (mer)  ngrps(object)   else NULL,
    nevents       = if (surv) object$nevents  else NULL,
    ncensor       = if (surv) object$ncensor  else NULL,
    ndelayed      = if (surv) object$ndelayed else NULL,
    print.digits  = digits,
    priors        = object$prior.info,
    class         = "summary.stanreg"
  )
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
  
  if (is.surv(atts)) { # survival models
    cat("\n function:       ", atts$stan_function)
    cat("\n baseline hazard:", atts$basehaz)
    cat("\n formula:        ", formula_string(atts$formula))
    cat("\n algorithm:      ", atts$algorithm)
    cat("\n priors:         ", "see help('prior_summary')")
    if (!is.null(atts$posterior_sample_size) && atts$algorithm == "sampling")
      cat("\n sample:         ", atts$posterior_sample_size, "(posterior sample size)")
    cat("\n observations:   ", atts$nobs)
    cat("\n events:         ", atts$nevents, percent_string(atts$nevents, atts$nobs))
    cat("\n censored:       ", atts$ncensor, percent_string(atts$ncensor, atts$nobs))
    cat("\n delayed entry:  ", yes_no_string(atts$ndelayed))
  } else { # anything except survival models
    cat("\n function:    ", atts$stan_function)
    cat("\n family:      ", atts$family)
    cat("\n formula:     ", formula_string(atts$formula))
    cat("\n algorithm:   ", atts$algorithm)
    cat("\n priors:      ", "see help('prior_summary')")
    if (!is.null(atts$posterior_sample_size) && atts$algorithm == "sampling")
      cat("\n sample:      ", atts$posterior_sample_size, "(posterior sample size)")
    cat("\n observations:", atts$nobs)
    if (!is.null(atts$npreds))
      cat("\n predictors:  ", atts$npreds)
    if (!is.null(atts$ngrps))
      cat("\n groups:      ", paste0(names(atts$ngrps), " (", 
                                     unname(atts$ngrps), ")", 
                                     collapse = ", "))
  }
  
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
      pars2 <- c(pars2, b_names(rownames(object$stan_summary), value = TRUE))
      pars[pars == "b"] <- NA
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
# @param ... Optionally include m to specify which submodel for stanmvreg models
family_plus_link <- function(x, ...) {
  if (is.stansurv(x)) {
    return(NULL)
  }
  fam <- family(x, ...)
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

# get name of aux parameter based on family
.aux_name <- function(object) {
  aux <- character()
  if (!is_polr(object)) {
    aux <- .rename_aux(family(object))
    if (is.na(aux)) {
      aux <- character()
    }
  }
  return(aux)
}

# print anova table for stan_aov models
# @param x stanreg object created by stan_aov()
print_anova_table <- function(x, digits, ...) {
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
  anova_table <- .median_and_madsd(effects)
  cat("\nANOVA-like table:\n")
  .printfr(anova_table, digits, ...)
}

# @param basehaz A list with info about the baseline hazard
basehaz_string <- function(basehaz, break_and_indent = TRUE) {
  nm <- get_basehaz_name(basehaz)
  switch(nm,
         exp      = "exponential",
         weibull  = "weibull",
         gompertz = "gompertz",
         ms       = "M-splines on hazard scale",
         bs       = "B-splines on log hazard scale",
         piecewise= "piecewise constant on log hazard scale",
         NULL)
}

# @param x A logical (or a scalar to be evaluated as a logical).
yes_no_string <- function(x) {
  if (x) "yes" else "no"
}

# @param numer,denom The numerator and denominator with which to evaluate a %.
percent_string <- function(numer, denom) {
  val <- round(100 * numer / denom, 1)
  paste0("(", val, "%)")
}

# equivalent to isFALSE(object$compute_mean_PPD)
no_mean_PPD <- function(object) {
  x <- object$compute_mean_PPD
  is.logical(x) && length(x) == 1L && !is.na(x) && !x
}
