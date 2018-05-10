# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2018 Sam Brilleman
# Copyright (C) 2018 Trustees of Columbia University
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

#' Bayesian proportional hazards regression
#'
#' \if{html}{\figure{stanlogo.png}{options: width="25px" alt="http://mc-stan.org/about/logo/"}}
#' Bayesian inference for proportional hazards regression models. The user
#' can specify a variety of standard parametric distributions for the
#' baseline hazard, or a flexible parametric model using M-splines for 
#' approximating the baseline hazard.
#'
#' @export
#'
#' @template args-prior_intercept
#' @template args-priors
#' @template args-prior_PD
#' @template args-algorithm
#' @template args-adapt_delta
#' @template args-max_treedepth
#' 
#' @param formula A two-sided formula object describing the model. 
#'   The left hand side of the formula should be a \code{Surv()} 
#'   object. See \code{\link[survival]{Surv}}.
#' @param data A data frame containing the variables specified in 
#'   \code{formula}.
#' @param basehaz A character string indicating which baseline hazard to use
#'   for the event submodel. Current options are: 
#'   \itemize{
#'     \item \code{"ms"}: a flexible parametric model using M-splines to 
#'     model the baseline hazard (or, equivalently, I-splines to model the 
#'     cumulative hazard). The default locations for the internal knots, as 
#'     well as the basis terms for the splines, are calculated with respect
#'     to log time. A closed form solution is available for both the hazard
#'     and cumulative hazard so this approach should be relatively fast.
#'     \item \code{"bs"}: a flexible parametric model using B-splines to model
#'     the \emph{log} baseline hazard. The default locations for the internal 
#'     knots, as well as the basis terms for the splines, are calculated with 
#'     respect to log time. A closed form solution for the cumulative hazard 
#'     is \strong{not} available and instead quadrature is used to evaluate 
#'     the cumulative hazard at each MCMC iteration, so this approach is much 
#'     slower than specifying \code{basehaz = "ms"}. The advantage is that it 
#'     can accomodate time-dependent effects (ie. non-proportional hazards).
#'     \item \code{"exp"}: an exponential distribution for the event times. 
#'     (i.e. a constant baseline hazard)
#'     \item \code{"weibull"}: a Weibull distribution for the event times.
#'     \item \code{"gompertz"}: a Gompertz distribution for the event times.
#'   }
#'   Note that all spline-based models used splines of degree 3 (ie. cubics).
#' @param basehaz_ops A named list specifying options related to the baseline
#'   hazard. Currently this can include: \cr
#'   \describe{
#'     \item{\code{df}}{A positive integer specifying the degrees of freedom 
#'     for the M-splines / I-splines. The default is 6.}
#'     \item{\code{knots}}{An optional numeric vector specifying the internal 
#'     knot locations for the splines if \code{basehaz = "ms"}. Knots cannot be
#'     specified if \code{df} is specified. If not specified, then the 
#'     default is to use \code{df - 4} knots, which are
#'     placed at equally spaced percentiles of the distribution of
#'     uncensored event times.}
#'     \item{\code{bknots}}{An optional numeric vector specifying the boundary 
#'     knot locations for the splines if \code{basehaz = "ms"}. 
#'     If not specified, then the default is to place the boundary knots at the
#'     minimum and maximum of the event times (including both censored and
#'     uncensored events).}     
#'   }
#' @param qnodes The number of nodes to use for the Gauss-Kronrod quadrature
#'   that is used to evaluate the cumulative hazard when \code{basehaz = "bs"}. 
#'   Options are 15 (the default), 11 or 7.
#' @param prior_aux The prior distribution for the parameters related to the
#'   baseline hazard. The "auxiliary" parameters refers to different  
#'   parameters depending on the type of baseline hazard specified in the 
#'   \code{basehaz} argument. All auxiliary parameters have a lower bound at
#'   zero.
#'   When \code{basehaz = "exp"} there is no auxiliary 
#'   parameter, since the log scale parameter is incorporated as an intercept 
#'   in the linear predictor.
#'   When \code{basehaz = "weibull"} the auxiliary parameter is the Weibull 
#'   shape parameter, whilst the log scale parameter is incorporated as an 
#'   intercept in the linear predictor. 
#'   When \code{basehaz = "ms"} or \code{basehaz = "bs"} the auxiliary 
#'   parameters are the coefficients for the spline basis terms.
#'   
#'   \code{prior_aux} can be a call to \code{exponential} to 
#'   use an exponential distribution, or \code{normal}, \code{student_t} or 
#'   \code{cauchy}, which results in a half-normal, half-t, or half-Cauchy 
#'   prior. See \code{\link{priors}} for details on these functions. To omit a 
#'   prior ---i.e., to use a flat (improper) uniform prior--- set 
#'   \code{prior_aux} to \code{NULL}.
#' 
#' @examples
#' options(mc.cores = parallel::detectCores())
#' 
#' #--- Simulated data
#' library(simsurv)
#' covs <- data.frame(id = 1:2000, trt = stats::rbinom(2000, 1L, 0.5))
#' dat1 <- simsurv(lambdas = 0.1, gammas = 1.5, 
#'                 betas = c(trt = -0.5),
#'                 x = covs, maxt = 5)
#' dat1 <- merge(dat1, covs)
#' mod1a <- stan_surv(Surv(eventtime, status) ~ trt, dat1, basehaz = "ms")
#' #mod1b <- stan_surv(Surv(eventtime, status) ~ trt, dat1, basehaz = "exp")
#' #mod1c <- stan_surv(Surv(eventtime, status) ~ trt, dat1, basehaz = "weibull")
#' mod1a$stanfit
#' #mod1b$stanfit
#' #mod1c$stanfit
#' 
#' #--- Breast cancer data
#' library(flexsurv)
#' dat2 <- flexsurv::bc
#' mod2a <- stan_surv(Surv(rectime, censrec) ~ group, dat2, basehaz = "ms")
#' #mod2b <- stan_surv(Surv(rectime, censrec) ~ group, dat2, basehaz = "exp")
#' #mod2c <- stan_surv(Surv(rectime, censrec) ~ group, dat2, basehaz = "weibull")
#' mod2d <- stan_surv(Surv(rectime, censrec) ~ group, dat2, basehaz = "bs")
#' mod2z <- flexsurvspline(Surv(rectime, censrec) ~ group, dat2, k = 3)
#' mod2a$stanfit
#' #mod2b$stanfit
#' #mod2c$stanfit
#' mod2d$stanfit
#' mod2z
#'                 
#' #--- PBC data
#' dat3 <- rstanarm::pbcSurv
#' mod3a <- stan_surv(Surv(futimeYears, death) ~ sex + trt, dat3)
#' mod3z <- flexsurvspline(Surv(futimeYears, death) ~ sex + trt, dat3, k = 3)
#' mod3a$stanfit
#' mod3z
#'
stan_surv <- function(formula, data, basehaz = "ms", basehaz_ops, qnodes = 15,
                      prior = normal(), prior_intercept = normal(),
                      prior_aux = cauchy(), prior_PD = FALSE,
                      algorithm = c("sampling", "meanfield", "fullrank"),
                      adapt_delta = 0.95, ...) {

  #-----------------------------
  # Pre-processing of arguments
  #-----------------------------
  
  if (!requireNamespace("survival"))
    stop("the 'survival' package must be installed to use this function.")
  
  if (missing(basehaz_ops)) basehaz_ops <- NULL
  
  dots <- list(...)
  algorithm <- match.arg(algorithm)

  # Formula
  formula <- parse_formula(formula, data)

  # Data
  data <- as.data.frame(data)

  #----------------
  # Construct data
  #----------------

  # standardised weights and nodes for quadrature
  qq <- get_quadpoints(nodes = qnodes)
  qp <- qq$points
  qw <- qq$weights
  
  # model data frame
  mf <- data
  
  #----- dimensions, response, predictor matrix

  # entry and exit times for each row of data
  t_end <- make_t(formula, mf, type = "end") # end time
  t_beg <- make_t(formula, mf, type = "beg") # beg time
  
  # event indicator for each row of data
  d <- make_d(formula, mf)
  event <- as.logical(d)
  
  # delayed entry indicator for each row of data
  delayed <- (!t_beg == 0)

  # time variables for stan
  t_events  <- t_end[event]
  t_censor  <- t_end[!event]
  t_delayed <- t_beg[delayed]

  # predictor matrix
  x_stuff <- make_x(formula, mf) 
  x <- x_stuff$x
  x_events  <- keep_rows(x, event)
  x_censor  <- keep_rows(x, !event)
  x_delayed <- keep_rows(x, delayed)
  
  # dimensions
  npats    <- 0L # not currently used
  K        <- ncol(x)
  nrows    <- nrow(x)
  nevents  <- nrow(x_events)
  ncensor  <- nrow(x_censor)
  ndelayed <- nrow(x_delayed)
  
  #----- time-dependent effects (i.e. non-proportional hazards)

  # degrees of freedom for time-dependent effects
  nvars_tde <- aa(rep(0L, K)) # not currently used
  has_tde <- any(nvars_tde > 0)

  #----- baseline hazard

  ok_basehaz <- c("exp", "weibull", "gompertz", "ms", "bs")
  ok_basehaz_ops <- get_ok_basehaz_ops(basehaz)
  basehaz <- handle_basehaz(basehaz = basehaz, 
                            basehaz_ops = basehaz_ops, 
                            ok_basehaz = ok_basehaz,
                            ok_basehaz_ops = ok_basehaz_ops,
                            times = t_end, status = event)
  type_name <- basehaz$type_name
  type      <- basehaz$type
  nvars     <- basehaz$nvars
  has_closed_form <- check_for_closed_form(type_name)
  
  # basis terms
  basis_events  <- make_basis(t_events,  basehaz)
  basis_censor  <- make_basis(t_censor,  basehaz)
  basis_delayed <- make_basis(t_delayed, basehaz)
  deriv_events  <- make_basis(t_events,  basehaz, deriv = TRUE)

  # flag if intercept is required for baseline hazard
  has_intercept  <- ai(has_intercept(type_name))
  
  # flag for quadrature
  has_quadrature <- has_tde || !has_closed_form
  
  if (has_quadrature) { # model uses quadrature
    
    # quadrature points & weights, evaluated for each row of data
    qpts <- uapply(qp, unstandardise_qpts, t_beg, t_end) # qpts for exit time
    qwts <- uapply(qw, unstandardise_qwts, t_beg, t_end) # qwts for exit time

    # quadrature points & weights, evaluated for rows with delayed entry
    if (length(t_delayed)) {
      qpts_delayed <- uapply(qp, unstandardise_qpts, 0, t_delayed) # qpts for entry time
      qwts_delayed <- uapply(qw, unstandardise_wpts, 0, t_delayed) # qwts for entry time
    } else {
      qpts_delayed <- rep(0,0)
      qwts_delayed <- rep(0,0)
    }
    
    # predictor matrix at quadrature points
    x_qpts         <- rep_rows(x, times = qnodes)         # tde not yet implemented
    x_qpts_delayed <- rep_rows(x_delayed, times = qnodes) # tde not yet implemented
    
    # basis terms at quadrature points
    basis_qpts         <- make_basis(qpts, basehaz)
    basis_qpts_delayed <- make_basis(qpts_delayed, basehaz)
    
    # dimensions
    qrows    <- nrow(x_qpts)
    qdelayed <- nrow(x_qpts_delayed)
    
  } else { # model does not use quadrature
    
    # dud entries for stan
    qpts               <- rep(0,0)
    qwts               <- rep(0,0)
    qpts_delayed       <- rep(0,0)
    qwts_delayed       <- rep(0,0)
    x_qpts             <- matrix(0,0,K)
    x_qpts_delayed     <- matrix(0,0,K)
    basis_qpts         <- matrix(0,0,nvars)
    basis_qpts_delayed <- matrix(0,0,nvars)
    qrows              <- 0L
    qdelayed           <- 0L
    
  }
  
  #----- stan data
  
  standata <- nlist(
    K,
    npats,
    nrows,
    nevents,
    ncensor,
    ndelayed,
    qnodes,
    qrows,
    qdelayed,
    nvars,
    nvars_tde,
    t_events,
    t_censor,
    t_delayed,
    x_events,
    x_censor,
    x_delayed,
    x_qpts,
    x_qpts_delayed,
    basis_events,
    basis_censor,
    basis_delayed,
    basis_qpts,
    basis_qpts_delayed,
    deriv_events,
    type,
    qwts,
    qwts_delayed,
    has_intercept,
    has_quadrature
  )
  
  #----- priors and hyperparameters

  # valid priors
  ok_dists <- nlist("normal", 
                    student_t = "t", 
                    "cauchy", 
                    "hs", 
                    "hs_plus", 
                    "laplace", 
                    "lasso") # disallow product normal
  ok_intercept_dists <- ok_dists[1:3]
  ok_aux_dists       <- ok_dists[1:3]
  
  # priors
  user_prior_stuff <- prior_stuff <-
    handle_glm_prior(prior, 
                     nvars = K,
                     default_scale = 2,
                     link = NULL,
                     ok_dists = ok_dists)
  
  user_prior_intercept_stuff <- prior_intercept_stuff <-
    handle_glm_prior(prior_intercept, 
                     nvars = 1,
                     default_scale = 20,
                     link = NULL,
                     ok_dists = ok_intercept_dists)
  
  user_prior_aux_stuff <- prior_aux_stuff <-
    handle_glm_prior(prior_aux, 
                     nvars = basehaz$nvars,
                     default_scale = get_default_aux_scale(basehaz$type_name),
                     link = NULL,
                     ok_dists = ok_aux_dists)
  
  # autoscaling of priors
  prior_stuff           <- autoscale_prior(prior_stuff, predictors = x)
  prior_intercept_stuff <- autoscale_prior(prior_intercept_stuff)
  prior_aux_stuff       <- autoscale_prior(prior_aux_stuff)

  # priors
  standata$prior_dist              <- prior_stuff$prior_dist
  standata$prior_dist_for_intercept<- prior_intercept_stuff$prior_dist
  standata$prior_dist_for_aux      <- prior_aux_stuff$prior_dist

  # hyperparameters
  standata$prior_mean               <- prior_stuff$prior_mean
  standata$prior_scale              <- prior_stuff$prior_scale
  standata$prior_df                 <- prior_stuff$prior_df
  standata$prior_mean_for_intercept <- c(prior_intercept_stuff$prior_mean)
  standata$prior_scale_for_intercept<- c(prior_intercept_stuff$prior_scale)
  standata$prior_df_for_intercept   <- c(prior_intercept_stuff$prior_df)
  standata$prior_scale_for_aux      <- prior_aux_stuff$prior_scale
  standata$prior_df_for_aux         <- prior_aux_stuff$prior_df
  standata$global_prior_scale       <- prior_stuff$global_prior_scale
  standata$global_prior_df          <- prior_stuff$global_prior_df
  standata$slab_df                  <- prior_stuff$slab_df
  standata$slab_scale               <- prior_stuff$slab_scale

  # any additional flags
  standata$prior_PD <- ai(prior_PD)

  #-----------
  # Fit model
  #-----------

  stanfit  <- stanmodels$surv
  
  stanpars <- c(if (standata$has_intercept) "gamma",
                if (standata$K)             "beta",
                if (standata$nvars > 0)     "coefs")
  
  if (algorithm == "sampling") {
    args <- set_sampling_args(
      object = stanfit, 
      data   = standata, 
      pars   = stanpars, 
      prior  = prior, 
      user_dots = list(...), 
      user_adapt_delta = adapt_delta, 
      show_messages = FALSE)
    stanfit <- do.call(rstan::sampling, args)
  } else {
    args <- nlist(
      object = stanfit,
      data   = standata,
      pars   = stanpars,
      algorithm # meanfield or fullrank vb
    )
    args[names(dots)] <- dots
    stanfit <- do.call(rstan::vb, args)
  }
  #check_stanfit(stanfit)

  fit <- nlist(stanfit, formula, data, basehaz, algorithm,
               stan_function = "stan_surv", call = match.call(expand.dots = TRUE))

  #out <- stansurv(fit)
  return(fit)
}


#---------- internal

make_basis <- function(times, basehaz, timescale = "log", deriv = FALSE) {

  if (is.null(times) || !length(times)) {
   
    X <- matrix(0, 0, basehaz$nvars)
    
  } else if (basehaz$type_name %in% c("exp", "weibull", "gompertz")) {
    
    X <- matrix(0, length(times), basehaz$nvars) # dud matrix for Stan
    
  } else if (basehaz$type_name == "bs") {
    
    if (is.null(basehaz$basis))
      stop2("Bug found: could not find info on spline basis.")
    
    X <- as.array(predict(basehaz$basis, times))
    
  } else if (basehaz$type_name == "piecewise") {
    
    if (is.null(basehaz$knots))
      stop2("Bug found: could not find info on knot locations.")
    
    X <- dummy_matrix(times, knots = basehaz$knots)
    
  } else if (basehaz$type_name == "ms") {
    
    timescale <- basehaz$timescale
    if (is.null(timescale)) {
      tt <- times
    } else if (timescale == "log") {
      tt <- log(times)
    }
    
    if (is.null(basehaz$basis))
      stop2("Bug found: could not find info on spline basis.")
    
    if (deriv) { # M-splines, i.e. derivative of I-spline basis
      X <- as.array(deriv(predict(basehaz$basis, tt)))
    } else { # I-splines
      X <- as.array(predict(basehaz$basis, tt))
    }
    
  } else {
    
    stop2("Bug found: type of baseline hazard unknown.") 
    
  }
  X
}

# Parse the model formula
#
# @param formula The user input to the formula argument.
# @param data The user input to the data argument (i.e. a data frame).
parse_formula <- function(formula, data) {
  
  formula <- validate_formula(formula, needs_response = TRUE)
  
  lhs <- lhs(formula) # full LHS of formula
  rhs <- rhs(formula) # full RHS of formula
  
  lhs_form <- reformulate_lhs(lhs)
  rhs_form <- reformulate_rhs(rhs)
  
  allvars <- all.vars(formula)
  allvars_form <- reformulate(allvars)
  
  surv <- eval(lhs, envir = data) # Surv object
  surv <- validate_surv(surv)
  type <- attr(surv, "type")
  
  if (type == "right") {
    tvar_beg <- NULL
    tvar_end <- as.character(lhs[[2L]])
    dvar     <- as.character(lhs[[3L]])
  } else if (type == "counting") {
    tvar_beg <- as.character(lhs[[2L]])
    tvar_end <- as.character(lhs[[3L]])
    dvar     <- as.character(lhs[[4L]])
  }
  
  nlist(lhs = lhs,
        rhs = rhs,
        lhs_form = lhs_form,
        rhs_form = rhs_form,
        fe_form = rhs_form, # no re terms accommodated yet
        re_form = NULL,     # no re terms accommodated yet
        allvars = allvars,
        allvars_form = allvars_form,
        tvar_beg = tvar_beg,
        tvar_end = tvar_end,
        dvar = dvar,
        surv_type = attr(surv, "type"))
}

# Check formula object
#
# @param formula The user input to the formula argument.
# @param needs_response A logical; if TRUE then formula must contain a LHS.
validate_formula <- function(formula, needs_response = TRUE) {
  
  if (!inherits(formula, "formula")) {
    stop2("'formula' must be a formula.")
  }
  
  if (needs_response) {
    len <- length(formula)
    if (len < 3) {
      stop2("'formula' must contain a response.")
    }
  }
  as.formula(formula)
}

# Check object is a Surv object with a valid type
#
# @param x A Surv object; the LHS of a formula evaluated in a data frame environment.
# @param ok_types A character vector giving the allowed types of Surv object.
validate_surv <- function(x, ok_types = c("right", "counting")) {
  
  if (!inherits(x, "Surv")) {
    stop2("LHS of 'formula' must be a 'Surv' object.")
  }
  
  if (!attr(x, "type") %in% ok_types) {
    stop2("Surv object type must be one of: ", comma(ok_types))
  }
  x
}


# Extract LHS of a formula
#
# @param x A formula object
# @param as_formula Logical. If TRUE then the result is reformulated.
lhs <- function(x, as_formula = FALSE) {
  len <- length(x)
  if (len == 3L) {
    out <- x[[2L]]
  } else {
    out <- NULL
  }
  out
}

# Extract RHS of a formula
#
# @param x A formula object
# @param as_formula Logical. If TRUE then the result is reformulated.
rhs <- function(x, as_formula = FALSE) {
  len <- length(x)
  if (len == 3L) {
    out <- x[[3L]]
  } else {
    out <- x[[2L]]
  }
  out
}

# Reformulate as LHS of a formula
#
# @param x A character string or expression object
# @param as_formula Logical. If TRUE then the result is reformulated.
reformulate_lhs <- function(x) {
  #x <- deparse(x, 500L)
  x <- formula(substitute(LHS ~ 1, list(LHS = x)))
  x
}

# Reformulate as RHS of a formula
#
# @param x A formula object
# @param as_formula Logical. If TRUE then the result is reformulated.
reformulate_rhs <- function(x) {
  #x <- deparse(x, 500L)
  x <- formula(substitute(~ RHS, list(RHS = x)))
  x
}

# Identify whether the type of baseline hazard requires an intercept in
# the linear predictor (NB splines incorporate the intercept into the basis)
#
# @param basehaz_name Character string, the type of baseline hazard
has_intercept <- function(basehaz_name) {
  (basehaz_name %in% c("exp", "weibull", "gompertz"))
}

# Return the response vector (time)
#
# @param formula The parsed model formula.
# @param data The model frame.
# @param type The type of time variable to return.
# @return A numeric vector
make_t <- function(formula, data, type = c("beg", "end", "gap")) {
  
  type <- match.arg(type)
  
  if (formula$surv_type == "right") {
    t_beg <- rep(0, nrow(data))
    t_end <- data[[formula$tvar_end]]
  } else if (formula$surv_type == "counting") {
    t_beg <- data[[formula$tvar_beg]]
    t_end <- data[[formula$tvar_end]]
  } else {
    stop2("Cannot yet handle '", formula$surv_type, "' type Surv objects.")
  }
  
  if (type == "beg") {
    out <- t_beg
  } else if (type == "end") {
    out <- t_end
  } else if (type == "gap") {
    out <- t_end - t_beg
  }
  out
}

# Return the response vector (status indicator)
#
# @param formula The parsed model formula.
# @param data The model frame.
# @return A numeric vector
make_d <- function(formula, data) {
  
  if (formula$surv_type == "right") {
    out <- data[[formula$dvar]]
  } else if (formula$surv_type == "counting") {
    out <- data[[formula$dvar]]
  } else {
    stop2("Bug found: cannot handle '", formula$surv_type, "' Surv objects.")
  }
  out
}

# Return the fe predictor matrix
#
# @param formula The parsed model formula.
# @param model_frame The model frame.
# @return A named list with the following elements:
#   x: the fe model matrix, not centred and may have intercept.
#   xtemp: fe model matrix, centred and no intercept.
#   x_form: the formula for the fe model matrix.
#   x_bar: the column means of the model matrix.
#   has_intercept: logical for whether the submodel has an intercept
#   N,K: number of rows (observations) and columns (predictors) in the
#     fixed effects model matrix
make_x <- function(formula, data) {
  
  x <- model.matrix(formula$fe_form, data)
  x <- drop_intercept(x)
  xbar <- colMeans(x)
  
  # identify any column of x with < 2 unique values (empty interaction levels)
  sel <- (apply(x, 2L, n_distinct) < 2)
  if (any(sel)) {
    cols <- paste(colnames(x)[sel], collapse = ", ")
    stop2("Cannot deal with empty interaction levels found in columns: ", cols)
  }
  
  nlist(x, xbar, N = NROW(x), K = NCOL(x))
}

# Replace an NA object, or NA entries in a vector
#
# @param x The vector with elements to potentially replace.
# @param replace_with The replacement value.
replace_na <- function(x, replace_with = "0") {
  if (is.na(x)) {
    x <- replace_with
  } else {
    x[is.na(x)] <- replace_with
  }
  x
}

# Replace an NULL object, or NULL entries in a vector
#
# @param x The vector with elements to potentially replace.
# @param replace_with The replacement value.
replace_null <- function(x, replace_with = "0") {
  if (is.null(x)) {
    x <- replace_with
  } else {
    x[is.null(x)] <- replace_with
  }
  x
}

# Replace named elements of 'x' with 'y'
replace_named_elements <- function(x, y) {
  x[names(y)] <- y
  x
}

# Check if all elements of a vector are zeros
all_zero <- function(x) {
  all(x == 0)
}

# Shorthand for as.integer, as.double, as.matrix, as.array
ai <- function(x, ...) as.integer(x, ...)
ad <- function(x, ...) as.double(x, ...)
am <- function(x, ...) as.matrix(x, ...)
aa <- function(x, ...) as.array(x, ...)

# Return a vector of 0's
zeros <- function(n) {
  rep(0, times = n)
}

# Return a vector of 1's
ones <- function(n) {
  rep(1, times = n)
}

# Return the maximum integer
max_integer <- function() {
  .Machine$integer.max
}

# Return the maximum double
max_double <- function() {
  .Machine$double.xmax
}

# Return the default scale parameter for 'prior_aux'
#
# @param basehaz Character string, the distribution for the baseline hazard.
get_default_aux_scale <- function(basehaz) {
  if (basehaz %in% c("weibull", "gompertz")) 2 else 20
}
