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

#' Bayesian survival models via Stan
#'
#' \if{html}{\figure{stanlogo.png}{options: width="25px" alt="http://mc-stan.org/about/logo/"}}
#' Bayesian inference for proportional or non-proportional hazards regression 
#' models. The user can specify a variety of standard parametric distributions 
#' for the baseline hazard, or a flexible parametric model (using either 
#' M-splines for modelling the baseline hazard, or B-splines for modelling 
#' the log baseline hazard). Covariate effects can be accommodated under
#' proportional hazards or non-proportional hazards (i.e. time-dependent 
#' effects).
#'
#' @export
#' @importFrom splines bs
#' 
#' @template args-prior_intercept
#' @template args-priors
#' @template args-prior_PD
#' @template args-algorithm
#' @template args-adapt_delta
#' 
#' @param formula A two-sided formula object describing the model. 
#'   The left hand side of the formula should be a \code{Surv()} 
#'   object. See \code{\link[survival]{Surv}}. If you wish to include
#'   time-dependent effect (i.e. time-dependent coefficients) in the model
#'   then the covariate(s) that you wish to estimate a time-dependent for 
#'   should be specified as \code{tde(varname)}  where \code{varname} is the 
#'   name of the covariate. See the \strong{Details} section for more 
#'   information on how the time-dependent effects are formulated, as well
#'   as the \strong{Examples} section.
#' @param data A data frame containing the variables specified in 
#'   \code{formula}.
#' @param basehaz A character string indicating which baseline hazard to use
#'   for the event submodel. Current options are: 
#'   \itemize{
#'     \item \code{"ms"}: a flexible parametric model using M-splines to 
#'     model the baseline hazard. The default locations for the internal knots, 
#'     as well as the basis terms for the splines, are calculated with respect
#'     to time. If the model does \emph{not} include any time-dependendent 
#'     effects then a closed form solution is available for both the hazard
#'     and cumulative hazard and so this approach should be relatively fast.
#'     On the other hand, if the model does include time-dependent effects then
#'     quadrature is used to evaluate the cumulative hazard at each MCMC
#'     iteration and, therefore, estimation of the model will be slower.
#'     \item \code{"bs"}: a flexible parametric model using B-splines to model
#'     the \emph{log} baseline hazard. The default locations for the internal 
#'     knots, as well as the basis terms for the splines, are calculated with 
#'     respect to time. A closed form solution for the cumulative hazard 
#'     is \strong{not} available (regardless of whether or not the model includes
#'     time-dependent effects) and therefore quadrature is used to evaluate 
#'     the cumulative hazard at each MCMC iteration. Therefore, if the model
#'     does not include any time-dependent effects, then estimation using the 
#'     \code{"ms"} baseline will be faster.
#'     \item \code{"exp"}: an exponential distribution for the event times. 
#'     (i.e. a constant baseline hazard)
#'     \item \code{"weibull"}: a Weibull distribution for the event times.
#'     \item \code{"gompertz"}: a Gompertz distribution for the event times.
#'   }
#'   Note that all spline-based models use splines of degree 3 (i.e. cubic
#'   splines).
#' @param basehaz_ops a named list specifying options related to the baseline
#'   hazard. Currently this can include: \cr
#'   \itemize{
#'     \item \code{df}: a positive integer specifying the degrees of freedom 
#'     for the M-splines / I-splines. The default is 6.
#'     \item \code{knots}: An optional numeric vector specifying the internal 
#'     knot locations for the splines if \code{basehaz = "ms"}. Knots cannot be
#'     specified if \code{df} is specified. If not specified, then the 
#'     default is to use \code{df - 4} knots, which are
#'     placed at equally spaced percentiles of the distribution of
#'     uncensored event times.
#'     \item \code{bknots}: an optional numeric vector specifying the boundary 
#'     knot locations for the splines if \code{basehaz = "ms"}. 
#'     If not specified, then the default is to place the boundary knots at the
#'     minimum and maximum of the event times (including both censored and
#'     uncensored events).
#'   }
#' @param qnodes The number of nodes to use for the Gauss-Kronrod quadrature
#'   that is used to evaluate the cumulative hazard when \code{basehaz = "bs"}
#'   or when time-dependent effects (i.e. non-proportional hazards) are 
#'   specified. Options are 15 (the default), 11 or 7.
#' @param prior_aux The prior distribution for the parameters related to the
#'   baseline hazard. The relevant "auxiliary" parameters differ depending on  
#'   on the type of baseline hazard specified in the \code{basehaz} 
#'   argument. The following applies:
#'   \itemize{
#'     \item \code{basehaz = "exp"}: there is \strong{no} auxiliary parameter, 
#'     since the log scale parameter for the exponential distribution is 
#'     incorporated as an intercept in the linear predictor.
#'     \item \code{basehaz = "weibull"}: the auxiliary parameter is the Weibull 
#'     shape parameter, while the log scale parameter for the Weibull 
#'     distribution is incorporated as an intercept in the linear predictor.
#'     The auxiliary parameter has a lower bound at zero. 
#'     \item \code{basehaz = "gompertz"}: the auxiliary parameter is the Gompertz 
#'     scale parameter, while the log shape parameter for the Gompertz 
#'     distribution is incorporated as an intercept in the linear predictor.
#'     The auxiliary parameter has a lower bound at zero. 
#'     \item \code{basehaz = "ms"}: the auxiliary parameters are the coefficients
#'     for the M-spline basis terms on the baseline hazard. These parameters
#'     have a lower bound at zero.
#'     \item \code{basehaz = "bs"}: the auxiliary parameters are the coefficients
#'     for the B-spline basis terms on the log baseline hazard. These parameters
#'     are unbounded.
#'   }
#'   Currently, \code{prior_aux} can be a call to \code{normal}, \code{student_t} 
#'   or \code{cauchy}. See \code{\link{priors}} for details on these functions. 
#'   To omit a prior ---i.e., to use a flat (improper) uniform prior--- set 
#'   \code{prior_aux} to \code{NULL}. 
#' @param prior_smooth This is only relevant when time-dependent effects are 
#'   specified in the model (i.e. the \code{tde()} function is used in the 
#'   model formula. When that is the case, \code{prior_smooth} determines the
#'   prior distribution for the hyperparameter of the smoothing function
#'   for the time-dependent coefficients (specifically the standard deviation 
#'   of the cubic B-spline coefficients). Lower values for the hyperparameter
#'   yield a less flexible smooth function. \code{prior_smooth} can be a call 
#'   to \code{exponential} to 
#'   use an exponential distribution, or \code{normal}, \code{student_t} or 
#'   \code{cauchy}, which results in a half-normal, half-t, or half-Cauchy 
#'   prior. See \code{\link{priors}} for details on these functions. To omit a 
#'   prior ---i.e., to use a flat (improper) uniform prior--- set 
#'   \code{prior_smooth} to \code{NULL}. The number of hyperparameters depends
#'   on the model specification but a scalar prior will be recylced as necessary
#'   to the appropriate length.
#'  
#' @details
#'   By default, any covariate effects specified in the \code{formula} are
#'   included in the model under a proportional hazards assumption. To relax
#'   this assumption, it is possible to estimate a time-dependent coefficient
#'   for a given covariate. This can be specified in the model \code{formula}
#'   by wrapping the covariate name in the \code{tde()} function (note that
#'   this function is not an exported function, rather it is an internal function
#'   that can only be evaluated within the formula of a \code{stan_surv} call).
#'   For example, if we wish to estimate a time-dependent effect for the 
#'   covariate \code{sex} then we can specify \code{tde(sex)} in the 
#'   \code{formula}, e.g. \code{Surv(time, status) ~ tde(sex) + age + trt}. 
#'   The coefficient for \code{sex} will then be modelled 
#'   using a flexible smooth function based on a cubic B-spline expansion of 
#'   time. The flexibility of the smooth function can be controlled through 
#'   the hyperparameters related the B-spline coefficients; see the 
#'   \code{prior_smooth} argument. Also, by default the cubic B-spline basis is 
#'   evaluated with 3 degrees of freedom (that is a cubic spline basis with  
#'   boundary knots at the limits of the time range, but no internal knots). If  
#'   you wish to increase the flexibility of the smooth function by using a 
#'   greater number of degrees of freedom, then you can specify this as part
#'   of the \code{tde} function call. For example, to use cubic B-splines with 
#'   7 degrees of freedom we could specify \code{tde(sex, df = 7)} in the
#'   model formula. See the \strong{Examples} section below for more details.
#'              
#' @examples
#' 
#' #---------- Proportional hazards
#' 
#' # Simulated data
#' library(simsurv)
#' covs <- data.frame(id  = 1:1000, 
#'                    trt = stats::rbinom(1000, 1L, 0.5))
#' dat1 <- simsurv(lambdas = 0.1, 
#'                 gammas  = 1.5, 
#'                 betas   = c(trt = -0.5),
#'                 x       = covs, 
#'                 maxt    = 5)
#' dat1 <- merge(dat1, covs)
#' fm1  <- Surv(eventtime, status) ~ trt
#' mod1a <- stan_surv(fm1, dat1, chains = 1, iter = 1000, basehaz = "ms")
#' mod1b <- stan_surv(fm1, dat1, chains = 1, iter = 1000, basehaz = "bs")
#' mod1c <- stan_surv(fm1, dat1, chains = 1, iter = 1000, basehaz = "exp")
#' mod1d <- stan_surv(fm1, dat1, chains = 1, iter = 1000, basehaz = "weibull")
#' mod1e <- stan_surv(fm1, dat1, chains = 1, iter = 1000, basehaz = "gompertz")
#' do.call(cbind, lapply(list(mod1a, mod1b, mod1c, mod1d, mod1e), fixef))
#' bayesplot::bayesplot_grid(plot(mod1a), plot(mod1b), 
#'                           plot(mod1c), plot(mod1d), 
#'                           ylim = c(0, 0.8))
#' 
#' # Breast cancer data
#' library(flexsurv)
#' dat2 <- flexsurv::bc
#' fm2  <- Surv(rectime, censrec) ~ group
#' mod2a <- stan_surv(fm2, dat2, chains = 1, iter = 1000)
#' mod2z <- flexsurv::flexsurvspline(fm2, dat2, k = 3)
#' print(mod2a, 4)
#' mod2z
#'                 
#' # PBC data
#' dat3 <- survival::pbc
#' dat3$timeYears <- dat3$time / 365.25
#' dat3$death     <- (dat3$status == 2)
#' fm3 <- Surv(timeYears, death) ~ sex + trt
#' mod3a <- stan_surv(fm3, dat3, chains = 1, iter = 1000)
#' mod3z <- flexsurv::flexsurvspline(fm3, dat3, k = 3)
#' print(mod3a, 4)
#' mod3z
#'
#' #---------- Non-proportional hazards
#' 
#' # Simulated data
#' library(simsurv)
#' library(rstpm2)
#' covs <- data.frame(id  = 1:1000, 
#'                    trt = stats::rbinom(1000, 1L, 0.5))
#' dat4 <- simsurv(lambdas = 0.1, 
#'                 gammas  = 1.5, 
#'                 betas   = c(trt = -0.5),
#'                 tde     = c(trt = 0.2),
#'                 x       = covs, 
#'                 maxt    = 5)
#' dat4 <- merge(dat4, covs)
#' fm4  <- Surv(eventtime, status) ~ tde(trt)
#' mod4a <- stan_surv(Surv(eventtime, status) ~ tde(trt), 
#'                    dat4, chains = 1, iter = 1000)
#' mod4z <- rstpm2::stpm2(Surv(eventtime, status) ~ trt, 
#'                        dat4, tvc = list(trt = 5))
#' print(mod4a, 4)
#' mod4z
#' plot(mod4a, "tde")
#' plot(mod4z, newdata = data.frame(trt = 0), type = "hr", var = "trt")
#' 
stan_surv <- function(formula, 
                      data, 
                      basehaz         = "ms", 
                      basehaz_ops, 
                      qnodes          = 15, 
                      prior           = normal(), 
                      prior_intercept = normal(),
                      prior_aux       = cauchy(), 
                      prior_smooth    = exponential(autoscale = FALSE), 
                      prior_PD        = FALSE,
                      algorithm       = c("sampling", "meanfield", "fullrank"),
                      adapt_delta     = 0.95, ...) {

  #-----------------------------
  # Pre-processing of arguments
  #-----------------------------
  
  if (!requireNamespace("survival"))
    stop("the 'survival' package must be installed to use this function.")
  
  if (missing(basehaz_ops)) 
    basehaz_ops <- NULL
  
  dots      <- list(...)
  algorithm <- match.arg(algorithm)
  
  formula   <- parse_formula(formula, data)
  data      <- make_model_data(formula$tf_form, data) # row subsetting etc.
  
  #----------------
  # Construct data
  #----------------
 
  #----- model frame stuff
  
  mf_stuff <- make_model_frame(formula$tf_form, data)
  
  mf <- mf_stuff$mf # model frame
  mt <- mf_stuff$mt # model terms
  
  #----- dimensions and response vectors

  # entry and exit times for each row of data
  t_end <- make_t(mf, type = "end") # end time
  t_beg <- make_t(mf, type = "beg") # beg time
  
  # event indicator for each row of data
  d <- make_d(mf)
  event <- as.logical(d == 1)

  # interval censoring indicator for each row of data
  interval <- as.logical(d == 3)
    
  # delayed entry indicator for each row of data
  delayed <- (!t_beg == 0)
  
  # time variables for stan
  t_events  <- t_end[event]
  t_censor  <- t_end[!event]
  t_delayed <- t_beg[delayed]

  # dimensions
  nevents  <- length(t_events)
  ncensor  <- length(t_censor)
  ndelayed <- length(t_delayed)
  
  #----- baseline hazard

  ok_basehaz <- c("exp", "weibull", "gompertz", "ms", "bs")
  ok_basehaz_ops <- get_ok_basehaz_ops(basehaz)
  basehaz <- handle_basehaz(basehaz = basehaz, 
                            basehaz_ops = basehaz_ops, 
                            ok_basehaz = ok_basehaz,
                            ok_basehaz_ops = ok_basehaz_ops,
                            times = t_end, status = event)
  nvars <- basehaz$nvars # number of basehaz aux parameters
  
  # flag if intercept is required for baseline hazard
  has_intercept   <- ai(has_intercept(basehaz))

  #----- define dimensions and times for quadrature

  # flag if formula uses time-dependent effects
  has_tde <- !is.null(formula$td_form)

  # flag if closed form available for cumulative baseline hazard
  has_closed_form <- check_for_closed_form(basehaz)

  # flag for quadrature
  has_quadrature <- has_tde || !has_closed_form
  
  if (has_quadrature) { # model uses quadrature
    
    # standardised weights and nodes for quadrature
    qq <- get_quadpoints(nodes = qnodes)
    qp <- qq$points
    qw <- qq$weights
    
    # quadrature points & weights, evaluated for each row of data
    qpts <- uapply(qp, unstandardise_qpts, t_beg, t_end) # qpts for exit time
    qwts <- uapply(qw, unstandardise_qwts, t_beg, t_end) # qwts for exit time

    # quadrature points & weights, evaluated for rows with delayed entry
    if (ndelayed) {
      qpts_delayed <- uapply(qp, unstandardise_qpts, 0, t_delayed) # qpts for entry time
      qwts_delayed <- uapply(qw, unstandardise_qwts, 0, t_delayed) # qwts for entry time
    } else {
      qpts_delayed <- rep(0,0)
      qwts_delayed <- rep(0,0)
    }
    
    # dimensions
    qrows    <- length(qpts)
    qdelayed <- length(qpts_delayed)
    
  } else { # model does not use quadrature
    
    qpts         <- rep(0,0) # dud entries for stan
    qwts         <- rep(0,0)
    qpts_delayed <- rep(0,0)
    qwts_delayed <- rep(0,0)
    qrows        <- 0L
    qdelayed     <- 0L
    
  }
  
  #----- basis terms for baseline hazard

  # basis terms at event times; used regardless of quadrature
  basis_events  <- make_basis(t_events, basehaz)
  ibasis_events <- make_basis(t_events, basehaz, integrate = TRUE)
  
  # basis terms at censoring times; used only without quadrature
  if (has_quadrature) {
    basis_censor  <- matrix(0,0,nvars) # dud entries for stan
    ibasis_censor <- matrix(0,0,nvars)
  } else {
    basis_censor  <- make_basis(t_censor, basehaz)
    ibasis_censor <- make_basis(t_censor, basehaz, integrate = TRUE)
  }

  # basis terms at delayed entry times; used only without quadrature
  if (has_quadrature) {
    basis_delayed  <- matrix(0,0,nvars) # dud entries for stan
    ibasis_delayed <- matrix(0,0,nvars)
  } else {
    basis_delayed  <- make_basis(t_delayed, basehaz)
    ibasis_delayed <- make_basis(t_delayed, basehaz, integrate = TRUE)
  }
  
  # basis terms at quadrature points; used only with quadrature
  if (has_quadrature) {
    basis_qpts         <- make_basis(qpts,         basehaz)
    basis_qpts_delayed <- make_basis(qpts_delayed, basehaz)
  } else {
    basis_qpts         <- matrix(0,0,nvars) # dud entries for stan
    basis_qpts_delayed <- matrix(0,0,nvars)
  }
    
  #----- predictor matrices
  
  # time-fixed predictor matrix
  x <- make_x(formula$tf_form, mf)$x
  x_events  <- keep_rows(x, event)
  x_censor  <- keep_rows(x, !event)
  x_delayed <- keep_rows(x, delayed)
  K <- ncol(x)
  if (has_quadrature) {
    x_qpts         <- rep_rows(x, times = qnodes)
    x_qpts_delayed <- rep_rows(x_delayed, times = qnodes)
  } else {
    x_qpts         <- matrix(0,0,K) # dud entries for stan
    x_qpts_delayed <- matrix(0,0,K)
  }
  
  # time-varying predictor matrix
  if (has_tde) { 
    tdfm  <- formula$td_form 
    xlevs <- .getXlevels(mt, mf)
    data_events <- keep_rows(data, event)
    data_qpts   <- rep_rows(data, times = qnodes)
    s_events    <- make_s(tdfm, data_events, times = t_events, xlevs = xlevs)
    s_qpts      <- make_s(tdfm, data_qpts,   times = qpts,     xlevs = xlevs)
    if (ndelayed) {
      data_qpts_delayed <- rep_rows(keep_rows(data, delayed), times = qnodes)
      s_qpts_delayed <- make_s(tdfm, data_qpts_delayed, qpts_delayed, xlevs = xlevs)
    } else {
      s_qpts_delayed <- matrix(0,0,ncol(s_events)) # dud entry for stan
    }
    smooth_map <- get_smooth_name(s_events, type = "smooth_map")
    S <- ncol(s_events) # num. of tde spline coefficients
  } else { # model does not have tde
    s_events       <- matrix(0,nevents ,0)
    s_qpts         <- matrix(0,qrows   ,0)
    s_qpts_delayed <- matrix(0,qdelayed,0)
    smooth_map     <- integer(0)
    S              <- 0L
  } 

  #----- stan data
  
  standata <- nlist(
    K, S,
    nevents,
    ncensor  = if (has_quadrature) 0L else ncensor,
    ndelayed = if (has_quadrature) 0L else ndelayed,
    qnodes,
    qrows,
    qdelayed,
    nvars,
    t_events,
    t_censor  = if (has_quadrature) rep(0,0) else t_censor,
    t_delayed = if (has_quadrature) rep(0,0) else t_delayed,
    x_events,
    x_censor  = if (has_quadrature) matrix(0,0,K) else x_censor,
    x_delayed = if (has_quadrature) matrix(0,0,K) else x_delayed,
    x_qpts,
    x_qpts_delayed,
    s_events,
    s_qpts,
    s_qpts_delayed,
    smooth_map,
    basis_events,
    basis_censor,
    basis_delayed,
    basis_qpts,
    basis_qpts_delayed,
    ibasis_events,
    ibasis_censor,
    ibasis_delayed,
    type = basehaz$type,
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
  ok_smooth_dists    <- c(ok_dists[1:3], "exponential")
  
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
                     default_scale = get_default_aux_scale(basehaz),
                     link = NULL,
                     ok_dists = ok_aux_dists)
  if (prior_PD && is.null(prior_aux))
    stop("'prior_aux' cannot be NULL if 'prior_PD' is TRUE")
  
  user_prior_smooth_stuff <- prior_smooth_stuff <-
    handle_glm_prior(prior_smooth, 
                     nvars = if (S) max(smooth_map) else 0,
                     default_scale = 1,
                     link = NULL,
                     ok_dists = ok_smooth_dists)
  
  # stop null priors if prior_PD is TRUE
  if (prior_PD) {
    if (is.null(prior))
      stop("'prior' cannot be NULL if 'prior_PD' is TRUE")
    if (is.null(prior_intercept) && has_intercept)
      stop("'prior_intercept' cannot be NULL if 'prior_PD' is TRUE")
    if (is.null(prior_aux))
      stop("'prior_aux' cannot be NULL if 'prior_PD' is TRUE")    
    if (is.null(prior_smooth) && (S > 0))
      stop("'prior_smooth' cannot be NULL if 'prior_PD' is TRUE")    
  }
  
  # autoscaling of priors
  prior_stuff           <- autoscale_prior(prior_stuff, predictors = x)
  prior_intercept_stuff <- autoscale_prior(prior_intercept_stuff)
  prior_aux_stuff       <- autoscale_prior(prior_aux_stuff)
  prior_smooth_stuff    <- autoscale_prior(prior_smooth_stuff)
  
  # priors
  standata$prior_dist              <- prior_stuff$prior_dist
  standata$prior_dist_for_intercept<- prior_intercept_stuff$prior_dist
  standata$prior_dist_for_aux      <- prior_aux_stuff$prior_dist
  standata$prior_dist_for_smooth   <- prior_smooth_stuff$prior_dist
  
  # hyperparameters
  standata$prior_mean               <- prior_stuff$prior_mean
  standata$prior_scale              <- prior_stuff$prior_scale
  standata$prior_df                 <- prior_stuff$prior_df
  standata$prior_mean_for_intercept <- c(prior_intercept_stuff$prior_mean)
  standata$prior_scale_for_intercept<- c(prior_intercept_stuff$prior_scale)
  standata$prior_df_for_intercept   <- c(prior_intercept_stuff$prior_df)
  standata$prior_scale_for_aux      <- prior_aux_stuff$prior_scale
  standata$prior_df_for_aux         <- prior_aux_stuff$prior_df
  standata$prior_mean_for_smooth    <- prior_smooth_stuff$prior_mean
  standata$prior_scale_for_smooth   <- prior_smooth_stuff$prior_scale
  standata$prior_df_for_smooth      <- prior_smooth_stuff$prior_df
  standata$global_prior_scale       <- prior_stuff$global_prior_scale
  standata$global_prior_df          <- prior_stuff$global_prior_df
  standata$slab_df                  <- prior_stuff$slab_df
  standata$slab_scale               <- prior_stuff$slab_scale

  # any additional flags
  standata$prior_PD <- ai(prior_PD)
  
  #---------------
  # Prior summary
  #---------------
  
  prior_info <- summarize_jm_prior(
    user_priorEvent           = user_prior_stuff,
    user_priorEvent_intercept = user_prior_intercept_stuff,
    user_priorEvent_aux       = user_prior_aux_stuff,
    adjusted_priorEvent_scale           = prior_stuff$prior_scale,
    adjusted_priorEvent_intercept_scale = prior_intercept_stuff$prior_scale,
    adjusted_priorEvent_aux_scale       = prior_aux_stuff$prior_scale,
    e_has_intercept  = has_intercept,
    e_has_predictors = K > 0,
    basehaz = basehaz
  )
  
  #-----------
  # Fit model
  #-----------

  # obtain stan model code
  stanfit  <- stanmodels$surv
  
  # specify parameters for stan to monitor
  stanpars <- c(if (standata$has_intercept) "gamma",
                if (standata$K)             "beta",
                if (standata$S)             "beta_tde",
                if (standata$S)             "smooth_sd",
                if (standata$nvars)         "coefs")
  
  # fit model using stan
  if (algorithm == "sampling") { # mcmc
    args <- set_sampling_args(
      object = stanfit, 
      data   = standata, 
      pars   = stanpars, 
      prior  = prior, 
      user_dots = list(...), 
      user_adapt_delta = adapt_delta, 
      show_messages = FALSE)
    stanfit <- do.call(rstan::sampling, args)
  } else { # meanfield or fullrank vb
    args <- nlist(
      object = stanfit,
      data   = standata,
      pars   = stanpars,
      algorithm
    )
    args[names(dots)] <- dots
    stanfit <- do.call(rstan::vb, args)
  }
  check_stanfit(stanfit)
  
  # define new parameter names
  nms_beta   <- colnames(x) # may be NULL
  nms_tde    <- get_smooth_name(s_events, type = "smooth_coefs") # may be NULL
  nms_smooth <- get_smooth_name(s_events, type = "smooth_sd")    # may be NULL
  nms_int    <- get_int_name_basehaz(basehaz)
  nms_aux    <- get_aux_name_basehaz(basehaz)
  nms_all    <- c(nms_int,
                  nms_beta,
                  nms_tde,
                  nms_smooth,
                  nms_aux,
                  "log-posterior")

  # substitute new parameter names into 'stanfit' object
  stanfit <- replace_stanfit_nms(stanfit, nms_all)

  # return an object of class 'stansurv'
  fit <- nlist(stanfit, 
               formula,
               has_tde,
               has_quadrature,
               model_data       = data,
               model_frame      = mf,
               terms            = mt,
               xlevels          = .getXlevels(mt, mf),
               x,
               s_events         = if (has_tde) s_events else NULL,
               t_beg, 
               t_end, 
               event, 
               delayed,
               basehaz,
               nobs             = nrow(mf),
               nevents          = length(t_events),
               ncensor          = length(t_censor),
               ndelayed         = length(t_delayed),
               prior_info,
               qnodes           = if (has_quadrature) qnodes else NULL,
               algorithm,
               stan_function    = "stan_surv",
               rstanarm_version = packageVersion("rstanarm"),
               call             = match.call(expand.dots = TRUE))
  stansurv(fit)
}


#---------- internal

# Identify whether the type of baseline hazard requires an intercept in
# the linear predictor (NB splines incorporate the intercept into the basis).
#
# @param basehaz A list with info about the baseline hazard; see 'handle_basehaz'.
# @return A Logical.
has_intercept <- function(basehaz) {
  nm <- get_basehaz_name(basehaz)
  (nm %in% c("exp", "weibull", "gompertz"))
}

# Return the name of the tde spline coefs or smoothing parameters.
#
# @param x The predictor matrix for the time-dependent effects, with column names.
# @param type The type of information about the smoothing parameters to return.
# @return A character or numeric vector, depending on 'type'.
get_smooth_name <- function(x, type = "smooth_coefs") {
  
  if (is.null(x) || !ncol(x))
    return(NULL)  

  nms <- gsub(":bs\\(times__.*\\)[0-9]*$", "", colnames(x))
  tally   <- table(nms)
  indices <- uapply(tally, seq_len)
  suffix  <- paste0(":tde-spline-coef", indices)
  
  switch(type,
         smooth_coefs = paste0(nms, suffix),
         smooth_sd    = paste0("smooth_sd[", unique(nms), "]"),
         smooth_map   = rep(seq_along(tally), tally),
         smooth_vars  = unique(nms),
         stop2("Bug found: invalid input to 'type' argument."))
}

# Return the default scale parameter for 'prior_aux'.
#
# @param basehaz A list with info about the baseline hazard; see 'handle_basehaz'.
# @return A scalar.
get_default_aux_scale <- function(basehaz) {
  nm <- get_basehaz_name(basehaz)
  if (nm %in% c("weibull", "gompertz")) 2 else 20
}

# Extract the name of baseline hazard from the list returned by 'handle_basehaz'.
#
# @param basehaz A list with info about the baseline hazard; see 'handle_basehaz'.
# @return A character string.
get_basehaz_name <- function(basehaz) {
  if (is.character(basehaz)) basehaz else basehaz$type_name
}

# Check if the type of baseline hazard has a closed form
#
# @param basehaz A list with info about the baseline hazard; see 'handle_basehaz'.
# @return A logical.
check_for_closed_form <- function(basehaz) {
  nm <- get_basehaz_name(basehaz)
  nm %in% c("exp",
            "weibull",
            "gompertz",
            "ms")
}

# Replace the parameter names slot of an object of class 'stanfit'.
#
# @param stanfit An object of class 'stanfit'.
# @param new_nms A character vector of new parameter names.
# @return A 'stanfit' object.
replace_stanfit_nms <- function(stanfit, new_nms) {
  stanfit@sim$fnames_oi <- new_nms
  stanfit
}

# Return the spline basis for the given type of baseline hazard.
# 
# @param times A numeric vector of times at which to evaluate the basis.
# @param basehaz A list with info about the baseline hazard, returned by a 
#   call to 'handle_basehaz'.
# @param integrate A logical, specifying whether to calculate the integral of
#   the specified basis.
# @return A matrix.
make_basis <- function(times, basehaz, integrate = FALSE) {
  N <- length(times)
  K <- basehaz$nvars
  if (!N) { # times is NULL or empty vector
    return(matrix(0, 0, K))
  } 
  switch(basehaz$type_name,
         "exp"       = matrix(0, N, K), # dud matrix for Stan
         "weibull"   = matrix(0, N, K), # dud matrix for Stan
         "gompertz"  = matrix(0, N, K), # dud matrix for Stan
         "ms"        = basis_matrix(times, basis = basehaz$basis, integrate = integrate),
         "bs"        = basis_matrix(times, basis = basehaz$basis),
         "piecewise" = dummy_matrix(times, knots = basehaz$knots),
         stop2("Bug found: type of baseline hazard unknown."))
}

# Evaluate a spline basis matrix at the specified times
#
# @param time A numeric vector.
# @param basis Info on the spline basis.
# @param integrate A logical, should the integral of the basis be returned?
# @return A two-dimensional array.
basis_matrix <- function(times, basis, integrate = FALSE) {
  out <- predict(basis, times)
  if (integrate) {
    stopifnot(inherits(basis, "mSpline"))
    out <- splines2:::predict.iSpline(basis, times)
  }
  aa(out)
}

# Parse the model formula
#
# @param formula The user input to the formula argument.
# @param data The user input to the data argument (i.e. a data frame).
parse_formula <- function(formula, data) {
  
  formula <- validate_formula(formula, needs_response = TRUE)
  
  lhs      <- lhs(formula) # full LHS of formula
  lhs_form <- reformulate_lhs(lhs)
  
  rhs        <- rhs(formula)         # RHS as expression
  rhs_form   <- reformulate_rhs(rhs) # RHS as formula
  rhs_terms  <- terms(rhs_form, specials = "tde")
  rhs_vars   <- rownames(attr(rhs_terms, "factors"))
  
  allvars <- all.vars(formula)
  allvars_form <- reformulate(allvars)
  
  surv <- eval(lhs, envir = data) # Surv object
  surv <- validate_surv(surv)
  type <- attr(surv, "type")
  
  if (type == "right") {
    tvar_beg <- NULL
    tvar_end <- as.character(lhs[[2L]])
    dvar     <- as.character(lhs[[3L]])
    min_t    <- 0
    max_t    <- max(surv[, "time"])
  } else if (type == "counting") {
    tvar_beg <- as.character(lhs[[2L]])
    tvar_end <- as.character(lhs[[3L]])
    dvar     <- as.character(lhs[[4L]])
    min_t    <- min(surv[, "start"])
    max_t    <- max(surv[, "stop"])
  }

  sel <- attr(rhs_terms, "specials")$tde
  
  if (!is.null(sel)) { # model has tde
    
    # replace 'tde(x, ...)' in formula with 'x'
    tde_oldvars <- rhs_vars
    tde_newvars <- sapply(tde_oldvars, function(oldvar) {
      if (oldvar %in% rhs_vars[sel]) {
        tde <- function(newvar, ...) { # define tde function locally
          safe_deparse(substitute(newvar)) 
        }
        eval(parse(text = oldvar))
      } else oldvar
    }, USE.NAMES = FALSE)
    term_labels <- attr(rhs_terms, "term.labels")
    for (i in sel) {
      sel_terms <- which(attr(rhs_terms, "factors")[i, ] > 0)
      for (j in sel_terms) {
        term_labels[j] <- gsub(tde_oldvars[i], 
                               tde_newvars[i], 
                               term_labels[j], 
                               fixed = TRUE)
      }
    }
    tf_form <- reformulate(term_labels, response = lhs)
    
    # extract 'tde(x, ...)' from formula and construct 'bs(times, ...)'
    tde_terms <- lapply(rhs_vars[sel], function(x) {
      tde <- function(vn, ...) { # define tde function locally
        dots <- list(...)
        ok_args <- c("df")
        if (!isTRUE(all(names(dots) %in% ok_args)))
          stop2("Invalid argument to 'tde' function. ",
                "Valid arguments are: ", comma(ok_args))
        df <- if (is.null(dots$df)) 3 else dots$df
        degree <- 3
        if (df == 3) {
          dots[["knots"]] <- numeric(0)
        } else {
          dx <- (max_t - min_t) / (df - degree + 1)
          dots[["knots"]] <- seq(min_t + dx, max_t - dx, dx)
        }
        dots[["Boundary.knots"]] <- c(min_t, max_t) 
        sub("^list\\(", "bs\\(times__, ", deparse(dots))
      }
      tde_calls <- eval(parse(text = x))
      sel_terms <- which(attr(rhs_terms, "factors")[x, ] > 0)
      new_calls <- sapply(seq_along(sel_terms), function(j) {
        paste0(term_labels[sel_terms[j]], ":", tde_calls)
      })
      nlist(tde_calls, new_calls)
    })
    td_basis <- fetch(tde_terms, "tde_calls")
    new_calls <- fetch_(tde_terms, "new_calls")
    td_form <- reformulate(new_calls, response = NULL, intercept = FALSE)
    
  } else { # model doesn't have tde
    tf_form  <- formula
    td_form  <- NULL
    td_basis <- NULL
  }

  nlist(formula,
        lhs,
        rhs,
        lhs_form,
        rhs_form,
        tf_form,
        td_form,
        td_basis,
        fe_form = rhs_form, # no re terms accommodated yet
        re_form = NULL,     # no re terms accommodated yet
        allvars,
        allvars_form,
        tvar_beg,
        tvar_end,
        dvar,
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


# Return the response vector (time) for estimation
#
# @param model_frame The model frame.
# @param type The type of time variable to return:
#   "beg": the entry time for the row in the survival data,
#   "end": the exit  time for the row in the survival data,
#   "gap": the difference between entry and exit times,
#   "upp": if the row involved interval censoring, then the exit time
#          would have been the lower limit of the interval, and "upp" 
#          is the upper limit of the interval.
# @return A numeric vector.
make_t <- function(model_frame, type = c("beg", "end", "gap", "upp")) {
  
  type <- match.arg(type)
  resp <- model.response(model_frame)
  surv <- attr(resp, "type")
  err  <- paste0("Bug found: cannot handle '", surv, "' Surv objects.")
  
  t_beg <- switch(surv,
                  "right"     = rep(0, nrow(model_frame)),
                  "interval"  = rep(0, nrow(model_frame)),
                  "interval2" = rep(0, nrow(model_frame)),
                  "counting"  = as.vector(resp[, "start"]),
                  stop(err))
  
  t_end <- switch(surv,
                  "right"     = as.vector(resp[, "time"]),
                  "interval"  = as.vector(resp[, "time1"]),
                  "interval2" = as.vector(resp[, "time1"]),
                  "counting"  = as.vector(resp[, "stop"]),
                  stop(err))

  t_upp <- switch(surv,
                  "right"     = rep(NA, nrow(model_frame)),
                  "counting"  = rep(NA, nrow(model_frame)),
                  "interval"  = as.vector(resp[, "time2"]),
                  "interval2" = as.vector(resp[, "time2"]),
                  stop(err))
    
  switch(type,
         "beg" = t_beg,
         "end" = t_end,
         "gap" = t_end - t_beg,
         "upp" = t_upp,
         stop("Bug found: cannot handle specified 'type'."))
}


# Return the response vector (status indicator)
#
# @param model_frame The model frame.
# @return A numeric vector.
make_d <- function(model_frame) {
  
  resp <- model.response(model_frame)
  surv <- attr(resp, "type")
  err  <- paste0("Bug found: cannot handle '", surv, "' Surv objects.")
  
  switch(surv,
         "right"     = as.vector(resp[, "status"]),
         "interval"  = as.vector(resp[, "status"]),
         "interval2" = as.vector(resp[, "status"]),
         "counting"  = as.vector(resp[, "status"]),
         stop(err))
}

# Return a data frame with NAs excluded
#
# @param formula The parsed model formula.
# @param data The user specified data frame.
make_model_data <- function(formula, data) {
  mf <- model.frame(formula, data, na.action = na.pass)
  include <- apply(mf, 1L, function(row) !any(is.na(row)))
  data[include, , drop = FALSE]
}

# Return the model frame
#
# @param formula The parsed model formula.
# @param data The model data frame.
make_model_frame <- function(formula, data, times = NULL) {

  # construct terms object from formula 
  Terms <- terms(formula)
  
  # construct model frame
  mf <- model.frame(Terms, data)
  mf <- check_constant_vars(mf)
  
  # check terms
  mt <- attr(mf, "terms")
  if (is.empty.model(mt)) 
    stop2("No intercept or predictors specified.")
  
  nlist(mf, mt)
}

# Return the fe predictor matrix for estimation
#
# @param formula The parsed model formula.
# @param model_frame The model frame.
# @return A named list with the following elements:
#   x: the fe model matrix, not centred and without intercept.
#   xbar: the column means of the model matrix.
#   N,K: number of rows (observations) and columns (predictors) in the
#     fixed effects model matrix
make_x <- function(formula, model_frame, xlevs = NULL) {

  # uncentred predictor matrix, without intercept
  x <- model.matrix(formula, model_frame, xlevs = xlevs)
  x <- drop_intercept(x)
  
  # column means of predictor matrix
  xbar <- colMeans(x)
  
  # identify any column of x with < 2 unique values (empty interaction levels)
  sel <- (apply(x, 2L, n_distinct) < 2)
  if (any(sel)) {
    cols <- paste(colnames(x)[sel], collapse = ", ")
    stop2("Cannot deal with empty interaction levels found in columns: ", cols)
  }
  
  nlist(x, xbar, N = NROW(x), K = NCOL(x))
}

# Return a predictor for the tde spline terms
#
# @param formula The formula for the time-dependent effects part of the model.
# @param data A data frame.
# @param times The vector of times at which the predictor matrix should be 
#   evaluated.
# @param xlevs The factor levels to use for the predictor matrix.
# @return A matrix.
make_s <- function(formula, data, times, xlevs = NULL) {
  
  # add times (as a new variable) to the model data
  if (!length(times) == nrow(data))
    stop("Bug found: 'times' is the incorrect length.")
  data <- data.frame(data, times__ = times)

  # make model frame and predictor matrix
  mf <- make_model_frame(formula, data, times)$mf
  x  <- make_x(formula, mf, xlevs = xlevs)$x
  return(x)
}

# Return the fe predictor matrix for prediction
#
# @param object A stansurv object.
# @param model_frame The model frame.
# @return A named list with the following elements:
#   x: the fe model matrix, not centred and may have intercept depending on
#     the requirement of the baseline hazard.
#   N,K: number of rows (observations) and columns (predictors) in the
#     fixed effects model matrix
make_pp_x <- function(object, model_frame) {
  
  # formula for fe predictor matrix
  tt <- delete.response(terms(object))
  
  # check data classes in the model frame match those used in model fitting
  if (!is.null(cl <- attr(tt, "dataClasses"))) 
    .checkMFClasses(cl, model_frame)
  
  # uncentered predictor matrix
  x <- model.matrix(tt, model_frame, contrasts.arg = object$contrasts)
  
  # drop intercept if baseline hazard doesn't require one
  if (!has_intercept(object$basehaz))
    x <- drop_intercept(x)
  
  nlist(x, N = NROW(x), K = NCOL(x))  
}

# apply b-spline time-dependent effect
apply_tde_fun <- function(model_terms, model_frame, times, bknots = NULL) {
  
  tde_stuff <- survival::untangle.specials(model_terms, "tde")

  if (!length(tde_stuff$terms)) 
    return(model_frame) # no time-dependent effects
  
  if (!nrow(model_frame))
    return(model_frame) # no rows in model frame (e.g. no delayed entry)
  
  vars  <- attr(model_terms, 'variables')
  pvars <- attr(model_terms, 'predvars')
  
  # loop over time-dependent terms in formula
  K <- length(tde_stuff$terms)
  for (i in 1:K) { 
    indx_i <- tde_stuff$terms[i] + 2 # index in call; +2 for 'list' & 'Surv()'
    var_i  <- vars [[indx_i]]        # var     in formula
    pvar_i <- pvars[[indx_i]]        # predvar in formula
    var_i  <- safe_deparse(var_i)    # treat call as a string
    pvar_i <- safe_deparse(pvar_i)   # treat call as a string
    # get the possible prefixes for the predvar (i.e. 'tde(x' or 'bs(x')
    prefix <- "^bs\\([^,]+,[[:blank:]]*|^tde\\([^,]+,[[:blank:]]*"
    # returns dots from 'tde(x, ...)' as a list
    chck <- grepl(prefix, pvar_i)
    if (chck) {
      args_i <- eval_string(sub(prefix, "list\\(", pvar_i)) 
    } else {
      args_i <- list()
    }
    # combine the dots with the times at which to evaluate the b-spline basis
    args_i$intercept <- TRUE
    if (!is.null(bknots))
      args_i$Boundary.knots <- bknots
    args_i <- c(list(x = times), args_i)
    # extract the variable from the model frame
    oldx_i  <- model_frame[[var_i]]
    # apply interaction with the b-spline basis evaluated at specified times
    newx_i <- oldx_i * do.call(splines::bs, args_i)
    # substitute back into the model frame
    model_frame[[var_i]] <- newx_i
  }
  
  return(model_frame)
}

update_tde_terms <- function(model_terms, model_frame) {
  tde_terms <- survival::untangle.specials(model_terms, "tde")$terms
  if (!length(tde_terms))
    return(model_frame) # no time-dependent effects
  vars  <- attr(model_terms, 'variables')
  pvars <- attr(model_terms, 'predvars')
  dclss <- attr(model_terms, "dataClasses")
  K <- length(tde_terms)
  for (i in 1:K) {
    indx_i <- tde_terms[i] + 2       # index in call; +2 for 'list' & 'Surv()'
    var_i  <- vars [[indx_i]]        # var     in formula
    pvar_i <- pvars[[indx_i]]        # predvar in formula
    var_i  <- safe_deparse(var_i)    # treat call as a string
    pvar_i <- safe_deparse(pvar_i)   # treat call as a string
    oldx_i <- model_frame[[var_i]]   # extract transformed variable from model frame
    dummy  <- as.call(list(as.name(class(oldx_i)[[1L]]), vars[[indx_i]][[2]]))
    ptemp  <- makepredictcall(oldx_i, dummy) # predvars call
    pvars[[indx_i]] <- ptemp
    dclss[[var_i]] <- class(oldx_i)[[1L]]
  }
  attr(model_terms, "predvars") <- pvars
  #attr(model_terms, "dataClasses") <- dclss
  return(model_terms)
}


#--------- not used; based on tt approach instead of tde approach

# # Validate the user input to the 'tt' argument. This draws on the 
# # code for the coxph modelling function in the survival package.
# #
# # Copyright (C) 2018 Sam Brilleman
# # Copyright (C) 2018 Terry Therneau, Thomas Lumley
# #
# # @param tt The user input to the 'tt' argument.
# # @param validate_length Integer specifying the required length of the 
# #   returned list.
# # @return A list of functions.
# validate_tt_fun <- function(tt, validate_length) {
#   
#   if (is.null(tt))
#     stop2("'tt' must be specified.")
#   
#   if (is.function(tt)) 
#     tt <- list(tt) # convert since function to a one element list
#   
#   if (!is.list(tt))
#     stop2("The 'tt' argument must contain a function or list of functions.")  
#   
#   if (!all(sapply(tt, is.function)))
#     stop2("The 'tt' argument must contain function or list of functions.")
#   
#   if (!length(tt) %in% c(1, validate_length)) 
#     stop2("The 'tt' argument contains a list of the incorrect length.")
#   
#   if (length(tt) == 1)
#     tt <- rep(tt, validate_length)
#   
#   return(tt)
# }
# 
# # apply time transform to the model frame; method based on survival package 
# apply_tt_fun <- function(model_frame, tt_funs, tt_vars, tt_terms, times) {
#   if (!length(tt_terms))
#     return(model_frame)
#   
#   for (i in 1:length(tt_terms)) { # loop over time transform terms
#     
#     # extract quantities used in time transform
#     varnm_i <- tt_vars[[i]] # varname in model frame
#     ttfun_i <- tt_funs[[i]] # user defined tt function
#     
#     # time transform at event times
#     oldx_i <- model_frame[[varnm_i]]   # extract var from model frame
#     newx_i <- (ttfun_i)(oldx_i, times) # evaluate tt function at times
#     model_frame[[varnm_i]] <- newx_i   # substitute back into model frame
#   }
#   
#   return(model_frame)
# }
#
# # update the predvars attribute for time transformed terms
# update_predvars <- function(model_terms, model_frame, tt_vars, tt_terms) {
#   tcall <- attr(model_terms, 'variables')[tt_terms + 2]
#   pvars <- attr(model_terms, 'predvars')
#   pmethod <- sub("makepredictcall.", "", as.vector(methods("makepredictcall")))
#   for (i in 1:length(tt_terms)) {
#     # update predvars if necessary
#     varnm_i <- tt_vars[[i]]       # varname in model frame
#     terms_i <- tt_terms[i] + 2    # index in terms object
#     x_i <- model_frame[[varnm_i]] # extract transformed variable from model frame
#     nclass <- class(x_i)          # check class of transformed variable
#     if (any(nclass %in% pmethod)) { # it has a makepredictcall method...
#       dummy <- as.call(list(as.name(class(x_i)[1]), tcall[[i]][[2]]))
#       ptemp <- makepredictcall(x_i, dummy)
#       pvars[[terms_i]] <- ptemp
#     }
#   }
#   attr(model_terms, "predvars") <- pvars
#   return(model_terms)
# }

