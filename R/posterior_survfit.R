# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2015, 2016, 2017 Trustees of Columbia University
# Copyright (C) 2016, 2017 Sam Brilleman
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

#' Posterior predictions for survival models
#' 
#' This function allows us to generate predicted quantities for survival
#' models at specified times. These quantities include the hazard rate, 
#' cumulative hazard, survival probability, or failure probability (i.e. CDF).
#' Note that the cumulative hazard, survival probability, or failure 
#' probability may be conditional on a last known survival time (see the 
#' \code{condition} argument discussed below). Predictions are obtained 
#' using unique draws from the posterior distribution of each of the model 
#' parameters and then summarised into a median and posterior uncertainty 
#' interval. For \code{stan_jm} models "dynamic" predictions are allowed and
#' are in fact the default when new data is provided (see the \code{dynamic}
#' argument discussed below).
#' 
#' 
#' @export
#' @import splines2
#' 
#' @templateVar stanregArg object
#' @template args-stansurv-stanjm-object
#' 
#' @param newdata Optionally, a data frame in which to look for variables with
#'   which to predict. If omitted, the model matrix is used. If \code{newdata}
#'   is provided and any variables were transformed (e.g. rescaled) in the data
#'   used to fit the model, then these variables must also be transformed in
#'   \code{newdata}. This only applies if variables were transformed before
#'   passing the data to one of the modeling functions and \emph{not} if
#'   transformations were specified inside the model formula. Also, 
#'   \code{newdata} can optionally include a variable with information 
#'   about the last known survival time for the new individuals -- 
#'   see the description for the \code{last_time} argument below
#'   -- however also note that when generating the survival probabilities it 
#'   is of course assumed that all individuals in \code{newdata} have not 
#'   yet experienced the event (that is, any variable in \code{newdataEvent} 
#'   that corresponds to the event indicator will be ignored).
#' @param newdataLong,newdataEvent An optional data frame (or in the case of 
#'   \code{newdataLong} this can be a list of data frames) in which to look 
#'   for variables with which to predict. If omitted, the model matrices are 
#'   used. If new data is provided, then it should also contain the longitudinal 
#'   outcome data on which to condition when drawing the new group-specific 
#'   coefficients for individuals in the new data unless the \code{dynamic}
#'   argument is set to \code{FALSE}. Note that there is only
#'   allowed to be one row of data for each individual in \code{newdataEvent}, 
#'   that is, time-varying covariates are not allowed in the prediction data for
#'   the event submodel. Also, \code{newdataEvent} can optionally include a 
#'   variable with information about the last known survival time for the new
#'   individuals -- see the description for the \code{last_time} argument below
#'   -- however also note that when generating the survival probabilities it 
#'   is of course assumed that all individuals in \code{newdataEvent} have not 
#'   yet experienced the event (that is, any variable in \code{newdataEvent} 
#'   that corresponds to the event indicator will be ignored).
#' @param type The type of prediction to return. The following are currently
#'   allowed:
#'   \itemize{
#'     \item \code{"surv"}: the estimated survival probability.
#'     \item \code{"cumhaz"}: the estimated cumulative hazard.
#'     \item \code{"haz"}: the estimated hazard rate.
#'     \item \code{"cdf"}: the estimated failure probability.
#'     \item \code{"logsurv"}: the estimated log survival probability.
#'     \item \code{"logcumhaz"}: the estimated log cumulative hazard.
#'     \item \code{"loghaz"}: the estimated log hazard rate.
#'     \item \code{"logcdf"}: the estimated log failure probability.
#'   }
#' @param extrapolate A logical specifying whether to extrapolate the estimated 
#'   survival probabilities beyond the times specified in the \code{times} argument.
#'   If \code{TRUE} then the extrapolation can be further controlled using
#'   the \code{control} argument.
#' @param control A named list with parameters controlling extrapolation 
#'   of the estimated survival function when \code{extrapolate = TRUE}. The list
#'   can contain one or more of the following named elements: \cr
#'   \itemize{
#'     \item \code{epoints}: a positive integer specifying the number of 
#'     discrete time points at which to calculate the forecasted survival 
#'     probabilities. The default is 100.
#'     \item \code{edist}: a positive scalar specifying the amount of time 
#'     across which to forecast the estimated survival function, represented 
#'     in the same units of time as were used for the event times in the fitted 
#'     model. The default is to extrapolate between the times specified in the 
#'     \code{times} argument and the maximum event or censoring time found in 
#'     the original data used to fit the model. If \code{edist} leads to times 
#'     that are beyond the maximum event or censoring time in the original data 
#'     then the estimated survival probabilities will be truncated at that 
#'     point, since an estimate for the baseline hazard is not available 
#'     beyond that time.
#'   }
#' @param condition A logical specifying whether the estimated 
#'     subject-specific survival probabilities at time \code{t} should be 
#'     conditioned on survival up to a fixed time point \code{u}. The default 
#'     is for \code{condition} to be set to \code{TRUE}, unless standardised survival
#'     probabilities have been requested (by specifying \code{standardise = TRUE}),
#'     in which case \code{condition} must (and will) be set to \code{FALSE}.
#'     When conditional survival probabilities are requested, the fixed
#'     time point \code{u} will be either: (i) the value specified via the 
#'     \code{last_time} argument; or if the \code{last_time} argument is 
#'     \code{NULL} then the latest observation time for each individual 
#'     (taken to be the value in the \code{times} argument if \code{newdataEvent} 
#'     is specified, or the observed event or censoring time if \code{newdataEvent} 
#'     is \code{NULL}.
#' @param last_time A scalar, character string, or \code{NULL}. This argument 
#'     specifies the last known survival time for each individual when
#'     conditional predictions are being obtained. If 
#'     \code{newdataEvent} is provided and conditional survival predictions are being
#'     obtained, then the \code{last_time} argument can be one of the following:
#'     (i) a scalar, this will use the same last time for each individual in 
#'     \code{newdataEvent}; (ii) a character string, naming a column in 
#'     \code{newdataEvent} in which to look for the last time for each individual;
#'     (iii) \code{NULL}, in which case the default is to use the time of the latest 
#'     longitudinal observation in \code{newdataLong}. If \code{newdataEvent} is
#'     \code{NULL} then the \code{last_time} argument cannot be specified 
#'     directly; instead it will be set equal to the event or censoring time for
#'     each individual in the dataset that was used to estimate the model. 
#'     If standardised survival probabilities are requested (i.e. 
#'     \code{standardise = TRUE}) then conditional survival probabilities are
#'     not allowed and therefore the \code{last_time} argument is ignored.
#' @param ids For \code{stan_jm} models. An optional vector specifying 
#'   a subset of IDs for whom the predictions should be obtained. 
#'   The default is to predict for all individuals
#'   who were used in estimating the model or, if \code{newdataLong} and 
#'   \code{newdataEvent} are specified, then all individuals contained in 
#'   the new data.
#' @param prob A scalar between 0 and 1 specifying the width to use for the 
#'   uncertainty interval (sometimes called credible interval) for the predictions. 
#'   For example \code{prob = 0.95} (the default) means that the 2.5th and 97.5th  
#'   percentiles will be provided.
#' @param times A scalar, a character string, or \code{NULL}. Specifies the  
#'   times at which the estimated survival probabilities should be calculated. 
#'   It can be either: (i) \code{NULL}, in which case it will default to the last known 
#'   survival time for each individual, as determined by the \code{last_time}
#'   argument; (ii) a scalar, specifying a time to estimate the survival probability
#'   for each of the individuals; or (iii) if \code{newdataEvent} is  
#'   provided, it can be the name of a variable in \code{newdataEvent} that 
#'   indicates the time at which the survival probabilities should be calculated  
#'   for each individual. 
#' @param standardise A logical specifying whether the estimated 
#'   subject-specific survival probabilities should be averaged
#'   across all individuals for whom the subject-specific predictions are 
#'   being obtained. This can be used to average over the covariate and random effects
#'   distributions of the individuals used in estimating the model, or the individuals 
#'   included in the \code{newdata} arguments. This approach of
#'   averaging across the observed distribution of the covariates is sometimes
#'   referred to as a "standardised" survival curve. If \code{standardise = TRUE}, 
#'   then the \code{times} argument must be specified and it must be constant across 
#'   individuals, that is, the survival probabilities must be calculated at the 
#'   same time for all individuals.
#' @param dynamic A logical that is only relevant for \code{stan_jm} models 
#'   when new data is provided via the \code{newdataLong} and \code{newdataEvent} 
#'   arguments. If \code{dynamic = TRUE}, then new group-specific parameters are 
#'   drawn for the individuals in the new data, conditional on their longitudinal 
#'   biomarker data contained in \code{newdataLong}. These group-specific
#'   parameters are then used to generate individual-specific survival probabilities
#'   for these individuals. These are often referred to as "dynamic predictions"
#'   in the joint modelling context, because the predictions can be updated
#'   each time additional longitudinal biomarker data is collected on the individual.
#'   On the other hand, if \code{dynamic = FALSE} then the survival probabilities
#'   will be obtained by marginalising over the distribution of the group-specific
#'   coefficients; this has the benefit that the user does not need to provide
#'   longitudinal outcome data for the new individuals, but it will also 
#'   mean that the survival predictions will incorporate all uncertainty due 
#'   to between-individual variation in the longitudinal trajectories and so 
#'   there is likely to be very wide credible intervals on the predicted 
#'   survival probabilities.
#' @param scale Only relevant for \code{stan_jm} models when new data
#'   is supplied and \code{dynamic = TRUE}, in which case new random effects 
#'   are simulated for the individuals in the new data using a 
#'   Metropolis-Hastings algorithm. The \code{scale} argument should be a 
#'   scalar. It specifies how much to multiply the asymptotic 
#'   variance-covariance matrix for the random effects by, which is then
#'   used as the "width" (ie. variance-covariance matrix) of the multivariate
#'   Student-t proposal distribution in the Metropolis-Hastings algorithm.
#' @param draws An integer specifying the number of MCMC draws to use when
#'   evaluating the predicted quantities. For \code{stan_surv} models, the
#'   default number of draws is the size of the posterior sample.
#'   For \code{stan_jm} models, the default number of draws is 200 (or the 
#'   size of the posterior sample if that is less than 200). The smaller
#'   default number of draws for \code{stan_jm} models is because dynamic
#'   predictions (when \code{dynamic = TRUE}) can be slow.
#' @param seed An optional \code{\link[=set.seed]{seed}} to use.
#' @param return_matrix A logical. If \code{TRUE} then a list of \code{draws} by 
#'   \code{nrow(newdata)} matrices is returned. Each matrix contains the actual
#'   simulations or draws from the posterior predictive distribution. Otherwise
#'   if \code{return_matrix} is set to \code{FALSE} (the default) then a 
#'   data frame is returned. See the \strong{Value} section below for more 
#'   detail.
#' @param ... Currently unused.
#' 
#' @details
#'   By default, the predicted quantities are evaluated conditional on observed 
#'   values of the fixed effect covariates. That is, predictions will be 
#'   obtained using either: 
#'   \itemize{
#'     \item the design matrices used in the original \code{\link{stan_surv}}
#'     or \code{\link{stan_jm}} model call, or
#'     \item the covariate values provided in the \code{newdata} argument
#'     (or \code{newdataLong} and \code{newdataEvent} arugments for the
#'     \code{stanjm} method). 
#'   }
#'   However, if you wish to average over the observed distribution 
#'   of the fixed effect covariates then this is possible -- such predictions
#'   are sometimes referred to as standardised survival probabilties -- see the 
#'   \code{standardise} argument.
#'   
#'   For \code{stansurv} objects, the predicted quantities are calculated for  
#'   \emph{each row of the prediction data}, at the specified \code{times} as
#'   well as any times generated through extrapolation (when 
#'   \code{extrapolate = TRUE}).
#'   
#'   For \code{stanjm} objects, the predicted quantities are calculated for 
#'   \emph{each individual}, at the specified \code{times} as well as any times 
#'   generated through extrapolation (when \code{extrapolate = TRUE}).
#'   
#'   \subsection{Dynamic versus marginalised predictions}{
#'   The following also applies for \code{stanjm} objects.
#'   By default the survival probabilities are conditional on an individual's 
#'   group-specific coefficients (i.e. their individual-level random
#'   effects). If prediction data is provided via the \code{newdataLong}  
#'   and \code{newdataEvent} arguments, then the default behaviour is to
#'   sample new group-specific coefficients for the individuals in the  
#'   new data using a Monte Carlo scheme that conditions on their 
#'   longitudinal outcome data provided in \code{newdataLong} 
#'   (sometimes referred to as "dynamic predictions", see Rizopoulos
#'   (2011)). This default behaviour can be stopped by specifying 
#'   \code{dynamic = FALSE}, in which case the predicted survival
#'   probabilities will be marginalised over the distribution of the 
#'   group-specific coefficients. This has the benefit that the user does
#'   not need to provide longitudinal outcome measurements for the new 
#'   individuals, however, it does mean that the predictions will incorporate
#'   all the uncertainty associated with between-individual variation in the
#'   biomarker (longitudinal outcome) values since the predictions aren't 
#'   conditional on any observed biomarker (longitudinal outcome) data for 
#'   the individual.
#'   }
#'   
#' @note 
#'   Note that if any variables were transformed (e.g. rescaled) in the data 
#'   used to fit the model, then these variables must also be transformed in 
#'   \code{newdataLong} and \code{newdataEvent}. This only applies if variables  
#'   were transformed before passing the data to one of the modeling functions 
#'   and \emph{not} if transformations were specified inside the model formula.
#'    
#' @return When \code{return_matrix = FALSE} (the default), a data frame of 
#'   class \code{survfit.stansurv} or \code{survfit.stanjm}. The data frame
#'   includes columns for each of the following: 
#'   (i) the median of the posterior predictions (\code{median});
#'   (ii) each of the lower and upper limits of the corresponding uncertainty 
#'   interval for the posterior predictions (\code{ci_lb} and \code{ci_ub});
#'   (iii) an observation identifier (for \code{stan_surv} models) or an 
#'   individual identifier (for \code{stan_jm} models), unless standardised 
#'   predictions were requested;
#'   (iv) the time that the prediction corresponds to (\code{time}).
#'   (v) the last known survival time on which the prediction is conditional 
#'   (\code{cond_time}); this will be set to NA if not relevant.
#'   The returned object also includes a number of additional attributes.
#'   
#'   When \code{return_matrix = TRUE} a list of matrices is returned. Each 
#'   matrix contains the predictions evaluated at one step of the 
#'   extrapolation time sequence (note that if \code{extrapolate = FALSE} 
#'   then the list will be of length one, i.e. the predictions are only 
#'   evaluated at \code{times} which corresponds to just one time point 
#'   for each individual). Each matrix will have \code{draws} rows and 
#'   \code{nrow(newdata)} columns, such that each row contains a 
#'   vector of predictions generated using a single draw of the model 
#'   parameters from the posterior distribution. The returned 
#'   list also includes a number of additional attributes.
#'    
#' @seealso 
#'   \code{\link{plot.survfit.stanjm}} for plotting the estimated survival  
#'   probabilities \cr
#'   \code{\link{ps_check}} for for graphical checks of the estimated 
#'   survival function \cr
#'   \code{\link{posterior_traj}} for estimating the
#'   marginal or subject-specific longitudinal trajectories \cr
#'   \code{\link{plot_stack_jm}} for combining plots of the estimated 
#'   subject-specific longitudinal trajectory and survival function
#'   
#' @references 
#'   Rizopoulos, D. (2011). Dynamic predictions and prospective accuracy in 
#'   joint models for longitudinal and time-to-event data. \emph{Biometrics}
#'   \strong{67}, 819.
#'      
#' @examples
#' \donttest{
#'   # Run example model if not already loaded
#'   if (!exists("example_jm")) example(example_jm)
#'   
#'   # Obtain subject-specific survival probabilities for a few
#'   # selected individuals in the estimation dataset who were  
#'   # known to survive up until their censoring time. By default
#'   # the posterior_survfit function will estimate the conditional
#'   # survival probabilities, that is, conditional on having survived
#'   # until the event or censoring time, and then by default will
#'   # extrapolate the survival predictions forward from there.  
#'   ps1 <- posterior_survfit(example_jm, ids = c(7,13,15))
#'   # We can plot the estimated survival probabilities using the
#'   # associated plot function
#'   plot(ps1)
#'   
#'   # If we wanted to estimate the survival probabilities for the
#'   # same three individuals as the previous example, but this time
#'   # we won't condition on them having survived up until their 
#'   # censoring time. Instead, we will estimate their probability
#'   # of having survived between 0 and 5 years given their covariates
#'   # and their estimated random effects.
#'   # The easiest way to achieve the time scale we want (ie, 0 to 5 years)
#'   # is to specify that we want the survival time estimated at time 0
#'   # and then extrapolated forward 5 years. We also specify that we
#'   # do not want to condition on their last known survival time.
#'   ps2 <- posterior_survfit(example_jm, ids = c(7,13,15), times = 0,
#'     extrapolate = TRUE, condition = FALSE, control = list(edist = 5))
#'     
#'   # Instead we may want to estimate subject-specific survival probabilities 
#'   # for a set of new individuals. To demonstrate this, we will simply take
#'   # the first two individuals in the estimation dataset, but pass their data
#'   # via the newdata arguments so that posterior_survfit will assume we are 
#'   # predicting survival for new individuals and draw new random effects 
#'   # under a Monte Carlo scheme (see Rizopoulos (2011)).
#'   ndL <- pbcLong[pbcLong$id %in% c(1,2),]
#'   ndE <- pbcSurv[pbcSurv$id %in% c(1,2),]
#'   ps3 <- posterior_survfit(example_jm,
#'     newdataLong = ndL, newdataEvent = ndE,
#'     last_time = "futimeYears", seed = 12345)
#'   head(ps3)
#'   # We can then compare the estimated random effects for these 
#'   # individuals based on the fitted model and the Monte Carlo scheme
#'   ranef(example_jm)$Long1$id[1:2,,drop=FALSE] # from fitted model
#'   colMeans(attr(ps3, "b_new"))                # from Monte Carlo scheme
#'   
#'   # Lastly, if we wanted to obtain "standardised" survival probabilities, 
#'   # (by averaging over the observed distribution of the fixed effect 
#'   # covariates, as well as averaging over the estimated random effects
#'   # for individuals in our estimation sample or new data) then we can 
#'   # specify 'standardise = TRUE'. We can then plot the resulting 
#'   # standardised survival curve.
#'   ps4 <- posterior_survfit(example_jm, standardise = TRUE, 
#'                            times = 0, extrapolate = TRUE)
#'   plot(ps4)
#' }
#'  
posterior_survfit <- function(object, ...) UseMethod("posterior_survfit")

#' @rdname posterior_survfit
#' @method posterior_survfit stansurv
#' @export
#'
posterior_survfit.stansurv <- function(object, 
                                       newdata       = NULL, 
                                       type          = "surv", 
                                       extrapolate   = TRUE, 
                                       control       = list(), 
                                       condition     = FALSE, 
                                       last_time     = NULL, 
                                       prob          = 0.95,
                                       times         = NULL,
                                       standardise   = FALSE, 
                                       draws         = NULL, 
                                       seed          = NULL,
                                       return_matrix = FALSE,
                                       ...) {
  
  validate_stansurv_object(object)
  
  basehaz  <- object$basehaz
  
  if (!is.null(seed)) 
    set.seed(seed)
  
  if (is.null(newdata) && object$ndelayed)
    stop("Prediction data for 'posterior_survfit' cannot include delayed ",
         "entry. If you estimated a model with delayed entry, you will ",
         "not be able to obtain predictions using the estimation data frame. ",
         "You must provide prediction data via the 'newdata' argument, and ",
         "indicate delayed entry via the 'last_time' argument.")
  
  dots <- list(...)
  
  newdata <- validate_newdata(newdata)
  has_newdata <- not.null(newdata)
  
  # Obtain a vector of unique subject ids 
  if (is.null(newdata)) {
    id_list <- seq(nrow(get_model_data(object)))
  } else {
    id_list <- seq(nrow(newdata))
  }
  
  # Error checks for conditional predictions
  if (condition) {
    if (standardise) 
      stop("'condition' cannot be TRUE for standardised predictions.")
    if (type %in% c("haz", "loghaz"))
      stop("'condition' cannot be TRUE when 'type = \"", type, "\"'.")
  }
  
  # Last known survival time for each individual
  if (is.null(newdata)) { # user did not specify newdata
    if (!is.null(last_time))
      stop("'last_time' cannot be provided when newdata is NULL, since times ",
           "are taken to be the event or censoring time for each individual.")
    last_time <- object$eventtime
  } else { # user specified newdata
    if (is.null(last_time)) { # assume at risk from time zero
      last_time <- rep(0, length(id_list))
    } else if (is.string(last_time)) {
      if (!last_time %in% colnames(newdata))
        stop("Cannot find 'last_time' column named in newdata")
      last_time <- newdata[[last_time]]      
    } else if (is.scalar(last_time)) {
      last_time <- rep(last_time, nrow(newdata)) 
    } else if (any(!is.numeric(last_time), !length(last_time) == nrow(newdata))) {
      stop("Bug found: could not reconcile 'last_time' argument.")
    }
    names(last_time) <- as.character(id_list)
  }
  
  # Prediction times
  if (standardise) { # standardised survival probs
    times <- 
      if (is.null(times)) {
        stop("'times' cannot be NULL for obtaining standardised survival probabilities.")
      } else if (is.scalar(times)) {
        rep(times, length(id_list))
      } else {
        stop("'times' should be a numeric vector of length 1 in order to obtain ",
             "standardised survival probabilities (the subject-specific survival ",
             "probabilities will be calculated at the specified time point, and ",
             "then averaged).")      
      }    
  } else if (is.null(newdata)) { # subject-specific survival probs without newdata
    times <- 
      if (is.null(times)) {
        object$eventtime
      } else if (is.scalar(times)) {
        rep(times, length(id_list))
      } else {
        stop("If newdata is NULL then 'times' must be NULL or a single number.")     
      }
  } else { # subject-specific survival probs with newdata
    times <- 
      if (is.null(times)) {
        times <- last_time
      } else if (is.scalar(times)) {
        rep(times, length(id_list))
      } else if (is.string(times)) {
        if (!times %in% colnames(newdata))
          stop("Variable specified in 'times' argument could not be found in newdata.")
        times <- newdata[[times]]
      } else {
        stop("If newdata is specified then 'times' can only be the name of a ",
             "variable in newdata, or a single number.")
      }
  } 
  
  maxtime <- max(object$eventtime)
  if (any(times > maxtime))
    stop("'times' are not allowed to be greater than the last event or ",
         "censoring time (since unable to extrapolate the baseline hazard).")
  
  # User specified extrapolation
  if (extrapolate) {
    control <- extrapolation_control(control, ok_args = c("epoints", "edist"))
    if (not.null(control$edist)) {
      endtime <- times + control$edist
    } else {
      endtime <- maxtime
    }
    endtime <- truncate(endtime, upper = maxtime)
    time_seq <- get_time_seq(control$epoints, times, endtime, simplify = FALSE)
  } else {
    time_seq <- list(times) # no extrapolation
  }
  
  # Get stanmat parameter matrix for specified number of draws
  stanmat <- sample_stanmat(object, draws = draws, default_draws = NA)
  pars    <- extract_pars(object, stanmat)
  
  # Calculate survival probability at each increment of extrapolation sequence
  surv <- lapply(time_seq, .pp_calculate_surv, 
                 object      = object,
                 newdata     = newdata,
                 pars        = pars,
                 type        = type,
                 standardise = standardise)
  
  # Calculate survival probability at last known survival time and then
  # use that to calculate conditional survival probabilities
  if (condition) {
    cond_surv <- .pp_calculate_surv(last_time,
                                    object  = object,
                                    newdata = newdata,
                                    pars    = pars,
                                    type    = type)
    surv <- lapply(surv, function(x) truncate(x / cond_surv, upper = 1))
    attr(surv, "last_time") <- last_time
  }
  
  # Optionally return draws rather than summarising into median and CI
  if (return_matrix) {
    return(structure(surv,
                     type        = type,
                     extrapolate = extrapolate, 
                     control     = control,
                     condition   = condition,
                     standardise = standardise, 
                     last_time   = if (condition) last_time else NULL,
                     ids         = id_list,
                     draws       = NROW(stanmat),
                     seed        = seed))
  }
  
  # Summarise posterior draws to get median and CI
  out <- .pp_summarise_surv(surv        = surv,
                            prob        = prob,
                            standardise = standardise)
  
  # Add attributes
  structure(out,
            id_var      = attr(out, "id_var"),
            time_var    = attr(out, "time_var"),
            type        = type,
            extrapolate = extrapolate, 
            control     = control, 
            condition   = condition, 
            standardise = standardise, 
            last_time   = if (condition) last_time else NULL, 
            ids         = id_list, 
            draws       = NROW(stanmat), 
            seed        = seed, 
            class       = c("survfit.stansurv", "data.frame"))
}

#' @rdname posterior_survfit
#' @method posterior_survfit stanjm
#' @export
#'
posterior_survfit.stanjm <- function(object, 
                                     newdataLong  = NULL, 
                                     newdataEvent = NULL,
                                     type         = "surv", 
                                     extrapolate  = TRUE, 
                                     control      = list(),
                                     condition    = NULL, 
                                     last_time    = NULL, 
                                     prob         = 0.95,
                                     ids, 
                                     times        = NULL, 
                                     standardise  = FALSE,
                                     dynamic      = TRUE, 
                                     scale        = 1.5,
                                     draws        = NULL, 
                                     seed         = NULL, 
                                     return_matrix = FALSE, 
                                     ...) {
  
  validate_stanjm_object(object)
  
  M        <- object$n_markers
  id_var   <- object$id_var
  time_var <- object$time_var
  basehaz  <- object$basehaz
  assoc    <- object$assoc
  family   <- family(object)
  
  if (!is.null(seed)) 
    set.seed(seed)
  if (missing(ids)) 
    ids <- NULL
  
  dots <- list(...)
  
  # Temporarily only allow survprob for stan_jm until refactoring is done
  if (!type == "surv")
    stop("Currently only 'type = \"surv\"' is allowed for stanjm models.")
  
  # Temporary stop, until make_assoc_terms can handle it
  sel_stop <- grep("^shared", rownames(object$assoc))
  if (any(unlist(object$assoc[sel_stop,])))
    stop("'posterior_survfit' cannot yet be used with shared_b or shared_coef ",
         "association structures.") 

  # Construct prediction data
  # ndL: dataLong to be used in predictions
  # ndE: dataEvent to be used in predictions
  if (!identical(is.null(newdataLong), is.null(newdataEvent)))
    stop("Both newdataLong and newdataEvent must be supplied together.")
  has_newdata <- not.null(newdataEvent)
  if (is.null(newdataLong)) { # user did not specify newdata
    dats <- get_model_data(object)
    ndL <- dats[1:M]
    ndE <- dats[["Event"]]
  } else { # user specified newdata
    newdatas <- validate_newdatas(object, newdataLong, newdataEvent, 
                                  response = dynamic, needs_time_var = dynamic)
    ndL <- newdatas[1:M]
    ndE <- newdatas[["Event"]]   
  }
  if (!is.null(ids)) { # user specified a subset of ids
    ndL <- subset_ids(ndL, ids, id_var)
    ndE <- subset_ids(ndE, ids, id_var)
  }  
  id_list <- factor(unique(ndE[[id_var]])) # order of ids from data, not ids arg

  # Last known survival time for each individual
  if (is.null(newdataLong)) { # user did not specify newdata
    if (!is.null(last_time))
      stop("'last_time' cannot be provided when newdata is NULL, since times ",
           "are taken to be the event or censoring time for each individual.")
    last_time <- object$eventtime[as.character(id_list)]
  } else { # user specified newdata
    if (is.null(last_time)) { # use latest longitudinal observation
      max_ytimes <- do.call("cbind", lapply(ndL, function(x) 
        tapply(x[[time_var]], x[[id_var]], FUN = max)))
      last_time <- apply(max_ytimes, 1L, max)
      # re-order last-time according to id_list
      last_time <- last_time[as.character(id_list)]
    } else if (is.character(last_time) && (length(last_time) == 1L)) {
      if (!last_time %in% colnames(ndE))
        stop("Cannot find 'last_time' column named in newdataEvent.")
      last_time <- ndE[[last_time]]      
    } else if (is.numeric(last_time) && (length(last_time) == 1L)) {
      last_time <- rep(last_time, length(id_list)) 
    } else if (is.numeric(last_time) && (length(last_time) > 1L)) {
      last_time <- last_time[as.character(id_list)]
    } else {
      stop("Bug found: could not reconcile 'last_time' argument.")
    }
    names(last_time) <- as.character(id_list)
  }   
  
  # Prediction times
  if (standardise) { # standardised survival probs
    times <- 
      if (is.null(times)) {
        stop("'times' cannot be NULL for obtaining standardised survival probabilities.")
      } else if (is.numeric(times) && (length(times) == 1L)) {
        rep(times, length(id_list))
      } else {
        stop("'times' should be a numeric vector of length 1 in order to obtain ",
             "standardised survival probabilities (the subject-specific survival ",
             "probabilities will be calculated at the specified time point, and ",
             "then averaged).")      
      }    
  } else if (is.null(newdataLong)) { # subject-specific survival probs without newdata
    times <- 
      if (is.null(times)) {
        object$eventtime[as.character(id_list)]
      } else if (is.numeric(times) && (length(times) == 1L)) {
        rep(times, length(id_list))
      } else {
        stop("If newdata is NULL then 'times' must be NULL or a single number.")     
      }
  } else { # subject-specific survival probs with newdata
    times <- 
      if (is.null(times)) {
        times <- last_time
      } else if (is.character(times) && (length(times) == 1L)) {
        if (!times %in% colnames(ndE))
          stop("Variable specified in 'times' argument could not be found in newdata.")
        tapply(ndE[[times]], ndE[[id_var]], FUN = max)
      } else if (is.numeric(times) && (length(times) == 1L)) {
        rep(times, length(id_list))
      } else {
        stop("If newdata is specified then 'times' can only be the name of a ",
             "variable in newdata, or a single number.")      
      }
  }
  if (!identical(length(times), length(id_list)))
    stop(paste0("length of the 'times' vector should be equal to the number of individuals ",
                "for whom predictions are being obtained (", length(id_list), ")."))     
  maxtime <- max(object$eventtime)
  if (any(times > maxtime))
    stop("'times' are not allowed to be greater than the last event or censoring ",
         "time (since unable to extrapolate the baseline hazard).")
  
  # User specified extrapolation
  if (extrapolate) {
    ok_args <- c("epoints", "edist")
    control <- extrapolation_control(control, ok_args = ok_args)
    endtime <- if (!is.null(control$edist)) times + control$edist else maxtime
    endtime[endtime > maxtime] <- maxtime # nothing beyond end of baseline hazard 
    time_seq <- get_time_seq(control$epoints, times, endtime, simplify = FALSE)
  } else time_seq <- list(times) # no extrapolation

  # Conditional survival times
  if (is.null(condition)) {
    condition <- ifelse(type == "surv", !standardise, FALSE)
  } else if (condition && standardise) {
    stop("'condition' cannot be set to TRUE if standardised survival ",
         "probabilities are requested.")
  }
  
  # Get stanmat parameter matrix for specified number of draws
  stanmat <- sample_stanmat(object, draws = draws, default_draws = 200)

  # Draw b pars for new individuals
  if (dynamic && has_newdata) {
    stanmat <- simulate_b_pars(object, 
                               stanmat = stanmat, 
                               ndL     = ndL, 
                               ndE     = ndE,
                               ids     = id_list, 
                               times   = last_time, 
                               scale   = scale)
    b_new           <- attr(stanmat, "b_new")
    acceptance_rate <- attr(stanmat, "acceptance_rate")
  }
  
  pars <- extract_pars(object, stanmat) # list of stanmat arrays
  
  # Matrix of surv probs at each increment of the extrapolation sequence
  # NB If no extrapolation then length(time_seq) == 1L
  surv_t <- lapply(time_seq, .pp_calculate_surv, 
                   object       = object,
                   newdataLong  = ndL, 
                   newdataEvent = ndE,
                   pars         = pars,
                   type         = type,
                   id_list      = id_list,
                   standardise  = standardise,
                   dynamic      = dynamic)
  
  # Calculate survival probability at last known survival time and then
  # use that to calculate conditional survival probabilities   
  if (condition) {
    if (!type == "surv")
      stop("'condition' can only be set to TRUE for survival probabilities.")
    cond_surv <- .pp_calculate_surv(last_time,
                                    object       = object,
                                    newdataLong  = ndL, 
                                    newdataEvent = ndE, 
                                    pars         = pars,
                                    type         = type,
                                    id_list      = id_list,
                                    dynamic      = dynamic)
    surv <- lapply(surv_t, function(x) truncate(x / cond_surv, upper = 1))
    attr(surv, "last_time") <- last_time
  } else {
    surv <- surv_t
  }
 
  # Optionally return draws rather than summarising into median and CI
  if (return_matrix) {
    return(structure(surv,
                     type        = type,
                     extrapolate = extrapolate, 
                     control     = control,
                     condition   = condition,
                     standardise = standardise, 
                     last_time   = if (condition) last_time else NULL,
                     ids         = id_list,
                     draws       = NROW(stanmat),
                     seed        = seed))
  }  
  
  # Summarise posterior draws to get median and CI
  out <- .pp_summarise_surv(surv        = surv,
                            prob        = prob,
                            id_var      = id_var,
                            time_var    = time_var,
                            standardise = standardise)
  
  # Temporary hack so that 'predictive_error' can call 'posterior_survfit'
  # with two separate conditioning times...
  fun_check <- isTRUE(grepl("predictive_error", get_calling_fun(), fixed = TRUE))
  dot_check <- isTRUE("last_time2" %in% names(dots))
  if (fun_check && dot_check) {
    if (!type == "surv")
      stop("'last_time2' can only be specified for survival probabilities.")
    cond_surv2 <- .pp_calculate_surv(ndE[[dots$last_time2]],
                                     object       = object,
                                     newdataLong  = ndL, 
                                     newdataEvent = ndE, 
                                     pars         = pars,
                                     type         = type,
                                     id_list      = id_list,
                                     dynamic      = dynamic)
    surv2 <- lapply(surv_t, function(x) truncate(x / cond_surv2, upper = 1))       
    out2 <- .pp_summarise_surv(surv        = surv2,
                               prob        = prob,
                               id_var      = id_var,
                               time_var    = time_var,
                               standardise = standardise,
                               colnames    = "survprob_eventtime")
    out <- merge(out, out2)
  }
  
  # Return object
  out <- structure(out, 
                   id_var      = id_var, 
                   time_var    = time_var, 
                   type        = type,
                   extrapolate = extrapolate,
                   control     = control, 
                   standardise = standardise, 
                   condition   = condition, 
                   last_time   = if (condition) last_time else NULL, 
                   ids         = id_list, 
                   draws       = NROW(stanmat), 
                   seed        = seed,
                   offset      = offset,
                   class       = c("survfit.stanjm", "data.frame"))
  
  if (dynamic && has_newdata) {
    out <- structure(out, b_new = b_new, acceptance_rate = acceptance_rate)
  }
  
  out
}


# -----------------  internal  ------------------------------------------------

# Calculate the desired prediction (e.g. hazard, cumulative hazard, survival
# probability) at the specified times
.pp_calculate_surv <- function(times, 
                               object,
                               newdata      = NULL,
                               newdataLong  = NULL,
                               newdataEvent = NULL, 
                               pars, 
                               type         = "surv", 
                               id_list      = NULL,
                               standardise  = FALSE,
                               dynamic      = TRUE) {
  
  if (is.stanjm(object) && !identical(length(times), length(id_list)))
    stop("Bug found: vector of ids should be same length as vector of times.")
  
  # Determine whether prediction type requires quadrature
  needs_quadrature <- type %in% c("cumhaz", 
                                  "surv", 
                                  "cdf", 
                                  "logcumhaz",
                                  "logsurv", 
                                  "logcdf")
  
  # Evaluate hazard, cumulative hazard, survival or failure probability
  if (is.stansurv(object)) {
    ppdat  <- .pp_data_surv(object, 
                            newdata       = newdata, 
                            times         = times,
                            at_quadpoints = needs_quadrature)
    out <- .pp_predict_surv(object, 
                            data = ppdat, 
                            pars = pars, 
                            type = type)
  } else if (is.stanjm(object)) {
    ppdat <- .pp_data_jm(object, 
                         newdataLong  = newdataLong, 
                         newdataEvent = newdataEvent, 
                         ids          = id_list, 
                         etimes       = times, 
                         long_parts   = FALSE,
                         response     = dynamic,
                         needs_time_var = dynamic)
    out  <- .ll_survival(object, # refactoring for stanjm not yet finished
                         data     = ppdat, 
                         pars     = pars, 
                         survprob = TRUE)
  }
  
  # Transform if only one individual 
  out <- transpose_vector(out)
  
  # Set survival probability == 1 if time == 0 (avoids possible NaN)
  if (type == "surv")
    out <- replace_where(out, times == 0, replacement = 1, margin = 2L)
  
  # Standardisation: within each iteration, calculate mean across individuals 
  if (standardise) {
    out   <- row_means(out)
    ids   <- "standardised_survprob"
    times <- unique(times)
  } else {
    ids   <- if (is.null(id_list)) seq(ncol(out)) else id_list
  }
  dimnames(out) <- list(iterations = NULL, ids = ids)
  
  # Add subject ids and prediction times as an attribute
  structure(out, ids = ids, times = times)
}


# Evaluate hazard, cumulative hazard, survival or failure probability
#
# @param object A stansurv or stanjm object.
# @param data Output from .pp_data_surv or .pp_data_jm.
# @param pars Output from extract_pars.
# @param type The type of prediction quantity to return.
.pp_predict_surv <- function(object, ...) UseMethod(".pp_predict_surv")

.pp_predict_surv.stansurv <- function(object,
                                      data, 
                                      pars,
                                      type = "surv") {
  
  args <- nlist(basehaz   = get_basehaz(object),
                intercept = pars$alpha,
                betas     = pars$beta,
                betas_tve = pars$beta_tve,
                b         = pars$b,
                aux       = pars$aux,
                times     = data$pts,
                x         = data$x,
                s         = data$s,
                z         = data$z)  
  
  if (type %in% c("loghaz", "haz")) { 
    # evaluate hazard; quadrature not relevant
    lhaz <- do.call(evaluate_log_haz, args)
  } else if (!data$has_quadrature){ 
    # evaluate survival; without quadrature
    lsurv <- do.call(evaluate_log_surv, args)
  } else { 
    # evaluate survival; with quadrature
    lhaz  <- do.call(evaluate_log_haz, args)
    lsurv <- -quadrature_sum(exp(lhaz), qnodes = data$qnodes, qwts = data$wts)
  }
  
  switch(type,
         loghaz    = lhaz,
         logcumhaz = log(-lsurv),
         logsurv   = lsurv,
         logcdf    = log(1 - exp(lsurv)),
         haz       = exp(lhaz),
         surv      = exp(lsurv),
         cumhaz    = -lsurv,
         cdf       = 1 - exp(lsurv),
         stop("Invalid input to the 'type' argument."))
}


# Summarise predictions into median, lower CI, upper CI
#
# @details Convert a list of matrices (with each element being a S by N matrix, 
#   where S is the number of MCMC draws and N the number of individuals)
#   and collapse it across the MCMC iterations summarising it into median 
#   and CI. The result is a data frame with K times N rows, where K was 
#   the length of the original list.
.pp_summarise_surv <- function(surv, 
                               prob        = NULL, 
                               id_var      = NULL, 
                               time_var    = NULL,
                               standardise = FALSE,
                               colnames    = NULL) {
  
  # Default variable names if not provided by the user
  if (is.null(id_var)) 
    id_var <- "id"
  if (is.null(time_var))
    time_var <- "time"
  
  # Define variable name for conditioning time
  cond_var <- paste0("cond_", time_var)
  
  # Extract ids and times for the predictions
  ids       <- uapply(surv, attr, "ids")
  times     <- uapply(surv, attr, "times")
  
  # Extract conditioning time that was used for predictions
  last_time <- attr(surv, "last_time")
  if (is.null(last_time)) { # if not using conditional survival
    last_time <- rep(NA, length(ids)) 
  }

  # Determine the quantiles corresponding to the median and CI limits
  if (is.null(prob)) {
    probs <- 0.5 # median only
    nms   <- c(id_var, cond_var, time_var, "median")
  } else {
    probs <- c(0.5, (1 - prob)/2, (1 + prob)/2) # median and CI
    nms   <- c(id_var, cond_var, time_var, "median", "ci_lb", "ci_ub")
  }
  
  # Possibly overide default variable names for the returned data frame
  if (!is.null(colnames)) {
    nms <- c(id_var, cond_var, time_var, colnames)
  }
  
  # Calculate mean and CI at each prediction time
  out <- data.frame(do.call("rbind", lapply(surv, col_quantiles_, probs)))
  out <- mutate_(out, id_var = ids, cond_var = last_time, time_var = times)
  out <- row_sort(out, id_var, time_var)
  out <- col_sort(out, id_var, cond_var, time_var)
  out <- set_rownames(out, NULL)
  out <- set_colnames(out, nms)
  
  # Drop excess info if standardised predictions were calculated
  if (standardise) { 
    out[[cond_var]] <- NULL
    out[[id_var]]   <- NULL
    id_var          <- NULL
  }
  
  structure(out, 
            id_var   = id_var, 
            time_var = time_var)
}


# ------------  print methods  ------------------------------------------------

#' Generic print method for \code{survfit.stansurv} and \code{survfit.stanjm} 
#' objects
#' 
#' @rdname print.survfit.stansurv
#' @method print survfit.stansurv
#' @keywords internal
#' @export
#' @param x An object of class \code{survfit.stansurv} or \code{survfit.stanjm}, 
#'   returned by a call to \code{\link{posterior_survfit}}.
#' @param digits Number of digits to use for formatting the time variable and 
#'   the survival probabilities.
#' @param ... Ignored.
#' 
print.survfit.stansurv <- function(x, digits = 4, ...) {
  
  x <- as.data.frame(x)
  sel <- c(attr(x, "time_var"), "median", "ci_lb", "ci_ub")
  for (i in sel) 
    x[[i]] <- format(round(x[[i]], digits), nsmall = digits)
  
  cat("stan_surv predictions\n")
  cat(" num. individuals:", length(attr(x, "ids")), "\n")
  cat(" prediction type: ", tolower(get_survpred_name(attr(x, "type"))), "\n")
  cat(" standardised?:   ", yes_no_string(attr(x, "standardise")), "\n")
  cat(" conditional?:    ", yes_no_string(attr(x, "condition")), "\n\n")
  print(x, quote = FALSE)
  invisible(x)
}

#' @rdname print.survfit.stansurv
#' @method print survfit.stanjm
#' @export
#' 
print.survfit.stanjm <- function(x, digits = 4, ...) {
  
  x <- as.data.frame(x)
  sel <- c(attr(x, "time_var"), "median", "ci_lb", "ci_ub")
  for (i in sel) 
    x[[i]] <- format(round(x[[i]], digits), nsmall = digits)
  
  cat("stan_jm predictions\n")
  cat(" num. individuals:", length(attr(x, "ids")), "\n")
  cat(" prediction type: ", tolower(get_survpred_name(attr(x, "type"))), "\n")
  cat(" standardised?:   ", yes_no_string(attr(x, "standardise")), "\n")
  cat(" conditional?:    ", yes_no_string(attr(x, "condition")), "\n\n")
  print(x, quote = FALSE)
  invisible(x)
}


# -----------------  plot methods  --------------------------------------------

#' Plot the estimated subject-specific or marginal survival function
#' 
#' This generic \code{plot} method for \code{survfit.stanjm} objects will
#' plot the estimated subject-specific or marginal survival function
#' using the data frame returned by a call to \code{\link{posterior_survfit}}.
#' The call to \code{posterior_survfit} should ideally have included an
#' "extrapolation" of the survival function, obtained by setting the 
#' \code{extrapolate} argument to \code{TRUE}.
#'    
#' @method plot survfit.stanjm
#' @export
#' @importFrom ggplot2 ggplot aes_string geom_line geom_ribbon 
#'   facet_wrap labs coord_cartesian
#'   
#' @templateVar idsArg ids
#' @templateVar labsArg xlab,ylab
#' @templateVar scalesArg facet_scales
#' @templateVar cigeomArg ci_geom_args
#' @template args-ids
#' @template args-labs
#' @template args-scales
#' @template args-ci-geom-args
#'  
#' @param x A data frame and object of class \code{survfit.stanjm}
#'   returned by a call to the function \code{\link{posterior_survfit}}.
#'   The object contains point estimates and uncertainty interval limits
#'   for estimated values of the survival function.
#' @param limits A quoted character string specifying the type of limits to
#'   include in the plot. Can be one of: \code{"ci"} for the Bayesian
#'   posterior uncertainty interval for the estimated survival probability
#'   (often known as a credible interval); or \code{"none"} for no interval 
#'   limits.
#' @param ... Optional arguments passed to 
#'   \code{\link[ggplot2]{geom_line}} and used to control features
#'   of the plotted survival function.
#'      
#' @return The plot method returns a \code{ggplot} object, also of class
#'   \code{plot.survfit.stanjm}. This object can be further customised using the
#'   \pkg{ggplot2} package. It can also be passed to the function
#'   \code{plot_stack_jm}.
#'   
#' @seealso \code{\link{posterior_survfit}}, \code{\link{plot_stack_jm}},
#'   \code{\link{posterior_traj}}, \code{\link{plot.predict.stanjm}}      
#'   
#' @examples 
#' \donttest{
#'   # Run example model if not already loaded
#'   if (!exists("example_jm")) example(example_jm)
#'   
#'   # Obtain subject-specific conditional survival probabilities
#'   # for all individuals in the estimation dataset.
#'   ps1 <- posterior_survfit(example_jm, extrapolate = TRUE)
#'   
#'   # We then plot the conditional survival probabilities for
#'   # a subset of individuals
#'   plot(ps1, ids = c(7,13,15))
#'   # We can change or add attributes to the plot
#'   plot(ps1, ids = c(7,13,15), limits = "none")
#'   plot(ps1, ids = c(7,13,15), xlab = "Follow up time")
#'   plot(ps1, ids = c(7,13,15), ci_geom_args = list(fill = "red"),
#'        color = "blue", linetype = 2)
#'   plot(ps1, ids = c(7,13,15), facet_scales = "fixed")
#'   
#'   # Since the returned plot is also a ggplot object, we can
#'   # modify some of its attributes after it has been returned
#'   plot1 <- plot(ps1, ids = c(7,13,15))
#'   plot1 + 
#'     ggplot2::theme(strip.background = ggplot2::element_blank()) +
#'     ggplot2::coord_cartesian(xlim = c(0, 15)) +
#'     ggplot2::labs(title = "Some plotted survival functions")
#'     
#'   # We can also combine the plot(s) of the estimated 
#'   # subject-specific survival functions, with plot(s) 
#'   # of the estimated longitudinal trajectories for the
#'   # same individuals
#'   ps1 <- posterior_survfit(example_jm, ids = c(7,13,15))
#'   pt1 <- posterior_traj(example_jm, , ids = c(7,13,15))
#'   plot_surv <- plot(ps1) 
#'   plot_traj <- plot(pt1, vline = TRUE, plot_observed = TRUE)
#'   plot_stack_jm(plot_traj, plot_surv)
#'    
#'   # Lastly, let us plot the standardised survival function
#'   # based on all individuals in our estimation dataset
#'   ps2 <- posterior_survfit(example_jm, standardise = TRUE, times = 0,
#'                           control = list(epoints = 20))
#'   plot(ps2)   
#' }
#'    
plot.survfit.stanjm <- function(x, 
                                ids    = NULL, 
                                limits = c("ci", "none"),  
                                xlab   = NULL, 
                                ylab   = NULL, 
                                facet_scales = "free", 
                                ci_geom_args = NULL, ...) {
  
  limits <- match.arg (limits)
  ci     <- as.logical(limits == "ci")
  
  type        <- attr(x, "type")
  standardise <- attr(x, "standardise")
  id_var      <- attr(x, "id_var")
  time_var    <- attr(x, "time_var")
  
  if (is.null(xlab)) xlab <- paste0("Time (", time_var, ")")
  if (is.null(ylab)) ylab <- get_survpred_name(type)
  
  if (!is.null(ids)) {
    if (standardise) 
      stop("'ids' argument cannot be specified when plotting standardised ",
           "survival probabilities.")
    x <- subset_ids(x, ids, id_var)
  } else {
    ids <- if (!standardise) attr(x, "ids") else NULL
  }
  if (!standardise) x$id <- factor(x[[id_var]])
  x$time <- x[[time_var]]
  
  geom_defaults <- list(color = "black")
  geom_mapp     <- list(mapping = aes_string(x = "time",
                                             y = "median"))
  geom_args     <- do.call("set_geom_args", 
                           c(defaults = list(geom_defaults), list(...)))  
  
  lim_defaults  <- list(alpha = 0.3)
  lim_mapp      <- list(mapping = aes_string(x = "time", 
                                             ymin = "ci_lb", 
                                             ymax = "ci_ub"))
  lim_args      <- do.call("set_geom_args", 
                           c(defaults = list(lim_defaults), ci_geom_args))
  
  if ((!standardise) && (length(ids) > 60L))
    stop("Too many individuals to plot for. Perhaps consider limiting ",
         "the number of individuals by specifying the 'ids' argument.")
  
  graph_base <-
    ggplot(x) + 
    theme_bw() + 
    coord_cartesian(ylim = get_survpred_ylim(type)) +
    do.call("geom_line", c(geom_mapp, geom_args))
  
  graph_facet <- 
    if ((!standardise) && (length(ids) > 1L)) {
      facet_wrap(~ id, scales = facet_scales)
    } else NULL
  
  graph_limits <- 
    if (ci) {
      do.call("geom_ribbon", c(lim_mapp, lim_args)) 
    } else NULL  
  
  graph_labels <- labs(x = xlab, y = ylab)
  
  gg        <- graph_base + graph_facet + graph_limits + graph_labels
  class_gg  <- class(gg)
  class(gg) <- c("plot.survfit.stanjm", class_gg)
  gg
}


#' @rdname plot.survfit.stanjm
#' @method plot survfit.stansurv
#' @export
#' 
plot.survfit.stansurv <- function(x, 
                                  ids = NULL, 
                                  limits = c("ci", "none"),  
                                  xlab = NULL, 
                                  ylab = NULL, 
                                  facet_scales = "free", 
                                  ci_geom_args = NULL, ...) {
  mc <- match.call(expand.dots = FALSE)
  mc[[1L]] <- quote(plot.survfit.stanjm)
  ret <- eval(mc)
  class(ret)[[1L]] <- "plot.survfit.stansurv"
  ret
}


#' @rdname plot.survfit.stanjm
#' @export
#' @importFrom ggplot2 ggplot_build facet_wrap aes_string expand_limits
#' 
#' @description The \code{plot_stack_jm} function takes arguments containing the plots of the estimated  
#' subject-specific longitudinal trajectory (or trajectories if a multivariate  
#' joint model was estimated) and the plot of the estimated subject-specific 
#' survival function and combines them into a single figure. This is most
#' easily understood by running the \strong{Examples} below.
#' 
#' @param yplot An object of class \code{plot.predict.stanjm}, returned by a 
#'   call to the generic \code{\link[=plot.predict.stanjm]{plot}} method for 
#'   objects of class \code{predict.stanjm}. If there is more than one 
#'   longitudinal outcome, then a list of such objects can be provided.
#' @param survplot An object of class \code{plot.survfit.stanjm}, returned by a 
#'   call to the generic \code{\link[=plot.survfit.stanjm]{plot}} method for 
#'   objects of class \code{survfit.stanjm}. 
#'   
#' @return \code{plot_stack_jm} returns an object of class
#'   \code{\link[bayesplot]{bayesplot_grid}} that includes plots of the
#'   estimated subject-specific longitudinal trajectories stacked on top of the 
#'   associated subject-specific survival curve.
#'   
#' @seealso \code{\link{plot.predict.stanjm}}, \code{\link{plot.survfit.stanjm}},
#'   \code{\link{posterior_predict}}, \code{\link{posterior_survfit}}
#'    
#' @examples
#' \donttest{
#'   if (!exists("example_jm")) example(example_jm)
#'   ps1 <- posterior_survfit(example_jm, ids = c(7,13,15))
#'   pt1 <- posterior_traj(example_jm, ids = c(7,13,15), extrapolate = TRUE)
#'   plot_surv <- plot(ps1) 
#'   plot_traj <- plot(pt1, vline = TRUE, plot_observed = TRUE)
#'   plot_stack_jm(plot_traj, plot_surv)
#' }
#'  
plot_stack_jm <- function(yplot, survplot) {
  
  if (!is(yplot, "list")) yplot <- list(yplot)
  
  lapply(yplot, function(x) {
    if (!is(x, "plot.predict.stanjm"))
      stop("'yplot' should be an object of class 'plot.predict.stanjm', ",
           "or a list of such objects.", call. = FALSE)
  })
  if (!is(survplot, "plot.survfit.stanjm"))
    stop("'survplot' should be an object of class 'plot.survfit.stanjm'.",
         call. = FALSE)   
  
  y_build <- lapply(yplot, ggplot_build)
  y_layout <- lapply(y_build, function(x) x$layout$panel_layout)
  y_ids <- lapply(y_layout, function(x)
    if (!"id" %in% colnames(x)) NULL else x[["id"]])
  
  e_build <- ggplot_build(survplot)
  e_layout <- e_build$layout$panel_layout    
  e_ids <- if (!"id" %in% colnames(e_layout)) NULL else e_layout[["id"]]
  
  if (!is.null(e_ids)) {
    lapply(y_ids, function(x, e_ids) {
      if (!all(sort(x) == sort(e_ids))) {
        stop("The individuals in the 'yplot' and 'survplot' appear to differ. Please ",
             "reestimate the plots using a common 'ids' argument.", call. = FALSE)
      }
    }, e_ids = e_ids)    
  }
  
  vline <- lapply(seq_along(y_build), function(m) {
    L <- length(y_build[[m]]$data)
    dat <- y_build[[m]]$data[[L]]
    if (!"xintercept" %in% colnames(dat)) {
      found <- FALSE
    } else {
      found <- TRUE
      dat <- dat[, c("PANEL", "xintercept"), drop = FALSE] 
      if (NROW(y_layout[[m]]) > 1) {
        panel_id_map <- y_layout[[m]][, c("PANEL", "id"), drop = FALSE]
        dat <- merge(dat, panel_id_map, by = "PANEL")
      }
      dat <- dat[, grep("PANEL", colnames(dat), invert = TRUE), drop = FALSE]
      colnames(dat) <- gsub("xintercept", paste0("xintercept", m), colnames(dat), fixed = TRUE)
    }
    list(dat = dat, found = found)
  })
  vline_found <- any(sapply(vline, function(x) x$found))
  if (!vline_found)
    cat("Could not find vertical line indicating last observation time in the",
        "plot of the longitudinal trajectory; you may wish to plot the longitudinal",
        "trajectories again with 'vline = TRUE' to aid interpretation.")
  vline_dat <- lapply(vline, function(x) x$dat)
  vline_alldat <- Reduce(function(...) merge(..., all = TRUE), vline_dat)
  vline_alldat$xintercept_max <- 
    apply(vline_alldat[, grep("id", colnames(vline_alldat), invert = TRUE), drop = FALSE], 1, max) 
  
  xmax <- max(sapply(c(y_build, list(e_build)), function(i) max(i$data[[1]]$x)))
  
  if ((!is.null(e_ids)) && (length(e_ids) > 20L)) {
    stop("Unable to generate 'plot_stack_jm' for this many individuals.", call. = FALSE)      
  } else if ((!is.null(e_ids)) && (length(e_ids) > 3L)) {
    warning("'plot_stack_jm' is unlikely to be legible with more than a few individuals.",
            immediate. = TRUE, call. = FALSE)
  }

  if (!is.null(e_ids)) {
    graph_facet <- facet_wrap(~ id, scales = "free", nrow = 1) 
  } else {
    graph_facet <- NULL
  }
  
  if (vline_found) {
    graph_vline <- geom_vline(aes_string(xintercept = "xintercept_max"), 
                              vline_alldat, linetype = 2)
  } else {
    graph_vline <- NULL
  }
  
  graph_xlims <- expand_limits(x = c(0, xmax))
  
  survplot_updated <- survplot + graph_xlims + graph_facet + graph_vline
 
  yplot_updated <- lapply(yplot, function(x) x + graph_xlims + graph_facet)
  
  bayesplot::bayesplot_grid(
    plots = c(yplot_updated, list(survplot_updated)), 
    grid_args = list(ncol = 1)
  )
}


# -----------------  helpers  ------------------------------------------------

# Return a user-friendly name for the prediction type
get_survpred_name <- function(x) {
  switch(x, 
         haz       = "Hazard rate",
         cumhaz    = "Cumulative hazard rate",
         surv      = "Event free probability",
         cdf       = "Failure probability",
         loghaz    = "log(Hazard rate)",
         logcumhaz = "log(Cumulative hazard rate)",
         logsurv   = "log(Event free probability)",
         logcdf    = "log(Failure probability)",
         stop("Bug found: invalid input to 'type' argument."))
}

# Return appropriate y-axis limits for the prediction type
get_survpred_ylim <- function(x) {
  switch(x, 
         surv = c(0,1),
         cdf  = c(0,1),
         NULL)
}

# Default plotting attributes
.PP_FILL <- "skyblue"
.PP_DARK <- "skyblue4"
.PP_VLINE_CLR <- "#222222"
.PP_YREP_CLR <- "#487575"
.PP_YREP_FILL <- "#222222"
