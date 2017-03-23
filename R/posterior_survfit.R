# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2015, 2016 Trustees of Columbia University
# Copyright (C) 2016 Sam Brilleman
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

#' Estimate subject-specific or marginal survival probabilities
#' 
#' This function allows us to generate estimated survival probabilities (either 
#' subject-specific, or by marginalising over the distribution of the random effects) 
#' based on draws from the posterior predictive distribution. In both the 
#' "subject-specifc" and "marginal" situations, the predicted survival probabilities 
#' will still be \emph{conditional} on observed values of the fixed effect covariates 
#' in the longitudinal and event submodels (ie, the predictions will be obtained 
#' using either the design matrices used in the original \code{\link{stan_jm}} model
#' call, or using the covariate values provided in the \code{newdata} argument). However, 
#' if you wish to also average over the observed distribution of the fixed effect 
#' covariates then this is possible -- however we refer to these
#' as standardised survival probabilties -- see the \code{standardise} 
#' argument below.
#' 
#' @export
#' @templateVar stanjmArg object
#' @template args-stanjm-object
#' 
#' @param newdata Optionally, a new data frame in which to look 
#'   for variables with which to predict. If omitted, the model matrices are used. 
#'   If new data is provided, then it should contain the covariate values needed
#'   for all longitudinal submodel(s) and the event submodel. There is only
#'   allowed to be one row of data for each individual in \code{newdata}, that
#'   is, time-varying covariates are not allowed in the prediction dataset. Also 
#'   note that if \code{newdata} is provided, then the \code{times} argument
#'   must also be specified. See the \strong{Details} section for further 
#'   important details that need to be considered when specifying \code{newdata}.
#' @param extrapolate A logical specifying whether to extrapolate the estimated 
#'   survival probabilities beyond the times specified in the \code{times} argument.
#'   If \code{TRUE} then the extrapolation can be further controlled using
#'   the \code{control} argument.
#' @param control A named list with parameters controlling extrapolation 
#'   of the estimated survival function when \code{extrapolate = TRUE}. The list
#'   can contain one or more of the following named elements: \cr
#'   \describe{
#'     \item{\code{epoints}}{a positive integer specifying the number of  
#'     discrete time points at which to calculate the forecasted survival 
#'     probabilities. The default is 10.}
#'     \item{\code{edist}}{a positive scalar specifying the amount of time 
#'     across which to forecast the estimated survival function, represented 
#'     in units of the time variable \code{time_var} (from fitting the model). 
#'     The default is to extrapolate between the times specified in the 
#'     \code{times} argument and the maximum event or censoring time in the 
#'     original data. If \code{ext_distance} leads to times that are beyond
#'     the maximum event or censoring time (in the original data) then the 
#'     estimated survival probabilities will be truncated at that point, since
#'     the estimate for the baseline hazard is not available beyond that time.}
#'     \item{\code{condition}}{a logical specifying whether the estimated 
#'     subject-specific survival probabilities at time \code{t} should be 
#'     conditioned on survival up to a fixed time point \code{u}. The default 
#'     is to condition on the latest observation time for each individual 
#'     (taken to be the event or censoring time if \code{newdata} is not 
#'     specified, or the value of the \code{times} argument if \code{newdata} is 
#'     specified but no \code{last_time} is provided in the \code{control} 
#'     list, or otherwise the times provided in the \code{last_time} element
#'     of the \code{control} list).}
#'     \item{\code{last_time}}{a scalar or a character string 
#'     specifying the last known survival time for each individual for whom
#'     conditional predictions are being obtained. Should only be specified if
#'     \code{newdata} is provided, and conditional survival predictions are being
#'     obtained. A scalar will use the same last time for each individual in 
#'     \code{newdata}. A character string will name a column in \code{newdata}
#'     in which to look for the last times. If \code{last_time} is not provided then
#'     the default is to use the value provided in the \code{times} argument.} 
#'   }
#' @param ids An optional vector specifying a subset of IDs for whom the 
#'   predictions should be obtained. The default is to predict for all individuals
#'   who were used in estimating the model or, if \code{newdata} is specified,
#'   then all individuals contained in \code{newdata}.
#' @param prob A scalar between 0 and 1 specifying the width to use for the 
#'   uncertainty interval (sometimes called credible interval) for the predictions. 
#'   For example \code{prob = 0.95} (the default) means that the 2.5th and 97.5th  
#'   percentiles will be provided.
#' @param times A numeric vector of length 1 or a character string. Specifies the  
#'   times at which to obtain the estimated survival probabilities. 
#'   If \code{newdata} is \code{NULL}, then the 
#'   \code{times} argument is optional; if it is not provided then \code{times} 
#'   will default to the last known event or censoring time for each individual,
#'   whereas if it is provided then it must be a numeric vector of length 1, and 
#'   the survival probabilities will be calculated at the same \code{times} for 
#'   all individuals in the estimation dataset.
#'   If \code{newdata} is provided, then the \code{times} argument cannot be
#'   \code{NULL}, rather, the user must provide a numeric vector of length 1 or
#'   the name of a variable in \code{newdata}, indicating the times at which 
#'   the survival probabilities should be calculated for the individuals in 
#'   \code{newdata}. 
#' @param standardise A logical specifying whether the estimated 
#'   subject-specific survival probabilities should be averaged
#'   across all individuals for whom the subject-specific predictions are 
#'   being obtained. This can be used to average over the covariate distribution
#'   of the individuals used in estimating the model, or the individuals 
#'   included in \code{newdata}. This approach of
#'   averaging across the observed distribution of the covariates is sometimes
#'   referred to as a "standardised" survival curve. If \code{standardise = TRUE}, 
#'   then the \code{times} argument must be specified and it must be constant across 
#'   individuals, that is, the survival probabilities must be calculated at the 
#'   same time for all individuals.
#' @param draws An integer indicating the number of MCMC draws to return. The default
#'   and maximum number of draws is the size of the posterior sample.
#' @param seed An optional \code{\link[=set.seed]{seed}} to use.
#' @param ... Currently unused.
#'
#' @details 
#'   If the user wishes to obtain predictions using the \code{newdata}
#'   argument then several things need to be considered: \cr 
#'   \cr
#'   First, if you wish to obtain survival probabilities for "new" individuals, 
#'   meaning those who were \strong{not} part of the dataset used to estimate
#'   the model, then you will likely want to marginalise over the distribution 
#'   of the individual-level random effects.
#'   To ensure that this happens, you must ensure that the IDs provided in the 
#'   \code{id_var} column of \code{newdata} do 
#'   \strong{not} coincide with the IDs of individuals who were used in estimating  
#'   the model. Otherwise the predictions will be obtained using draws of the random  
#'   effects for the specific individual in the estimation data with the matching ID. 
#'   (Note that in the situation where you do want to obtain predictions for a given
#'   individual who was used in the estimation data but using new values for their 
#'   covariates, for example changing their treatment code or predicting at times 
#'   other than their actual observation times, then you could do this by specifying   
#'   the relevant individuals ID in the \code{id_var} columns of \code{newdata}. \cr
#'   \cr
#'   Second, if any variables were transformed (e.g. rescaled) in the data 
#'   used to fit the model, then these variables must also be transformed in 
#'   \code{newdata}. This only applies if variables  
#'   were transformed before passing the data to one of the modeling functions and  
#'   \emph{not} if transformations were specified inside the model formula. Also  
#'   see the \strong{Note} section in \code{\link{posterior_predict}} for a note  
#'   about using the \code{newdata} argument with binomial models.
#'    
#' @return A data frame of class \code{survfit.stanjm}. The data frame includes 
#'   columns for each of the following: 
#'   (i) the median of the posterior predictions of the estimated survival
#'   probabilities (\code{survfit});
#'   (ii) each of the lower and upper limits of the corresponding uncertainty 
#'   interval for the estimated survival probabilities (\code{ci_lb} and 
#'   \code{ci_ub});
#'   (iii) a subject identifier (\code{id_var}), unless standardised survival
#'   probabilities were estimated;
#'   (iv) the time that the estimated survival probability is calculated for 
#'   (\code{time_var}).
#'   The returned object also includes a number of additional attributes.
#' 
#' @seealso \code{\link{plot.survfit.stanjm}} for plotting the estimated survival  
#'   probabilities, \code{\link{ps_check}} for for graphical checks of the estimated 
#'   survival function, and \code{\link{posterior_traj}} for estimating the
#'   marginal or subject-specific longitudinal trajectories.
#'   
#' @examples
#' 
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
#'   head(pbcSurv[pbcSurv$status == 0,])
#'   ps1 <- posterior_survfit(example_jm, ids = c(7,13,16))
#'   head(ps1)
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
#'   ps2 <- posterior_survfit(example_jm, ids = c(7,13,16), times == 0,
#'     extrapolate = TRUE, control = list(edist = 5, condition = FALSE))
#'   ps2
#'   
#'   # Instead of estimating survival probabilities for a specific individual 
#'   # in the estimation dataset, we may want to estimate the marginal 
#'   # survival probability, that is, marginalising over the individual-level
#'   # random effects. 
#'   # Here we will estimate survival between baseline and 5 years, for a 
#'   # female who received either (i) D-penicillamine or (ii) placebo. 
#'   # To do this we will need to provide the necessary values  
#'   # of the predictors via the 'newdata' argument. However, it is important
#'   # to realise that by marginalising over the random effects 
#'   # distribution we will introduce a large amount of uncertainty into
#'   # the estimated survival probabilities. This is because we have no 
#'   # longitudinal measurements for these "new" individuals and therefore do
#'   # not have any specific information with which to estimate their random
#'   # effects. As such, there is a very wide 95% uncertainty interval 
#'   # associated with the estimated survival probabilities.
#'   nd <- data.frame(id = c("new1", "new2"),
#'                    sex = c("f", "f"), 
#'                    trt = c(1, 0))
#'   ps3 <- posterior_survfit(example_jm, newdata = nd, times = 0,
#'     extrapolate = TRUE, control = list(edist = 5, condition = FALSE))
#'   ps3
#'   
#'   # We can then plot the estimated survival functions to compare
#'   # them. To do this, we use the generic plot function.
#'   plot(ps3, limits = "none")                          
#'   
#'   # Lastly, if we wanted to obtain "standardised" survival probabilities, 
#'   # (by averaging over the observed distribution of the fixed effect 
#'   # covariates, as well as averaging over the estimated random effects
#'   # for individuals in our estimation sample) then we can specify
#'   # 'standardise = TRUE'. We can then plot the resulting standardised
#'   # survival curve.
#'   ps4 <- posterior_survfit(example_jm, standardise = TRUE, 
#'                            times = 0, extrapolate = TRUE)
#'   plot(ps4)                         
#' 
#'  
posterior_survfit <- function(object, newdata = NULL, extrapolate = TRUE, 
                              control = list(), prob = 0.95, ids,
                              times = NULL, standardise = FALSE, 
                              draws = NULL, seed = NULL, ...) {
  
  validate_stanjm_object(object)
  M        <- object$n_markers
  id_var   <- object$id_var
  time_var <- object$time_var
  if (!is.null(seed)) 
    set.seed(seed)
  if (missing(ids)) 
    ids <- NULL
  
  # Construct prediction data
  # ndL: dataLong to be used in predictions
  # ndE: dataEvent to be used in predictions
  newdata <- validate_newdata(newdata)
  if (is.null(newdata)) { # user did not specify newdata
    ndL <- model.frame(object)[1:M]
    ndE <- model.frame(object)$Event
  } else { # user specified newdata
    if (!id_var %in% colnames(newdata))
      stop("id_var from the original model call must appear 'newdata'.")
    if (any(duplicated(newdata[[id_var]])))
      stop("'newdata' should only contain one row per individual, since ",
           "time varying covariates are not allowed in the prediction data.")
    ndL <- rep(list(newdata), M)
    ndE <- newdata
  }
  
  # User specified a subset of ids
  if (!is.null(ids)) {
    check_for_missing_ids(ndE, id_var, ids)
    ndE <- ndE[ndE[[id_var]] %in% ids, , drop = FALSE]
    ndL <- lapply(ndL, function(x) x[x[[id_var]] %in% ids, , drop = FALSE])
  }  
  id_list <- unique(ndE[[id_var]]) # order of ids from data, not ids arg
  if (!is.null(newdata))
    check_for_estimation_ids(object, id_list) # warn if ids not in newdata
  
  # Prediction times
  if (standardise) { # standardised survival probs
    if (is.null(times)) {
      stop("'times' cannot be NULL for obtaining standardised survival probabilities.")
    } else if (is.numeric(times) && (length(times) == 1L)) {
      times <- rep(times, length(id_list))
    } else {
      stop("'times' should be a numeric vector of length 1 in order to obtain ",
           "standardised survival probabilities (the subject-specific survival ",
           "probabilities will be calculated at the specified time point, and ",
           "then averaged).")      
    }    
  } else if (is.null(newdata)) { # subject-specific survival probs without newdata
    if (is.null(times)) {
      times <- object$eventtime[as.character(id_list)]
    } else if (is.numeric(times) && (length(times) == 1L)) {
      times <- rep(times, length(id_list))
    } else {
      stop("If 'newdata' is NULL then 'times' can only be NULL or a ",
           "numeric vector of length 1.")     
    }
  } else { # subject-specific survival probs with newdata
    if (is.null(times)) {
      stop("'times' cannot be NULL if newdata is specified.")
    } else if (is.character(times) && (length(times) == 1L)) {
      if (!times %in% colnames(ndE))
        stop("Variable specified in 'times' argument could not be found in 'newdata'.")
      times <- tapply(ndE[[times]], ndE[[id_var]], FUN = max)
    } else if (is.numeric(times) && (length(times) == 1L)) {
      times <- rep(times, length(id_list))
    } else {
      stop("If 'newdata' is specified then 'times' can only be the name of a ",
           "variable in newdata, or a numeric vector of length 1.")      
    }
  }
  if (!identical(length(times), length(id_list)))
    stop(paste0("length of the 'times' vector should be equal to the number of individuals ",
                "for whom predictions are being obtained (", length(id_list), ")."))     
  maxtime <- max(object$eventtime)
  if (any(times > maxtime))
    stop("'times' are not allowed to be greater than the last event or censoring ",
         "time (since unable to extrapolate the baseline hazard).")
  
  # Last known survival time for each individual
  if (is.null(newdata)) { # user did not specify newdata
    if (!is.null(control$last_time))
      stop("'last_time' cannot be provided when 'newdata' is NULL, since times ",
           "are taken to be the event or censoring time for each individual.")
    last_time <- object$eventtime[as.character(id_list)]
  } else { # user specified newdata
    user_arg <- control$last_time # can be NULL
    if (is.null(user_arg)) {
      last_time <- times
    } else if (is.character(user_arg) && (length(user_arg) == 1L)) {
      if (!user_arg %in% colnames(ndE))
        stop("Cannot find 'last_time' column named in 'newdata'")
      last_time <- ndE[[user_arg]]      
    } else if (is.numeric(user_arg) && (length(user_arg) == 1L)) {
      last_time <- rep(user_arg, length(id_list)) 
    } else if (is.numeric(user_arg) && (length(user_arg) > 1L)) {
      last_time <- user_arg[as.character(id_list)]
    } else {
      stop("Bug found: could not reconcile last_time argument.")
    }
  }  
  
  # User specified extrapolation
  if (extrapolate) {
    ok_control_args <- c("epoints", "edist", "condition", "last_time")
    control <- get_extrapolation_control(control, ok_control_args = ok_control_args,
                                         standardise = standardise)
    endtime <- if (!is.null(control$edist)) times + control$edist else maxtime
    endtime[endtime > maxtime] <- maxtime # nothing beyond end of baseline hazard 
    time_seq <- get_time_seq(control$epoints, times, endtime, simplify = FALSE)
  } else time_seq <- list(times) # no extrapolation
  
  # Obtain survival prob matrix at each increment of time_seq
  # NB: length(time_seq) == 1 if no extrapolation
  ndL <- lapply(ndL, prepare_data_table, id_var = id_var, time_var = time_var)
  ndE <- prepare_data_table(ndE, id_var = id_var, time_var = time_var)
  surv <- lapply(time_seq, function(t, quadnodes) {  
    if (!identical(length(t), length(id_list))) # check length of times vector
      stop("Bug found: the vector of prediction times is not the same length ",
           "as the number of individuals.")
    # evaluate quadrature points and weights for t        
    quadpoints <- lapply(get_quadpoints(quadnodes)$points, unstandardise_quadpoints, 0, t)
    quadweights <- lapply(get_quadpoints(quadnodes)$weights, unstandardise_quadweights, 0, t)
    # evaluate design matrices for various longitudinal submodel contributions 
    # to the association structure (e.g. x_eta, Zt_eta, x_eps, Zt_eps, etc)
    y_X <- lapply(1:M, function(m)
      make_assoc_parts(ndL[[m]], assoc = object$assoc, id_var = id_var, 
                       time_var = time_var, id_list = id_list, times = quadpoints, 
                       use_function = pp_data, object = object,
                       m = m, re.form = NULL, offset = offset))
    # design matrix for event submodel at current increment  
    e_mf <- rolling_merge(ndE, ids = id_list, times = quadpoints)
    e_X <- .pp_data_mer_x(object, newdata = e_mf, re.form = NULL, 
                          offset = offset, m = "Event")       
    # matrix of survival probs at current increment  
    surv <- pp_survcalc(object, y_X = y_X, e_X = e_X, eventtime = t, 
                        quadpoints = quadpoints, quadweights = quadweights, 
                        draws = draws)
    # standardised survival probs
    if (standardise) { 
      surv <- matrix(rowMeans(surv), ncol = 1)
      dimnames(surv) <- list(iterations = NULL, "standardised_survprob") 
    } else {
      dimnames(surv) <- list(iterations = NULL, ids = id_list)
    }
    surv
  }, quadnodes = object$quadnodes)
  
  # If conditioning then also need to obtain surv probs at last known survival time
  if (extrapolate && control$condition) {
    # evaluate quadrature points and weights for last_time        
    quadpoints <- lapply(get_quadpoints(object$quadnodes)$points, unstandardise_quadpoints, 0, last_time)
    quadweights <- lapply(get_quadpoints(object$quadnodes)$weights, unstandardise_quadweights, 0, last_time)  
    # evaluate design matrices for various longitudinal submodel contributions 
    # to the association structure (e.g. x_eta, Zt_eta, x_eps, Zt_eps, etc)
    y_X <- lapply(1:M, function(m)
      make_assoc_parts(ndL[[m]], assoc = object$assoc, id_var = id_var, 
                       time_var = time_var, id_list = id_list, times = quadpoints, 
                       use_function = pp_data, object = object,
                       m = m, re.form = NULL, offset = offset))
    # design matrix for event submodel at last_time
    e_mf <- rolling_merge(ndE, ids = id_list, times = quadpoints)
    e_X <- .pp_data_mer_x(object, newdata = e_mf, re.form = NULL, 
                          offset = offset, m = "Event")       
    # matrix of survival probs at last_time 
    cond_surv <- pp_survcalc(object, y_X = y_X, e_X = e_X, eventtime = last_time, 
                             quadpoints = quadpoints, quadweights = quadweights, 
                             draws = draws)
    # conditional survival probs
    surv <- lapply(surv, function(x) {
      vec <- x / cond_surv
      vec[is.na(vec)] <- 1
      # If last_time was after the time of the estimated probability, 
      # leading to a survival probability greater than 1
      vec[vec > 1] <- 1  
      vec
    })        
  }
  
  # Summarise posterior draws to get median and ci
  out <- do.call("rbind", 
                 lapply(seq_along(surv), function(x, standardise, id_list, time_seq, prob) {
                   val <- median_and_bounds(surv[[x]], prob) 
                   cbind(IDVAR   = if (!standardise) id_list, 
                         TIMEVAR = if (!standardise) time_seq[[x]] else unique(time_seq[[x]]),
                         val$med, val$lb, val$ub)
                 }, standardise = standardise, id_list = id_list, time_seq = time_seq, prob = prob))
  rownames(out) <- NULL
  colnames(out) <- c(if ("IDVAR" %in% colnames(out)) id_var,
                     time_var, "survpred", "ci_lb", "ci_ub")
  out <- data.frame(out)
  if (id_var %in% colnames(out)) { # data has id column -- sort by id and time
    out <- out[order(out[, id_var, drop = F], out[, time_var, drop = F]), , drop = F]
  } else { # data does not have id column -- sort by time only
    out <- out[order(out[, time_var, drop = F]), , drop = F]
  }
  class(out) <- c("survfit.stanjm", "data.frame")
  structure(out, id_var = id_var, time_var = time_var, extrapolate = extrapolate, 
            control = control, standardise = standardise, ids = id_list, 
            draws = draws, seed = seed, offset = offset)
}

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
#' @return A \code{ggplot} object, also of class \code{plot.survfit.stanjm}.
#'   This object can be further customised using the \pkg{ggplot2} package.
#'   It can also be passed to the function \code{\link{plot_stack}}.
#'   
#' @seealso \code{\link{posterior_survfit}}, \code{\link{plot_stack}},
#'   \code{\link{posterior_traj}}, \code{\link{plot.predict.stanjm}}      
#'   
#' @examples 
#' 
#'   # Run example model if not already loaded
#'   if (!exists("example_jm")) example(example_jm)
#'   
#'   # Obtain subject-specific conditional survival probabilities
#'   # for all individuals in the estimation dataset.
#'   ps1 <- posterior_survfit(example_jm, extrapolate = TRUE)
#'   
#'   # We then plot the conditional survival probabilities for
#'   # a subset of individuals
#'   plot(ps1, ids = c(7,13,16))
#'   
#'   # We can change or add attributes to the plot
#'   plot(ps1, ids = c(7,13,16), limits = "none")
#'   plot(ps1, ids = c(7,13,16), xlab = "Follow up time")
#'   plot(ps1, ids = c(7,13,16), ci_geom_args = list(fill = "red"),
#'        color = "blue", linetype = 2)
#'   plot(ps1, ids = c(7,13,16), facet_scales = "fixed")
#'   
#'   # Since the returned plot is also a ggplot object, we can
#'   # modify some of its attributes after it has been returned
#'   plot1 <- plot(ps1, ids = c(7,13,16))
#'   plot1 + 
#'     ggplot2::theme(strip.background = ggplot2::element_blank()) +
#'     ggplot2::coord_cartesian(xlim = c(0, 15)) +
#'     ggplot2::labs(title = "Some plotted survival functions")
#'     
#'   # We can also combine the plot(s) of the estimated 
#'   # subject-specific survival functions, with plot(s) 
#'   # of the estimated longitudinal trajectories for the
#'   # same individuals
#'   ps1 <- posterior_survfit(example_jm, ids = c(7,13,16))
#'   pt1 <- posterior_traj(example_jm, , ids = c(7,13,16))
#'   plot_surv <- plot(ps1) 
#'   plot_traj <- plot(pt1, vline = TRUE, plot_observed = TRUE)
#'   plot_stack(plot_traj, plot_surv)
#'    
#'   # Lastly, let us plot the standardised survival function
#'   # based on all individuals in our estimation dataset
#'   ps2 <- posterior_survfit(example_jm, standardise = TRUE, times = 0,
#'                           control = list(epoints = 20))
#'   plot(ps2)   
#' 
#'    
plot.survfit.stanjm <- function(x, ids = NULL, 
                                limits = c("ci", "none"),  
                                xlab = NULL, ylab = NULL, facet_scales = "free", 
                                ci_geom_args = NULL, ...) {
  
  limits <- match.arg(limits)
  ci <- (limits == "ci")
  standardise <- attr(x, "standardise")
  id_var <- attr(x, "id_var")
  time_var <- attr(x, "time_var")
  if (is.null(xlab)) xlab <- paste0("Time (", time_var, ")")
  if (is.null(ylab)) ylab <- "Event free probability"
  if (!is.null(ids)) {
    if (standardise) 
      stop("'ids' argument cannot be specified when plotting standardised ",
           "survival probabilities.")
    if (!id_var %in% colnames(x))
      stop("Bug found: could not find 'id_var' column in the data frame.")
    ids_missing <- which(!ids %in% x[[id_var]])
    if (length(ids_missing))
      stop("The following 'ids' are not present in the survfit.stanjm object: ",
           paste(ids[[ids_missing]], collapse = ", "), call. = FALSE)
    x <- x[(x[[id_var]] %in% ids), , drop = FALSE]
  } else {
    ids <- if (!standardise) attr(x, "ids") else NULL
  }
  if (!standardise) x$id <- x[[id_var]]
  x$time <- x[[time_var]]
  
  geom_defaults <- list(color = "black")
  geom_args <- set_geom_args(geom_defaults, ...)  
  
  lim_defaults <- list(alpha = 0.3)
  lim_args <- do.call("set_geom_args", c(defaults = list(lim_defaults), ci_geom_args))
  
  if ((!standardise) && (length(ids) > 60L)) {
    stop("Too many individuals to plot for. Perhaps consider limiting ",
         "the number of individuals by specifying the 'ids' argument.")
  } else if ((!standardise) && (length(ids) > 1L)) {
    graph <- ggplot(x, aes_string(x = "time", y = "survpred")) +
      theme_bw() +
      do.call("geom_line", geom_args) +
      coord_cartesian(ylim = c(0, 1)) +      
      facet_wrap(~ id, scales = facet_scales)
    if (ci) {
      lim_mapp <- list(mapping = aes_string(ymin = "ci_lb", ymax = "ci_ub"))
      graph_limits <- do.call("geom_ribbon", c(lim_mapp, lim_args))
    } else graph_limits <- NULL
  } else {
    graph <- ggplot(x, aes_string(x = "time", y = "survpred")) + 
      theme_bw() +
      do.call("geom_line", geom_args) + 
      coord_cartesian(ylim = c(0, 1))
    if (ci) {
      lim_mapp <- list(mapping = aes_string(ymin = "ci_lb", ymax = "ci_ub"))
      graph_limits <- do.call("geom_ribbon", c(lim_mapp, lim_args))
    } else graph_limits <- NULL
    
  }    
  
  ret <- graph + graph_limits + labs(x = xlab, y = ylab) 
  class_ret <- class(ret)
  class(ret) <- c("plot.survfit.stanjm", class_ret)
  ret
}


# internal ----------------------------------------------------------------

# Return a data.table with the key set using the appropriate time variable
# 
# @param data A data frame
# @param id_var The name of the ID variable 
# @param time_var The name of the time variable
# @return A data.table (which will be used in a rolling merge against the
#   event times and/or quadrature times)
prepare_data_table <- function(data, id_var, time_var) {
  if (survival::is.Surv(data[[1]])) { # event submodel model.frame from fitted JM
    # If the design matrix is for the event submodel and data is model.frame obtained 
    # from the fitted model (ie, not supplied by the user) then the time point
    # for merging on covariate values is taken to be either:
    # (i) the unique observation time (single row per individual surv data), or 
    # (ii) "start" of the start/stop interval (multiple row per individual surv data)
    resp_type <- attr(data[[1]], "type")
    data <- cbind(unclass(data[[1]]), data[,-1])
    if (resp_type == "right") { # single row data
      data <- data.table::data.table(data, key = c(id_var, "time"))
      data[["time"]] <- as.numeric(data[["time"]])
    } else if (resp_type == "counting") { # start/stop multiple row data
      data <- data.table::data.table(data, key = c(id_var, "start"))
      data[["start"]] <- as.numeric(data[["start"]])
    } else {
      stop("Bug found: 'data' arg appears to be the model.frame from the fitted model, ",
           "but cannot find an appropriate time variable in the Surv(.) response.")
    }
  } else { # user provided new data, or the data is for the long. submodel
    # Alternatively, the user provided the new data for the event submodel 
    # which must be single row per individual (since multiple row per individual
    # data is not allowed by posterior_survfit), or, the data is for a 
    # longitudinal submodel. In either case, the time_var variable is used
    # for merging on covariate values -- if time_var doesn't already exist in
    # the single row per individual data then we create a dud time_var variable.
    if (!time_var %in% colnames(data)) 
      data[[time_var]] <- rep(0.0, nrow(data))  
    data <- data.table::data.table(data, key = c(id_var, time_var))
  }
  return(data)
}

# Create matrix of posterior survival probabilities at time t
# 
# @param object A stanjm object
# @param pp_dataLong A list equal in length to the number of markers.
#   Each element contains a named list with elements $mod_eta, $mod_eps,
#   $mod_lag, $mod_auc, which each contain output returned by pp_data()
#   namely the X and Zt matrices evaluated at the relevant quadpoints.
#   The list also includes elements $xmat_data and $K_data.
# @param pp_dataEvent A list returned by a call to .pp_data_mer_X containing
#   the design matrix for the event submodel evaluated at the quadpoints.
# @param draws Integer specifying the number of draws
# @return A matrix of survival probabilities at time t, with the S
#   rows corresponding to different MCMC draws of the parameters 
#   from the posterior, and each column corresponding to one of the Npat
#   individuals in the new data
pp_survcalc <- function(object, y_X, e_X, eventtime, quadpoints, 
                        quadweights, draws = NULL) {

  # Get stanmat for specified number of draws
  S <- posterior_sample_size(object)
  if (is.null(draws)) 
    draws <- S
  if (draws > S) {
    err <- paste0("'draws' should be <= posterior sample size (", S, ").")
    stop(err)
  }
  stanmat <- as.matrix(object$stanfit)
  some_draws <- isTRUE(draws < S)
  if (some_draws) {
    samp <- sample(S, draws)
    stanmat <- stanmat[samp, , drop = FALSE]
  }

  # Extract parameters for constructing linear predictor  
  M <- object$n_markers
  nms <- collect_nms(colnames(stanmat), M)
  b_ord  <- lapply(1:M, function(m) y_X[[m]]$mod_eta$Z_names)
  y_b    <- lapply(1:M, function(m) pp_b_ord(stanmat[, nms$y_b[[m]], drop = FALSE], b_ord[[m]]))
  y_beta <- lapply(1:M, function(m) stanmat[, nms$y[[m]], drop = FALSE])
  e_beta <- stanmat[, nms$e, drop = FALSE]
  a_beta <- stanmat[, nms$a, drop = FALSE]
  
  # Linear predictor for the event submodel
  e_eta <- linear_predictor.matrix(e_beta, e_X, offset = NULL)
  assoc <- object$assoc
  sel <- grep("which|null", rownames(assoc), invert = TRUE)
  if (any(unlist(assoc[sel,]))) { # has association structure
    if (any(unlist(assoc["shared_b|shared_coef",])))
      stop("posterior_survfit cannot yet be used with shared_b or shared_coef ",
           "association structures.") # until make_assoc_terms can handle it
    for (s in seq(draws)) {
      beta_s <- lapply(y_beta, function(x) x[s,])
      b_s    <- lapply(y_b,    function(x) x[s,])
      a_Xs <- make_assoc_terms(parts = y_X, assoc = assoc, family = object$family, 
                               beta = beta_s, b = b_s)
      e_eta[s,] <- e_eta[s,] + linear_predictor.default(a_beta[s,], a_Xs)
    }
  }
  
  # Baseline hazard
  if (object$basehaz$type_name == "weibull") {
    shape <- stanmat[, nms$e_extra, drop = FALSE]
    log_basehaz <- as.vector(log(shape)) + 
      (shape - 1) %*% matrix(log(unlist(quadpoints)), nrow = 1)
  } else if (object$basehaz$type_name == "bs") {
    coefs <- stanmat[, nms$e_extra, drop = FALSE]
    log_basehaz <- 
      coefs %*% t(predict(object$basehaz$bs_basis, unlist(quadpoints)))
  } else {
    stop("posterior_survfit not yet implemented for basehaz = ", 
         object$basehaz$type_name)
  }
  
  # Evaluate cumulative hazard up to time t, and survival prob at time t
  haz <- exp(log_basehaz + e_eta)
  qweighted_haz <- t(apply(haz, 1L, function(row) unlist(quadweights) * row))
  cum_haz <- Reduce('+', array2list(qweighted_haz, nsplits = length(quadweights)))  
  surv_t <- exp(-cum_haz) # should be S * Npat matrix
  if (is.vector(surv_t) == 1L) # transform if only one individual
    surv_t <- t(surv_t)
  # set survprob matrix at time 0 to S(t) = 1 
  # (otherwise some NaN possible due to numerical inaccuracies)
  surv_t[, (eventtime == 0)] <- 1
  return(surv_t) # returns S x Npat matrix of survival probabilities at t
}


# default plotting attributes
.PP_FILL <- "skyblue"
.PP_DARK <- "skyblue4"
.PP_VLINE_CLR <- "#222222"
.PP_YREP_CLR <- "#487575"
.PP_YREP_FILL <- "#222222"




