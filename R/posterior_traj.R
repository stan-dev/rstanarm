# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2013, 2014, 2015, 2016, 2017 Trustees of Columbia University
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

#' Estimate the subject-specific or marginal longitudinal trajectory 
#' 
#' This function allows us to generate an estimated longitudinal trajectory 
#' (either subject-specific, or by marginalising over the distribution of the 
#' group-specific parameters) based on draws from the posterior predictive 
#' distribution. 
#' 
#' @export
#' 
#' @templateVar stanjmArg object
#' @templateVar mArg m
#' @template args-stanjm-object
#' @template args-m
#' @param newdata \strong{Deprecated}: please use \code{newdataLong} instead.
#'   Optionally, a data frame in which to look for variables with
#'   which to predict. If omitted, the model matrix is used. If \code{newdata}
#'   is provided and any variables were transformed (e.g. rescaled) in the data
#'   used to fit the model, then these variables must also be transformed in
#'   \code{newdata}. This only applies if variables were transformed before
#'   passing the data to one of the modeling functions and \emph{not} if
#'   transformations were specified inside the model formula.
#' @param newdataLong,newdataEvent Optionally, a data frame (or in the case of 
#'   \code{newdataLong} this can be a list of data frames) in which to look 
#'   for variables with which to predict. If omitted, the model matrices are used. 
#'   If new data is provided, then two options are available. Either one can
#'   provide observed covariate and outcome data, collected up to some time 
#'   \emph{t}, and use this data to draw new individual-specific coefficients
#'   (i.e. individual-level random effects). This is the default behaviour when 
#'   new data is provided, determined by the argument \code{dynamic = TRUE}, and
#'   requiring both \code{newdataLong} and \code{newdataEvent} to be specified.
#'   Alternatively, one can specify \code{dynamic = FALSE}, and then predict 
#'   using just covariate data, by marginalising over the distribution
#'   of the group-specific coefficients; in this case, only \code{newdataLong}
#'   needs to be specified and it only needs to be a single data frame with
#'   the covariate data for the predictions for the one longitudinal submodel.
#' @param interpolate A logical specifying whether to interpolate the estimated 
#'   longitudinal trajectory in between the observation times. This can be used
#'   to achieve a smooth estimate of the longitudinal trajectory across the 
#'   entire follow up time. If \code{TRUE} then the interpolation can be further 
#'   controlled using the \code{control} argument.
#' @param extrapolate A logical specifying whether to extrapolate the estimated 
#'   longitudinal trajectory beyond the time of the last known observation time.
#'   If \code{TRUE} then the extrapolation can be further controlled using
#'   the \code{control} argument.
#' @param control A named list with parameters controlling the interpolation or
#'   extrapolation of the estimated longitudinal trajectory when either 
#'   \code{interpolate = TRUE} or \code{extrapolate = TRUE}. The 
#'   list can contain one or more of the following named elements: \cr
#'   \describe{
#'     \item{\code{ipoints}}{a positive integer specifying the number of discrete 
#'     time points at which to calculate the estimated longitudinal response for
#'     \code{interpolate = TRUE}. These time points are evenly spaced starting at 
#'     0 and ending at the last known observation time for each individual. The
#'     last observation time for each individual is taken to be either: the
#'     event or censoring time if no new data is provided; the time specified
#'     in the "last_time" column if provided in the new data (see \strong{Details}
#'     section below); or the time of the last longitudinal measurement if new
#'     data is provided but no "last_time" column is included. The default is 15.}
#'     \item{\code{epoints}}{a positive integer specifying the number of discrete 
#'     time points at which to calculate the estimated longitudinal response for
#'     \code{extrapolate = TRUE}. These time points are evenly spaced between the 
#'     last known observation time for each individual and the extrapolation 
#'     distance specifed using either \code{edist} or \code{eprop}.
#'     The default is 15.}
#'     \item{\code{eprop}}{a positive scalar between 0 and 1 specifying the 
#'     amount of time across which to extrapolate the longitudinal trajectory,
#'     represented as a proportion of the total observed follow up time for each
#'     individual. For example specifying \code{eprop = 0.2} means that for an
#'     individual for whom the latest of their measurement, event or censoring times
#'     was 10 years, their estimated longitudinal trajectory will be extrapolated 
#'     out to 12 years (i.e. 10 + (0.2 * 10)). The default value is 0.2.}
#'     \item{\code{edist}}{a positive scalar specifying the amount of time 
#'     across which to extrapolate the longitudinal trajectory for each individual,
#'     represented in units of the time variable \code{time_var} (from fitting the
#'     model). This cannot be specified if \code{eprop} is specified.} 
#'   }
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
#' @param ids An optional vector specifying a subset of subject IDs for whom the 
#'   predictions should be obtained. The default is to predict for all individuals
#'   who were used in estimating the model or, if \code{newdata} is specified,
#'   then all individuals contained in \code{newdata}.
#' @param prob A scalar between 0 and 1 specifying the width to use for the 
#'   uncertainty interval (sometimes called credible interval) for the predicted
#'   mean response and the prediction interval for the predicted (raw) response. 
#'   For example \code{prob = 0.95} (the default) means that the 2.5th and 97.5th  
#'   percentiles will be provided. Only relevant when \code{return_matrix} is 
#'   \code{FALSE}. 
#' @param dynamic A logical that is only relevant if new data is provided
#'   via the \code{newdata} argument. If 
#'   \code{dynamic = TRUE}, then new group-specific parameters are drawn for 
#'   the individuals in the new data, conditional on their longitudinal 
#'   biomarker data contained in \code{newdata}. These group-specific
#'   parameters are then used to generate individual-specific survival probabilities
#'   for these individuals. These are often referred to as "dynamic predictions"
#'   in the joint modelling context, because the predictions can be updated
#'   each time additional longitudinal biomarker data is collected on the individual.
#'   On the other hand, if \code{dynamic = FALSE} then the survival probabilities
#'   will just be marginalised over the distribution of the group-specific
#'   coefficients; this will mean that the predictions will incorporate all
#'   uncertainty due to between-individual variation so there will likely be
#'   very wide credible intervals on the predicted survival probabilities.
#' @param scale A scalar, specifying how much to multiply the asymptotic 
#'   variance-covariance matrix for the random effects by, which is then
#'   used as the "width" (ie. variance-covariance matrix) of the multivariate
#'   Student-t proposal distribution in the Metropolis-Hastings algorithm. This
#'   is only relevant when \code{newdataEvent} is supplied and 
#'   \code{dynamic = TRUE}, in which case new random effects are simulated
#'   for the individuals in the new data using the Metropolis-Hastings algorithm.
#' @param draws An integer indicating the number of MCMC draws to return. 
#'   The default is to set the number of draws equal to 200, or equal to the 
#'   size of the posterior sample if that is less than 200. 
#' @param seed An optional \code{\link[=set.seed]{seed}} to use.
#' @param return_matrix A logical. If \code{TRUE} then a \code{draws} by 
#'   \code{nrow(newdata)} matrix is returned which contains all the actual
#'   simulations or draws from the posterior predictive distribution. Otherwise
#'   if \code{return_matrix} is set to \code{FALSE} (the default) then a 
#'   data frame is returned, as described in the \strong{Value} section below.
#' @param ... Other arguments passed to \code{\link{posterior_predict}}, for
#'   example \code{draws}, \code{re.form}, \code{seed}, etc.
#'   
#' @details The \code{posterior_traj} function acts as a wrapper to the 
#' \code{\link{posterior_predict}} function, but allows predictions to be 
#' easily generated at time points that are interpolated and/or extrapolated 
#' between time zero (baseline) and the last known survival time for the 
#' individual, thereby providing predictions that correspond to a smooth estimate
#' of the longitudinal trajectory (useful for the plotting via the associated
#' \code{\link{plot.predict.stanjm}} method). In addition it returns a data 
#' frame by default, whereas the \code{\link{posterior_predict}} function 
#' returns a matrix; see the \strong{Value} section below for details. Also,
#' \code{posterior_traj} allows predictions to only be generated for a subset
#' of individuals, via the \code{ids} argument. 
#' 
#' @return When \code{return_matrix = FALSE}, a data frame 
#'   of class \code{predict.stanjm}. The data frame includes a column for the median 
#'   of the posterior predictions of the mean longitudinal response (\code{yfit}),
#'   a column for each of the lower and upper limits of the uncertainty interval
#'   corresponding to the posterior predictions of the mean longitudinal response 
#'   (\code{ci_lb} and \code{ci_ub}), and a column for each of the lower and upper
#'   limits of the prediction interval corresponding to the posterior predictions
#'   of the (raw) longitudinal response. The data frame also includes columns for
#'   the subject ID variable, and each of the predictor variables. The returned
#'   object also includes a number of attributes.
#'   
#'   When \code{return_matrix = TRUE}, the returned object is the same as that 
#'   described for \code{\link{posterior_predict}}.
#'   
#' @seealso \code{\link{plot.predict.stanjm}}, \code{\link{posterior_predict}},
#'   \code{\link{posterior_survfit}}
#' 
#' @examples
#' \donttest{
#'   # Run example model if not already loaded
#'   if (!exists("example_jm")) example(example_jm)
#'   
#'   # Obtain subject-specific predictions for all individuals 
#'   # in the estimation dataset
#'   pt1 <- posterior_traj(example_jm, interpolate = FALSE, extrapolate = FALSE)
#'   head(pt1)
#'   
#'   # Obtain subject-specific predictions only for a few selected individuals
#'   pt2 <- posterior_traj(example_jm, ids = c(1,3,8))
#'   
#'   # If we wanted to obtain subject-specific predictions in order to plot the 
#'   # longitudinal trajectories, then we might want to ensure a full trajectory 
#'   # is obtained by interpolating and extrapolating time. We can then use the 
#'   # generic plot function to plot the subject-specific predicted trajectories
#'   # for the first three individuals. Interpolation and extrapolation is 
#'   # carried out by default.
#'   pt3 <- posterior_traj(example_jm)
#'   head(pt3) # predictions at additional time points compared with pt1 
#'   plot(pt3, ids = 1:3)
#'   
#'   # If we wanted to extrapolate further in time, but decrease the number of 
#'   # discrete time points at which we obtain predictions for each individual, 
#'   # then we could specify a named list in the 'control' argument
#'   pt4 <- posterior_traj(example_jm, control = list(ipoints = 10, epoints = 10, eprop = 0.5))
#'   
#'   # If we have prediction data for a new individual, and we want to
#'   # estimate the longitudinal trajectory for that individual conditional
#'   # on this new data (perhaps extrapolating forward from our last
#'   # longitudinal measurement) then we can do that. It requires drawing
#'   # new individual-specific parameters, based on the full likelihood,
#'   # so we must supply new data for both the longitudinal and event 
#'   # submodels. These are sometimes known as dynamic predictions.
#'   ndL <- pbcLong[pbcLong$id == 8, , drop = FALSE]
#'   ndE <- pbcSurv[pbcSurv$id == 8, , drop = FALSE]
#'   ndL$id <- "new_subject" # new id can't match one used in training data
#'   ndE$id <- "new_subject"
#'   pt5 <- posterior_traj(example_jm, 
#'                         newdataLong = ndL,
#'                         newdataEvent = ndE)
#'                         
#'   # By default it is assumed that the last known survival time for 
#'   # the individual is the time of their last biomarker measurement,
#'   # but if we know they survived to some later time then we can
#'   # condition on that information using the last_time argument
#'   pt6 <- posterior_traj(example_jm, 
#'                         newdataLong = ndL,
#'                         newdataEvent = ndE, 
#'                         last_time = "futimeYears")
#'   
#'   # Alternatively we may want to estimate the marginal longitudinal
#'   # trajectory for a given set of covariates. To do this, we can pass
#'   # the desired covariate values in a new data frame (however the only
#'   # covariate in our fitted model was the time variable, year). To make sure  
#'   # that we marginalise over the random effects, we need to specify an ID value
#'   # which does not correspond to any of the individuals who were used in the
#'   # model estimation and specify the argument dynamic=FALSE.
#'   # The marginal prediction is obtained by generating subject-specific 
#'   # predictions using a series of random draws from the random 
#'   # effects distribution, and then integrating (ie, averaging) over these. 
#'   # Our marginal prediction will therefore capture the between-individual 
#'   # variation associated with the random effects.
#'   nd <- data.frame(id = rep("new1", 11), year = (0:10 / 2))
#'   pt7 <- posterior_traj(example_jm, newdataLong = nd, dynamic = FALSE)
#'   head(pt7)  # note the greater width of the uncertainty interval compared 
#'              # with the subject-specific predictions in pt1, pt2, etc
#'   
#'   # Alternatively, we could have estimated the "marginal" trajectory by 
#'   # ignoring the random effects (ie, assuming the random effects were set 
#'   # to zero). This will generate a predicted longitudinal trajectory only
#'   # based on the fixed effect component of the model. In essence, for a 
#'   # linear mixed effects model (ie, a model that uses an identity link 
#'   # function), we should obtain a similar point estimate ("yfit") to the
#'   # estimates obtained in pt5 (since the mean of the estimated random effects
#'   # distribution will be approximately 0). However, it is important to note that
#'   # the uncertainty interval will be much more narrow, since it completely
#'   # ignores the between-individual variability captured by the random effects.
#'   # Further, if the model uses a non-identity link function, then the point
#'   # estimate ("yfit") obtained only using the fixed effect component of the
#'   # model will actually provide a biased estimate of the marginal prediction.
#'   # Nonetheless, to demonstrate how we can obtain the predictions only using 
#'   # the fixed effect component of the model, we simply specify 're.form = NA'. 
#'   # (We will use the same covariate values as used in the prediction for 
#'   # example for pt5).
#'   pt8 <- posterior_traj(example_jm, newdataLong = nd, dynamic = FALSE, 
#'                         re.form = NA)
#'   head(pt8)  # note the much narrower ci, compared with pt5
#' }
#' 
posterior_traj <- function(object, m = 1, newdata = NULL, newdataLong = NULL,
                           newdataEvent = NULL, interpolate = TRUE, extrapolate = FALSE,
                           control = list(), last_time = NULL, prob = 0.95, ids, 
                           dynamic = TRUE, scale = 1.5, draws = NULL, seed = NULL, 
                           return_matrix = FALSE, ...) {
  if (!requireNamespace("data.table"))
    stop("the 'data.table' package must be installed to use this function")
  validate_stanjm_object(object)
  M         <- object$n_markers; validate_positive_scalar(m, M)
  id_var    <- object$id_var
  time_var  <- object$time_var
  grp_stuff <- object$grp_stuff[[m]]
  if (!is.null(seed)) 
    set.seed(seed)
  if (missing(ids)) 
    ids <- NULL
  dots <- list(...)
  
  # Deal with deprecate newdata argument
  if (!is.null(newdata)) {
    warning("The 'newdata' argument is deprecated. Use 'newdataLong' instead.")
    if (!is.null(newdataLong))
      stop2("'newdata' and 'newdataLong' cannot both be specified.")
    newdataLong <- newdata
  }
  
  # Construct prediction data, NB dats == observed data to return to user
  if (is.null(newdataLong)) { # user did not specify newdata
    if (!is.null(newdataEvent))
      stop2("'newdataEvent' can only be specified when 'newdataLong' is provided.")
    dats <- get_model_data(object)
    ndL <- dats[1:M]
    ndE <- dats[["Event"]]
  } else { # user specified newdataLong
    if (dynamic && is.null(newdataEvent))
      stop2("Dynamic predictions require both 'newdataLong' and 'newdataEvent' ",
            "to be specified. Either specify data for both the longitudinal and ",
            "event submodels or, alternatively, specify argument 'dynamic = FALSE' ",
            "to marginalise over the distribution of group-specific parameters.")
    dats <- validate_newdatas(object, newdataLong, newdataEvent, 
                              response = isTRUE(dynamic))
    ndL <- dats[1:M]
    ndE <- dats[["Event"]]
  }
  if (!is.null(ids)) { # user specified a subset of ids
    ndL <- subset_ids(object, ndL, ids)
    ndE <- subset_ids(object, ndE, ids)
  }  
  id_list <- factor(unique(ndL[[m]][[id_var]])) # order of ids from data, not ids arg
  
  # Last known survival time for each individual
  if (is.null(newdataLong)) { # user did not provide newdata
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
      if (!is.null(ndE)) { # user provided newdataEvent for dynamic predictions
        if (!last_time %in% colnames(ndE))
          stop("Cannot find 'last_time' column named in newdataEvent.")
        last_time <- ndE[[last_time]]      
      }
    } else if (is.numeric(last_time) && (length(last_time) == 1L)) {
      last_time <- rep(last_time, length(id_list)) 
    } else if (is.numeric(last_time) && (length(last_time) > 1L)) {
      last_time <- last_time[as.character(id_list)]
    } else {
      stop("Bug found: could not reconcile 'last_time' argument.")
    }
    names(last_time) <- as.character(id_list)
  }
  
  # Get stanmat parameter matrix for specified number of draws
  S <- posterior_sample_size(object)
  if (is.null(draws)) 
    draws <- if (S > 200L) 200L else S 
  if (draws > S)
    stop("'draws' should be <= posterior sample size (", S, ").")
  stanmat <- as.matrix(object$stanfit)
  some_draws <- isTRUE(draws < S)
  if (some_draws) {
    samp <- sample(S, draws)
    stanmat <- stanmat[samp, , drop = FALSE]
  }
  
  # Draw b pars for new individuals
  if (dynamic && !is.null(newdataEvent)) {
    stanmat <- simulate_b_pars(object, stanmat = stanmat, ndL = ndL, ndE = ndE,
                               ids = id_list, times = last_time, scale = scale)
    b_new <- attr(stanmat, "b_new")
    acceptance_rate <- attr(stanmat, "acceptance_rate")
  }
  
  newX <- ndL[[m]] # design matrix used for predictions
  if (interpolate || extrapolate) { # user specified interpolation or extrapolation
    if (return_matrix) 
      stop("'return_matrix' cannot be TRUE if 'interpolate' or 'extrapolate' is TRUE.")
    ok_control_args <- c("ipoints", "epoints", "edist", "eprop")
    control <- get_extrapolation_control(control, ok_control_args = ok_control_args)
    dist <- if (!is.null(control$eprop)) control$eprop * (last_time - 0) else control$edist
    iseq <- if (interpolate) get_time_seq(control$ipoints, 0, last_time) else NULL
    eseq <- if (extrapolate) get_time_seq(control$epoints, last_time, last_time + dist) else NULL
    time_seq <- as.data.frame(cbind(iseq, eseq))
    colnames(time_seq) <- paste0("V", 1:NCOL(time_seq))
    time_seq <- reshape(time_seq, direction = "long", varying = colnames(time_seq), 
                        v.names = time_var, timevar = "obs", ids = id_list, idvar = id_var)
    newX[[time_var]] <- as.numeric(newX[[time_var]]) # ensures no rounding during data.table merge
    if (grp_stuff$has_grp) {
      grp_var <- grp_stuff$grp_var
      time_seq <- merge(time_seq, unique(newX[, c(id_var, grp_var)]), by = id_var)
      time_seq <- time_seq[order(time_seq[["obs"]], time_seq[[id_var]], time_seq[[grp_var]]), ]
      newX <- prepare_data_table(newX, id_var = id_var, time_var = time_var, grp_var = grp_var)
      newX <- rolling_merge(newX, time_seq[[id_var]], time_seq[[time_var]], time_seq[[grp_var]])
    } else {
      newX <- prepare_data_table(newX, id_var = id_var, time_var = time_var)
      newX <- rolling_merge(newX, time_seq[[id_var]], time_seq[[time_var]])
    }
  }
    
  ytilde <- posterior_predict(object, newdata = newX, m = m, stanmat = stanmat, ...)
  if (return_matrix) {
    attr(ytilde, "mu") <- NULL # remove attribute mu
    return(ytilde) # return S * N matrix, instead of data frame 
  } 
  mutilde <- attr(ytilde, "mu")
  if (!is.null(newX) && nrow(newX) == 1L) 
    mutilde <- t(mutilde)
  ytilde_bounds  <- median_and_bounds(ytilde,  prob) # median and prob% CrI limits
  mutilde_bounds <- median_and_bounds(mutilde, prob) # median and prob% CrI limits
  out <- data.frame(IDVAR = newX[[id_var]], 
                    TIMEVAR = newX[[time_var]], 
                    yfit = mutilde_bounds$med,
                    ci_lb = mutilde_bounds$lb, ci_ub = mutilde_bounds$ub,
                    pi_lb = ytilde_bounds$lb,  pi_ub = ytilde_bounds$ub)
  if (grp_stuff$has_grp) {
    out$GRPVAR = newX[[grp_var]] # add grp_var and reorder cols
    out <- out[, c("IDVAR", "GRPVAR", "TIMEVAR", 
                   "yfit", "ci_lb", "ci_ub", "pi_lb", "pi_ub")]
  }
  colnames(out) <- c(id_var, if (grp_stuff$has_grp) grp_var, time_var, 
                     "yfit", "ci_lb", "ci_ub", "pi_lb", "pi_ub")
  class(out) <- c("predict.stanjm", "data.frame")
  Terms <- terms(formula(object, m = m))
  vars  <- rownames(attr(Terms, "factors"))
  y_var <- vars[[attr(Terms, "response")]]
  out <- structure(out, observed_data = ndL[[m]], last_time = last_time,
                   y_var = y_var, id_var = id_var, time_var = time_var,
                   grp_var = if (grp_stuff$has_grp) grp_var else NULL, 
                   interpolate = interpolate, extrapolate = extrapolate, 
                   control = control, call = match.call())  
  if (dynamic && !is.null(newdataEvent)) {
    out <- structure(out, b_new = b_new, acceptance_rate = acceptance_rate)
  }
  out  
}


#' Plot the estimated subject-specific or marginal longitudinal trajectory
#' 
#' This generic \code{plot} method for \code{predict.stanjm} objects will
#' plot the estimated subject-specific or marginal longitudinal trajectory
#' using the data frame returned by a call to \code{\link{posterior_traj}}.
#' To ensure that enough data points are available to plot the longitudinal
#' trajectory, it is assumed that the call to \code{\link{posterior_traj}}
#' would have used the default \code{interpolate = TRUE}, and perhaps also 
#' \code{extrapolate = TRUE} (the latter being optional, depending on 
#' whether or not the user wants to see extrapolation of the longitudinal 
#' trajectory beyond the last observation time).
#' 
#' @method plot predict.stanjm
#' @export
#' @importFrom ggplot2 ggplot aes aes_string geom_line geom_smooth geom_ribbon 
#'   geom_point facet_wrap geom_vline labs ggplot_build theme_bw
#'    
#' @templateVar labsArg xlab,ylab
#' @templateVar scalesArg facet_scales
#' @templateVar cigeomArg ci_geom_args
#' @template args-ids
#' @template args-labs
#' @template args-scales
#' @template args-ci-geom-args
#' 
#' @param x A data frame and object of class \code{predict.stanjm}
#'   returned by a call to the function \code{\link{posterior_traj}}.
#'   The object contains point estimates and uncertainty interval limits
#'   for the fitted values of the longitudinal response.
#' @param limits A quoted character string specifying the type of limits to
#'   include in the plot. Can be one of: \code{"ci"} for the Bayesian
#'   posterior uncertainty interval for the estimated mean longitudinal
#'   response (often known as a credible interval);
#'   \code{"pi"} for the prediction interval for the estimated (raw)
#'   longitudinal response; or \code{"none"} for no interval limits.
#' @param vline A logical. If \code{TRUE} then a vertical dashed line
#'   is added to the plot indicating the event or censoring time for
#'   the individual. Can only be used if each plot within the figure
#'   is for a single individual.
#' @param plot_observed A logical. If \code{TRUE} then the observed
#'   longitudinal measurements are overlaid on the plot.
#' @param grp_overlay Only relevant if the model had lower level units 
#'   clustered within an individual. If \code{TRUE}, then the fitted trajectories 
#'   for the lower level units will be overlaid in the same plot region (that 
#'   is, all lower level units for a single individual will be shown within a 
#'   single facet). If \code{FALSE}, then the fitted trajectories for each lower
#'   level unit will be shown in a separate facet.
#' @param ... Optional arguments passed to 
#'   \code{\link[ggplot2]{geom_smooth}} and used to control features
#'   of the plotted longitudinal trajectory.
#'   
#' @return A \code{ggplot} object, also of class \code{plot.predict.stanjm}.
#'   This object can be further customised using the \pkg{ggplot2} package.
#'   It can also be passed to the function \code{\link{plot_stack_jm}}.
#'   
#' @seealso \code{\link{posterior_traj}}, \code{\link{plot_stack_jm}},
#'   \code{\link{posterior_survfit}}, \code{\link{plot.survfit.stanjm}}   
#'     
#' @examples 
#' \donttest{
#'   # Run example model if not already loaded
#'   if (!exists("example_jm")) example(example_jm)
#'   
#'   # For a subset of individuals in the estimation dataset we will
#'   # obtain subject-specific predictions for the longitudinal submodel 
#'   # at evenly spaced times between 0 and their event or censoring time.
#'   pt1 <- posterior_traj(example_jm, ids = c(7,13,15), interpolate = TRUE)
#'   plot(pt1)                  # credible interval for mean response
#'   plot(pt1, limits = "pi")   # prediction interval for raw response
#'   plot(pt1, limits = "none") # no uncertainty interval
#'   
#'   # We can also extrapolate the longitudinal trajectories.
#'   pt2 <- posterior_traj(example_jm, ids = c(7,13,15), interpolate = TRUE,
#'                         extrapolate = TRUE)
#'   plot(pt2)
#'   plot(pt2, vline = TRUE)    # add line indicating event or censoring time
#'   plot(pt2, vline = TRUE, plot_observed = TRUE)  # overlay observed longitudinal data
#'  
#'   # We can change or add attributes to the plot
#'   plot1 <- plot(pt2, ids = c(7,13,15), xlab = "Follow up time",
#'                      vline = TRUE, plot_observed = TRUE, 
#'                      facet_scales = "fixed", color = "blue", linetype = 2,
#'                      ci_geom_args = list(fill = "red"))
#'   plot1
#'        
#'   # Since the returned plot is also a ggplot object, we can
#'   # modify some of its attributes after it has been returned
#'   plot1 + 
#'     ggplot2::theme(strip.background = ggplot2::element_blank()) +
#'     ggplot2::labs(title = "Some plotted longitudinal trajectories")
#' }
plot.predict.stanjm <- function(x, ids = NULL, limits = c("ci", "pi", "none"), 
                                xlab = NULL, ylab = NULL, vline = FALSE, 
                                plot_observed = FALSE, facet_scales = "free_x", 
                                ci_geom_args = NULL, grp_overlay = FALSE, ...) {
  
  limits <- match.arg(limits)
  if (!(limits == "none")) ci <- (limits == "ci")
  y_var <- attr(x, "y_var")
  id_var <- attr(x, "id_var")
  time_var <- attr(x, "time_var")
  grp_var <- attr(x, "grp_var")
  obs_dat <- attr(x, "observed_data")
  if (is.null(ylab)) ylab <- paste0("Long. response (", y_var, ")")
  if (is.null(xlab)) xlab <- paste0("Time (", time_var, ")")
  if (!id_var %in% colnames(x))
    stop("Bug found: could not find 'id_var' column in the data frame.")
  if (!is.null(grp_var) && (!grp_var %in% colnames(x)))
    stop("Bug found: could not find 'grp_var' column in the data frame.")
  if (!is.null(ids)) {
    ids_missing <- which(!ids %in% x[[id_var]])
    if (length(ids_missing))
      stop("The following 'ids' are not present in the predict.stanjm object: ",
           paste(ids[ids_missing], collapse = ", "), call. = FALSE)
    plot_dat <- x[x[[id_var]] %in% ids, , drop = FALSE]
    obs_dat <- obs_dat[obs_dat[[id_var]] %in% ids, , drop = FALSE]
  } else {
    plot_dat <- x
  }
  
  # 'id_list' provides unique IDs sorted in the same order as plotting data
  id_list <- unique(plot_dat[[id_var]])
  if (!is.null(grp_var))
    grp_list <- unique(plot_dat[[grp_var]])
    
  plot_dat$id <- factor(plot_dat[[id_var]])
  plot_dat$time <- plot_dat[[time_var]]
  if (!is.null(grp_var))
    plot_dat$grp <- plot_dat[[grp_var]]
  
  geom_defaults <- list(color = "black", method = "loess", se = FALSE)
  geom_args <- set_geom_args(geom_defaults, ...)
  
  lim_defaults <- list(alpha = 0.3)
  lim_args <- do.call("set_geom_args", c(defaults = list(lim_defaults), ci_geom_args))
  
  obs_defaults <- list()
  obs_args <- set_geom_args(obs_defaults)

  if (is.null(grp_var)) { # no lower level clusters
    group_var <- NULL
    facet_var <- "id"
  } else if (grp_overlay) { # overlay lower level clusters
    group_var <- "grp"
    facet_var <- "id"
  } else { # separate facets for lower level clusters
    group_var <- NULL
    facet_var <- "grp"
  }  
  n_facets <- if (facet_var == "id") length(id_list) else length(grp_list)
  
  if (n_facets > 60L) {
    stop("Too many facets (ie. individuals) to plot. Perhaps limit the ",
         "number of individuals by specifying the 'ids' argument.")
  } else if (n_facets > 1L) {
    geom_mapp <- list(
      mapping = aes_string(x = "time", y = "yfit", group = group_var), 
      data = plot_dat)
    graph <- ggplot() + theme_bw() +
      do.call("geom_smooth", c(geom_mapp, geom_args)) +
      facet_wrap(facet_var, scales = facet_scales)
    if (!limits == "none") {
      graph_smoothlim <- ggplot(plot_dat) + 
        geom_smooth(
          aes_string(x = "time", y = if (ci) "ci_lb" else "pi_lb", group = group_var),
          method = "loess", se = FALSE) +
        geom_smooth(
          aes_string(x = "time", y = if (ci) "ci_ub" else "pi_ub", group = group_var), 
          method = "loess", se = FALSE) +
        facet_wrap(facet_var, scales = facet_scales)
      build_smoothlim <- ggplot_build(graph_smoothlim)
      df_smoothlim <- data.frame(PANEL = build_smoothlim$data[[1]]$PANEL,
                                 time = build_smoothlim$data[[1]]$x,
                                 lb = build_smoothlim$data[[1]]$y,
                                 ub = build_smoothlim$data[[2]]$y,
                                 group = build_smoothlim$data[[1]]$group)
      panel_id_map <- build_smoothlim$layout$layout[, c("PANEL", facet_var), drop = FALSE]
      df_smoothlim <- merge(df_smoothlim, panel_id_map)
      lim_mapp <- list(
        mapping = aes_string(x = "time", ymin = "lb", ymax = "ub", group = "group"), 
        data = df_smoothlim)
      graph_limits <- do.call("geom_ribbon", c(lim_mapp, lim_args))
    } else graph_limits <- NULL
  } else {
    geom_mapp <- list(mapping = aes_string(x = "time", y = "yfit", group = group_var), 
                      data = plot_dat)
    graph <- ggplot() + theme_bw() + 
      do.call("geom_smooth", c(geom_mapp, geom_args))
    if (!(limits == "none")) {
      graph_smoothlim <- ggplot(plot_dat) + 
        geom_smooth(aes_string(x = "time", y = if (ci) "ci_lb" else "pi_lb"), 
                    method = "loess", se = FALSE) +
        geom_smooth(aes_string(x = "time", y = if (ci) "ci_ub" else "pi_ub"), 
                    method = "loess", se = FALSE)
      build_smoothlim <- ggplot_build(graph_smoothlim)
      df_smoothlim <- data.frame(time = build_smoothlim$data[[1]]$x,
                                 lb = build_smoothlim$data[[1]]$y,
                                 ub = build_smoothlim$data[[2]]$y,
                                 group = build_smoothlim$data[[1]]$group) 
      lim_mapp <- list(
        mapping = aes_string(x = "time", ymin = "lb", ymax = "ub", group = "group"), 
        data = df_smoothlim)
      graph_limits <- do.call("geom_ribbon", c(lim_mapp, lim_args))
    } else graph_limits <- NULL
  }    
  if (plot_observed) {
    if (y_var %in% colnames(obs_dat)) {
      obs_dat$y <- obs_dat[[y_var]]
    } else {
      obs_dat$y <- try(eval(parse(text = y_var), obs_dat))
      if (inherits(obs_dat$y, "try-error"))
        stop("Could not find ", y_var, "in observed data, nor able to parse ",
             y_var, "as an expression.")
    }
    obs_dat$id <- factor(obs_dat[[id_var]])
    obs_dat$time <- obs_dat[[time_var]]
    if (!is.null(grp_var))
      obs_dat$grp <- obs_dat[[grp_var]]
    if (is.null(obs_dat[["y"]]))
      stop("Cannot find observed outcome data to add to plot.")
    obs_mapp <- list(
      mapping = aes_string(x = "time", y = "y", group = group_var), 
      data = obs_dat)
    graph_obs <- do.call("geom_point", c(obs_mapp, obs_args)) 
  } else graph_obs <- NULL
  if (vline) {
    if (facet_var == "id") {
      facet_list <- unique(plot_dat[, id_var])
      last_time <- attr(x, "last_time")[as.character(facet_list)] # potentially reorder last_time to match plot_dat
    } else {
      facet_list <- unique(plot_dat[, c(id_var, grp_var)])
      last_time <- attr(x, "last_time")[as.character(facet_list[[id_var]])] # potentially reorder last_time to match plot_dat
      facet_list <- facet_list[[grp_var]]
    }
    vline_dat <- data.frame(FACETVAR = facet_list, last_time = last_time)
    colnames(vline_dat) <- c(facet_var, "last_time")
    graph_vline <- geom_vline(
      mapping = aes_string(xintercept = "last_time"), 
      data = vline_dat, linetype = 2)
  } else graph_vline <- NULL
  
  ret <- graph + graph_limits + graph_obs + graph_vline + labs(x = xlab, y = ylab) 
  class_ret <- class(ret)
  class(ret) <- c("plot.predict.stanjm", class_ret)
  ret
  
}


# internal ----------------------------------------------------------------

# Return a list with the control arguments for interpolation and/or
# extrapolation in posterior_predict.stanmvreg and posterior_survfit.stanjm
#
# @param control A named list, being the user input to the control argument
#   in the posterior_predict.stanmvreg or posterior_survfit.stanjm call
# @param ok_control_args A character vector of allowed control arguments
# @return A named list
get_extrapolation_control <- 
  function(control = list(), ok_control_args = c("epoints", "edist", "eprop")) {
  defaults <- list(ipoints = 15, epoints = 15, edist = NULL, eprop = 0.2, last_time = NULL)
  if (!is.list(control)) {
    stop("'control' should be a named list.")
  } else if (!length(control)) {
    control <- defaults[ok_control_args] 
  } else {  # user specified control list
    nms <- names(control)
    if (!length(nms))
      stop("'control' should be a named list.")
    if (any(!nms %in% ok_control_args))
      stop(paste0("'control' list can only contain the following named arguments: ",
                  paste(ok_control_args, collapse = ", ")))
    if (all(c("edist", "eprop") %in% nms))
      stop("'control' list cannot include both 'edist' and 'eprop'.")        
    if (("ipoints" %in% ok_control_args) && is.null(control$ipoints)) 
      control$ipoints <- defaults$ipoints   
    if (("epoints" %in% ok_control_args) && is.null(control$epoints)) 
      control$epoints <- defaults$epoints  
    if (is.null(control$edist) && is.null(control$eprop)) 
      control$eprop <- defaults$eprop
  }
  return(control)
}

# Set plotting defaults
set_geom_args <- function(defaults, ...) {
  dots <- list(...)
  if (!length(dots)) 
    return(defaults)
  dot_names <- names(dots)
  def_names <- names(defaults)
  for (j in seq_along(def_names)) {
    if (def_names[j] %in% dot_names)
      defaults[[j]] <- dots[[def_names[j]]]
  }
  extras <- setdiff(dot_names, def_names)
  if (length(extras)) {
    for (j in seq_along(extras))
      defaults[[extras[j]]] <- dots[[extras[j]]]
  }
  return(defaults)
}


