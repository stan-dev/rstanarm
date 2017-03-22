# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2013, 2014, 2015, 2016 Trustees of Columbia University
# Copyright (C) 2016 Monash University
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
#' @templateVar mArg m
#' @template args-m
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
#'   if (!exists("examplejm")) example(examplejm)
#'   
#'   # Obtain subject-specific predictions for all individuals 
#'   # in the estimation dataset
#'   pt1 <- posterior_traj(examplejm, interpolate = FALSE, extrapolate = FALSE)
#'   head(pt1)
#'   
#'   # Obtain subject-specific predictions only for a few selected individuals
#'   pt2 <- posterior_traj(examplejm, ids = c(1,3,8))
#'   
#'   # If we wanted to obtain subject-specific predictions in order to plot the 
#'   # longitudinal trajectories, then we might want to ensure a full trajectory 
#'   # is obtained by interpolating and extrapolating time. We can then use the 
#'   # generic plot function to plot the subject-specific predicted trajectories
#'   # for the first three individuals. Interpolation and extrapolation is 
#'   # carried out by default.
#'   pt3 <- posterior_traj(examplejm)
#'   head(pt3) # predictions at additional time points compared with pt1 
#'   plot(pt3, ids = 1:3)
#'   
#'   # If we wanted to extrapolate further in time, but decrease the number of 
#'   # discrete time points at which we obtain predictions for each individual, 
#'   # then we could specify a named list in the 'control' argument
#'   pt4 <- posterior_traj(examplejm, control = list(ipoints = 10, epoints = 10, eprop = 0.5))
#'   
#'   # Alternatively we may want to estimate the marginal longitudinal
#'   # trajectory for a given set of covariates. To do this, we can pass
#'   # the desired covariate values in a new data frame (however the only
#'   # covariate in our fitted model was the time variable, year). To make sure  
#'   # that we marginalise over the random effects, we need to specify an ID value
#'   # which does not correspond to any of the individuals who were used in the
#'   # model estimation. (The marginal prediction is obtained by generating 
#'   # subject-specific predictions using a series of random draws from the random 
#'   # effects distribution, and then integrating (ie, averaging) over these. 
#'   # Our marginal prediction will therefore capture the between-individual 
#'   # variation associated with the random effects.)
#'   nd <- data.frame(id = rep("new1", 11), year = (0:10 / 2))
#'   pt5 <- posterior_traj(examplejm, newdata = nd)
#'   head(pt5)  # note the greater width of the uncertainty interval compared 
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
#'   pt6 <- posterior_traj(examplejm, newdata = nd, re.form = NA)
#'   head(pt6)  # note the much narrower ci, compared with pt5
#' }
#' 
posterior_traj <- function(object, m = 1, newdata = NULL, 
                           interpolate = TRUE, extrapolate = FALSE,
                           prob = 0.95, ids, control = list(), 
                           return_matrix = FALSE, ...) {
  validate_stanjm_object(object)
  M <- get_M(object)
  m <- validate_positive_scalar(m, not_greater_than = M)
  if (missing(ids)) 
    ids <- NULL
  
  # data: observed data to return to user
  # newX: design matrix used for predictions
  newdata <- validate_newdata(newdata)
  if (is.null(newdata)) {      
    data <- newX <- as.data.frame(model.frame(object)[[m]])
  } else {
    data <- newX <- as.data.frame(newdata)
  }  
  
  id_var <- object$id_var
  if (!id_var %in% names(data))
    stop("The ID variable (", id_var, ") from the original fitted model ",
         "must appear in 'newdata'. To predict for new individuals ",
         "you should enter values for the ID variable that do not coincide ",
         "with those individuals used to fit the model. See the help file ",
         "for more details.", call. = FALSE)
  time_var <- object$time_var
  if (!time_var %in% names(data))
    stop("The time variable from the original fitted model must appear ",
         "in 'newdata'.", call. = FALSE) 
  
  # User specified a subset of ids
  if (!is.null(ids)) {
    check_for_missing_ids(data, id_var, ids)
    data <- data[data[[id_var]] %in% ids, , drop = FALSE] 
    newX <- newX[newX[[id_var]] %in% ids, , drop = FALSE]     
  }
  
  # Ordering of IDs taken from data, not from ids argument  
  id_list <- unique(data[[id_var]])
  
  # Issue warning if IDs in newdata match individuals in the estimation data
  if (!is.null(newdata))
    check_for_estimation_ids(object, id_list, m = m)
  
  # Last known survival time for each individual
  if (is.null(newdata)) { # user did not provide newdata
    last_time <- object$eventtime[as.character(id_list)]
  } else {
    if ("last_time" %in% names(data)) { # user provided newdata with last_time column
      if (!all(tapply(data[["last_time"]], data[[id_var]], FUN = sd) == 0))
        stop("'last_time' column in 'newdata' should be constant within individuals")
      last_time <- tapply(data[["last_time"]], data[[id_var]], FUN = max)
    } else { # user provided newdata but no last_time column, last_time inferred from time_var
      last_time <- tapply(data[[time_var]], data[[id_var]], FUN = max)
    }
  }
  
  # User specified interpolation or extrapolation
  if (interpolate || extrapolate) {  
    if (return_matrix) 
      stop("'return_matrix' cannot be TRUE if either 'interpolate' or ",
           "'extrapolate' is set to TRUE.")
    ok_control_args <- c("ipoints", "epoints", "edist", "eprop")
    control <- get_extrapolation_control(control, ok_control_args = ok_control_args)
    dist <- if (!is.null(control$eprop)) control$eprop * (last_time - 0) else control$edist
    iseq <- if (interpolate) get_time_seq(control$ipoints, 0, last_time) else NULL
    eseq <- if (extrapolate) get_time_seq(control$epoints, last_time, last_time + dist) else NULL
    time_seq <- as.data.frame(cbind(iseq, eseq))
    colnames(time_seq) <- paste0("V", 1:NCOL(time_seq))
    time_seq <- reshape(data.frame(time_seq, id = id_list), 
                        direction = "long", varying = colnames(time_seq), 
                        v.names = time_var, timevar = "obs", idvar = id_var)
    newX[[time_var]] <- as.numeric(newX[[time_var]]) # ensures no rounding during data.table merge
    newX <- data.table::data.table(newX, key = c(id_var, time_var))
    newX <- rolling_merge(newX, time_seq[[id_var]], time_seq[[time_var]])
  }
  
  ytilde <- posterior_predict(object, newdata = newX, ...)
  if (return_matrix) {
    attr(ytilde, "mu") <- NULL # remove attribute mu
    return(ytilde) # return S * N matrix, instead of data frame 
  } 
  mutilde <- attr(ytilde, "mu")
  if (!is.null(newX) && nrow(newX) == 1L) 
    mutilde <- t(mutilde)
  ytilde_bounds  <- median_and_bounds(ytilde,  prob) # median and prob% CrI limits
  mutilde_bounds <- median_and_bounds(mutilde, prob) # median and prob% CrI limits
  
  Terms <- terms(model.frame(object)[[m]])
  y_var <- rownames(attr(Terms, "factors"))[attr(Terms, "response")]
  xvars <- rownames(attr(Terms, "factors"))[-attr(Terms, "response")]
  ret_dat <- if (length(xvars)) as.data.frame(newX)[, xvars, drop = FALSE] else NULL
  out <- data.frame(cbind(ret_dat, yfit = mutilde_bounds$med, 
                          ci_lb = mutilde_bounds$lb, ci_ub = mutilde_bounds$ub, 
                          pi_lb = ytilde_bounds$lb,  pi_ub = ytilde_bounds$ub))
  class(out) <- c("predict.stanjm", "data.frame")
  structure(out, observed_data = data, last_time = last_time,
            y_var = y_var, id_var = id_var, time_var = time_var,
            interpolate = interpolate, extrapolate = extrapolate, 
            control = control, call = match.call())  
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
#' @param ... Optional arguments passed to 
#'   \code{\link[ggplot2]{geom_smooth}} and used to control features
#'   of the plotted longitudinal trajectory.
#'   
#' @return A \code{ggplot} object, also of class \code{plot.predict.stanjm}.
#'   This object can be further customised using the \pkg{ggplot2} package.
#'   It can also be passed to the function \code{\link{plot_stack}}.
#'   
#' @seealso \code{\link{posterior_traj}}, \code{\link{plot_stack}},
#'   \code{\link{posterior_survfit}}, \code{\link{plot.survfit.stanjm}}   
#'     
#' @examples 
#' 
#'   # Run example model if not already loaded
#'   if (!exists("examplejm")) example(examplejm)
#'   
#'   # For a subset of individuals in the estimation dataset we will
#'   # obtain subject-specific predictions for the longitudinal submodel 
#'   # at evenly spaced times between 0 and their event or censoring time.
#'   pt1 <- posterior_traj(examplejm, ids = c(7,13,16), interpolate = TRUE)
#'   plot(pt1)                  # credible interval for mean response
#'   plot(pt1, limits = "pi")   # prediction interval for raw response
#'   plot(pt1, limits = "none") # no uncertainty interval
#'   
#'   # We can also extrapolate the longitudinal trajectories.
#'   pt2 <- posterior_traj(examplejm, ids = c(7,13,16), interpolate = TRUE,
#'                            extrapolate = TRUE)
#'   plot(pt2)
#'   plot(pt2, vline = TRUE)    # add line indicating event or censoring time
#'   plot(pt2, vline = TRUE, plot_observed = TRUE)  # overlay observed longitudinal data
#'  
#'   # We can change or add attributes to the plot
#'   plot1 <- plot(pt2, ids = c(7,13,16), xlab = "Follow up time",
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
#' 
#' 
plot.predict.stanjm <- function(x, ids = NULL, limits = c("ci", "pi", "none"), 
                                xlab = NULL, ylab = NULL, vline = FALSE, 
                                plot_observed = FALSE, facet_scales = "free_x", 
                                ci_geom_args = NULL, ...) {
  
  limits <- match.arg(limits)
  if (!(limits == "none")) ci <- (limits == "ci")
  y_var <- attr(x, "y_var")
  id_var <- attr(x, "id_var")
  time_var <- attr(x, "time_var")
  obs_dat <- attr(x, "observed_data")
  if (is.null(ylab)) ylab <- paste0("Long. response (", y_var, ")")
  if (is.null(xlab)) xlab <- paste0("Time (", time_var, ")")
  if (!id_var %in% colnames(x))
    stop("Bug found: could not find 'id_var' column in the data frame.")
  if (!is.null(ids)) {
    ids_missing <- which(!ids %in% x[[id_var]])
    if (length(ids_missing))
      stop("The following 'ids' are not present in the predict.stanjm object: ",
           paste(ids[[ids_missing]], collapse = ", "), call. = FALSE)
    plot_dat <- x[x[[id_var]] %in% ids, , drop = FALSE]
    obs_dat <- obs_dat[obs_dat[[id_var]] %in% ids, , drop = FALSE]
  } else {
    plot_dat <- x
  }
  
  # 'id_list' provides unique IDs sorted in the same order as plotting data
  id_list <- unique(plot_dat[[id_var]])
  last_time <- attr(x, "last_time")[as.character(id_list)]  # potentially reorder last_time to match plot_dat
  
  plot_dat$time <- plot_dat[[time_var]]
  plot_dat$id <- plot_dat[[id_var]]
  
  obs_dat$y <- obs_dat[[y_var]]
  obs_dat$time <- obs_dat[[time_var]]
  obs_dat$id <- obs_dat[[id_var]]
  
  geom_defaults <- list(color = "black", method = "loess", se = FALSE)
  geom_args <- set_geom_args(geom_defaults, ...)
  
  lim_defaults <- list(alpha = 0.3)
  lim_args <- do.call("set_geom_args", c(defaults = list(lim_defaults), ci_geom_args))
  
  obs_defaults <- list()
  obs_args <- set_geom_args(obs_defaults)
  
  if (length(id_list) > 60L) {
    stop("Too many individuals to plot for. Perhaps limit the number of ",
         "individuals by specifying the 'ids' argument.")
  } else if (length(id_list) > 1L) {
    geom_mapp <- list(mapping = aes_string(x = "time", y = "yfit"), 
                      data = plot_dat)
    graph <- ggplot() + theme_bw() +
      do.call("geom_smooth", c(geom_mapp, geom_args)) +
      facet_wrap(~ id, scales = facet_scales)
    if (!(limits == "none")) {
      graph_smoothlim <- ggplot(plot_dat) + 
        geom_smooth(aes_string(x = "time", y = if (ci) "ci_lb" else "pi_lb"), 
                    method = "loess", se = FALSE) +
        geom_smooth(aes_string(x = "time", y = if (ci) "ci_ub" else "pi_ub"), 
                    method = "loess", se = FALSE) +
        facet_wrap(~ id, scales = facet_scales)
      build_smoothlim <- ggplot_build(graph_smoothlim)
      df_smoothlim <- data.frame(PANEL = build_smoothlim$data[[1]]$PANEL,
                                 time = build_smoothlim$data[[1]]$x,
                                 lb = build_smoothlim$data[[1]]$y,
                                 ub = build_smoothlim$data[[2]]$y)
      panel_id_map <- build_smoothlim$layout$panel_layout[, c("PANEL", "id"), drop = FALSE]
      df_smoothlim <- merge(df_smoothlim, panel_id_map)
      lim_mapp <- list(mapping = aes_string(x = "time", ymin = "lb", ymax = "ub"), 
                       data = df_smoothlim)
      graph_limits <- do.call("geom_ribbon", c(lim_mapp, lim_args))
    } else graph_limits <- NULL
  } else {
    geom_mapp <- list(mapping = aes_string(x = "time", y = "yfit"), 
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
                                 ub = build_smoothlim$data[[2]]$y) 
      lim_mapp <- list(mapping = aes_string(x = "time", ymin = "lb", ymax = "ub"), 
                       data = df_smoothlim)
      graph_limits <- do.call("geom_ribbon", c(lim_mapp, lim_args))
    } else graph_limits <- NULL
  }    
  if (plot_observed) {
    if (is.null(obs_dat[["y"]]))
      stop("Cannot find observed outcome data to add to plot.")
    obs_mapp <- list(mapping = aes_string(x = "time", y = "y"), 
                     data = obs_dat)
    graph_obs <- do.call("geom_point", c(obs_mapp, obs_args)) 
  } else graph_obs <- NULL
  if (vline) {
    graph_vline <- geom_vline(aes_string(xintercept = "last_time"), 
                              data.frame(id = id_list, last_time = last_time), 
                              linetype = 2)
  } else graph_vline <- NULL
  
  ret <- graph + graph_limits + graph_obs + graph_vline + labs(x = xlab, y = ylab) 
  class_ret <- class(ret)
  class(ret) <- c("plot.predict.stanjm", class_ret)
  ret
  
}


# internal ----------------------------------------------------------------

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

# Check whether individuals in the ids argument are present in the data
#
# @param data A data frame
# @param id_var The name of the ID variable in data
# @param ids A vector of ids
check_for_missing_ids <- function(data, id_var, ids) {
  if (!is.data.frame(data))
    stop("'data' should be a data frame.")
  if (!all(is.character(id_var), length(id_var == 1L)))
    stop("'id_var' should be a variable name.")
  if (!id_var %in% colnames(data))
    stop(paste(substitute(deparse(id_var)), "not found in 'data'."))
  sel <- which(!ids %in% data[[id_var]])
  if (length(sel))
    stop("The following 'ids' do not appear in the data: ", 
         paste(ids[[sel]], collapse = ", "))
}

# Check whether individuals in the ids argument were used in the 
# model estimation
#
# @param object A stanjm object
# @param ids A vector of ids
# @param m Integer specifying which submodel to get the estimation IDs from
check_for_estimation_ids <- function(object, ids, m = 1) {
  ids2 <- unique(model.frame(object, m = m)[[object$id_var]])
  if (any(ids %in% ids2))
    warning("Some of the IDs in the 'newdata' correspond to individuals in the ",
            "estimation dataset. Please be sure you want to obtain subject-",
            "specific predictions using the estimated random effects for those ",
            "individuals. If you instead meant to marginalise over the distribution ",
            "of the random effects, then please make sure the ID values do not ",
            "correspond to individuals in the estimation dataset.", immediate. = TRUE)
}

# Return a list with the control arguments for interpolation and/or
# extrapolation in posterior_predict.stanjm and posterior_survfit.stanjm
#
# @param control A named list, being the user input to the control argument
#   in the posterior_predict.stanjm or posterior_survfit.stanjm call
# @param ok_control_args A character vector of allowed control arguments
# @param standardise A logical, being the user input to the standardise
#   argument in a posterior_survfit.stanjm call.
# @return A named list
get_extrapolation_control <- function(control = list(), 
                                      ok_control_args = c("epoints", "edist", "eprop"), 
                                      standardise = FALSE) {
  defaults <- list(ipoints = 15, epoints = 15, edist = NULL, eprop = 0.2,
                   condition = TRUE, last_time = NULL)
  if (!is.list(control)) {
    stop("'control' should be a named list.")
  } else if (!length(control)) {
    control <- defaults[ok_control_args] 
    if (("condition" %in% ok_control_args) && standardise) 
      control$condition <- FALSE
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
    if (("condition" %in% ok_control_args) && is.null(control$condition)) {
      control$condition <- if (!standardise) defaults$condition else FALSE
    } else if (("condition" %in% ok_control_args) && control$condition && standardise) {
      stop("'condition' cannot be set to TRUE if standardised survival ",
           "probabilities are requested.")
    }
  }
  return(control)
}

# Return an array or list with the time sequence used for posterior predictions
#
# @param increments An integer with the number of increments (time points) at
#   which to predict the outcome for each individual
# @param t0,t1 Numeric vectors giving the start and end times across which to
#   generate prediction times
# @param simplify Logical specifying whether to return each increment as a 
#   column of an array (TRUE) or as an element of a list (FALSE) 
get_time_seq <- function(increments, t0, t1, simplify = TRUE) {
  val <- sapply(0:(increments - 1), function(x, t0, t1) {
    t0 + (t1 - t0) * (x / (increments - 1))
  }, t0 = t0, t1 = t1, simplify = simplify)
  if (simplify && is.vector(val)) {
    # need to transform if there is only one individual
    val <- t(val)
    rownames(val) <- if (!is.null(names(t0))) names(t0) else 
      if (!is.null(names(t1))) names(t1) else NULL
  }
  return(val)
}
