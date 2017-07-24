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

#' Fitted multivariate model object
#' 
#' The \pkg{rstanarm} package contains several model-fitting functions which fit a 
#' model that is multivariate, in the sense that it includes one or more component
#' \emph{submodels}. These functions currently include: \code{\link{stan_mvmer}},
#' \code{\link{stan_jm}}. These functions return an object of class 
#' 'stanmvreg', which is a list containing some or all of the components listed 
#' below. In most cases each of the components return a named list, with each 
#' element of the list containing the described item (for example, coefficients, 
#' standard errors, etc) for one of the submodels. \cr
#' \cr
#' Note that each stanmvreg object also inherits the class 
#' \code{\link[=stanreg-objects]{stanreg}}.
#' 
#' @name stanmvreg-objects 
#'   
#' @section stanmvreg objects:   
#' \describe{
#'   \item{\code{coefficients}}{
#'   Point estimates for each of the submodels.
#'   As described in \code{\link{print.stanreg}}.
#'   }
#'   \item{\code{ses}}{
#'   Standard errors for each of the submodels.
#'   Based on \code{\link[stats]{mad}}, as described in
#'   \code{\link{print.stanreg}}.
#'   }
#'   \item{\code{covmat}}{
#'   Variance-covariance matrix for the coefficients based on draws from the
#'   posterior distribution, the variational approximation, or the asymptotic 
#'   sampling distribution, depending on the estimation algorithm.
#'   }
#'   \item{\code{formula}}{
#'   The model \code{\link[stats]{formula}} for each of the submodels.
#'   }
#'   \item{\code{prior.weights}}{
#'   For models estimated using \code{\link{stan_jm}} and with weights 
#'   specified, the data frame with the prior weights for each individual.
#'   }
#'   \item{\code{prior.info}}{
#'   A list with information about the prior distributions used.
#'   }    
#'   \item{\code{algorithm}}{
#'   The estimation method used.
#'   }    
#'   \item{\code{stan_function}}{
#'   The model fitting function used.
#'   }
#'   \item{\code{cnms}}{
#'   A list of column names of the random effects according to the grouping 
#'   factors, collapsed across all longitudinal orGLM submodels. 
#'   See \code{\link[lme4]{mkReTrms}}. These can be obtained separately
#'   for each submodel through the components 
#'   \code{object$glmod[[1]]@cnms}, \code{object$glmod[[2]]@cnms}, etc.
#'   }
#'   \item{\code{n_markers}}{
#'   The number of longitudinal or GLM submodels.
#'   }
#'   \item{\code{n_yobs}}{
#'   The number of observations for each longitudinal or GLM submodel.
#'   }
#'   \item{\code{n_grps}}{
#'   The number of levels for each grouping factor. For models estimated using 
#'   \code{\link{stan_jm}},will be equal to \code{n_subjects} if the 
#'   individual is the only grouping factor.
#'   }
#'   \item{\code{n_subjects}}{
#'   For models estimated using \code{\link{stan_jm}}, the number of 
#'   individuals.
#'   }
#'   \item{\code{n_events}}{
#'   For models estimated using \code{\link{stan_jm}}, the number of 
#'   non-censored event times.
#'   }   
#'   \item{\code{id_var,time_var}}{
#'   For models estimated using \code{\link{stan_jm}}, the names of 
#'   the variables distinguishing between individuals, and 
#'   representing time, in the longitudinal submodel.
#'   } 
#'   \item{\code{y}}{
#'   The response vector for each of the GLM submodels.
#'   }  
#'   \item{\code{eventtime,status}}{
#'   For models estimated using \code{\link{stan_jm}}, the event (or censoring) 
#'   times and the status indicator for each individual.
#'   }   
#'   \item{\code{data,dataLong,dataEvent}}{
#'   The data frames that were passed to the data arguments.
#'   }
#'   \item{\code{call}}{
#'   The matched call.
#'   }
#'   \item{\code{runtime}}{
#'   The time taken to fit the model (in minutes).
#'   }
#'   \item{\code{family}}{
#'   The \code{\link[stats]{family}} object used for each of the GLM submodels.
#'   }
#'   \item{\code{basehaz}}{
#'   For models estimated using \code{\link{stan_jm}}, a list containing 
#'   information about the baseline hazard.
#'   }
#'   \item{\code{assoc}}{
#'   For models estimated using \code{\link{stan_jm}}, a list containing 
#'   information about the association structure.
#'   }
#'   \item{\code{epsilon}}{
#'   For models estimated using \code{\link{stan_jm}}, the width of the 
#'   one-sided difference used to numerically evaluate the slope of the
#'   longitudinal trajectory; this is only relevant if a slope-based 
#'   association structure was specified (e.g. etaslope, muslope, etc).
#'   }   
#'   \item{\code{quadnodes}}{
#'   For models estimated using \code{\link{stan_jm}}, the number of 
#'   Gauss-Kronrod quadrature nodes used to evaluate the 
#'   cumulative hazard in the joint likelihood function.
#'   }
#'   \item{\code{runtime}}{
#'   The time taken (in mins) for model estimation, presented separately for  
#'   each MCMC chain.
#'   }
#'   \item{\code{stanfit,stan_summary}}{
#'   The object of \code{\link[rstan]{stanfit-class}} returned by RStan and a
#'   matrix of various summary statistics from the stanfit object.
#'   }
#'   \item{\code{glmod}}{
#'   The fitted model object(s) returned by a call(s) to \code{\link[lme4]{lmer}}
#'   or \code{\link[lme4]{glmer}} when estimating the separate GLM model(s)
#'   prior to fitting the multivariate model.
#'   }
#'   \item{\code{coxmod}}{
#'   For models estimated using \code{\link{stan_jm}}, the fitted model object 
#'   returned by a call to \code{\link[survival]{coxph}}
#'   when estimating the separate time-to-event model prior to fitting the joint model.
#'   This time-to-event model is evaluated using the observed data in \code{dataEvent};
#'   it is not evaluated at the quadrature points necessary for fitting the joint model.
#'   }   
#'}
#'
#' @seealso Objects of this class have a number of generic 
#'   methods described in \code{\link{stanmvreg-methods}},
#'   as well as \code{\link{print.stanmvreg}}, \code{\link{summary.stanmvreg}},
#'   as well as the non-generic functions
#'   \code{\link{posterior_traj}}, \code{\link{posterior_survfit}}, 
#'   \code{\link{posterior_predict}}, \code{\link{posterior_interval}},
#'   \code{\link{pp_check}} and \code{\link{ps_check}}.
#'   
NULL