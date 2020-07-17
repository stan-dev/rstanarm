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

#' Fitted model objects 
#' 
#' The \pkg{rstanarm} model-fitting functions return an object of class 
#' \code{'stanreg'}, which is a list containing at a minimum the components listed 
#' below. Each \code{stanreg} object will also have additional classes (e.g. 'aov', 
#' 'betareg', 'glm', 'polr', etc.) and several additional components depending
#' on the model and estimation algorithm. \cr
#' \cr
#' Some additional details apply to models estimated using the \code{\link{stan_mvmer}}
#' or \code{\link{stan_jm}} modelling functions. The \code{\link{stan_mvmer}} modelling 
#' function returns an object of class \code{'stanmvreg'}, which inherits the 
#' \code{'stanreg'} class, but has a number of additional elements described in the 
#' subsection below. The \code{\link{stan_jm}} modelling function returns an object of class
#' \code{'stanjm'}, which inherits both the \code{'stanmvreg'} and \code{'stanreg'} 
#' classes, but has a number of additional elements described in the subsection below. 
#' Both the \code{'stanjm'} and \code{'stanmvreg'} classes have several of their own 
#' methods for situations in which the default \code{'stanreg'} methods are not 
#' suitable; see the \strong{See Also} section below.
#' 
#' @name stanreg-objects 
#'   
#' @section Elements for \code{stanreg} objects:   
#' \describe{
#'   \item{\code{coefficients}}{
#'   Point estimates, as described in \code{\link{print.stanreg}}.
#'   }
#'   \item{\code{ses}}{
#'   Standard errors based on \code{\link[stats]{mad}}, as described in
#'   \code{\link{print.stanreg}}.
#'   }
#'   \item{\code{residuals}}{
#'   Residuals of type \code{'response'}.
#'   }
#'   \item{\code{fitted.values}}{
#'   Fitted mean values. For GLMs the linear predictors are transformed by the
#'   inverse link function.
#'   }
#'   \item{\code{linear.predictors}}{
#'   Linear fit on the link scale. For linear models this is the same as
#'   \code{fitted.values}.
#'   }
#'   \item{\code{covmat}}{
#'   Variance-covariance matrix for the coefficients based on draws from the
#'   posterior distribution, the variational approximation, or the asymptotic 
#'   sampling distribution, depending on the estimation algorithm.
#'   }
#'   \item{\code{model,x,y}}{
#'   If requested, the the model frame, model matrix and response variable used, 
#'   respectively.
#'   }
#'   \item{\code{family}}{
#'   The \code{\link[stats]{family}} object used.
#'   }
#'   \item{\code{call}}{
#'   The matched call.
#'   }
#'   \item{\code{formula}}{
#'   The model \code{\link[stats]{formula}}.
#'   }
#'   \item{\code{data,offset,weights}}{
#'   The \code{data}, \code{offset}, and \code{weights} arguments.
#'   }
#'   \item{\code{algorithm}}{
#'   The estimation method used.
#'   }
#'   \item{\code{prior.info}}{
#'   A list with information about the prior distributions used.
#'   }
#'   \item{\code{stanfit,stan_summary}}{
#'   The object of \code{\link[rstan]{stanfit-class}} returned by RStan and a
#'   matrix of various summary statistics from the stanfit object.
#'   }
#'   \item{\code{rstan_version}}{
#'   The version of the \pkg{rstan} package that was used to fit the model.
#'   }
#' }
#'   
#' @section Elements for \code{stanmvreg} objects:   
#' \describe{
#'   The \code{stanmvreg} objects contain the majority of the elements described
#'   above for \code{stanreg} objects, but in most cases these will be a list with each
#'   elements of the list correponding to one of the submodels (for example, 
#'   the \code{family} element of a \code{stanmvreg} object will be a list with each 
#'   element of the list containing the \code{\link[stats]{family}} object for one
#'   submodel). In addition, \code{stanmvreg} objects contain the following additional 
#'   elements:
#'   \item{\code{cnms}}{
#'   The names of the grouping factors and group specific parameters, collapsed 
#'   across the longitudinal or glmer submodels.
#'   }
#'   \item{\code{flevels}}{
#'   The unique factor levels for each grouping factor, collapsed across the 
#'   longitudinal or glmer submodels.
#'   } 
#'   \item{\code{n_markers}}{
#'   The number of longitudinal or glmer submodels.
#'   } 
#'   \item{\code{n_yobs}}{
#'   The number of observations for each longitudinal or glmer submodel.
#'   }
#'   \item{\code{n_grps}}{
#'   The number of levels for each grouping factor (for models estimated using 
#'   \code{\link{stan_jm}}, this will be equal to \code{n_subjects} if the 
#'   individual is the only grouping factor).
#'   }
#'   \item{\code{runtime}}{
#'   The time taken to fit the model (in minutes).
#'   }
#' }
#' 
#' @section Additional elements for \code{stanjm} objects:   
#' \describe{
#'   The \code{stanjm} objects contain the elements described above for 
#'   \code{stanmvreg} objects, but also contain the following additional 
#'   elements:
#'   \item{\code{id_var,time_var}}{
#'   The names of the variables distinguishing between individuals, and 
#'   representing time in the longitudinal submodel.
#'   } 
#'   \item{\code{n_subjects}}{
#'   The number of individuals.
#'   }
#'   \item{\code{n_events}}{
#'   The number of non-censored events.
#'   }  
#'   \item{\code{eventtime,status}}{
#'   The event (or censoring) time and status indicator for each individual.
#'   }
#'   \item{\code{basehaz}}{
#'   A list containing information about the baseline hazard.
#'   }
#'   \item{\code{assoc}}{
#'   An array containing information about the association structure.
#'   }
#'   \item{\code{epsilon}}{
#'   The width of the one-sided difference used to numerically evaluate the 
#'   slope of the longitudinal trajectory; only relevant if a slope-based 
#'   association structure was specified (e.g. etaslope, muslope, etc).
#'   }
#'   \item{\code{qnodes}}{
#'   The number of Gauss-Kronrod quadrature nodes used to evaluate the 
#'   cumulative hazard in the joint likelihood function.
#'   }        
#' }
#' 
#' @note The \code{\link{stan_biglm}} function is an exception. It returns a 
#'   \link[rstan:stanfit-class]{stanfit} object rather than a stanreg object.
#'
#' @seealso \code{\link{stanreg-methods}}, \code{\link{stanmvreg-methods}}
#'   
NULL
