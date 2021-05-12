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

#' \code{adapt_delta}: Target average acceptance probability
#' 
#' Details about the \code{adapt_delta} argument to \pkg{rstanarm}'s modeling
#' functions.
#' 
#' @name adapt_delta
#' @template reference-stan-manual
#' @references Brief Guide to Stan's Warnings:
#'   \url{https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup}
#'   
#'   
#' @details For the No-U-Turn Sampler (NUTS), the variant of Hamiltonian Monte
#'   Carlo used used by \pkg{rstanarm}, \code{adapt_delta} is the target average
#'   proposal acceptance probability during Stan's adaptation period.
#'   \code{adapt_delta} is ignored by \pkg{rstanarm} if the \code{algorithm} argument 
#'   is not set to \code{"sampling"}.
#'   
#'   The default value of \code{adapt_delta} is 0.95, except when the prior for 
#'   the regression coefficients is \code{\link{R2}}, \code{\link{hs}}, or 
#'   \code{\link{hs_plus}}, in which case the default is 0.99.
#'   
#'   These defaults are higher (more conservative) than the default of
#'   \code{adapt_delta=0.8} used in the \pkg{rstan} package, which may result in
#'   slower sampling speeds but will be more robust to posterior distributions
#'   with high curvature.
#'   
#'   In general you should not need to change \code{adapt_delta} unless you see
#'   a warning message about divergent transitions, in which case you can
#'   increase \code{adapt_delta} from the default to a value \emph{closer} to 1
#'   (e.g. from 0.95 to 0.99, or from 0.99 to 0.999, etc). The step size used by
#'   the numerical integrator is a function of \code{adapt_delta} in that
#'   increasing \code{adapt_delta} will result in a smaller step size and fewer
#'   divergences. Increasing \code{adapt_delta} will typically result in a
#'   slower sampler, but it will always lead to a more robust sampler.
NULL
