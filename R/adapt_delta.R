# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2015, 2016 Trustees of Columbia University
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
#' Target average acceptance probability
#' 
#' @name adapt_delta
#' @template reference-stan-manual
#'   
#' @details For the No-U-Turn Sampler (NUTS), the variant of Hamiltonian Monte
#'   Carlo used used by \pkg{rstanarm}, \code{adapt_delta} is the target average
#'   proposal acceptance probability for adaptation. \code{adapt_delta} is
#'   ignored if \code{algorithm} is not \code{"sampling"}.
#' 
#' The default value of \code{adapt_delta} is 0.95, except when the prior for
#' the regression coefficients is \code{\link{R2}}, \code{\link{hs}},
#' \code{\link{hs_plus}}, or \code{\link{student_t}} with \code{df <= 2}, in
#' which case the default is 0.99.
#' 
#' In general you should not need to change \code{adapt_delta} unless you see a 
#' warning message about divergent transitions, in which case you can increase 
#' \code{adapt_delta} from the default to a value \emph{closer} to 1 (e.g. from 
#' 0.95 to 0.99, or from 0.99 to 0.999, etc). The step size used by the 
#' numerical integrator is a function of \code{adapt_delta} in that increasing 
#' \code{adapt_delta} will result in a smaller step size and fewer divergences. 
#' Increasing \code{adapt_delta} will typically result in a slower sampler, but
#' it will always lead to a more robust sampler.
NULL