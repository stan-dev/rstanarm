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

#' Using the ShinyStan GUI with stanreg objects
#' 
#' The \code{\link[shinystan]{launch_shinystan}} function will accept a 
#' \link[=stanreg-objects]{'stanreg' object} as input. Currently, almost any
#' model fit using one of \pkg{rstanarm}'s model-fitting functions can be
#' used with ShinyStan. The only exception is that ShinyStan does not currently
#' support \pkg{rstanarm} models fit using \code{algorithm='optimizing'}.
#' 
#' See the \pkg{\link[=shinystan-package]{shinystan}} package documentation for
#' more information.
#' 
#' @name shinystan
#' @aliases launch_shinystan
#' @importFrom shinystan launch_shinystan
#'   
#' @examples
#' \dontrun{
#' # Launch the ShinyStan app (saving resulting shinystan object as sso)
#' sso <- launch_shinystan(example_model)
#' }
#' 
#' 
NULL