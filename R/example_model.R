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
#' Example model
#' 
#' A pre-fit model for use in \pkg{rstanarm} examples. 
#' 
#' @name example_model
#' @format A \code{\link[=stanreg-objects]{stanreg}} object containing the
#'   output from fitting the model in the Examples section, below.
#'   The \code{chains} and \code{iter} arguments are specified to make this example be 
#'   small in size. In practice, we recommend that they be left unspecified in
#'   order to use the default values (4 and 2000 respectively) or increased if
#'   there are convergence problems. The \code{cores} argument is optional and
#'   on a multicore system, the user may well want to set that equal to the
#'   number of chains being executed.
#' 
#' @seealso \code{\link[lme4]{cbpp}} for a description of the data.
#' @examples
#' \dontrun{
#' example_model <- 
#'   stan_glmer(cbind(incidence, size - incidence) ~ size + period + (1|herd),
#'              data = lme4::cbpp, family = binomial,
#'              # this next line is only to keep the example small in size!
#'              chains = 2, cores = 1, seed = 12345, iter = 500)
#' }
#' example_model
NULL
