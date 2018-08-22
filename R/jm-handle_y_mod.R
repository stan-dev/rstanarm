# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2016, 2017, 2018 Sam Brilleman
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

# Construct a list with information about the longitudinal submodel
#
# @param formula The model formula for the glmer submodel.
# @param data The data for the glmer submodel.
# @param family The family object for the glmer submodel.
# @return A named list with the following elements:
#   y: named list with the reponse vector and related info.
#   x: named list with the fe design matrix and related info.
#   z: named list with the re design matrices and related info.
#   terms: the model.frame terms object with bars "|" replaced by "+".
#   model_frame: The model frame with all variables used in the 
#     model formula.
#   formula: The model formula.
#   reTrms: returned by lme4::glFormula$reTrms.
#   family: the (modified) family object for the glmer submodel.
#   intercept_type: named list with info about the type of 
#     intercept required for the glmer submodel.
#   has_aux: logical specifying whether the glmer submodel 
#     requires an auxiliary parameter.
handle_y_mod <- function(formula, data, family, stub) {
  mf <- stats::model.frame(lme4::subbars(formula), data)
  if (!length(formula) == 3L)
    stop2("An outcome variable must be specified.")
  
  # lme4 parts
  lme4_parts <- lme4::glFormula(formula, data)
  reTrms <- lme4_parts$reTrms
  
  # Response vector, design matrices
  y <- make_y_for_stan(formula, mf, family) 
  x <- make_x_for_stan(formula, mf)
  z <- make_z_for_stan(formula, mf) 
  
  # Terms
  terms <- attr(mf, "terms")
  terms <- append_predvars_attribute(terms, formula, data)
  
  # Binomial with >1 trials not allowed by stan_{mvmver,jm}
  is_binomial <- is.binomial(family$family)
  is_bernoulli <- is_binomial && NCOL(y$y) == 1L && all(y$y %in% 0:1)
  if (is_binomial && !is_bernoulli)
    STOP_binomial()
  
  # Various flags
  intercept_type <- check_intercept_type(x, family)
  has_aux <- check_for_aux(family)
  family <- append_mvmer_famlink(family, is_bernoulli)
  
  nlist(y, x, z, reTrms, model_frame = mf, formula, terms, 
        family, intercept_type, has_aux, stub)
}
