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

# Return design matrices for evaluating longitudinal submodel quantities 
# at specified quadrature points/times
#
# @param data A data frame, the data for the longitudinal submodel.
# @param assoc A list with information about the association structure for 
#   the one longitudinal submodel. 
# @param y_mod A named list returned by a call to handle_y_mod (the
#   fit for a single longitudinal submodel)
# @param grp_stuff A list with information about any lower level grouping
#   factors that are clustered within patients and how to handle them in 
#   the association structure.
# @param ids,times The subject IDs and times vectors that correspond to the
#   event and quadrature times at which the design matrices will
#   need to be evaluated for the association structure.
# @param id_var The name on the ID variable.
# @param time_var The name of the time variable.
# @param epsilon The half-width of the central difference used for 
#   numerically calculating the derivative of the design matrix for slope
#   based association structures.
# @param auc_qnodes Integer specifying the number of GK quadrature nodes to
#   use in the integral/AUC based association structures.
# @return The list returned by make_assoc_parts.
handle_assocmod <- function(data, assoc, y_mod, e_mod, grp_stuff, meta) {
  
  if (!requireNamespace("data.table"))
    stop2("the 'data.table' package must be installed to use this function.")
  
  id_var   <- meta$id_var
  time_var <- meta$time_var
  
  # before turning data into a data.table (for a rolling merge against 
  # the quadrature points) we want to make sure that the data does not
  # include any NAs for the predictors or assoc formula variables
  tt <- attr(y_mod$terms, "term.labels")
  tt <- c(tt, uapply(assoc[["which_formulas"]], all.vars))
  fm <- reformulate(tt, response = NULL)
  df <- get_all_vars(fm, data)
  df <- df[complete.cases(df), , drop = FALSE]
  
  # declare df as a data.table for merging with quadrature points
  dt <- prepare_data_table(df, 
                           id_var   = id_var, 
                           time_var = time_var, 
                           grp_var  = grp_stuff$grp_var) # grp_var may be NULL
  
  # design matrices for calculating association structure based on 
  # (possibly lagged) eta, slope, auc and any interactions with data
  args <- list(use_function = make_assoc_parts_for_stan,
               newdata      = dt, 
               y_mod        = y_mod, 
               grp_stuff    = grp_stuff, 
               meta         = meta, 
               assoc        = assoc,
               ids          = e_mod$cids,
               times        = e_mod$cpts)
  
  do.call(make_assoc_parts, args)
}
