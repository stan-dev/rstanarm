# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2013, 2014, 2015, 2016 Trustees of Columbia University
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

# Center a matrix x and return extra stuff
#
# @param x A design matrix
# @param sparse A flag indicating whether x is to be treated as sparse
center_x <- function(x, sparse) {
  x <- as.matrix(x)
  has_intercept <- if (ncol(x) == 0) 
    FALSE else grepl("(Intercept", colnames(x)[1L], fixed = TRUE)
  
  xtemp <- if (has_intercept) x[, -1L, drop=FALSE] else x
  if (!sparse) {
    xbar <- colMeans(xtemp)
    xtemp <- sweep(xtemp, 2, xbar, FUN = "-")
  }
  else xbar <- rep(0, ncol(xtemp))
  
  sel <- (2 > apply(xtemp, 2L, function(x) length(unique(x))))
  if (any(sel)) {
    # drop any column of x with < 2 unique values (empty interaction levels)
    warning("Dropped empty interaction levels: ",
            paste(colnames(xtemp)[sel], collapse = ", "))
    xtemp <- xtemp[, !sel, drop = FALSE]
    xbar <- xbar[!sel]
  }
  
  return(nlist(xtemp, xbar, has_intercept))
}

# Deal with priors
#
# @param prior A list
# @param nvars An integer indicating the number of variables
# @param ok_dists A character vector of admissible distributions
handle_glm_prior <- function(prior, nvars, default_scale, link,
                             ok_dists = nlist("normal", student_t = "t", "cauchy", 
                                              "hs", "hs_plus")) {
  if (!length(prior))
    return(list(prior_dist = 0L, prior_mean = as.array(rep(0, nvars)),
                prior_scale = as.array(rep(1, nvars)),
                prior_df = as.array(rep(1, nvars)), 
                global_prior_scale = 0, global_prior_df = 0))
  if (!is.list(prior)) 
    stop(sQuote(deparse(substitute(prior))), " should be a named list")
  
  prior_dist <- prior$dist
  prior_scale <- prior$scale
  prior_mean <- prior$location
  prior_df <- prior$df
  prior_df[is.na(prior_df)] <- 1
  if (!prior_dist %in% unlist(ok_dists)) {
    stop("The prior distribution should be one of ",
         paste(names(ok_dists), collapse = ", "))
  } else if (prior_dist %in% c("normal", "t", "cauchy")) {
    prior_dist <- ifelse(prior_dist == "normal", 1L, 2L)
    prior_scale <- set_prior_scale(prior_scale, default = default_scale, 
                                   link = link)
    global_prior_scale <- 0
    global_prior_df <- 0
  } else {
    prior_dist <- ifelse(prior_dist == "hs", 3L, 4L)
    global_prior_scale <- prior$global_scale
    global_prior_df <- prior$global_df
  }
  
  prior_df <- maybe_broadcast(prior_df, nvars)
  prior_df <- as.array(pmin(.Machine$double.xmax, prior_df))
  prior_mean <- maybe_broadcast(prior_mean, nvars)
  prior_mean <- as.array(prior_mean)
  prior_scale <- maybe_broadcast(prior_scale, nvars)
  return(nlist(prior_dist, prior_mean, prior_scale, prior_df, 
               global_prior_scale, global_prior_df)) 
}
