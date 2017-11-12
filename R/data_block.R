# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2013, 2014, 2015, 2016, 2017 Trustees of Columbia University
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
  if (has_intercept && !sparse) {
    xbar <- colMeans(xtemp)
    xtemp <- sweep(xtemp, 2, xbar, FUN = "-")
  }
  else xbar <- rep(0, ncol(xtemp))
  
  sel <- apply(xtemp, 2L, function(x) !all(x == 1) && length(unique(x)) < 2)
  if (any(sel)) {
    # drop any column of x with < 2 unique values (empty interaction levels)
    # exception is column of 1s isn't dropped 
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
# @param default_scale Default value to use to scale if not specified by user
# @param link String naming the link function.
# @param ok_dists A list of admissible distributions.
handle_glm_prior <- function(prior, nvars, default_scale, link,
                             ok_dists = nlist("normal", student_t = "t", 
                                              "cauchy", "hs", "hs_plus", 
                                              "laplace", "lasso", "product_normal")) {
  if (!length(prior))
    return(list(prior_dist = 0L, prior_mean = as.array(rep(0, nvars)),
                prior_scale = as.array(rep(1, nvars)),
                prior_df = as.array(rep(1, nvars)), prior_dist_name = NA,
                global_prior_scale = 0, global_prior_df = 0,
                slab_df = 0, slab_scale = 0,
                prior_autoscale = FALSE))

  if (!is.list(prior)) 
    stop(sQuote(deparse(substitute(prior))), " should be a named list")
  
  prior_dist_name <- prior$dist
  prior_scale <- prior$scale
  prior_mean <- prior$location
  prior_df <- prior$df
  prior_mean[is.na(prior_mean)] <- 0
  prior_df[is.na(prior_df)] <- 1
  global_prior_scale <- 0
  global_prior_df <- 0
  slab_df <- 0
  slab_scale <- 0
  if (!prior_dist_name %in% unlist(ok_dists)) {
    stop("The prior distribution should be one of ",
         paste(names(ok_dists), collapse = ", "))
  } else if (prior_dist_name %in% 
             c("normal", "t", "cauchy", "laplace", "lasso", "product_normal")) {
    if (prior_dist_name == "normal") prior_dist <- 1L
    else if (prior_dist_name == "t") prior_dist <- 2L
    else if (prior_dist_name == "laplace") prior_dist <- 5L
    else if (prior_dist_name == "lasso") prior_dist <- 6L
    else if (prior_dist_name == "product_normal") prior_dist <- 7L
    prior_scale <- set_prior_scale(prior_scale, default = default_scale, 
                                   link = link)
  } else if (prior_dist_name %in% c("hs", "hs_plus")) {
    prior_dist <- ifelse(prior_dist_name == "hs", 3L, 4L)
    global_prior_scale <- prior$global_scale
    global_prior_df <- prior$global_df
    slab_df <- prior$slab_df
    slab_scale <- prior$slab_scale
  } else if (prior_dist_name %in% "exponential") {
    prior_dist <- 3L # only used for scale parameters so 3 not a conflict with 3 for hs
  }
  
  prior_df <- maybe_broadcast(prior_df, nvars)
  prior_df <- as.array(pmin(.Machine$double.xmax, prior_df))
  prior_mean <- maybe_broadcast(prior_mean, nvars)
  prior_mean <- as.array(prior_mean)
  prior_scale <- maybe_broadcast(prior_scale, nvars)

  nlist(prior_dist, 
        prior_mean, 
        prior_scale, 
        prior_df, 
        prior_dist_name, 
        global_prior_scale,
        global_prior_df,
        slab_df,
        slab_scale,
        prior_autoscale = isTRUE(prior$autoscale))
}
