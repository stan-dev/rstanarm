# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2013, 2014, 2015, 2016, 2017 Trustees of Columbia University
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

# Function to construct a design matrix for the association structure in
# the event submodel, to be multiplied by a vector of association parameters
#
# @param assoc An array with information about the desired association 
#   structure, returned by a call to validate_assoc.
# @param parts A list equal in length to the number of markers. Each element
#   parts[[m]] should contain a named list with components $mod_eta, $mod_eps,
#   $mod_auc, etc, which each contain either the linear predictor at quadtimes, 
#   quadtimes + eps, and auc quadtimes, or the design matrices
#   used for constructing the linear predictor. Each element parts[[m]] should 
#   also contain $X_data and $K_data.
# @param family A list of family objects, equal in length to the number of 
#   longitudinal submodels.
# @param ... If parts does not contain the linear predictors, then this should
#   include elements beta and b, each being a length M list of parameters for the
#   longitudinal submodels.
# @return A design matrix containing the association terms to be multiplied by
#   the association paramters.
make_assoc_terms <- function(parts, assoc, family, ...) {
  M <- length(parts)
  a_X <- list()
  mark <- 1
  for (m in 1:M) {
    times   <- attr(parts[[m]], "times")
    epsilon <- attr(parts[[m]], "epsilon")  
    qnodes  <- attr(parts[[m]], "auc_qnodes")
    qwts    <- attr(parts[[m]], "auc_qwts")
    
    eps_uses_derivative_of_x <- 
      attr(parts[[m]], "eps_uses_derivative_of_x") # experimental
    
    has_assoc <- !assoc["null",][[m]]
    
    if (has_assoc) {
      assoc_m   <- assoc[,m]
      invlink_m <- family[[m]]$linkinv    
      eta_m    <- get_element(parts, m = m, "eta", ...)
      eps_m    <- get_element(parts, m = m, "eps", ...)
      auc_m    <- get_element(parts, m = m, "auc", ...)
      X_data_m <- get_element(parts, m = m, "X_data", ...)
      K_data_m <- get_element(parts, m = m, "K_data", ...)
      grp_m    <- get_element(parts, m = m, "grp_stuff", ...)
      
      has_grp   <- grp_m$has_grp # TRUE/FALSE
      if (has_grp) {
        # method for collapsing information across clusters within patients 
        grp_assoc <- grp_m$grp_assoc
        # indexing for collapsing across grps (based on the ids and times
        # used to generate the design matrices in make_assoc_parts)
        grp_idx <- attr(parts[[m]], "grp_idx") 
      }
      
      #---  etavalue and any interactions  ---#
      
      # etavalue
      if (assoc_m[["etavalue"]]) { 
        if (has_grp) {
          a_X[[mark]] <- collapse_within_groups(eta_m, grp_idx, grp_assoc)
        } else {
          a_X[[mark]] <- eta_m
        }
        mark <- mark + 1
      }
      
      # etavalue * data interactions
      if (assoc_m[["etavalue_data"]]) { 
        X_temp <- X_data_m[["etavalue_data"]]
        K_temp <- K_data_m[["etavalue_data"]]
        for (i in 1:K_temp) {
          if (is.matrix(eta_m)) {
            val <- sweep(eta_m, 2L, X_temp[, i], `*`)
          } else {
            val <- as.vector(eta_m) * X_temp[, i]
          }
          if (has_grp) {
            a_X[[mark]] <- collapse_within_groups(val, grp_idx, grp_assoc)
          } else {
            a_X[[mark]] <- val
          }
          mark <- mark + 1
        }
      }
      
      # etavalue * etavalue interactions
      if (assoc_m[["etavalue_etavalue"]]) {
        sel <- assoc_m[["which_interactions"]][["etavalue_etavalue"]]
        for (j in sel) {
          eta_j <- get_element(parts, m = j, "eta", ...)
          val <- eta_m * eta_j 
          a_X[[mark]] <- val
          mark <- mark + 1
        }
      }
      
      # etavalue * muvalue interactions
      if (assoc_m[["etavalue_muvalue"]]) {
        sel <- assoc_m[["which_interactions"]][["etavalue_muvalue"]]
        for (j in sel) {
          eta_j <- get_element(parts, m = j, "eta", ...)
          invlink_j <- family[[j]]$linkinv
          val <- eta_m * invlink_j(eta_j) 
          a_X[[mark]] <- val
          mark <- mark + 1             
        }
      }
      
      #---  etaslope and any interactions  ---#
    
      if (assoc_m[["etaslope"]] || assoc_m[["etaslope_data"]]) {
        if (eps_uses_derivative_of_x) {
          deta_m <- eps_m
        } else {
          deta_m <- (eps_m - eta_m) / epsilon
        }
      }
        
      # etaslope
      if (assoc_m[["etaslope"]]) {
        if (has_grp) {
          a_X[[mark]] <- collapse_within_groups(deta_m, grp_idx, grp_assoc)
        } else {
          a_X[[mark]] <- deta_m
        } 
        mark <- mark + 1             
      }
      
      # etaslope * data interactions
      if (assoc_m[["etaslope_data"]]) { 
        X_temp <- X_data_m[["etaslope_data"]]
        K_temp <- K_data_m[["etaslope_data"]]
        for (i in 1:K_temp) {
          if (is.matrix(deta_m)) {
            val <- sweep(deta_m, 2L, X_temp[, i], `*`)
          } else {
            val <- as.vector(deta_m) * X_temp[, i]
          }
          if (has_grp) {
            a_X[[mark]] <- collapse_within_groups(val, grp_idx, grp_assoc)
          } else {
            a_X[[mark]] <- val
          }
          mark <- mark + 1
        }
      }
      
      #---  etaauc  ---#
      
      if (assoc_m[["etaauc"]]) {
        if (is.matrix(eta_m)) {
          nr <- nrow(eta_m)
          nc <- ncol(eta_m)
          val   <- matrix(NA, nrow = nr, ncol = nc) 
          for (j in 1:nc) {
            wgt_j <- qwts[((j-1) * qnodes + 1):(j * qnodes)]
            auc_j <- auc_m[, ((j-1) * qnodes + 1):(j * qnodes), drop = FALSE]
            tmp_j <- sweep(auc_j, 2L, wgt_j, `*`)
            val[,j] <- rowSums(tmp_j)
          }
        } else {
          val <- c()
          for (j in 1:length(eta_m)) {
            wgt_j <- qwts[((j-1) * qnodes + 1):(j * qnodes)]
            auc_j <- auc_m[((j-1) * qnodes + 1):(j * qnodes)]
            val[j] <- sum(wgt_j * auc_j)
          }          
        }
        a_X[[mark]] <- val
        mark <- mark + 1            
      }
      
      #---  muvalue and any interactions  ---#
      
      # muvalue
      if (assoc_m[["muvalue"]]) {
        mu_m <- invlink_m(eta_m) 
        a_X[[mark]] <- mu_m
        mark <- mark + 1            
      }
      
      # muvalue * data interactions
      if (assoc_m[["muvalue_data"]]) {
        mu_m <- invlink_m(eta_m) 
        X_temp <- X_data_m[["muvalue_data"]]
        K_temp <- K_data_m[["muvalue_data"]]
        for (i in 1:K_temp) {
          if (is.matrix(mu_m)) {
            val <- sweep(mu_m, 2L, X_temp[, i], `*`)
          } else {
            val <- as.vector(mu_m) * X_temp[, i]
          }
          if (has_grp) {
            a_X[[mark]] <- collapse_within_groups(val, grp_idx, grp_assoc)
          } else {
            a_X[[mark]] <- val
          }
          mark <- mark + 1
        }
      }
      
      # muvalue * etavalue interactions
      if (assoc_m[["muvalue_etavalue"]]) {
        sel <- assoc_m[["which_interactions"]][["muvalue_etavalue"]]
        for (j in sel) {
          eta_j <- get_element(parts, m = j, "eta", ...)
          val   <- invlink_m(eta_m) * eta_j 
          a_X[[mark]] <- val
          mark <- mark + 1           
        }
      } 
      
      # muvalue * muvalue interactions
      if (assoc_m[["muvalue_muvalue"]]) {
        sel <- assoc_m[["which_interactions"]][["muvalue_muvalue"]]
        for (j in sel) {
          eta_j <- get_element(parts, m = j, "eta", ...)
          invlink_j <- family[[j]]$linkinv
          val <- invlink_m(eta_m) * invlink_j(eta_j) 
          a_X[[mark]] <- val
          mark <- mark + 1                   
        }
      }
      
      #---  muslope and any interactions  ---#
      
      if (assoc_m[["muslope"]] || assoc_m[["muslope_data"]]) {
        if (eps_uses_derivative_of_x) {
          stop2("Cannot currently use muslope interaction structure.")
        } else {
          dmu_m <- (invlink_m(eps_m) - invlink_m(eta_m)) / epsilon
        }
      }
      
      # muslope
      if (assoc_m[["muslope"]]) {
        a_X[[mark]] <- dmu_m
        mark <- mark + 1                   
      }
      
      # muslope * data interactions
      if (assoc_m[["muslope_data"]]) {
        X_temp <- X_data_m[["muslope_data"]]
        K_temp <- K_data_m[["muslope_data"]]
        for (i in 1:K_temp) {
          if (is.matrix(dmu_m)) {
            val <- sweep(dmu_m, 2L, X_temp[, i], `*`)
          } else {
            val <- as.vector(dmu_m) * X_temp[, i]
          }
          if (has_grp) {
            a_X[[mark]] <- collapse_within_groups(val, grp_idx, grp_assoc)
          } else {
            a_X[[mark]] <- val
          }
          mark <- mark + 1
        }
      }
      
      #---  muauc  ---#
      
      if (assoc_m[["muauc"]]) {
        if (is.matrix(eta_m)) {
          nr <- nrow(eta_m)
          nc <- ncol(eta_m)
          val   <- matrix(NA, nrow = nr, ncol = nc) 
          for (j in 1:nc) {
            wgt_j <- qwts[((j-1) * qnodes + 1):(j * qnodes)]
            auc_j <- invlink_m(auc_m[, ((j-1) * qnodes + 1):(j * qnodes), drop = FALSE])
            tmp_j <- sweep(auc_j, 2L, wgt_j, `*`)
            val[,j] <- rowSums(tmp_j)
          }
        } else {
          val <- c()
          for (j in 1:length(eta_m)) {
            wgt_j <- qwts[((j-1) * qnodes + 1):(j * qnodes)]
            auc_j <- invlink_m(auc_m[((j-1) * qnodes + 1):(j * qnodes)])
            val[j] <- sum(wgt_j * auc_j)
          }          
        }
        a_X[[mark]] <- val
        mark <- mark + 1 
      }
      
    }
  }
  for (m in 1:M) {
    # shared_b
    if (assoc["shared_b",][[m]]) {
      sel <- assoc["which_b_zindex",][[m]]
      val <- get_element(parts, m = m, "b_mat", ...)[,sel]
      a_X[[mark]] <- val
      mark <- mark + 1                   
    }
  }    
  for (m in 1:M) {
    # shared_coef
    if (assoc["shared_coef",][[m]]) {
      sel <- assoc["which_coef_zindex",][[m]]
      val <- get_element(parts, m = m, "b_mat", ...)[,sel]
      a_X[[mark]] <- val
      mark <- mark + 1                   
    }
  }
  
  if (is.matrix(a_X[[1L]])) a_X else do.call("cbind", a_X)
}

# Function to get an "element" (e.g. a linear predictor, a linear predictor
# evaluated at epsilon shift, linear predictor evaluated at auc quadpoints, 
# etc) constructed from the "parts" (e.g. mod_eta, mod_eps, mod_auc, etc) 
# returned by a call to the function 'make_assoc_parts'.
#
# @param parts A named list containing the parts for constructing the association 
#   structure. It may contain elements $mod_eta, $mod_eps, $mod_auc, etc. as 
#   well as $X_data, $K_data, $grp_stuff. It is returned by a call to the
#   function 'make_assoc_parts'.
# @param m An integer specifying which submodel to get the element for.
# @param which A character string specifying which element to get.
get_element <- function(parts, m = 1, which = "eta", ...) {
  
  ok_which_args <- c("eta", "eps", "auc", "X_data", "K_data", 
                     "b_mat", "grp_stuff")
  if (!which %in% ok_which_args)
    stop("'which' must be one of: ", paste(ok_which_args, collapse = ", "))
  
  if (which %in% c("eta", "eps", "auc")) {
    part <- parts[[m]][[paste0("mod_", which)]]
    if (is.null(part)) { 
      # model doesn't include an assoc related to 'which'
      return(NULL)
    } else { 
      # construct linear predictor for the 'which' part
      x <- part$x
      Zt <- part$Zt
      Znames  <- part$Z_names
      if (is.null(x) || is.null(Zt))
        stop2("Bug found: cannot find x and Zt in 'parts'. They are ",
              "required to build the linear predictor for '", which, "'.")
      
      dots <- list(...)
      beta <- dots$beta[[m]]
      b    <- dots$b[[m]]
      if (is.null(beta) || is.null(b))
        stop2("Bug found: beta and b must be provided to build the ",
              "linear predictor for '", which, "'.")
      
      eta <- linear_predictor(beta, x)
      if (NCOL(b) == 1) {
        eta <- eta + as.vector(b %*% Zt)
      } else {
        eta <- eta + as.matrix(b %*% Zt)
      }
      return(eta)
    }
  } else if (which %in% c("X_data", "K_data", "b_mat", "grp_stuff")) {
    return(parts[[m]][[which]])
  } else {
    stop("'which' argument doesn't include a valid entry.")
  }
}

# Collapse the linear predictor across the lower level units
# clustered an individual, using the function specified in the
# 'grp_assoc' argument
#
# @param eta The linear predictor evaluated for all lower level groups
#   at the quadrature points.
# @param grp_idx An N*2 array providing the indices of the first (col 1)
#   and last (col 2) observations in eta that correspond to individuals
#   i = 1,...,N.
# @param grp_assoc Character string, the function to use to collapse
#   across the lower level units clustered within individuals.
# @return A vector or matrix, depending on the method called.
collapse_within_groups <- function(eta, grp_idx, grp_assoc = "sum") {
  UseMethod("collapse_within_groups")
}
collapse_within_groups.default <- function(eta, grp_idx, grp_assoc) {
  N <- nrow(grp_idx)
  val <- rep(NA, N)
  for (n in 1:N) {
    tmp <- eta[grp_idx[n,1]:grp_idx[n,2]]
    val[n] <- do.call(grp_assoc, list(tmp))
  }
  val
}
collapse_within_groups.matrix <- function(eta, grp_idx, grp_assoc) {
  N <- nrow(grp_idx)
  val <- matrix(NA, nrow = nrow(eta), ncol = N)
  for (n in 1:N) {
    tmp <- eta[, grp_idx[n,1]:grp_idx[n,2], drop = FALSE]
    val[,n] = apply(tmp, 1L, grp_assoc)
  }
  val
}
