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

# Construct a list with information about the event submodel
#
# @param formula The model formula for the event submodel
# @param data The data for the event submodel
# @param qnodes An integer specifying the number of GK quadrature nodes
# @param id_var The name of the ID variable
# @param y_id_list A character vector with a unique list of subject IDs 
#   (factor levels) that appeared in the longitudinal submodels
# @return A named list with the following elements:
#   mod: The fitted Cox model.
#   entrytime: Named vector of numeric entry times.
#   eventtime: Named vector of numeric event times.
#   status: Named vector of event/failure indicators.
#   Npat: Number of individuals.
#   Nevents: Total number of events/failures.
#   id_list: A vector of unique subject IDs, as a factor.
#   qnodes: The number of GK quadrature nodes.
#   qwts,qpts: Vector of unstandardised quadrature weights and points.
#     The vector is ordered such that the first Npat items are the
#     weights/locations of the first quadrature point, then the second
#     Npat items are the weights/locations for the second quadrature
#     point, and so on. 
#   qids: The subject IDs corresponding to each element of qwts/qpts.
#   epts: The event times, but only for individuals who were NOT censored
#     (i.e. those individual who had an event).
#   eids: The subject IDs corresponding to each element of epts.
#   cpts: Combined vector of failure and quadrature times: c(epts, qpts).
#   cids: Combined vector subject IDs: c(eids, qids).
#   Xq: The model matrix for the event submodel, centred and no intercept.
#   Xbar: Vector of column means for the event submodel model matrix.
#   K: Number of predictors for the event submodel.
#   norm_const: Scalar, the constant used to shift the event submodel
#     linear predictor (equal to the log of the mean incidence rate). 
#   model_frame: The model frame for the fitted Cox model, but with the
#     subject ID variable also included.
#   tvc: Logical, if TRUE then a counting type Surv() object was used
#     in the fitted Cox model (ie. time varying covariates). 
handle_e_mod <- function(formula, data, meta) {
  
  if (!requireNamespace("survival"))
    stop("the 'survival' package must be installed to use this function")
  if (!requireNamespace("data.table"))
    stop("the 'data.table' package must be installed to use this function")
  
  id_var      <- meta$id_var
  id_list     <- meta$id_list
  qnodes      <- meta$qnodes
  basehaz     <- meta$basehaz
  basehaz_ops <- meta$basehaz_ops
  
  # parse formula, create model data & frame
  formula   <- parse_formula(formula, data)
  formula2  <- addto_formula(formula$formula, id_var) # includes id_var
  data      <- make_model_data (formula2, data)       # row subsetting etc.
  mf_stuff  <- make_model_frame(formula2, data)       # returns Surv object
  
  mf <- mf_stuff$mf # model frame
  mt <- mf_stuff$mt # model terms
  
  mf[[id_var]] <- promote_to_factor(mf[[id_var]]) # same as lme4
  ids <- factor(mf[[id_var]])
  
  # error checks for the id variable
  validate_jm_ids(y_ids = id_list, e_ids = ids)
  
  # entry and exit times for each row of data
  t_beg <- make_t(mf, type = "beg") # entry time
  t_end <- make_t(mf, type = "end") # exit  time
  t_upp <- make_t(mf, type = "upp") # upper time for interval censoring
  
  # event indicator for each row of data
  d <- make_d(mf)
  
  event <- as.logical(d == 1)  
  rcens <- as.logical(d == 0)
  lcens <- as.logical(d == 2)
  icens <- as.logical(d == 3)
  
  if (any(d < 0 || d > 3))
    stop2("Invalid status indicator in Surv object.")
 
  # delayed entry indicator for each row of data
  delayed  <- as.logical(!t_beg == 0)
  
  # time variables for stan
  t_event <- t_end[event]   # exact event time
  t_lcens <- t_end[lcens]   # left  censoring time
  t_rcens <- t_end[rcens]   # right censoring time
  t_icenl <- t_end[icens]   # lower limit of interval censoring time
  t_icenu <- t_upp[icens]   # upper limit of interval censoring time
  t_delay <- t_beg[delayed] # delayed entry time
  
  # entry and exit times for each individual
  t_tmp <- t_end; t_tmp[icens] <- t_upp[icens]
  entrytime <- tapply(t_beg, ids, min)
  eventtime <- tapply(t_tmp, ids, max)
  status    <- tapply(d,     ids, max)
  
  # dimensions
  nevent <- sum(event)
  nrcens <- sum(rcens)
  nlcens <- sum(lcens)
  nicens <- sum(icens)
  ndelay <- sum(delayed)
  
  # baseline hazard
  ok_basehaz <- c("weibull", "bs", "piecewise")
  ok_basehaz_ops <- get_ok_basehaz_ops(basehaz)
  basehaz <- handle_basehaz(basehaz        = basehaz, 
                            basehaz_ops    = basehaz_ops, 
                            ok_basehaz     = ok_basehaz,
                            ok_basehaz_ops = ok_basehaz_ops,
                            times          = t_end, 
                            status         = event,
                            min_t          = min(t_beg),
                            max_t          = max(c(t_end,t_upp), na.rm = TRUE))
  nvars <- basehaz$nvars # number of basehaz aux parameters
  
  # flag if intercept is required for baseline hazard
  has_intercept <- ai(has_intercept(basehaz))
  
  # standardised weights and nodes for quadrature
  qq <- get_quadpoints(nodes = qnodes)
  qp <- qq$points
  qw <- qq$weights
  
  # quadrature points & weights, evaluated for each row of data
  qpts_event <- uapply(qp, unstandardise_qpts, 0, t_event)
  qpts_lcens <- uapply(qp, unstandardise_qpts, 0, t_lcens)
  qpts_rcens <- uapply(qp, unstandardise_qpts, 0, t_rcens)
  qpts_icenl <- uapply(qp, unstandardise_qpts, 0, t_icenl)
  qpts_icenu <- uapply(qp, unstandardise_qpts, 0, t_icenu)
  qpts_delay <- uapply(qp, unstandardise_qpts, 0, t_delay)
  
  qwts_event <- uapply(qw, unstandardise_qwts, 0, t_event)
  qwts_lcens <- uapply(qw, unstandardise_qwts, 0, t_lcens)
  qwts_rcens <- uapply(qw, unstandardise_qwts, 0, t_rcens)
  qwts_icenl <- uapply(qw, unstandardise_qwts, 0, t_icenl)
  qwts_icenu <- uapply(qw, unstandardise_qwts, 0, t_icenu)
  qwts_delay <- uapply(qw, unstandardise_qwts, 0, t_delay)
  
  eids_event <- ids[event]
  qids_event <- rep(ids[event],   times = qnodes)
  qids_lcens <- rep(ids[lcens],   times = qnodes)
  qids_rcens <- rep(ids[rcens],   times = qnodes)
  qids_icens <- rep(ids[icens],   times = qnodes)
  qids_delay <- rep(ids[delayed], times = qnodes)
  
  # times at events and all quadrature points
  cids <- c(eids_event,
            qids_event,
            qids_lcens,
            qids_rcens,
            qids_icens,
            qids_icens,
            qids_delay)
  cpts_list <- list(t_event,
                    qpts_event,
                    qpts_lcens,
                    qpts_rcens,
                    qpts_icenl,
                    qpts_icenu,
                    qpts_delay)
  idx_cpts <- get_idx_array(sapply(cpts_list, length))
  cpts     <- unlist(cpts_list) # as vector for stan
  len_cpts <- length(cpts)
  
  # number of quadrature points
  qevent <- length(qwts_event)
  qlcens <- length(qwts_lcens)
  qrcens <- length(qwts_rcens)
  qicens <- length(qwts_icenl)
  qdelay <- length(qwts_delay)

  # basis terms for baseline hazard
  basis_cpts <- make_basis(cpts, basehaz) 
  
  # predictor matrices
  x <- make_x(formula$tf_form, mf)$x
  x_event <- keep_rows(x, d == 1)
  x_lcens <- keep_rows(x, d == 2)
  x_rcens <- keep_rows(x, d == 0)
  x_icens <- keep_rows(x, d == 3)
  x_delay <- keep_rows(x, delayed)
  K <- ncol(x)
  x_cpts <- rbind(x_event,
                  rep_rows(x_event, times = qnodes),
                  rep_rows(x_lcens, times = qnodes),
                  rep_rows(x_rcens, times = qnodes),
                  rep_rows(x_icens, times = qnodes),
                  rep_rows(x_delay, times = qnodes))
  
  # fit a cox model
  if (formula$surv_type %in% c("right", "counting")) {
    mod <- survival::coxph(formula$formula, data = data, x = TRUE)
  } else if (formula$surv_type %in% c("interval", "interval2")) {
    mod <- survival::survreg(formula$formula, data = data, x = TRUE)
  } else {
    stop("Bug found: Invalid Surv type.")
  }
  
  # calculate mean log incidence, used as a shift in log baseline hazard
  norm_const <- log(nevent / sum(eventtime - entrytime))
  
  nlist(mod,
        surv_type = formula$surv_type,
        qnodes,
        basehaz,
        has_intercept,
        has_icens = as.logical(nicens),
        model_frame = mf,
        entrytime,
        eventtime,
        d, status,
        norm_const,
        t_beg, 
        t_end,
        t_upp,
        t_event,
        t_rcens,
        t_lcens,
        t_icenl,
        t_icenu,
        t_delay,
        nevent,
        nlcens,
        nrcens,
        nicens,
        ndelay,
        qevent,
        qlcens,
        qrcens,
        qicens,
        qdelay,
        cids,
        cpts,
        len_cpts,
        idx_cpts,
        qwts_event, 
        qwts_lcens, 
        qwts_rcens, 
        qwts_icenl, 
        qwts_icenu,
        qwts_delay, 
        eids_event,
        qids_event,
        qids_lcens,
        qids_rcens,
        qids_icens,
        qids_delay,
        x, 
        x_cpts,
        basis_cpts,
        x_bar = colMeans(x),
        K)
}
