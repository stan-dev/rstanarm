# Part of the rstanarm package for estimating model parameters
# Copyright 2015 Douglas Bates, Martin Maechler, Ben Bolker, Steve Walker
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

pp_data <-
  function(object,
           newdata = NULL,
           re.form = NULL,
           offset = NULL,
           m = NULL,
           ...) {
    validate_stanreg_object(object)
    if (is.mer(object)) {
      if (is.nlmer(object))
        out <- .pp_data_nlmer(object, newdata = newdata, re.form = re.form, m = m, ...)
      else
        out <- .pp_data_mer(object, newdata = newdata, re.form = re.form, m = m, ...)
      if (!is.null(offset)) out$offset <- offset
      return(out)
    }
    .pp_data(object, newdata = newdata, offset = offset, ...)
  }

# for models without lme4 structure
.pp_data <- function(object, newdata = NULL, offset = NULL, ...) {
  if (is(object, "gamm4")) {
    requireNamespace("mgcv", quietly = TRUE)
    if (is.null(newdata))   x <- predict(object$jam, type = "lpmatrix")
    else x <- predict(object$jam, newdata = newdata, type = "lpmatrix")
    if (is.null(offset)) 
      offset <- object$offset %ORifNULL% rep(0, nrow(x))
    return(nlist(x, offset))
  }
  if (is.null(newdata)) {
    x <- get_x(object)
    if (is.null(offset)) {
      offset <- object$offset %ORifNULL% rep(0, nrow(x))
    }
    if (inherits(object, "betareg")) {
      return(nlist(x, offset, z_betareg = object$z))
    }
    
    return(nlist(x, offset))
  }

  offset <- .pp_data_offset(object, newdata, offset)
  Terms <- delete.response(terms(object))
  m <- model.frame(Terms, newdata, xlev = object$xlevels)
  if (!is.null(cl <- attr(Terms, "dataClasses"))) 
    .checkMFClasses(cl, m)
  x <- model.matrix(Terms, m, contrasts.arg = object$contrasts)
  if (is(object, "polr") && !is_scobit(object)) 
    x <- x[,colnames(x) != "(Intercept)", drop = FALSE]
  
  if (inherits(object, "betareg")) {
    mf <- model.frame(delete.response(object$terms$precision), 
                      data = newdata, na.action = object$na.action, 
                      xlev = object$levels$precision)
    z_betareg <- model.matrix(object$terms$precision, mf, contrasts = object$contrasts$precision)
    return(nlist(x, offset, z_betareg))
  }
  
  return(nlist(x, offset))
}


# for models fit using stan_(g)lmer or stan_gamm4
.pp_data_mer <- function(object, newdata, re.form, m = NULL, ...) {
  if (is(object, "gamm4")) {
    requireNamespace("mgcv", quietly = TRUE)
    if (is.null(newdata))   x <- predict(object$jam, type = "lpmatrix")
    else x <- predict(object$jam, newdata = newdata, type = "lpmatrix")
    if (is.null(re.form)) {
      re.form <- as.formula(object$call$random)
      if (length(re.form) == 0) re.form <- NA
      z <- .pp_data_mer_z(object, newdata, re.form, ...)
    }
    else z <- .pp_data_mer_z(object, newdata, re.form, ...)
  } else {
    x <- .pp_data_mer_x(object, newdata, m = m, ...)
    z <- .pp_data_mer_z(object, newdata, re.form, m = m, ...)
  }
  offset <- model.offset(model.frame(object, m = m))
  if (!missing(newdata) && (!is.null(offset) || !is.null(object$call$offset))) {
    offset <- try(eval(object$call$offset, newdata), silent = TRUE)
    if (!is.numeric(offset)) offset <- NULL
  }
  return(nlist(x, offset = offset, Zt = z$Zt, Z_names = z$Z_names))
}

# for models fit using stan_nlmer
.pp_data_nlmer <- function(object, newdata, re.form, offset = NULL, m = NULL, ...) {
  inputs <- parse_nlf_inputs(object$glmod$respMod)
  if (is.null(newdata)) {
    arg1 <- arg2 <- NULL
  } else if (object$family$link == "inv_SSfol") {
    arg1 <- newdata[[inputs[2]]]
    arg2 <- newdata[[inputs[3]]]
  } else {
    arg1 <- newdata[[inputs[2]]]
    arg2 <- NULL
  }
  f <- formula(object, m = m)
  if (!is.null(re.form) && !is.na(re.form)) {
    f <- as.character(f)
    f[3] <- as.character(re.form)
    f <- as.formula(f[-1])
  }
  if (is.null(newdata)) newdata <- model.frame(object)
  else {
    yname <- names(model.frame(object))[1]
    newdata[[yname]] <- 0
  }
  mc <- match.call(expand.dots = FALSE)
  mc$re.form <- mc$offset <- mc$object <- mc$newdata <- NULL
  mc$data <- newdata
  mc$formula <- f
  mc$start <- fixef(object)
  nlf <- nlformula(mc)
  offset <- .pp_data_offset(object, newdata, offset)

  group <- with(nlf$reTrms, pad_reTrms(Ztlist, cnms, flist))
  if (!is.null(re.form) && !is(re.form, "formula") && is.na(re.form)) 
    group$Z@x <- 0
  return(nlist(x = nlf$X, offset = offset, Z = group$Z,
               Z_names = make_b_nms(group), arg1, arg2))
}

# the functions below are heavily based on a combination of 
# lme4:::predict.merMod and lme4:::mkNewReTrms, although they do also have 
# substantial modifications
.pp_data_mer_x <- function(object, newdata, m = NULL, ...) {
  x <- get_x(object, m = m)
  if (is.null(newdata)) return(x)
  form <- if (is.null(m)) attr(object$glmod$fr, "formula") else 
    formula(object, m = m)
  L <- length(form)
  form[[L]] <- lme4::nobars(form[[L]])
  RHS <- formula(substitute(~R, list(R = form[[L]])))
  Terms <- terms(object, m = m)
  mf <- model.frame(object, m = m)
  ff <- formula(form)
  vars <- rownames(attr(terms.formula(ff), "factors"))
  mf <- mf[vars]
  isFac <- vapply(mf, is.factor, FUN.VALUE = TRUE)
  isFac[attr(Terms, "response")] <- FALSE
  orig_levs <- if (length(isFac) == 0) 
    NULL else lapply(mf[isFac], levels)
  mfnew <- model.frame(delete.response(Terms), newdata, xlev = orig_levs)
  x <- model.matrix(RHS, data = mfnew, contrasts.arg = attr(x, "contrasts"))
  return(x)
}

.pp_data_mer_z <- function(object, newdata, re.form = NULL,
                           allow.new.levels = TRUE, na.action = na.pass, 
                           m = NULL, ...) {
  NAcheck <- !is.null(re.form) && !is(re.form, "formula") && is.na(re.form)
  fmla0check <- (is(re.form, "formula") && 
                   length(re.form) == 2 && 
                   identical(re.form[[2]], 0))
  if (NAcheck || fmla0check) return(list())
  if (is.null(newdata) && is.null(re.form)) {
    Z <- get_z(object, m = m) 
    if (!is.stanmvreg(object)) {
      # Z_names not needed for stanreg with no newdata
      return(list(Zt = t(Z)))
    } else {
      # must supply Z_names for stanmvreg since b pars
      # might be for multiple submodels and Zt will only
      # be for one submodel, so their elements may not 
      # correspond exactly
      ReTrms <- object$glmod[[m]]$reTrms
      Z_names <- make_b_nms(ReTrms, m = m, stub = get_stub(object))
      return(nlist(Zt = ReTrms$Zt, Z_names))
    }
  }
  else if (is.null(newdata)) {
    rfd <- mfnew <- model.frame(object, m = m)
  } 
  else if (inherits(object, "gamm4")) {
    requireNamespace("mgcv", quietly = TRUE)
    if (is.null(newdata))   x <- predict(object$jam, type = "lpmatrix")
    else x <- predict(object$jam, newdata = newdata, type = "lpmatrix")
    NAs <- apply(is.na(x), 1, any)
    rfd <- mfnew <- newdata[!NAs,, drop=FALSE]
    attr(rfd,"na.action") <- "na.omit"
  } else {
    terms_fixed <- delete.response(terms(object, fixed.only = TRUE, m = m))
    mfnew <- model.frame(terms_fixed, newdata, na.action = na.action)
    newdata.NA <- newdata
    if (!is.null(fixed.na.action <- attr(mfnew,"na.action"))) {
      newdata.NA <- newdata.NA[-fixed.na.action,]
    }
    tt <- delete.response(terms(object, random.only = TRUE, m = m))
    rfd <- model.frame(tt, newdata.NA, na.action = na.pass)
    if (!is.null(fixed.na.action))
      attr(rfd,"na.action") <- fixed.na.action
  }
  if (is.null(re.form)) 
    re.form <- justRE(formula(object, m = m))
  if (!inherits(re.form, "formula"))
    stop("'re.form' must be NULL, NA, or a formula.")
  if (length(fit.na.action <- attr(mfnew,"na.action")) > 0) {
    newdata <- newdata[-fit.na.action,]
  }
  ReTrms <- lme4::mkReTrms(lme4::findbars(re.form[[2]]), rfd)
  if (!allow.new.levels && any(vapply(ReTrms$flist, anyNA, NA)))
    stop("NAs are not allowed in prediction data",
         " for grouping variables unless 'allow.new.levels' is TRUE.")
  ns.re <- names(re <- ranef(object, m = m))
  nRnms <- names(Rcnms <- ReTrms$cnms)
  if (!all(nRnms %in% ns.re))
    stop("Grouping factors specified in re.form that were not present in original model.")
  new_levels <- lapply(ReTrms$flist, function(x) levels(factor(x)))
  Zt <- ReTrms$Zt
  Z_names <- make_b_nms(ReTrms, m = m, stub = get_stub(object))
  z <- nlist(Zt = ReTrms$Zt, Z_names)
  return(z)
}



# handle offsets ----------------------------------------------------------
null_or_zero <- function(x) {
  isTRUE(is.null(x) || all(x == 0))
}

.pp_data_offset <- function(object, newdata = NULL, offset = NULL) {
  if (is.null(newdata)) {
    # get offset from model object (should be null if no offset)
    if (is.null(offset)) 
      offset <- object$offset %ORifNULL% model.offset(model.frame(object))
  } else {
    if (!is.null(offset))
      stopifnot(length(offset) == nrow(newdata))
    else {
      # if newdata specified but not offset then confirm that model wasn't fit
      # with an offset (warning, not error)
      if (!is.null(object$call$offset) || 
          !null_or_zero(object$offset) || 
          !null_or_zero(model.offset(model.frame(object)))) {
        warning(
          "'offset' argument is NULL but it looks like you estimated ", 
          "the model using an offset term.", 
          call. = FALSE
        )
      }
      offset <- rep(0, nrow(newdata))
    }
  }
  return(offset)
}


#----------------------- pp_data for joint models --------------------------

# Return the design matrices required for evaluating the linear predictor or
# log-likelihood in post-estimation functions for a \code{stan_jm} model
#
# @param object A stanmvreg object
# @param newdataLong A data frame or list of data frames with the new 
#   covariate data for the longitudinal submodel
# @param newdataEvent A data frame with the new covariate data for the
#   event submodel
# @param ids An optional vector of subject IDs specifying which individuals
#   should be included in the returned design matrices.
# @param etimes An optional vector of times at which the event submodel
#   design matrices should be evaluated (also used to determine the 
#   quadrature times). If NULL then times are taken to be the eventimes in
#   the fitted object (if newdataEvent is NULL) or in newdataEvent.
# @param long_parts,event_parts A logical specifying whether to return the
#   design matrices for the longitudinal and/or event submodels.
# @return A named list (with components M, Npat, ndL, ndE, yX, tZt, 
#   yZnames, eXq, assoc_parts) 
.pp_data_jm <- function(object, newdataLong = NULL, newdataEvent = NULL, 
                        ids = NULL, etimes = NULL, long_parts = TRUE, 
                        event_parts = TRUE) {
  M <- get_M(object)
  id_var   <- object$id_var
  time_var <- object$time_var
  
  if (!is.null(newdataLong) || !is.null(newdataEvent))
    newdatas <- validate_newdatas(object, newdataLong, newdataEvent)
  
  # prediction data for longitudinal submodels
  ndL <- if (is.null(newdataLong)) 
    get_model_data(object)[1:M] else newdatas[1:M]
  
  # prediction data for event submodel
  ndE <- if (is.null(newdataEvent)) 
    get_model_data(object)[["Event"]] else newdatas[["Event"]]   
  
  # possibly subset
  if (!is.null(ids)) {
    ndL <- subset_ids(object, ndL, ids)
    ndE <- subset_ids(object, ndE, ids)
  }
  id_list <- unique(ndE[[id_var]]) # unique subject id list
  
  # evaluate the last known survival time and status
  if (!is.null(newdataEvent) && is.null(etimes)) {
    # prediction data for the event submodel was provided but  
    # no event times were explicitly specified by the user, so
    # they must be evaluated using the data frame
    surv <- eval(formula(object, m = "Event")[[2L]], ndE)
    etimes  <- unclass(surv)[,"time"]
    estatus <- unclass(surv)[,"status"]    
  } else if (is.null(etimes)) {
    # if no prediction data was provided then event times are 
    # taken from the fitted model
    etimes  <- object$eventtime[as.character(id_list)]
    estatus <- object$status[as.character(id_list)]
  } else { 
    # otherwise, event times ('etimes') are only directly specified for dynamic   
    # predictions via posterior_survfit in which case the 'etimes' correspond 
    # to the last known survival time and therefore we assume everyone has survived
    # up to that point (ie, set estatus = 0 for all individuals), this is true 
    # even if there is an event indicated in the data supplied by the user.
    estatus <- rep(0, length(etimes))
  }
  res <- nlist(M, Npat = length(id_list), ndL, ndE)
  
  if (long_parts && event_parts)
    lapply(ndL, function(x) {
      if (!time_var %in% colnames(x)) 
        STOP_no_var(time_var)
      if (!id_var %in% colnames(x)) 
        STOP_no_var(id_var)
      if (any(x[[time_var]] < 0))
        stop2("Values for the time variable (", time_var, ") should not be negative.")
      mt <- tapply(x[[time_var]], factor(x[[id_var]]), max)
      if (any(mt > etimes))
        stop2("There appears to be observation times in the longitudinal data that ",
              "are later than the event time specified in the 'etimes' argument.")      
    }) 
  
  # response and design matrices for longitudinal submodels
  if (long_parts) {
    y <- lapply(1:M, function(m) eval(formula(object, m = m)[[2L]], ndL[[m]]))
    ydat <- lapply(1:M, function(m) pp_data(object, ndL[[m]], m = m))
    yX <- fetch(ydat, "x")
    yZt <- fetch(ydat, "Zt")
    yZ_names <- fetch(ydat, "Z_names")
    flist <- lapply(ndL, function(x) factor(x[[id_var]]))
    res <- c(res, nlist(y, yX, yZt, yZ_names, flist))
  }
  
  # design matrices for event submodel and association structure
  if (event_parts) {
    qnodes <- object$qnodes
    qq <- get_quadpoints(qnodes)
    qtimes <- uapply(qq$points,  unstandardise_qpts, 0, etimes)
    qwts   <- uapply(qq$weights, unstandardise_qwts, 0, etimes)
    starttime <- deparse(formula(object, m = "Event")[[2L]][[2L]])
    edat <- prepare_data_table(ndE, id_var, time_var = starttime)
    id_rep <- rep(id_list, qnodes + 1)
    times <- c(etimes, qtimes) # times used to design event submodel matrices
    edat <- rolling_merge(edat, ids = id_rep, times = times)
    eXq  <- .pp_data_mer_x(object, newdata = edat, m = "Event")
    assoc_parts <- lapply(1:M, function(m) {
      ymf <- ndL[[m]]
      grp_stuff <- object$grp_stuff[[m]]
      if (grp_stuff$has_grp) {
        grp_stuff <- get_extra_grp_info( # update grp_info with new data
          grp_stuff, flist = ymf, id_var = id_var, 
          grp_assoc = grp_stuff$grp_assoc)
      }
      ymf <- prepare_data_table(ymf, id_var = id_var, time_var = time_var,
                                grp_var = grp_stuff$grp_var)
      make_assoc_parts(
        ymf, assoc = object$assoc[,m], id_var = id_var, time_var = time_var, 
        ids = id_rep, times = times, grp_stuff = grp_stuff,
        use_function = pp_data, object = object, m = m)
    })
    assoc_attr <- nlist(.Data = assoc_parts, qnodes, qtimes, qwts, etimes, estatus)
    assoc_parts <- do.call("structure", assoc_attr)
    res <- c(res, nlist(eXq, assoc_parts))
  }
  
  return(res)
}

# Return a data frame for each submodel that:
# (1) only includes variables used in the model formula
# (2) only includes rows contained in the glmod/coxmod model frames
# (3) ensures that additional variables that are required
#     such as the ID variable or variables used in the 
#     interaction-type association structures, are included.
#
# It is necessary to drop unneeded variables though so that 
# errors are not encountered if the original data contained 
# NA values for variables unrelated to the model formula.
# We generate a data frame here for in-sample predictions 
# rather than using a model frame, since some quantities will
# need to be recalculated at quadrature points etc, for example
# in posterior_survfit.
#
# @param object A stanmvreg object.
# @param m Integer specifying which submodel to get the
#   prediction data frame for.
# @return A data frame or list of data frames with all the
#   (unevaluated) variables required for predictions.
get_model_data <- function(object, m = NULL) {
  validate_stanmvreg_object(object)
  M <- get_M(object)
  terms <- terms(object, fixed.only = FALSE)
  
  # identify variables to add to the terms objects
  if (is.jm(object)) {
    extra_vars <- lapply(1:M, function(m) {
      # for each submodel loop over the four possible assoc  
      # interaction formulas and collect any variables used
      forms_m <- object$assoc["which_formulas",][[m]]
      uapply(forms_m, function(x) {
        if (length(x)) {
          rownames(attr(terms.formula(x), "factors")) 
        } else NULL
      })
    })
    # also ensure that id_var is in the event data
    extra_vars$Event <- object$id_var
    
    if (!identical(length(terms), length(extra_vars)))
      stop2("Bug found: terms and extra_vars should be same length.")
    
    # add the extra variables to the terms formula for each submodel
    terms <- xapply(terms, extra_vars, FUN = function(x, y) {
      lhs <- x[[2L]]
      rhs <- deparse(x[[3L]], 500L)
      if (!is.null(y))
        rhs <- c(rhs, y)
      reformulate(rhs, response = lhs)
    })
    
    datas <- c(object$dataLong, list(object$dataEvent))
  } else {
    datas <- object$data
  }
  
  # identify rows that were in the model frame
  row_nms <- lapply(model.frame(object), rownames)
  
  # drop rows and variables not required for predictions
  mfs <- xapply(w = terms, x = datas, y = row_nms,
                FUN = function(w, x, y) 
                  get_all_vars(w, x)[y, , drop = FALSE])
  
  mfs <- list_nms(mfs, M, stub = get_stub(object))
  if (is.null(m)) mfs else mfs[[m]]
}
