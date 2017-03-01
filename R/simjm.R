# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2013, 2014, 2015, 2016 Trustees of Columbia University
# Copyright (C) 2016 Monash University
# Copyright (C) 2016 Margarita Moreno-Betancur
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

#' Simulate data for a univariate or multivariate joint model
#' 
#' Returns a data frame containing data simulated from a joint model for
#' longitudinal and time-to-event data.
#' 
#' @export
#' @keywords internal
#' 
#' @param n Number of subjects
#' @param M Number of longitudinal markers
#' @param fixed_trajectory The desired type of trajectory in the fixed effects
#'   part of the longitudinal model
#' @param random_trajectory The desired type of trajectory in the random effects 
#'   part of the longitudinal model
#' @param assoc The desired type of association structure   
#' @param betaLong_intercept True intercept in the longitudinal submodel
#' @param betaLong_binary True coefficient for the binary covariate in the 
#'   longitudinal submodel
#' @param betaLong_continuous True coefficient for the continuous covariate in the 
#'   longitudinal submodel   
#' @param betaLong_slope True coefficient for the fixed effect slope in the 
#'   longitudinal submodel
#' @param betaEvent_intercept True intercept term (log hazard scale) in the 
#'   event submodel
#' @param betaEvent_binary True coefficient (log hazard ratio) for the binary 
#'   covariate in the event submodel
#' @param betaEvent_continuous True coefficient (log hazard ratio) for the 
#'   continuous covariate in the event submodel   
#' @param betaEvent_assoc True association parameter (log hazard ratio) in the 
#'   event submodel
#' @param b_sd Vector of standard deviations for the random effects
#' @param b_rho Correlation between the random effects
#' @param error_sd Standard deviation of residual error terms in the longitudinal
#'   submodel
#' @param max_yobs The maximum allowed number of longitudinal measurements. The 
#'   actual number of observed measurements will depend on the individuals event time.
#' @param max_fuptime The maximum follow up time in whatever the desired time 
#'   units are. This time will also be used as the censoring time (i.e. for subjects
#'   who have a simulated survival time that is after \code{max_fuptime})  
#' 
simjm <- function(n = 200, M = 3,
                  fixed_trajectory = c("linear", "none"),
                  random_trajectory = c("linear", "none"),
                  assoc = c("etavalue"),
                  betaLong_intercept = 90, 
                  betaLong_binary = -1.5, 
                  betaLong_continuous = 1, 
                  betaLong_slope = 2.5, 
                  betaEvent_intercept = -11.87,
                  betaEvent_binary = 0.6,
                  betaEvent_continuous = 0.08,
                  betaEvent_assoc = 0.01,
                  b_sd = c(20,3), b_rho = 0.5,
                  error_sd = 10,
                  max_yobs = 8, 
                  max_fuptime = 10)
{
  
  fixed_trajectory  <- match.arg(fixed_trajectory)
  random_trajectory <- match.arg(random_trajectory)
  assoc             <- match.arg(assoc)             # only etavalue currently implemented
  
  if (random_trajectory == "linear" && fixed_trajectory == "none")
    stop("Cannot use a linear random slope without a fixed linear slope.")
  
  betaLong_intercept  <- maybe_broadcast(betaLong_intercept,  M)
  betaLong_binary     <- maybe_broadcast(betaLong_binary,     M)
  betaLong_continuous <- maybe_broadcast(betaLong_continuous, M)
  betaLong_slope      <- maybe_broadcast(betaLong_slope,      M)
  betaEvent_assoc     <- maybe_broadcast(betaEvent_assoc,     M)
  error_sd            <- maybe_broadcast(error_sd,            M)

  weibull_shape <- 2  # must be 2, for correct closed-form solution to integral in cumhCox function
    
  # Generate baseline covariates - binary
  prob_Z1 <- 0.45                  # probability for binary covariate
  Z1 <- rbinom(n, 1, prob_Z1)      # covariate value for each subject

  # Generate baseline covariates - continuous
  mean_Z2 <- 44                    # mean for continuous covariate
  sd_Z2   <- 8.5                   # sd for continuous covariate
  Z2 <- rnorm(n, mean_Z2, sd_Z2)   # covariate value for each subject
  
  # Generate subject-specific random effects
  b_dim_perM <- switch(random_trajectory,
                       none   = 1L, # random intercepts model
                       linear = 2L, # random slopes model
                       poly   = 3L) # random poly (degree = 2) model
  b_dim <- M * b_dim_perM           # total num of random effects
  if (length(b_sd) == b_dim_perM) {
    b_sd <- rep(b_sd, times = M)
  } else if (length(b_sd) != b_dim) {
    stop("b_sd appears to be the wrong length.")
  }
  corr_mat <- matrix(rep(b_rho, b_dim ^ 2), ncol = b_dim)
  diag(corr_mat) <- 1
  dd <- MASS::mvrnorm(n = n, mu = rep(0, b_dim), Sigma = corr_mat)
  b <- sapply(1:length(b_sd), function(x) b_sd[x] * dd[,x])
  colnames(b) <- paste0("b", 1:b_dim)
  
  # Randomly generate observation times between 0 and max_fuptime
  tij <- runif(n * max_yobs, 0, max_fuptime)
  
  # Residual error terms
  eij <- sapply(error_sd, function(x) rnorm(n * max_yobs, 0, x))
  colnames(eij) <- paste0("eij_", 1:M)
  
  # Construct data frame
  dat.id <- data.frame(id = 1:n, b, Z1, Z2)         # single row per subject
  dat <- dat.id[rep(row.names(dat.id), max_yobs), ] # multiple row per subject
  dat <- data.frame(dat, tij, eij)                  # merge on measurement times and errors
  dat <- dat[order(dat$id, dat$tij), ]              # sort on ID and time
  
  # Calculate longitudinal outcomes
  for (m in 1:M) {
    nm_intercept <- paste0("b", (m - 1) * b_dim_perM + 1)
    dat[[paste0("Xij_", m)]] <- 
      betaLong_intercept[m] + 
      betaLong_binary[m]     * dat$Z1 + 
      betaLong_continuous[m] * dat$Z2 +
      dat[[nm_intercept]]
    if (fixed_trajectory == "linear")
      dat[[paste0("Xij_", m)]] <- dat[[paste0("Xij_", m)]] + (betaLong_slope[m] * dat$tij)
    if (random_trajectory == "linear") {
      nm_slope <- paste0("b", (m-1) * b_dim_perM + 2)
      dat[[paste0("Xij_", m)]] <- dat[[paste0("Xij_", m)]] + (dat[[nm_slope]] * dat$tij) 
    }  
    dat[[paste0("Yij_", m)]] <- dat[[paste0("Xij_", m)]] + dat[[paste0("eij_", m)]]
  }
  
  # Sum parameters for time-fixed part of each longitudinal submodel
  # for feeding into the calculation of the cumulative hazard
  time_fixed_part <- lapply(1:M, function(m) {
    nm_intercept <- paste0("b", (m - 1) * b_dim_perM + 1)
    betaLong_intercept[m] + 
    betaLong_binary[m]     * dat.id$Z1 + 
    betaLong_continuous[m] * dat.id$Z2 +
    dat.id[[nm_intercept]]
  })
  
  # Sum parameters for time-varying part of each longitudinal submodel 
  # for feeding into the calculation of the cumulative hazard
  time_varying_part <- lapply(1:M, function(m) {
    val <- 0
    if (fixed_trajectory == "linear")
      val <- val + betaLong_slope[m] 
    if (random_trajectory == "linear") {
      nm_slope <- paste0("b", (m-1) * b_dim_perM + 2)
      val <- val + dat.id[[nm_slope]]
    }
    val
  })  

  # Random uniform variable used for generating survival time
  u <- runif(n) 
  
  # Time-fixed and time-varying parts used in calculating cumulative hazard
  K1 <- weibull_shape * exp(betaEvent_intercept + 
                            betaEvent_binary     * dat.id$Z1 + 
                            betaEvent_continuous * dat.id$Z2)
  for (i in 1:4)
    if (M > (i-1)) K1 <- K1 * exp(betaEvent_assoc[i] * time_fixed_part[[i]]) # equivalent to adding terms onto the linear predictor
  K2 <- if (M > 0) (betaEvent_assoc[1] * time_varying_part[[1]])
  for (i in 2:4)
    if (M > (i-1)) K2 <- K2 + (betaEvent_assoc[i] * time_varying_part[[i]])
  parts <- cbind(u, K1, K2)
  
  # Calculate survival time under Weibull model
  true_eventtime <- apply(parts, 1, function(x) {
    if (cumhCox(500, x[1], x[2], x[3]) < 0) return(202) else
      stats::uniroot(cumhCox, x[1], x[2], x[3], interval = c(0,500))$root
  })
  
  # Observed event times and event indicator
  event     <- (true_eventtime <= max_fuptime)
  eventtime <- (true_eventtime * event) + (max_fuptime * (1 - event))
  eventdat  <- data.frame(id = 1:n, eventtime, event)
  
  # Final dataset
  ret <- merge(dat, eventdat, by = "id")
  ret <- ret[ret$tij <= ret$eventtime, ]  # only keep rows before event time
  sel <- grep("^id|^Z|^tij|^Yij|event", colnames(ret))
  ret <- ret[, sel, drop = FALSE]
  
  # Store 'true' parameter values
  long_params <- nlist(
    betaLong_intercept, 
    betaLong_binary, 
    betaLong_continuous
  )
  if (fixed_trajectory == "linear") 
    long_params$betaLong_slope <- betaLong_slope  
  event_params <- nlist( 
    betaEvent_intercept,
    betaEvent_binary,
    betaEvent_continuous,
    betaEvent_assoc,
    weibull_shape
  )
  re_params <- nlist(
    b_sd,
    b_corr = corr_mat
  )
  
  # Return object
  structure(ret, params = c(long_params, event_params, re_params), 
            n = n, M = M, max_yobs = max_yobs, max_fuptime = max_fuptime, assoc = assoc,
            fixed_trajectory = fixed_trajectory, random_trajectory = random_trajectory)
} 
  
  
#--------------------- internal  
  
#' Cumulative hazard for Weibull proportional hazards model
#'
#' Internal function used to generate cumulative hazard function based on a 
#' proportional hazards model with Weibull baseline hazard with shape parameter 
#' equal to 2.
#' 
#' @export
#' @keywords internal
#' 
#' @author Margarita Moreno-Betancur <margarita.moreno@mcri.edu.au>
#' 
#' @param tt Upper limit of integral
#' @param u Random uniform variate on range 0 to 1
#' @param K1 The time-fixed contribution to the hazard
#' @param K2 Sum of parameters contained in the time-varying contribution to the  
#'   hazard
#'   
cumhCox <- function(tt,u,K1,K2){
  if (K2 != 0) {
    return(K1 * ((tt*exp(K2*tt)/K2) - (exp(K2*tt)/(K2^2)) + (1/(K2^2))) + log(u))
  } else {
    return(K1 * (tt ^ 2) / 2 + log(u))
  }
}
