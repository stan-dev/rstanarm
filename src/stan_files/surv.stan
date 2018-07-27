#include /pre/Columbia_copyright.stan
#include /pre/Brilleman_copyright.stan
#include /pre/license.stan

functions {

  /**
   * Hierarchical shrinkage parameterization
   *
   * @param z_beta A vector of primitive coefficients
   * @param global A real array of positive numbers
   * @param local A vector array of positive numbers
   * @param global_prior_scale A positive real number
   * @param error_scale 1 or sigma in the Gaussian case
   * @param c2 A positive real number
   * @return A vector of coefficientes
   */
  vector hs_prior(vector z_beta, real[] global, vector[] local,
                  real global_prior_scale, real error_scale, real c2) {
    int K = rows(z_beta);
    vector[K] lambda = local[1] .* sqrt(local[2]);
    real tau = global[1] * sqrt(global[2]) * global_prior_scale * error_scale;
    vector[K] lambda2 = square(lambda);
    vector[K] lambda_tilde = sqrt( c2 * lambda2 ./ (c2 + square(tau) * lambda2) );
    return z_beta .* lambda_tilde * tau;
  }

  /**
   * Hierarchical shrinkage plus parameterization
   *
   * @param z_beta A vector of primitive coefficients
   * @param global A real array of positive numbers
   * @param local A vector array of positive numbers
   * @param global_prior_scale A positive real number
   * @param error_scale 1 or sigma in the Gaussian case
   * @param c2 A positive real number
   * @return A vector of coefficientes
   */
  vector hsplus_prior(vector z_beta, real[] global, vector[] local,
                      real global_prior_scale, real error_scale, real c2) {
    int K = rows(z_beta);
    vector[K] lambda = local[1] .* sqrt(local[2]);
    vector[K] eta = local[3] .* sqrt(local[4]);
    real tau = global[1] * sqrt(global[2]) * global_prior_scale * error_scale;
    vector[K] lambda_eta2 = square(lambda .* eta);
    vector[K] lambda_tilde = sqrt( c2 * lambda_eta2 ./
                                 ( c2 + square(tau) * lambda_eta2) );
    return z_beta .* lambda_tilde * tau;
  }

  /**
   * Cornish-Fisher expansion for standard normal to Student t
   *
   * See result 26.7.5 of
   * http://people.math.sfu.ca/~cbm/aands/page_949.htm
   *
   * @param z A scalar distributed standard normal
   * @param df A scalar degrees of freedom
   * @return An (approximate) Student t variate with df degrees of freedom
   */
  real CFt(real z, real df) {
    real z2 = square(z);
    real z3 = z2 * z;
    real z5 = z2 * z3;
    real z7 = z2 * z5;
    real z9 = z2 * z7;
    real df2 = square(df);
    real df3 = df2 * df;
    real df4 = df2 * df2;
    return z + (z3 + z) / (4 * df) + (5 * z5 + 16 * z3 + 3 * z) / (96 * df2)
           + (3 * z7 + 19 * z5 + 17 * z3 - 15 * z) / (384 * df3)
           + (79 * z9 + 776 * z7 + 1482 * z5 - 1920 * z3 - 945 * z) / (92160 * df4);
  }

  /**
  * Return the lower bound for the baseline hazard parameters
  *
  * @param type An integer indicating the type of baseline hazard
  * @return A real
  */
  real coefs_lb(int type) {
    real lb;
    if (type == 2) // B-splines, on log haz scale
      lb = negative_infinity();
    else if (type == 3) // piecewise constant, on log haz scale
      lb = negative_infinity();
    else
      lb = 0;
    return lb;
  }

  /**
  * Return the required number of local hs parameters
  *
  * @param prior_dist An integer indicating the prior distribution
  * @return An integer
  */
  int get_nvars_for_hs(int prior_dist) {
    int hs = 0;
    if (prior_dist == 3) hs = 2;
    else if (prior_dist == 4) hs = 4;
    return hs;
  }

  /**
  * Scale the primitive population level parameters based on prior information
  *
  * @param z_beta A vector of primitive parameters
  * @param prior_dist Integer, the type of prior distribution
  * @param prior_mean,prior_scale Vectors of mean and scale parameters
  *   for the prior distributions
  * @return A vector containing the population level parameters (coefficients)
  */
  vector make_beta(vector z_beta, int prior_dist, vector prior_mean,
                   vector prior_scale, vector prior_df, real global_prior_scale,
                   real[] global, vector[] local, real[] ool, vector[] mix,
                   real[] aux, int family, real slab_scale, real[] caux) {
    vector[rows(z_beta)] beta;
    if (prior_dist == 0) beta = z_beta;
    else if (prior_dist == 1) beta = z_beta .* prior_scale + prior_mean;
    else if (prior_dist == 2) for (k in 1:rows(prior_mean)) {
      beta[k] = CFt(z_beta[k], prior_df[k]) * prior_scale[k] + prior_mean[k];
    }
    else if (prior_dist == 3) {
      real c2 = square(slab_scale) * caux[1];
      if (family == 1) // don't need is_continuous since family == 1 is gaussian in mvmer
        beta = hs_prior(z_beta, global, local, global_prior_scale, aux[1], c2);
      else
        beta = hs_prior(z_beta, global, local, global_prior_scale, 1, c2);
    }
    else if (prior_dist == 4) {
      real c2 = square(slab_scale) * caux[1];
      if (family == 1) // don't need is_continuous since family == 1 is gaussian in mvmer
        beta = hsplus_prior(z_beta, global, local, global_prior_scale, aux[1], c2);
      else
        beta = hsplus_prior(z_beta, global, local, global_prior_scale, 1, c2);
    }
    else if (prior_dist == 5) // laplace
      beta = prior_mean + prior_scale .* sqrt(2 * mix[1]) .* z_beta;
    else if (prior_dist == 6) // lasso
      beta = prior_mean + ool[1] * prior_scale .* sqrt(2 * mix[1]) .* z_beta;
    return beta;
  }

  /**
  * Log-prior for coefficients
  *
  * @param z_beta Vector of primative coefficients
  * @param prior_dist Integer, the type of prior distribution
  * @param prior_scale Real, scale for the prior distribution
  * @param prior_df Real, df for the prior distribution
  * @param global_prior_df Real, df for the prior for the global hs parameter
  * @param local Vector of hs local parameters
  * @param global Real, the global parameter
  * @param mix Vector of shrinkage parameters
  * @param one_over_lambda Real
  * @return nothing
  */
  void beta_lp(vector z_beta, int prior_dist, vector prior_scale,
               vector prior_df, real global_prior_df, vector[] local,
               real[] global, vector[] mix, real[] one_over_lambda,
               real slab_df, real[] caux) {
    if      (prior_dist == 1) target += normal_lpdf(z_beta | 0, 1);
    else if (prior_dist == 2) target += normal_lpdf(z_beta | 0, 1); // Student t
    else if (prior_dist == 3) { // hs
      target += normal_lpdf(z_beta | 0, 1);
      target += normal_lpdf(local[1] | 0, 1);
      target += inv_gamma_lpdf(local[2] | 0.5 * prior_df, 0.5 * prior_df);
      target += normal_lpdf(global[1] | 0, 1);
      target += inv_gamma_lpdf(global[2] | 0.5 * global_prior_df, 0.5 * global_prior_df);
      target += inv_gamma_lpdf(caux | 0.5 * slab_df, 0.5 * slab_df);
    }
    else if (prior_dist == 4) { // hs+
      target += normal_lpdf(z_beta | 0, 1);
      target += normal_lpdf(local[1] | 0, 1);
      target += inv_gamma_lpdf(local[2] | 0.5 * prior_df, 0.5 * prior_df);
      target += normal_lpdf(local[3] | 0, 1);
      // unorthodox useage of prior_scale as another df hyperparameter
      target += inv_gamma_lpdf(local[4] | 0.5 * prior_scale, 0.5 * prior_scale);
      target += normal_lpdf(global[1] | 0, 1);
      target += inv_gamma_lpdf(global[2] | 0.5 * global_prior_df, 0.5 * global_prior_df);
      target += inv_gamma_lpdf(caux | 0.5 * slab_df, 0.5 * slab_df);
    }
    else if (prior_dist == 5) { // laplace
      target += normal_lpdf(z_beta | 0, 1);
      target += exponential_lpdf(mix[1] | 1);
    }
    else if (prior_dist == 6) { // lasso
      target += normal_lpdf(z_beta | 0, 1);
      target += exponential_lpdf(mix[1] | 1);
      target += chi_square_lpdf(one_over_lambda[1] | prior_df[1]);
    }
    else if (prior_dist == 7) { // product_normal
      target += normal_lpdf(z_beta | 0, 1);
    }
    /* else prior_dist is 0 and nothing is added */
  }

  /**
  * Log-prior for intercept parameters
  *
  * @param gamma Real, the intercept parameter
  * @param dist Integer, the type of prior distribution
  * @param mean Real, mean of prior distribution
  * @param scale Real, scale for the prior distribution
  * @param df Real, df for the prior distribution
  * @return nothing
  */
  void gamma_lp(real gamma, int dist, real mean, real scale, real df) {
    if (dist == 1)  // normal
      target += normal_lpdf(gamma | mean, scale);
    else if (dist == 2)  // student_t
      target += student_t_lpdf(gamma | df, mean, scale);
    /* else dist is 0 and nothing is added */
  }

  /**
  * Log-prior for baseline hazard parameters
  *
  * @param aux_unscaled Vector (potentially of length 1) of unscaled
  *   auxiliary parameter(s)
  * @param dist Integer specifying the type of prior distribution
  * @param df Real specifying the df for the prior distribution
  * @return nothing
  */
  void basehaz_lp(vector aux_unscaled, int dist, vector df) {
    if (dist > 0) {
      if (dist == 1)
        target += normal_lpdf(aux_unscaled | 0, 1);
      else if (dist == 2)
        target += student_t_lpdf(aux_unscaled | df, 0, 1);
      else
        target += exponential_lpdf(aux_unscaled | 1);
    }
  }

}

data {

  // dimensions
  int<lower=0> K;          // num. cols in predictor matrix
  int<lower=0> nevents;    // num. rows w/ an event (ie. not censored)
  int<lower=0> ncensor;    // num. rows w/ censoring
  int<lower=0> ndelayed;   // num. rows w/ delayed entry
  int<lower=0> qnodes;     // num. nodes for GK quadrature
  int<lower=0> qrows;      // num. rows used for quadrature
  int<lower=0> qdelayed;   // num. rows used for quadrature for delayed entry
  int<lower=0> nvars;      // num. aux parameters for baseline hazard

  // response variables
  vector[nevents]  t_events;  // time of events
  vector[ncensor]  t_censor;  // time of censoring
  vector[ndelayed] t_delayed; // time of entry for delayed entry

  // predictor matrices
  matrix[nevents,K]  x_events;       // for rows with events
  matrix[ncensor,K]  x_censor;       // for rows with censoring
  matrix[ndelayed,K] x_delayed;      // for rows with delayed entry
  matrix[qrows,K]    x_qpts;         // for rows at quadpoints
  matrix[qdelayed,K] x_qpts_delayed; // for rows at quadpoints for delayed entry

  // design matrices for baseline hazard
  matrix[nevents,nvars]  basis_events;       // spline basis for rows with events
  matrix[ncensor,nvars]  basis_censor;       // spline basis for rows with censoring
  matrix[ndelayed,nvars] basis_delayed;      // spline basis for rows with delayed entry
  matrix[qrows,nvars]    basis_qpts;         // spline basis for rows at quadpoints
  matrix[qdelayed,nvars] basis_qpts_delayed; // spline basis for rows at quadpoints
                                             //   for delayed entry
  matrix[nevents,nvars]  ibasis_events;      // integral of spline basis for rows with 
	                                           //   events; only used for M-splines
  matrix[ncensor,nvars]  ibasis_censor;      // integral of spline basis for rows with
	                                           //   censoring; only used for M-splines
  matrix[ndelayed,nvars] ibasis_delayed;     // integral of spline basis for rows with
	                                           //   delayed entry; only used for M-splines

  // baseline hazard type:
  //   1 = weibull
  //   2 = B-splines
  //   3 = piecewise
  //   4 = M-splines
  //   5 = exponential
  //   6 = gompertz
  int<lower=1,upper=7> type;

  // GK quadrature weights, with (b-a)/2 scaling already incorporated
  vector[qrows]    qwts;
  vector[qdelayed] qwts_delayed;

  // flags
  int<lower=0,upper=1> has_quadrature;// log surv is calculated using quadrature
  int<lower=0,upper=1> has_intercept; // basehaz requires intercept
  int<lower=0,upper=1> prior_PD;      // draw only from prior predictive dist.

  // prior family:
  //   0 = none
  //   1 = normal
  //   2 = student_t
  //   3 = hs
  //   4 = hs_plus
  //   5 = laplace
  //   6 = lasso
  int<lower=0,upper=6> prior_dist;

  // prior family:
  //   0 = none
  //   1 = normal
  //   2 = student_t
  int<lower=0,upper=2> prior_dist_for_intercept;

  // prior family:
  //   0 = none
  //   1 = normal
  //   2 = student_t
  //   3 = exponential
  int<lower=0,upper=3> prior_dist_for_aux;

  // hyperparameter (log hazard ratios), set to 0 if there is no prior
  vector[K]           prior_mean;
  vector<lower=0>[K]  prior_scale;
  vector<lower=0>[K]  prior_df;
  real<lower=0>       global_prior_scale; // for hs priors only
  real<lower=0>       global_prior_df;
  real<lower=0>       slab_scale;
  real<lower=0>       slab_df;

  // hyperparameters (intercept), set to 0 if there is no prior
  real                prior_mean_for_intercept;
  real<lower=0>       prior_scale_for_intercept;
  real<lower=0>       prior_df_for_intercept;

  // hyperparameters (basehaz pars), set to 0 if there is no prior
  vector<lower=0>[nvars] prior_scale_for_aux;
  vector<lower=0>[nvars] prior_df_for_aux;
}

transformed data {
  int<lower=0> hs = get_nvars_for_hs(prior_dist);

  vector[nevents]  log_t_events  = log(t_events);  // log time of events
  vector[ncensor]  log_t_censor  = log(t_censor);  // log time of censoring
  vector[ndelayed] log_t_delayed = log(t_delayed); // log time of entry for delayed entry

  real sum_t_events = sum(t_events);         // sum of time of events
  real sum_log_t_events = sum(log_t_events); // sum of log time of events
}

parameters {

  // primitive log hazard ratios
  vector[K] z_beta;

  // intercept
  real gamma[has_intercept == 1];

  // unscaled basehaz parameters
  //   exp model:      nvars = 0, ie. no aux parameter
  //   weibull model:  nvars = 1, ie. shape parameter
  //   gompertz model: nvars = 1, ie. scale parameter
  //   M-spline model: nvars = number of basis terms, ie. spline coefs
  //   B-spline model: nvars = number of basis terms, ie. spline coefs
  vector<lower=coefs_lb(type)>[nvars] z_coefs;

  // parameters for priors
  real<lower=0> global[hs];
  vector<lower=0>[hs > 0 ? K : 0] local[hs];
  real<lower=0> caux[hs > 0];
  vector<lower=0>[K] mix[prior_dist == 5 || prior_dist == 6];
  real<lower=0> ool[prior_dist == 6];
}

transformed parameters {

  // log hazard ratios
  vector[K] beta;

  // basehaz parameters
  vector[nvars] coefs;

  // define log hazard ratios
  if (K > 0) {
    beta = make_beta(z_beta, prior_dist, prior_mean,
                     prior_scale, prior_df, global_prior_scale,
                     global, local, ool, mix, rep_array(1.0, 0), 0,
                     slab_scale, caux);
  }

  // define basehaz parameters
  if (nvars > 0) {
    coefs = z_coefs .* prior_scale_for_aux;
  }
}

model {

  //-------- models without quadrature

  if (has_quadrature == 0) {

    vector[nevents]  eta_events;  // linear predictor at event times
    vector[ncensor]  eta_censor;  // linear predictor at censoring times
    vector[ndelayed] eta_delayed; // linear predictor at entry times

    real lsur_events  = 0;  // summation of log surv at event times
    real lsur_censor  = 0;  // summation of log surv at censoring times
    real lsur_delayed = 0;  // summation of log surv at entry times

    real lhaz = 0; // summation of log hazard at event times

    // linear predictor
    if (K > 0) {
      if (nevents > 0)
        eta_events = x_events * beta;
      if (ncensor > 0)
        eta_censor = x_censor * beta;
      if (ndelayed > 0)
        eta_delayed = x_delayed * beta;
    }
    else {
      if (nevents > 0)
        eta_events = rep_vector(0.0, nevents);
      if (ncensor > 0)
        eta_censor = rep_vector(0.0, ncensor);
      if (ndelayed > 0)
        eta_delayed = rep_vector(0.0, ndelayed);
    }

    // add intercept
    if (has_intercept == 1) {
      if (nevents > 0)
        eta_events = eta_events + gamma[1];
      if (ncensor > 0)
        eta_censor = eta_censor + gamma[1];
      if (ndelayed > 0)
        eta_delayed = eta_delayed + gamma[1];
    }

    // evaluate log hazard and log survival
    if (type == 5) { // exponential model
      if (nevents > 0)
        lsur_events = - dot_product(t_events, exp(eta_events));
      if (ncensor > 0)
        lsur_censor = - dot_product(t_censor, exp(eta_censor));
      if (ndelayed > 0)
        lsur_delayed = - dot_product(t_delayed, exp(eta_delayed));
      if (nevents > 0)
        lhaz = sum(eta_events);
    }
    else if (type == 1) { // weibull model
      real shape = coefs[1];
      real log_shape = log(shape);
      if (nevents > 0)
        lsur_events = - dot_product(exp(shape * log_t_events), exp(eta_events));
      if (ncensor > 0)
        lsur_censor = - dot_product(exp(shape * log_t_censor), exp(eta_censor));
      if (ndelayed > 0)
        lsur_delayed = - dot_product(exp(shape * log_t_delayed), exp(eta_delayed));
      if (nevents > 0)
        lhaz = (nevents * log_shape) + (shape - 1) * sum_log_t_events + sum(eta_events);
    }
    else if (type == 6) { // gompertz model
      real scale = coefs[1];
      if (nevents > 0) {
        vector[nevents] temp_events = (exp(scale * t_events) - 1) / scale;
        lsur_events = - dot_product(temp_events, exp(eta_events));
      }
      if (ncensor > 0) {
        vector[ncensor] temp_censor = (exp(scale * t_censor) - 1) / scale;
        lsur_censor = - dot_product(temp_censor, exp(eta_censor));
      }
      if (ndelayed > 0) {
        vector[ndelayed] temp_delayed = (exp(scale * t_delayed) - 1) / scale;
        lsur_delayed = - dot_product(temp_delayed, exp(eta_delayed));
      }
      if (nevents > 0)
        lhaz = scale * sum_t_events + sum(eta_events);
    }
    else if (type == 4) { // M-splines, on haz scale
      if (nevents > 0)
        lsur_events = - dot_product(ibasis_events * coefs, exp(eta_events));
      if (ncensor > 0)
        lsur_censor = - dot_product(ibasis_censor * coefs, exp(eta_censor));
      if (ndelayed > 0)
        lsur_delayed = - dot_product(ibasis_delayed * coefs, exp(eta_delayed));
      if (nevents > 0)
        lhaz = sum(log(basis_events * coefs) + eta_events);
    }
    else {
      reject("Bug found: invalid baseline hazard (without quadrature).");
    }

    // increment target
    if (prior_PD == 0) { // unweighted log likelihood
      target += lhaz + lsur_events + lsur_censor - lsur_delayed;
    }
  }

  //-------- models with quadrature

  else {

    vector[nevents]  eta_events;       // linear pred at event times
    vector[qrows]    eta_qpts;         // linear pred at quadpoints
    vector[qdelayed] eta_qpts_delayed; // linear pred at quadpoints for entry times

    real lsur  = 0;        // summation of log surv at event & censoring times
    real lsur_delayed = 0; // summation of log surv at entry times

    real lhaz = 0; // summation of log hazard at event times

    // linear predictor
    if (K > 0) {
      if (nevents > 0)
        eta_events = x_events * beta;
      if (qrows > 0)
        eta_qpts = x_qpts * beta;
      if (qdelayed > 0)
        eta_qpts_delayed = x_qpts_delayed * beta;
    }
    else {
      if (nevents > 0)
        eta_events = rep_vector(0.0, nevents);
      if (qrows > 0)
        eta_qpts = rep_vector(0.0, qrows);
      if (qdelayed > 0)
        eta_qpts_delayed = rep_vector(0.0, qdelayed);
    }

    // add intercept
    if (has_intercept == 1) {
      if (nevents > 0)
        eta_events = eta_events + gamma[1];
      if (qrows > 0)
        eta_qpts = eta_qpts + gamma[1];
      if (qdelayed > 0)
        eta_qpts_delayed = eta_qpts_delayed + gamma[1];
    }

    // evaluate log hazard and log survival
    if (type == 2) { // B-splines, on log haz scale
      if (qrows > 0) {
        vector[qrows] lhaz_qpts;
        lhaz_qpts = basis_qpts * coefs + eta_qpts;
        lsur = - dot_product(qwts, exp(lhaz_qpts));
      }
      if (qdelayed > 0) {
        vector[qdelayed] lhaz_qpts_delayed;
        lhaz_qpts_delayed = basis_qpts_delayed * coefs + eta_qpts_delayed;
        lsur_delayed = - dot_product(qwts_delayed, exp(lhaz_qpts_delayed));
      }
      if (nevents > 0) {
        lhaz = sum(basis_events * coefs + eta_events);
      }
    }
    else {
      reject("Bug found: invalid baseline hazard (with quadrature).");
    }

    // increment target
    if (prior_PD == 0) { // unweighted log likelihood
      target += lhaz + lsur - lsur_delayed;
    }
  }

  //-------- log priors

  // log priors for coefficients
  if (K > 0) {
    beta_lp(z_beta, prior_dist, prior_scale, prior_df, global_prior_df,
            local, global, mix, ool, slab_df, caux);
  }

  // log prior for intercept
  if (has_intercept == 1) {
    gamma_lp(gamma[1], prior_dist_for_intercept, prior_mean_for_intercept,
             prior_scale_for_intercept, prior_df_for_intercept);
  }

  // log priors for baseline hazard parameters
  if (nvars > 0) {
    basehaz_lp(z_coefs, prior_dist_for_aux, prior_df_for_aux);
  }

}
