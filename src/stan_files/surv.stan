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

  /**
  * Log-prior for tde spline coefficients and their smoothing parameters
  *
  * @param z_beta_tde Vector of unscaled spline coefficients
  * @param smooth_sd_raw Vector (potentially of length 1) of smoothing sds
  * @param dist Integer specifying the type of prior distribution for the
  *   smoothing sds
  * @param df Vector of reals specifying the df for the prior distribution
  *   for the smoothing sds
  * @return nothing
  */
  void smooth_lp(vector z_beta_tde, int[,] smooth_idx,
                 vector smooth_sd_raw, int dist, vector df) {
    int J = size(smooth_idx);
    for (j in 1:J) {
      int K = smooth_idx[j,2] - smooth_idx[j,1] + 1;
      vector[K] z_beta_tde_tmp = z_beta_tde[smooth_idx[j,1]:smooth_idx[j,2]];
      matrix[K,K] mat = crossprod(penalty_matrix(K)) + diag_matrix(rep_vector(1e-6,K));
      target += multi_normal_prec_lpdf(z_beta_tde_tmp | rep_vector(0.0,K), mat);
    }
    if (dist > 0) {
      real log_half = -0.693147180559945286;
      if (dist == 1)
        target += normal_lpdf(smooth_sd_raw | 0, 1) - log_half;
      else if (dist == 2)
        target += student_t_lpdf(smooth_sd_raw | df, 0, 1) - log_half;
      else if (dist == 3)
        target += exponential_lpdf(smooth_sd_raw | 1);
    }
  }

  /**
  * Raise each element of x to the power of y
  *
  * @param x Vector
  * @param y Real, the power to raise to
  * @return vector
  */
  vector pow_vec(vector x, real y) {
    int N = rows(x);
    vector[N] res;
    for (n in 1:N)
      res[n] = pow(x[n], y);
    return res;
  }

  /**
  * Log hazard, survival and CDF for exponential distribution
  *
  * @param eta Vector, linear predictor
  * @param t Vector, event or censoring times
  * @return nothing
  */
  vector exponential_log_haz(vector eta) {
    return eta;
  }

  vector exponential_log_cumhaz(vector eta, vector t) {
    vector[rows(eta)] res;
    res = log(t) + eta;
    return res;
  }

  vector exponential_log_surv(vector eta, vector t) {
    vector[rows(eta)] res;
    res = - t .* exp(eta);
    return res;
  }

  vector exponential_log_cdf(vector eta, vector t) {
    vector[rows(eta)] res;
    res = log(1 - exp(-t .* exp(eta)));
    return res;
  }

  vector exponential_log_cdf2(vector eta, vector t_lower, vector t_upper) {
    int N = rows(eta);
    vector[N] exp_eta = exp(eta);
    vector[N] surv_lower = exp(-t_lower .* exp_eta);
    vector[N] surv_upper = exp(-t_upper .* exp_eta);
    vector[N] res;
    res = log(surv_lower - surv_upper);
    return res;
  }

  /**
  * Log hazard, survival and CDF for Weibull distribution
  *
  * @param eta Vector, linear predictor
  * @param t Vector, event or censoring times
  * @param shape Real, Weibull shape
  * @return nothing
  */
  vector weibull_log_haz(vector eta, vector t, real shape) {
    vector[rows(eta)] res;
    res = log(shape) + (shape - 1) * log(t) + eta;
    return res;
  }

  vector weibull_log_cumhaz(vector eta, vector t, real shape) {
    vector[rows(eta)] res;
    res = shape * log(t) + eta;
    return res;
  }

  vector weibull_log_surv(vector eta, vector t, real shape) {
    vector[rows(eta)] res;
    res = - pow_vec(t, shape) .* exp(eta);
    return res;
  }

  vector weibull_log_cdf(vector eta, vector t, real shape) {
    vector[rows(eta)] res;
    res = log(1 - exp(- pow_vec(t, shape) .* exp(eta)));
    return res;
  }

  vector weibull_log_cdf2(vector eta, vector t_lower, vector t_upper, real shape) {
    int N = rows(eta);
    vector[N] exp_eta = exp(eta);
    vector[N] surv_lower = exp(- pow_vec(t_lower, shape) .* exp_eta);
    vector[N] surv_upper = exp(- pow_vec(t_upper, shape) .* exp_eta);
    vector[N] res;
    res = log(surv_lower - surv_upper);
    return res;
  }

  /**
  * Log hazard, survival and CDF for M-spline model
  *
  * @param eta Vector, linear predictor
  * @param t Vector, event or censoring times
  * @param coefs Vector, M-spline coefficients
  * @return nothing
  */
  vector mspline_log_haz(vector eta, matrix basis, vector coefs) {
    vector[rows(eta)] res;
    res = log(basis * coefs) + eta;
    return res;
  }

  vector mspline_log_cumhaz(vector eta, matrix ibasis, vector coefs) {
    vector[rows(eta)] res;
    res = log(ibasis * coefs) + eta;
    return res;
  }

  vector mspline_log_surv(vector eta, matrix ibasis, vector coefs) {
    vector[rows(eta)] res;
    res = - (ibasis * coefs) .* exp(eta);
    return res;
  }

  vector mspline_log_cdf(vector eta, matrix ibasis, vector coefs) {
    vector[rows(eta)] res;
    res = log(1 - exp(-(ibasis * coefs) .* exp(eta)));
    return res;
  }

  vector mspline_log_cdf2(vector eta, matrix ibasis_lower, matrix ibasis_upper, vector coefs) {
    int N = rows(eta);
    vector[N] exp_eta = exp(eta);
    vector[N] surv_lower = exp(-(ibasis_lower * coefs) .* exp_eta);
    vector[N] surv_upper = exp(-(ibasis_upper * coefs) .* exp_eta);
    vector[N] res;
    res = log(surv_lower - surv_upper);
    return res;
  }

}

data {

  // dimensions
  int<lower=0> K;          // num. cols in predictor matrix (time-fixed)
  int<lower=0> S;          // num. cols in predictor matrix (time-varying)
  int<lower=0> nevents;    // num. rows w/ an event (ie. not censored)
  int<lower=0> nlcens;     // num. rows w/ left censoring
  int<lower=0> nrcens;     // num. rows w/ right censoring
  int<lower=0> nicens;     // num. rows w/ interval censoring
  int<lower=0> ndelayed;   // num. rows w/ delayed entry
  int<lower=0> qnodes;     // num. nodes for GK quadrature
  int<lower=0> qrows;      // num. rows used for quadrature
  int<lower=0> qdelayed;   // num. rows used for quadrature for delayed entry
  int<lower=0> nvars;      // num. aux parameters for baseline hazard
  int<lower=1> smooth_map[S]; // indexing of smooth sds for tde spline coefs
  int<lower=0> smooth_idx[S > 0 ? max(smooth_map) : 0, 2];

  // response and time variables
  vector[nevents]  t_events;     // time of events
  vector[nlcens]   t_lcens;      // time of left censoring
  vector[nrcens]   t_rcens;      // time of right censoring
  vector[nicens]   t_icensl;     // time of lower limit for interval censoring
  vector[nicens]   t_icensu;     // time of upper limit for interval censoring
  vector[ndelayed] t_delayed;    // time of entry for delayed entry
  vector[qrows]    qpts;         // time at quadpoints
  vector[qdelayed] qpts_delayed; // time at quadpoints for delayed entry

  // predictor matrices (time-fixed)
  matrix[nevents, K] x_events;       // for rows with events
  matrix[nlcens,  K] x_lcens;        // for rows with left censoring
  matrix[nrcens,  K] x_rcens;        // for rows with right censoring
  matrix[nicens,  K] x_icens;        // for rows with interval censoring
  matrix[ndelayed,K] x_delayed;      // for rows with delayed entry
  matrix[qrows,K]    x_qpts;         // for rows at quadpoints
  matrix[qdelayed,K] x_qpts_delayed; // for rows at quadpoints for delayed entry

  // predictor matrices (time-varying)
  matrix[nevents,S]  s_events;       // for rows with events
  matrix[qrows,S]    s_qpts;         // for rows at quadpoints
  matrix[qdelayed,S] s_qpts_delayed; // for rows at quadpoints for delayed entry

  // design matrices for baseline hazard
  matrix[nevents, nvars] basis_events;       // spline basis for rows with events
  matrix[nlcens,  nvars] basis_lcens;        // spline basis for rows with left censoring
  matrix[nrcens,  nvars] basis_rcens;        // spline basis for rows with right censoring
  matrix[nicens,  nvars] basis_icensl;       // spline basis for at lower limit of interval censoring
  matrix[nicens,  nvars] basis_icensu;       // spline basis for at upper limit of interval censoring
  matrix[ndelayed,nvars] basis_delayed;      // spline basis for rows with delayed entry
  matrix[qrows,   nvars] basis_qpts;         // spline basis for rows at quadpoints
  matrix[qdelayed,nvars] basis_qpts_delayed; // spline basis for rows at quadpoints
                                             //   for delayed entry
  matrix[nevents, nvars] ibasis_events;      // integral of spline basis for rows with
                                             //   events; only used for M-splines
  matrix[nlcens,  nvars] ibasis_lcens;       // integral of spline basis for rows with
                                             //   left censoring; only used for M-splines
  matrix[nrcens,  nvars] ibasis_rcens;       // integral of spline basis for rows with
                                             //   left censoring; only used for M-splines
  matrix[nicens,  nvars] ibasis_icensl;      // integral of spline basis at lower limit of
                                             //   interval censoring; only used for M-splines
  matrix[nicens,  nvars] ibasis_icensu;      // integral of spline basis at upper limit of
                                             //   interval censoring; only used for M-splines
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

  // prior family:
  //   0 = none
  //   1 = normal
  //   2 = student_t
  //   3 = exponential
  int<lower=0,upper=3> prior_dist_for_smooth;

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

  // hyperparameters (tde smooths), set to 0 if there is no prior
  vector         [S > 0 ? max(smooth_map) : 0] prior_mean_for_smooth;
  vector<lower=0>[S > 0 ? max(smooth_map) : 0] prior_scale_for_smooth;
  vector<lower=0>[S > 0 ? max(smooth_map) : 0] prior_df_for_smooth;

}

transformed data {
  int<lower=0> hs = get_nvars_for_hs(prior_dist);

  vector[nevents]  log_t_events     = log(t_events);     // log time at event
  vector[qrows]    log_qpts         = log(qpts);         // log time at qpts
  vector[qdelayed] log_qpts_delayed = log(qpts_delayed); // log time at qpts
                                                         // for delayed entry

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

  // unscaled tde spline coefficients
  vector[S] z_beta_tde;

  // hyperparameter, the prior sd for the tde spline coefs
  vector<lower=0>[S > 0 ? max(smooth_map) : 0] smooth_sd_raw;

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

  // tde spline coefficients and their hyperparameters
  vector[S] beta_tde;
  vector[S > 0 ? max(smooth_map) : 0] smooth_sd; // sd for tde splines

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

  // define tde spline coefficients
  if (S > 0) {
    smooth_sd = smooth_sd_raw .* prior_scale_for_smooth + prior_mean_for_smooth;
    beta_tde = z_beta_tde .* smooth_sd[smooth_map];
  }
}

model {

  if (prior_PD == 0) {

    //-------- models without quadrature

    if (has_quadrature == 0) {

      vector[nevents]  eta_events;  // linear predictor for events
      vector[nlcens]   eta_lcens;   // linear predictor for left censored
      vector[nrcens]   eta_rcens;   // linear predictor for right censored
      vector[nicens]   eta_icens;   // linear predictor for interval censored
      vector[ndelayed] eta_delayed; // linear predictor for delayed entry

      // linear predictor
      if (K > 0) {
        if (nevents  > 0) eta_events  = x_events  * beta;
        if (nlcens   > 0) eta_lcens   = x_lcens   * beta;
        if (nrcens   > 0) eta_rcens   = x_rcens   * beta;
        if (nicens   > 0) eta_icens   = x_icens   * beta;
        if (ndelayed > 0) eta_delayed = x_delayed * beta;
      }
      else {
        if (nevents  > 0) eta_events  = rep_vector(0.0, nevents);
        if (nlcens   > 0) eta_lcens   = rep_vector(0.0, nlcens);
        if (nrcens   > 0) eta_rcens   = rep_vector(0.0, nrcens);
        if (nicens   > 0) eta_icens   = rep_vector(0.0, nicens);
        if (ndelayed > 0) eta_delayed = rep_vector(0.0, ndelayed);
      }

      // add intercept
      if (has_intercept == 1) {
        if (nevents  > 0) eta_events  = eta_events  + gamma[1];
        if (nlcens   > 0) eta_lcens   = eta_lcens   + gamma[1];
        if (nrcens   > 0) eta_rcens   = eta_rcens   + gamma[1];
        if (nicens   > 0) eta_icens   = eta_icens   + gamma[1];
        if (ndelayed > 0) eta_delayed = eta_delayed + gamma[1];
      }

      // evaluate log hazard and log survival
      if (type == 5) { // exponential model
        if (nevents  > 0) target +=  exponential_log_haz (eta_events);
        if (nevents  > 0) target +=  exponential_log_surv(eta_events,  t_events);
        if (nlcens   > 0) target +=  exponential_log_cdf (eta_lcens,   t_lcens);
        if (nrcens   > 0) target +=  exponential_log_surv(eta_rcens,   t_rcens);
        if (nicens   > 0) target +=  exponential_log_cdf2(eta_icens ,  t_icensl, t_icensu);
        if (ndelayed > 0) target += -exponential_log_surv(eta_delayed, t_delayed);
      }
      else if (type == 1) { // weibull model
        real shape = coefs[1];
        if (nevents  > 0) target +=  weibull_log_haz (eta_events,  t_events, shape);
        if (nevents  > 0) target +=  weibull_log_surv(eta_events,  t_events, shape);
        if (nlcens   > 0) target +=  weibull_log_cdf (eta_lcens,   t_lcens,  shape);
        if (nrcens   > 0) target +=  weibull_log_surv(eta_rcens,   t_rcens,  shape);
        if (nicens   > 0) target +=  weibull_log_cdf2(eta_icens,   t_icensl, t_icensu, shape);
        if (ndelayed > 0) target += -weibull_log_surv(eta_delayed, t_delayed, shape);
      }
      else if (type == 6) { // gompertz model
        real scale = coefs[1];
        reject("Gompertz temporarily removed.");
      }
      else if (type == 4) { // M-splines, on haz scale
        if (nevents  > 0) target +=  mspline_log_haz (eta_events,  basis_events, coefs);
        if (nevents  > 0) target +=  mspline_log_surv(eta_events,  ibasis_events, coefs);
        if (nlcens   > 0) target +=  mspline_log_cdf (eta_lcens,   ibasis_lcens,  coefs);
        if (nrcens   > 0) target +=  mspline_log_surv(eta_rcens,   ibasis_rcens,  coefs);
        if (nicens   > 0) target +=  mspline_log_cdf2(eta_icens,   ibasis_icensl, ibasis_icensu, coefs);
        if (ndelayed > 0) target += -mspline_log_surv(eta_delayed, ibasis_delayed, coefs);
      }
      else {
        reject("Bug found: invalid baseline hazard (without quadrature).");
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

      // linear predictor (time-fixed part)
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

      // add on time-varying part of linear predictor
      if (S > 0) {
        if (nevents > 0)
          eta_events = eta_events + s_events * beta_tde;
        if (qrows > 0)
          eta_qpts = eta_qpts + s_qpts * beta_tde;
        if (qdelayed > 0)
          eta_qpts_delayed = eta_qpts_delayed + s_qpts_delayed * beta_tde;
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
      if (type == 5) { // exponential model
        if (nevents > 0) {
          lhaz = sum(eta_events);
        }
        if (qrows > 0) {
          lsur = - dot_product(qwts, exp(eta_qpts));
        }
        if (qdelayed > 0) {
          lsur_delayed = - dot_product(qwts_delayed, exp(eta_qpts_delayed));
        }
      }
      else if (type == 1) { // weibull model
        real shape = coefs[1];
        real log_shape = log(shape);
        if (nevents > 0) {
          lhaz = (nevents * log_shape) + (shape - 1) * sum_log_t_events + sum(eta_events);
        }
        if (qrows > 0) {
          vector[qrows] lhaz_qpts;
          lhaz_qpts = log_shape + (shape - 1) * log_qpts + eta_qpts;
          lsur = - dot_product(qwts, exp(lhaz_qpts));
        }
        if (qdelayed > 0) {
          vector[qdelayed] lhaz_qpts_delayed;
          lhaz_qpts_delayed = log_shape + (shape - 1) * log_qpts_delayed + eta_qpts_delayed;
          lsur_delayed = - dot_product(qwts_delayed, exp(lhaz_qpts_delayed));
        }
      }
      else if (type == 6) { // gompertz model
        real scale = coefs[1];
        if (nevents > 0) {
          lhaz = scale * sum_t_events + sum(eta_events);
        }
        if (qrows > 0) {
          vector[qrows] lhaz_qpts;
          lhaz_qpts = scale * qpts + eta_qpts;
          lsur = - dot_product(qwts, exp(lhaz_qpts));
        }
        if (qdelayed > 0) {
          vector[qdelayed] lhaz_qpts_delayed;
          lhaz_qpts_delayed = scale * qpts_delayed + eta_qpts_delayed;
          lsur_delayed = - dot_product(qwts_delayed, exp(lhaz_qpts_delayed));
        }
      }
      else if (type == 4) { // M-splines, on haz scale
        if (nevents > 0)
          lhaz = sum(log(basis_events * coefs) + eta_events);
        if (qrows > 0) {
          vector[qrows] lhaz_qpts;
          lhaz_qpts = log(basis_qpts * coefs) + eta_qpts;
          lsur = - dot_product(qwts, exp(lhaz_qpts));
        }
        if (qdelayed > 0) {
          vector[qdelayed] lhaz_qpts_delayed;
          lhaz_qpts_delayed = log(basis_qpts_delayed * coefs) + eta_qpts_delayed;
          lsur_delayed = - dot_product(qwts_delayed, exp(lhaz_qpts_delayed));
        }
      }
      else if (type == 2) { // B-splines, on log haz scale
        if (nevents > 0) {
          lhaz = sum(basis_events * coefs + eta_events);
        }
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
      }
      else {
        reject("Bug found: invalid baseline hazard (with quadrature).");
      }

      // increment target
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

  // log priors for tde spline coefficients and their smoothing parameters
  if (S > 0) {
    smooth_lp(z_beta_tde, smooth_idx, smooth_sd_raw,
              prior_dist_for_smooth, prior_df_for_smooth);
  }

}
