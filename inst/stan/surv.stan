#include /include/Columbia_copyright.stan
#include /include/Brilleman_copyright.stan
#include /include/license.stan

functions {

  #include /functions/common_functions.stan
  #include /functions/hazard_functions.stan

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
                   array[] real global, array[] vector local, array[] real ool, array[] vector mix,
                   array[] real aux, int family, real slab_scale, array[] real caux) {
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
  * @return Real, the log probability.
  */
  real beta_custom_lpdf(vector z_beta, int prior_dist, vector prior_scale,
               vector prior_df, real global_prior_df, array[] vector local,
               array[] real global, array[] vector mix, array[] real one_over_lambda,
               real slab_df, array[] real caux) {
    real lp = 0;
    if      (prior_dist == 1) lp += normal_lpdf(z_beta | 0, 1);
    else if (prior_dist == 2) lp += normal_lpdf(z_beta | 0, 1); // Student t
    else if (prior_dist == 3) { // hs
      lp += normal_lpdf(z_beta | 0, 1);
      lp += normal_lpdf(local[1] | 0, 1);
      lp += inv_gamma_lpdf(local[2] | 0.5 * prior_df, 0.5 * prior_df);
      lp += normal_lpdf(global[1] | 0, 1);
      lp += inv_gamma_lpdf(global[2] | 0.5 * global_prior_df, 0.5 * global_prior_df);
      lp += inv_gamma_lpdf(caux | 0.5 * slab_df, 0.5 * slab_df);
    }
    else if (prior_dist == 4) { // hs+
      lp += normal_lpdf(z_beta | 0, 1);
      lp += normal_lpdf(local[1] | 0, 1);
      lp += inv_gamma_lpdf(local[2] | 0.5 * prior_df, 0.5 * prior_df);
      lp += normal_lpdf(local[3] | 0, 1);
      // unorthodox useage of prior_scale as another df hyperparameter
      lp += inv_gamma_lpdf(local[4] | 0.5 * prior_scale, 0.5 * prior_scale);
      lp += normal_lpdf(global[1] | 0, 1);
      lp += inv_gamma_lpdf(global[2] | 0.5 * global_prior_df, 0.5 * global_prior_df);
      lp += inv_gamma_lpdf(caux | 0.5 * slab_df, 0.5 * slab_df);
    }
    else if (prior_dist == 5) { // laplace
      lp += normal_lpdf(z_beta | 0, 1);
      lp += exponential_lpdf(mix[1] | 1);
    }
    else if (prior_dist == 6) { // lasso
      lp += normal_lpdf(z_beta | 0, 1);
      lp += exponential_lpdf(mix[1] | 1);
      lp += chi_square_lpdf(one_over_lambda[1] | prior_df[1]);
    }
    else if (prior_dist == 7) { // product_normal
      lp += normal_lpdf(z_beta | 0, 1);
    }
    /* else prior_dist is 0 and nothing is added */
    return lp;
  }

  /**
  * Log-prior for intercept parameters
  *
  * @param gamma Real, the intercept parameter
  * @param dist Integer, the type of prior distribution
  * @param mean Real, mean of prior distribution
  * @param scale Real, scale for the prior distribution
  * @param df Real, df for the prior distribution
  * @return Real, the log probability
  */
  real gamma_custom_lpdf(real gamma, int dist, real mean, real scale, real df) {
    real lp = 0;
    if (dist == 1)  // normal
      lp += normal_lpdf(gamma | mean, scale);
    else if (dist == 2)  // student_t
      lp += student_t_lpdf(gamma | df, mean, scale);
    /* else dist is 0 and nothing is added */
    return lp;
  }

  /**
  * Log-prior for baseline hazard parameters
  *
  * @param aux_unscaled Vector (potentially of length 1) of unscaled
  *   auxiliary parameter(s)
  * @param dist Integer specifying the type of prior distribution
  * @param df Real specifying the df for the prior distribution, or in the case
  *   of the dirichlet distribution it is the concentration parameter(s)
  * @return Real, the log probability
  */
  real basehaz_lpdf(vector aux_unscaled, int dist, vector df) {
    real lp = 0;
    if (dist > 0) {
      if (dist == 1)
        lp += normal_lpdf(aux_unscaled | 0, 1);
      else if (dist == 2)
        lp += student_t_lpdf(aux_unscaled | df, 0, 1);
      else if (dist == 3)
        lp += exponential_lpdf(aux_unscaled | 1);
      else
        lp += dirichlet_lpdf(aux_unscaled | df); // df is concentration here
    }
    return lp;
  }

  /**
  * Log-prior for tve spline coefficients and their smoothing parameters
  *
  * @param z_beta_tve Vector of unscaled spline coefficients
  * @param smooth_sd_raw Vector (potentially of length 1) of smoothing sds
  * @param dist Integer specifying the type of prior distribution for the
  *   smoothing sds
  * @param df Vector of reals specifying the df for the prior distribution
  *   for the smoothing sds
  * @return Real, the log probability
  */
  real smooth_lpdf(vector z_beta_tve, vector smooth_sd_raw, int dist, vector df) {
    real lp = 0;
    lp += normal_lpdf(z_beta_tve | 0, 1);
    if (dist > 0) {
      real log_half = -0.693147180559945286;
      if (dist == 1)
        lp += normal_lpdf(smooth_sd_raw | 0, 1) - log_half;
      else if (dist == 2)
        lp += student_t_lpdf(smooth_sd_raw | df, 0, 1) - log_half;
      else if (dist == 3)
        lp += exponential_lpdf(smooth_sd_raw | 1);
    }
    return lp;
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
  * Log survival and log CDF for exponential distribution
  *
  * @param eta Vector, linear predictor
  * @param t Vector, event or censoring times
  * @return A vector
  */
  vector exponential_log_surv(vector eta, vector t) {
    vector[rows(eta)] res;
    res = - t .* exp(eta);
    return res;
  }

  vector exponential_log_cdf1(vector eta, vector t) {
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
  * Log survival and log CDF for exponential distribution; AFT parameterisation
  *
  * @param caf Vector, cumulative acceleration factor
  * @return A vector
  */
  vector exponentialAFT_log_surv(vector caf) {
    vector[rows(caf)] res;
    res = - caf;
    return res;
  }

  vector exponentialAFT_log_cdf1(vector caf) {
    vector[rows(caf)] res;
    res = log(1 - exp(-caf));
    return res;
  }

  vector exponentialAFT_log_cdf2(vector caf_lower, vector caf_upper) {
    int N = rows(caf_lower);
    vector[N] surv_lower = exp(-caf_lower);
    vector[N] surv_upper = exp(-caf_upper);
    vector[N] res;
    res = log(surv_lower - surv_upper);
    return res;
  }

  /**
  * Log survival and log CDF for Weibull distribution
  *
  * @param eta Vector, linear predictor
  * @param t Vector, event or censoring times
  * @param shape Real, Weibull shape
  * @return A vector
  */
  vector weibull_log_surv(vector eta, vector t, real shape) {
    vector[rows(eta)] res;
    res = - pow_vec(t, shape) .* exp(eta);
    return res;
  }

  vector weibull_log_cdf1(vector eta, vector t, real shape) {
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
  * Log survival and log CDF for Weibull distribution; AFT parameterisation
  *
  * @param caf Vector, cumulative acceleration factor
  * @param shape Real, Weibull shape
  * @return A vector
  */
  vector weibullAFT_log_surv(vector caf, real shape) {
    vector[rows(caf)] res;
    res = - pow_vec(caf, shape);
    return res;
  }

  vector weibullAFT_log_cdf1(vector caf, real shape) {
    vector[rows(caf)] res;
    res = log(1 - exp(- pow_vec(caf, shape)));
    return res;
  }

  vector weibullAFT_log_cdf2(vector caf_lower, vector caf_upper, real shape) {
    int N = rows(caf_lower);
    vector[N] surv_lower = exp(- pow_vec(caf_lower, shape));
    vector[N] surv_upper = exp(- pow_vec(caf_upper, shape));
    vector[N] res;
    res = log(surv_lower - surv_upper);
    return res;
  }

  /**
  * Log survival and log CDF for Gompertz distribution
  *
  * @param eta Vector, linear predictor
  * @param t Vector, event or censoring times
  * @param scale Real, Gompertz scale
  * @return A vector
  */
  vector gompertz_log_surv(vector eta, vector t, real scale) {
    vector[rows(eta)] res;
    res = inv(scale) * -(exp(scale * t) - 1) .* exp(eta);
    return res;
  }

  vector gompertz_log_cdf1(vector eta, vector t, real scale) {
    vector[rows(eta)] res;
    res = log(1 - exp(inv(scale) * -(exp(scale * t) - 1) .* exp(eta)));
    return res;
  }

  vector gompertz_log_cdf2(vector eta, vector t_lower, vector t_upper, real scale) {
    int N = rows(eta);
    real inv_scale = inv(scale);
    vector[N] exp_eta = exp(eta);
    vector[N] surv_lower = exp(inv_scale * -(exp(scale * t_lower) - 1) .* exp_eta);
    vector[N] surv_upper = exp(inv_scale * -(exp(scale * t_upper) - 1) .* exp_eta);
    vector[N] res;
    res = log(surv_lower - surv_upper);
    return res;
  }

  /**
  * Log survival and log CDF for M-spline model
  *
  * @param eta Vector, linear predictor
  * @param t Vector, event or censoring times
  * @param coefs Vector, M-spline coefficients
  * @return A vector
  */
  vector mspline_log_surv(vector eta, matrix ibasis, vector coefs) {
    vector[rows(eta)] res;
    res = - (ibasis * coefs) .* exp(eta);
    return res;
  }

  vector mspline_log_cdf1(vector eta, matrix ibasis, vector coefs) {
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
  int<lower=0> nevent;     // num. rows w/ an event (ie. not censored)
  int<lower=0> nlcens;     // num. rows w/ left censoring
  int<lower=0> nrcens;     // num. rows w/ right censoring
  int<lower=0> nicens;     // num. rows w/ interval censoring
  int<lower=0> ndelay;     // num. rows w/ delayed entry
  int<lower=0> Nevent;     // num. rows w/ an event;      used only w/ quadrature
  int<lower=0> Nlcens;     // num. rows w/ left cens;     used only w/ quadrature
  int<lower=0> Nrcens;     // num. rows w/ right cens;    used only w/ quadrature
  int<lower=0> Nicens;     // num. rows w/ interval cens; used only w/ quadrature
  int<lower=0> Ndelay;     // num. rows w/ delayed entry; used only w/ quadrature
  int<lower=0> qnodes;     // num. nodes for GK quadrature
  int<lower=0> qevent;     // num. quadrature points for rows w/ an event
  int<lower=0> qlcens;     // num. quadrature points for rows w/ left censoring
  int<lower=0> qrcens;     // num. quadrature points for rows w/ right censoring
  int<lower=0> qicens;     // num. quadrature points for rows w/ interval censoring
  int<lower=0> qdelay;     // num. quadrature points for rows w/ delayed entry
  int<lower=0> nvars;      // num. aux parameters for baseline hazard
  array[S] int<lower=1> smooth_map; // indexing of smooth sds for tve spline coefs
  array[S > 0 ? max(smooth_map) : 0, 2] int<lower=0> smooth_idx;

  // dimensions for random efffects structure, see table 3 of
  // https://cran.r-project.org/web/packages/lme4/vignettes/lmer.pdf
  int<lower=0> t;          // num. terms (maybe 0) with a | in the glmer formula
  array[t] int<lower=1> p;       // num. variables on the LHS of each |
  array[t] int<lower=1> l;       // num. levels for the factor(s) on the RHS of each |
  int<lower=0> q;          // conceptually equals \sum_{i=1}^t p_i \times l_i

  // log crude event rate / time (for centering linear predictor)
  real log_crude_event_rate;

  // response and time variables
  vector[nevent] t_event;  // time of events
  vector[nlcens] t_lcens;  // time of left censoring
  vector[nrcens] t_rcens;  // time of right censoring
  vector[nicens] t_icenl;  // time of lower limit for interval censoring
  vector[nicens] t_icenu;  // time of upper limit for interval censoring
  vector[ndelay] t_delay;  // time of entry for delayed entry

  vector[Nevent] epts_event;  // time of events
  vector[qevent] qpts_event;  // qpts for time of events
  vector[qlcens] qpts_lcens;  // qpts for time of left censoring
  vector[qrcens] qpts_rcens;  // qpts for time of right censoring
  vector[qicens] qpts_icenl;  // qpts for time of lower limit for interval censoring
  vector[qicens] qpts_icenu;  // qpts for time of upper limit for interval censoring
  vector[qdelay] qpts_delay;  // qpts for time of entry for delayed entry

  // predictor matrices (time-fixed), without quadrature
  vector[K] x_bar;          // predictor means
  matrix[nevent,K] x_event; // for rows with events
  matrix[nlcens,K] x_lcens; // for rows with left censoring
  matrix[nrcens,K] x_rcens; // for rows with right censoring
  matrix[nicens,K] x_icens; // for rows with interval censoring
  matrix[ndelay,K] x_delay; // for rows with delayed entry

  // predictor matrices (time-fixed), with quadrature
  matrix[Nevent,K] x_epts_event; // for rows with events
  matrix[qevent,K] x_qpts_event; // for rows with events
  matrix[qlcens,K] x_qpts_lcens; // for rows with left censoring
  matrix[qrcens,K] x_qpts_rcens; // for rows with right censoring
  matrix[qicens,K] x_qpts_icens; // for rows with interval censoring
  matrix[qdelay,K] x_qpts_delay; // for rows with delayed entry

  // predictor matrices (time-varying)
  matrix[Nevent,S] s_epts_event; // for rows with events
  matrix[qevent,S] s_qpts_event; // for rows with events
  matrix[qlcens,S] s_qpts_lcens; // for rows with left censoring
  matrix[qrcens,S] s_qpts_rcens; // for rows with right censoring
  matrix[qicens,S] s_qpts_icenl; // for rows with interval censoring
  matrix[qicens,S] s_qpts_icenu; // for rows with interval censoring
  matrix[qdelay,S] s_qpts_delay; // for rows with delayed entry

  // random effects structure, without quadrature
  //   nnz: number of non-zero elements in the Z matrix
  //   w: non-zero elements in the implicit Z matrix
  //   v: column indices for w
  //   u: where the non-zeros start in each row
  int<lower=0> nnz_event;
  int<lower=0> nnz_lcens;
  int<lower=0> nnz_rcens;
  int<lower=0> nnz_icens;
  int<lower=0> nnz_delay;

  vector[nnz_event] w_event;
  vector[nnz_lcens] w_lcens;
  vector[nnz_rcens] w_rcens;
  vector[nnz_icens] w_icens;
  vector[nnz_delay] w_delay;

  array[nnz_event] int<lower=0,upper=q-1> v_event;
  array[nnz_lcens] int<lower=0,upper=q-1> v_lcens;
  array[nnz_rcens] int<lower=0,upper=q-1> v_rcens;
  array[nnz_icens] int<lower=0,upper=q-1> v_icens;
  array[nnz_delay] int<lower=0,upper=q-1> v_delay;

  array[(t > 0 && nevent > 0) ? nevent + 1 : 0] int<lower=0,upper=rows(w_event)+1> u_event;
  array[(t > 0 && nlcens > 0) ? nlcens + 1 : 0] int<lower=0,upper=rows(w_lcens)+1> u_lcens;
  array[(t > 0 && nrcens > 0) ? nrcens + 1 : 0] int<lower=0,upper=rows(w_rcens)+1> u_rcens;
  array[(t > 0 && nicens > 0) ? nicens + 1 : 0] int<lower=0,upper=rows(w_icens)+1> u_icens;
  array[(t > 0 && ndelay > 0) ? ndelay + 1 : 0] int<lower=0,upper=rows(w_delay)+1> u_delay;

  // random effects structure, with quadrature
  //   nnz: number of non-zero elements in the Z matrix
  //   w: non-zero elements in the implicit Z matrix
  //   v: column indices for w
  //   u: where the non-zeros start in each row
  int<lower=0> nnz_epts_event;
  int<lower=0> nnz_qpts_event;
  int<lower=0> nnz_qpts_lcens;
  int<lower=0> nnz_qpts_rcens;
  int<lower=0> nnz_qpts_icens;
  int<lower=0> nnz_qpts_delay;

  vector[nnz_epts_event] w_epts_event;
  vector[nnz_qpts_event] w_qpts_event;
  vector[nnz_qpts_lcens] w_qpts_lcens;
  vector[nnz_qpts_rcens] w_qpts_rcens;
  vector[nnz_qpts_icens] w_qpts_icens;
  vector[nnz_qpts_delay] w_qpts_delay;

  array[nnz_epts_event] int<lower=0,upper=q-1> v_epts_event;
  array[nnz_qpts_event] int<lower=0,upper=q-1> v_qpts_event;
  array[nnz_qpts_lcens] int<lower=0,upper=q-1> v_qpts_lcens;
  array[nnz_qpts_rcens] int<lower=0,upper=q-1> v_qpts_rcens;
  array[nnz_qpts_icens] int<lower=0,upper=q-1> v_qpts_icens;
  array[nnz_qpts_delay] int<lower=0,upper=q-1> v_qpts_delay;

  array[(t > 0 && Nevent > 0) ? Nevent + 1 : 0] int<lower=0,upper=rows(w_epts_event)+1> u_epts_event;
  array[(t > 0 && qevent > 0) ? qevent + 1 : 0] int<lower=0,upper=rows(w_qpts_event)+1> u_qpts_event;
  array[(t > 0 && qlcens > 0) ? qlcens + 1 : 0] int<lower=0,upper=rows(w_qpts_lcens)+1> u_qpts_lcens;
  array[(t > 0 && qrcens > 0) ? qrcens + 1 : 0] int<lower=0,upper=rows(w_qpts_rcens)+1> u_qpts_rcens;
  array[(t > 0 && qicens > 0) ? qicens + 1 : 0] int<lower=0,upper=rows(w_qpts_icens)+1> u_qpts_icens;
  array[(t > 0 && qdelay > 0) ? qdelay + 1 : 0] int<lower=0,upper=rows(w_qpts_delay)+1> u_qpts_delay;

  // basis matrices for M-splines / I-splines, without quadrature
  matrix[nevent,nvars] basis_event;  // at event time
  matrix[nevent,nvars] ibasis_event; // at event time
  matrix[nlcens,nvars] ibasis_lcens; // at left  censoring time
  matrix[nrcens,nvars] ibasis_rcens; // at right censoring time
  matrix[nicens,nvars] ibasis_icenl; // at lower limit of interval censoring
  matrix[nicens,nvars] ibasis_icenu; // at upper limit of interval censoring
  matrix[ndelay,nvars] ibasis_delay; // at delayed entry time

  // basis matrices for M-splines, with quadrature
  matrix[Nevent,nvars] basis_epts_event; // at event time
  matrix[qevent,nvars] basis_qpts_event; // at qpts for event time
  matrix[qlcens,nvars] basis_qpts_lcens; // at qpts for left  censoring time
  matrix[qrcens,nvars] basis_qpts_rcens; // at qpts for right censoring time
  matrix[qicens,nvars] basis_qpts_icenl; // at qpts for lower limit of icens time
  matrix[qicens,nvars] basis_qpts_icenu; // at qpts for upper limit of icens time
  matrix[qdelay,nvars] basis_qpts_delay; // at qpts for delayed entry time

  // baseline hazard type:
  //   1 = weibull
  //   2 = B-splines
  //   3 = piecewise
  //   4 = M-splines
  //   5 = exponential
  //   6 = gompertz
  //   7 = exponential AFT
  //   8 = weibull AFT
  int<lower=1,upper=8> type;

  // GK quadrature weights, with (b-a)/2 scaling already incorporated
  vector[qevent] qwts_event;
  vector[qlcens] qwts_lcens;
  vector[qrcens] qwts_rcens;
  vector[qicens] qwts_icenl;
  vector[qicens] qwts_icenu;
  vector[qdelay] qwts_delay;

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
  //   4 = dirichlet
  int<lower=0,upper=4> prior_dist_for_aux;

  // prior family:
  //   0 = none
  //   1 = normal
  //   2 = student_t
  //   3 = exponential
  int<lower=0,upper=3> prior_dist_for_smooth;

  // hyperparameters (log hazard ratios), set to 0 if there is no prior
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
  vector<lower=0>[nvars] prior_conc_for_aux; // dirichlet concentration pars

  // hyperparameters (tve smooths), set to 0 if there is no prior
  vector         [S > 0 ? max(smooth_map) : 0] prior_mean_for_smooth;
  vector<lower=0>[S > 0 ? max(smooth_map) : 0] prior_scale_for_smooth;
  vector<lower=0>[S > 0 ? max(smooth_map) : 0] prior_df_for_smooth;

  // hyperparameters (random effects structure), set to 0 if there is no prior
  vector<lower=0>[t]   b_prior_shape;
  vector<lower=0>[t]   b_prior_scale;
  int<lower=0>         len_theta_L; // length of the theta_L vector
  int<lower=0>         len_concentration;
  int<lower=0>         len_regularization;
  array[len_concentration] real<lower=0>        concentration;
  array[len_regularization] real<lower=0>        regularization;
  int<lower=0,upper=1> special_case; // is the only term (1|group)

}

transformed data {

  int<lower=0> hs = get_nvars_for_hs(prior_dist);

  int sc = special_case;

  array[sc ? t : 0, nevent] int<lower=1> V_event = make_V(nevent, sc ? t : 0, v_event);
  array[sc ? t : 0, nlcens] int<lower=1> V_lcens = make_V(nlcens, sc ? t : 0, v_lcens);
  array[sc ? t : 0, nrcens] int<lower=1> V_rcens = make_V(nrcens, sc ? t : 0, v_rcens);
  array[sc ? t : 0, nicens] int<lower=1> V_icens = make_V(nicens, sc ? t : 0, v_icens);
  array[sc ? t : 0, ndelay] int<lower=1> V_delay = make_V(ndelay, sc ? t : 0, v_delay);

  array[sc ? t : 0, Nevent] int<lower=1> V_epts_event = make_V(Nevent, sc ? t : 0, v_epts_event);
  array[sc ? t : 0, qevent] int<lower=1> V_qpts_event = make_V(qevent, sc ? t : 0, v_qpts_event);
  array[sc ? t : 0, qlcens] int<lower=1> V_qpts_lcens = make_V(qlcens, sc ? t : 0, v_qpts_lcens);
  array[sc ? t : 0, qrcens] int<lower=1> V_qpts_rcens = make_V(qrcens, sc ? t : 0, v_qpts_rcens);
  array[sc ? t : 0, qicens] int<lower=1> V_qpts_icens = make_V(qicens, sc ? t : 0, v_qpts_icens);
  array[sc ? t : 0, qdelay] int<lower=1> V_qpts_delay = make_V(qdelay, sc ? t : 0, v_qpts_delay);

  int<lower=1>  pos = 1;
  int<lower=0>  len_z_T = 0;
  int<lower=0>  len_rho = sum(p) - t;
  array[len_concentration] real<lower=0> delta;

  for (i in 1:t) {
    if (p[i] > 1) {
      for (j in 1:p[i]) {
        delta[pos] = concentration[j];
        pos += 1;
      }
    }
    for (j in 3:p[i]) len_z_T += p[i] - 1;
  }



}

parameters {

  // primitive log hazard ratios
  vector[K] z_beta;

  // intercept
  array[has_intercept == 1] real gamma;

  // unscaled basehaz parameters
  //   exp model:      nvars = 0, ie. no aux parameter
  //   weibull model:  nvars = 1, ie. shape parameter
  //   gompertz model: nvars = 1, ie. scale parameter
  //   M-spline model: nvars = number of basis terms, ie. spline coefs
  //   B-spline model: nvars = number of basis terms, ie. spline coefs
  vector<lower=coefs_lb(type)>[type == 4 ? 0 : nvars] z_coefs;
  simplex[type == 4 ? nvars : 1] ms_coefs; // constrained coefs for M-splines

  // unscaled tve spline coefficients
  vector[S] z_beta_tve;

  // hyperparameter, the prior sd for the tve spline coefs
  vector<lower=0>[S > 0 ? max(smooth_map) : 0] smooth_sd_raw;

  // parameters for random effects
  vector[q] z_b;
  vector[len_z_T] z_T;
  vector<lower=0,upper=1>[len_rho] rho;
  vector<lower=0>[len_concentration] zeta;
  vector<lower=0>[t] tau;

  // parameters for priors
  array[hs] real<lower=0> global;
  array[hs] vector<lower=0>[hs > 0 ? K : 0] local;
  array[hs > 0] real<lower=0> caux;
  array[prior_dist == 5 || prior_dist == 6] vector<lower=0>[K] mix;
  array[prior_dist == 6] real<lower=0> ool;
}

transformed parameters {

  // declare log hazard ratios
  vector[K] beta;

  // declare basehaz parameters
  vector[type == 4 ? 0 : nvars] coefs;

  // declare tve spline coefficients and their hyperparameters
  vector[S] beta_tve;
  vector[S > 0 ? max(smooth_map) : 0] smooth_sd; // sd for tve splines

  // declare random effects and var-cov parameters
  vector[q] b;
  vector[len_theta_L] theta_L;

  // define log hazard ratios
  if (K > 0) {
    beta = make_beta(z_beta,
                     prior_dist, prior_mean, prior_scale, prior_df,
                     global_prior_scale, global, local, ool, mix,
                     rep_array(1.0, 0), 0, slab_scale, caux);
  }

  // define basehaz parameters
  if (type != 4 && nvars > 0) {
    coefs = z_coefs .* prior_scale_for_aux;
  }

  // define tve spline coefficients using random walk
  if (S > 0) {
    smooth_sd = smooth_sd_raw .* prior_scale_for_smooth + prior_mean_for_smooth;
    for (i in 1:max(smooth_map)) {
      int beg = smooth_idx[i,1];        // index of first spline coef
      int end = smooth_idx[i,2];        // index of last  spline coef
      beta_tve[beg] = z_beta_tve[beg];  // define first spline coef
      if (end > beg) {                  // define subsequent spline coefs
        for (j in (beg+1):end) {
          real tmp = beta_tve[j-1];
          beta_tve[j] = tmp + z_beta_tve[j] * smooth_sd[smooth_map[j]];
        }
      }
    }
  }

  // define random effects and var-cov parameters
  if (t > 0) {
    if (special_case == 1) {
      int start = 1;
      theta_L = b_prior_scale .* tau * 1.0;
      if (t == 1) {
        b = theta_L[1] * z_b;
      }
      else for (i in 1:t) {
        int end = start + l[i] - 1;
        b[start:end] = theta_L[i] * z_b[start:end];
        start = end + 1;
      }
    }
    else {
      theta_L = make_theta_L(len_theta_L, p, 1.0, tau, b_prior_scale, zeta, rho, z_T);
      b = make_b(z_b, theta_L, p, l);
    }
  }

}

model {

  if (prior_PD == 0) {

    //-------- models without quadrature

    if (has_quadrature == 0) {

      vector[nevent] eta_event; // for events
      vector[nlcens] eta_lcens; // for left censored
      vector[nrcens] eta_rcens; // for right censored
      vector[nicens] eta_icens; // for interval censored
      vector[ndelay] eta_delay; // for delayed entry

      // linear predictor
      if (K > 0) {
        if (nevent > 0) eta_event = x_event * beta;
        if (nlcens > 0) eta_lcens = x_lcens * beta;
        if (nrcens > 0) eta_rcens = x_rcens * beta;
        if (nicens > 0) eta_icens = x_icens * beta;
        if (ndelay > 0) eta_delay = x_delay * beta;
      }
      else {
        if (nevent > 0) eta_event = rep_vector(0.0, nevent);
        if (nlcens > 0) eta_lcens = rep_vector(0.0, nlcens);
        if (nrcens > 0) eta_rcens = rep_vector(0.0, nrcens);
        if (nicens > 0) eta_icens = rep_vector(0.0, nicens);
        if (ndelay > 0) eta_delay = rep_vector(0.0, ndelay);
      }

      // add on log crude event rate / time (helps to center intercept)
      if (nevent > 0) eta_event += log_crude_event_rate;
      if (nlcens > 0) eta_lcens += log_crude_event_rate;
      if (nrcens > 0) eta_rcens += log_crude_event_rate;
      if (nicens > 0) eta_icens += log_crude_event_rate;
      if (ndelay > 0) eta_delay += log_crude_event_rate;

      // add on intercept to linear predictor
      if (has_intercept == 1) {
        if (nevent > 0) eta_event += gamma[1];
        if (nlcens > 0) eta_lcens += gamma[1];
        if (nrcens > 0) eta_rcens += gamma[1];
        if (nicens > 0) eta_icens += gamma[1];
        if (ndelay > 0) eta_delay += gamma[1];
      }

      // add on random effects terms to linear predictor
      if (t > 0) {
        if (special_case) for (i in 1:t) {
          if (nevent > 0) eta_event += b[V_event[i]];
          if (nlcens > 0) eta_lcens += b[V_lcens[i]];
          if (nrcens > 0) eta_rcens += b[V_rcens[i]];
          if (nicens > 0) eta_icens += b[V_icens[i]];
          if (ndelay > 0) eta_delay += b[V_delay[i]];
        }
        else {
          if (nevent > 0) eta_event +=
            csr_matrix_times_vector(nevent, q, w_event, v_event, u_event, b);
          if (nlcens > 0) eta_lcens +=
            csr_matrix_times_vector(nlcens, q, w_lcens, v_lcens, u_lcens, b);
          if (nrcens > 0) eta_rcens +=
            csr_matrix_times_vector(nrcens, q, w_rcens, v_rcens, u_rcens, b);
          if (nicens > 0) eta_icens +=
            csr_matrix_times_vector(nicens, q, w_icens, v_icens, u_icens, b);
          if (ndelay > 0) eta_delay +=
            csr_matrix_times_vector(ndelay, q, w_delay, v_delay, u_delay, b);
        }
      }

      // aft models
      if (type == 7 || type == 8) {

        // acceleration factor at event times
        vector[nevent] af_event = exp(-eta_event);

        // cumulative acceleration factors
        vector[nevent] caf_event = t_event .* exp(-eta_event);
        vector[nlcens] caf_lcens = t_lcens .* exp(-eta_lcens);
        vector[nrcens] caf_rcens = t_rcens .* exp(-eta_rcens);
        vector[nicens] caf_icenl = t_icenl .* exp(-eta_icens);
        vector[nicens] caf_icenu = t_icenu .* exp(-eta_icens);
        vector[ndelay] caf_delay = t_delay .* exp(-eta_delay);

        // increment target with log-lik contributions
        if (type == 7) { // exponential AFT model
          if (nevent > 0) target +=  exponentialAFT_log_haz (af_event);
          if (nevent > 0) target +=  exponentialAFT_log_surv(caf_event);
          if (nlcens > 0) target +=  exponentialAFT_log_cdf1(caf_lcens);
          if (nrcens > 0) target +=  exponentialAFT_log_surv(caf_rcens);
          if (nicens > 0) target +=  exponentialAFT_log_cdf2(caf_icenl, caf_icenu);
          if (ndelay > 0) target += -exponentialAFT_log_surv(caf_delay);
        } else if (type == 8) { // weibull AFT model
          real shape = coefs[1];
          if (nevent > 0) target +=  weibullAFT_log_haz (af_event, caf_event, shape);
          if (nevent > 0) target +=  weibullAFT_log_surv(caf_event, shape);
          if (nlcens > 0) target +=  weibullAFT_log_cdf1(caf_lcens, shape);
          if (nrcens > 0) target +=  weibullAFT_log_surv(caf_rcens, shape);
          if (nicens > 0) target +=  weibullAFT_log_cdf2(caf_icenl, caf_icenu, shape);
          if (ndelay > 0) target += -weibullAFT_log_surv(caf_delay, shape);
        }

      }

      // hazard models
      else {

        // evaluate log hazard and log survival
        if (type == 5) { // exponential model
          if (nevent > 0) target +=  exponential_log_haz (eta_event);
          if (nevent > 0) target +=  exponential_log_surv(eta_event, t_event);
          if (nlcens > 0) target +=  exponential_log_cdf1(eta_lcens, t_lcens);
          if (nrcens > 0) target +=  exponential_log_surv(eta_rcens, t_rcens);
          if (nicens > 0) target +=  exponential_log_cdf2(eta_icens, t_icenl, t_icenu);
          if (ndelay > 0) target += -exponential_log_surv(eta_delay, t_delay);
        }
        else if (type == 1) { // weibull model
          real shape = coefs[1];
          if (nevent > 0) target +=  weibull_log_haz (eta_event, t_event, shape);
          if (nevent > 0) target +=  weibull_log_surv(eta_event, t_event, shape);
          if (nlcens > 0) target +=  weibull_log_cdf1(eta_lcens, t_lcens, shape);
          if (nrcens > 0) target +=  weibull_log_surv(eta_rcens, t_rcens, shape);
          if (nicens > 0) target +=  weibull_log_cdf2(eta_icens, t_icenl, t_icenu, shape);
          if (ndelay > 0) target += -weibull_log_surv(eta_delay, t_delay, shape);
        }
        else if (type == 6) { // gompertz model
          real scale = coefs[1];
          if (nevent > 0) target +=  gompertz_log_haz (eta_event, t_event, scale);
          if (nevent > 0) target +=  gompertz_log_surv(eta_event, t_event, scale);
          if (nlcens > 0) target +=  gompertz_log_cdf1(eta_lcens, t_lcens, scale);
          if (nrcens > 0) target +=  gompertz_log_surv(eta_rcens, t_rcens, scale);
          if (nicens > 0) target +=  gompertz_log_cdf2(eta_icens, t_icenl, t_icenu, scale);
          if (ndelay > 0) target += -gompertz_log_surv(eta_delay, t_delay, scale);
        }
        else if (type == 4) { // M-splines, on haz scale
          if (nevent > 0) target +=  mspline_log_haz (eta_event,  basis_event, ms_coefs);
          if (nevent > 0) target +=  mspline_log_surv(eta_event, ibasis_event, ms_coefs);
          if (nlcens > 0) target +=  mspline_log_cdf1(eta_lcens, ibasis_lcens, ms_coefs);
          if (nrcens > 0) target +=  mspline_log_surv(eta_rcens, ibasis_rcens, ms_coefs);
          if (nicens > 0) target +=  mspline_log_cdf2(eta_icens, ibasis_icenl, ibasis_icenu, ms_coefs);
          if (ndelay > 0) target += -mspline_log_surv(eta_delay, ibasis_delay, ms_coefs);
        }
        else {
          reject("Bug found: invalid baseline hazard (without quadrature).");
        }

      }
    }

    //-------- models with quadrature

    else {

      vector[Nevent] eta_epts_event; // for event times
      vector[qevent] eta_qpts_event; // for qpts for event time
      vector[qlcens] eta_qpts_lcens; // for qpts for left  censoring time
      vector[qrcens] eta_qpts_rcens; // for qpts for right censoring time
      vector[qicens] eta_qpts_icenl; // for qpts for lower limit of icens time
      vector[qicens] eta_qpts_icenu; // for qpts for upper limit of icens time
      vector[qdelay] eta_qpts_delay; // for qpts for delayed entry time

      // linear predictor (time-fixed part)
      if (K > 0) {
        if (Nevent > 0) eta_epts_event = x_epts_event * beta;
        if (qevent > 0) eta_qpts_event = x_qpts_event * beta;
        if (qlcens > 0) eta_qpts_lcens = x_qpts_lcens * beta;
        if (qrcens > 0) eta_qpts_rcens = x_qpts_rcens * beta;
        if (qicens > 0) eta_qpts_icenl = x_qpts_icens * beta;
        if (qicens > 0) eta_qpts_icenu = x_qpts_icens * beta;
        if (qdelay > 0) eta_qpts_delay = x_qpts_delay * beta;
      }
      else {
        if (Nevent > 0) eta_epts_event = rep_vector(0.0, Nevent);
        if (qevent > 0) eta_qpts_event = rep_vector(0.0, qevent);
        if (qlcens > 0) eta_qpts_lcens = rep_vector(0.0, qlcens);
        if (qrcens > 0) eta_qpts_rcens = rep_vector(0.0, qrcens);
        if (qicens > 0) eta_qpts_icenl = rep_vector(0.0, qicens);
        if (qicens > 0) eta_qpts_icenu = rep_vector(0.0, qicens);
        if (qdelay > 0) eta_qpts_delay = rep_vector(0.0, qdelay);
      }

      // add on time-varying part to linear predictor
      if (S > 0) {
        if (Nevent > 0) eta_epts_event += s_epts_event * beta_tve;
        if (qevent > 0) eta_qpts_event += s_qpts_event * beta_tve;
        if (qlcens > 0) eta_qpts_lcens += s_qpts_lcens * beta_tve;
        if (qrcens > 0) eta_qpts_rcens += s_qpts_rcens * beta_tve;
        if (qicens > 0) eta_qpts_icenl += s_qpts_icenl * beta_tve;
        if (qicens > 0) eta_qpts_icenu += s_qpts_icenu * beta_tve;
        if (qdelay > 0) eta_qpts_delay += s_qpts_delay * beta_tve;
      }

      // add on log crude event rate / time (helps to center intercept)
      if (Nevent > 0) eta_epts_event += log_crude_event_rate;
      if (qevent > 0) eta_qpts_event += log_crude_event_rate;
      if (qlcens > 0) eta_qpts_lcens += log_crude_event_rate;
      if (qrcens > 0) eta_qpts_rcens += log_crude_event_rate;
      if (qicens > 0) eta_qpts_icenl += log_crude_event_rate;
      if (qicens > 0) eta_qpts_icenu += log_crude_event_rate;
      if (qdelay > 0) eta_qpts_delay += log_crude_event_rate;

      // add on intercept to linear predictor
      if (has_intercept == 1) {
        if (Nevent > 0) eta_epts_event += gamma[1];
        if (qevent > 0) eta_qpts_event += gamma[1];
        if (qlcens > 0) eta_qpts_lcens += gamma[1];
        if (qrcens > 0) eta_qpts_rcens += gamma[1];
        if (qicens > 0) eta_qpts_icenl += gamma[1];
        if (qicens > 0) eta_qpts_icenu += gamma[1];
        if (qdelay > 0) eta_qpts_delay += gamma[1];
      }

      // add on random effects terms to linear predictor
      if (t > 0) {
        if (special_case) for (i in 1:t) {
          if (Nevent > 0) eta_epts_event += b[V_epts_event[i]];
          if (qevent > 0) eta_qpts_event += b[V_qpts_event[i]];
          if (qlcens > 0) eta_qpts_lcens += b[V_qpts_lcens[i]];
          if (qrcens > 0) eta_qpts_rcens += b[V_qpts_rcens[i]];
          if (qicens > 0) eta_qpts_icenl += b[V_qpts_icens[i]];
          if (qicens > 0) eta_qpts_icenu += b[V_qpts_icens[i]];
          if (qdelay > 0) eta_qpts_delay += b[V_qpts_delay[i]];
        }
        else {
          if (Nevent > 0) eta_epts_event +=
            csr_matrix_times_vector(Nevent, q, w_epts_event, v_epts_event, u_epts_event, b);
          if (qevent > 0) eta_qpts_event +=
            csr_matrix_times_vector(qevent, q, w_qpts_event, v_qpts_event, u_qpts_event, b);
          if (qlcens > 0) eta_qpts_lcens +=
            csr_matrix_times_vector(qlcens, q, w_qpts_lcens, v_qpts_lcens, u_qpts_lcens, b);
          if (qrcens > 0) eta_qpts_rcens +=
            csr_matrix_times_vector(qrcens, q, w_qpts_rcens, v_qpts_rcens, u_qpts_rcens, b);
          if (qicens > 0) eta_qpts_icenl +=
            csr_matrix_times_vector(qicens, q, w_qpts_icens, v_qpts_icens, u_qpts_icens, b);
          if (qicens > 0) eta_qpts_icenu +=
            csr_matrix_times_vector(qicens, q, w_qpts_icens, v_qpts_icens, u_qpts_icens, b);
          if (qdelay > 0) eta_qpts_delay +=
            csr_matrix_times_vector(qdelay, q, w_qpts_delay, v_qpts_delay, u_qpts_delay, b);
        }

      }

      // aft models
      if (type == 7 || type == 8) {

        vector[Nevent] af_event;

        vector[Nevent] caf_event;
        vector[Nlcens] caf_lcens;
        vector[Nrcens] caf_rcens;
        vector[Nicens] caf_icenl;
        vector[Nicens] caf_icenu;
        vector[Ndelay] caf_delay;

        // acceleration factor at event time
        if (Nevent > 0) af_event = exp(-eta_epts_event);

        // evaluate cumulative acceleration factors
        if (Nevent > 0) caf_event = quadrature_aft(qwts_event, eta_qpts_event, qnodes, Nevent);
        if (Nlcens > 0) caf_lcens = quadrature_aft(qwts_lcens, eta_qpts_lcens, qnodes, Nlcens);
        if (Nrcens > 0) caf_rcens = quadrature_aft(qwts_rcens, eta_qpts_rcens, qnodes, Nrcens);
        if (Nicens > 0) caf_icenl = quadrature_aft(qwts_icenl, eta_qpts_icenl, qnodes, Nicens);
        if (Nicens > 0) caf_icenu = quadrature_aft(qwts_icenu, eta_qpts_icenu, qnodes, Nicens);
        if (Ndelay > 0) caf_delay = quadrature_aft(qwts_delay, eta_qpts_delay, qnodes, Ndelay);

        // increment target with log-lik contributions
        if (type == 7) { // exponential AFT model
          if (Nevent > 0) target +=  exponentialAFT_log_haz (af_event);
          if (Nevent > 0) target +=  exponentialAFT_log_surv(caf_event);
          if (Nlcens > 0) target +=  exponentialAFT_log_cdf1(caf_lcens);
          if (Nrcens > 0) target +=  exponentialAFT_log_surv(caf_rcens);
          if (Nicens > 0) target +=  exponentialAFT_log_cdf2(caf_icenl, caf_icenu);
          if (Ndelay > 0) target += -exponentialAFT_log_surv(caf_delay);
        } else if (type == 8) { // weibull AFT model
          real shape = coefs[1];
          if (Nevent > 0) target +=  weibullAFT_log_haz (af_event, caf_event, shape);
          if (Nevent > 0) target +=  weibullAFT_log_surv(caf_event, shape);
          if (Nlcens > 0) target +=  weibullAFT_log_cdf1(caf_lcens, shape);
          if (Nrcens > 0) target +=  weibullAFT_log_surv(caf_rcens, shape);
          if (Nicens > 0) target +=  weibullAFT_log_cdf2(caf_icenl, caf_icenu, shape);
          if (Ndelay > 0) target += -weibullAFT_log_surv(caf_delay, shape);
        }

      }

      // hazard models
      else {

        vector[Nevent] lhaz_epts_event;
        vector[qevent] lhaz_qpts_event;
        vector[qlcens] lhaz_qpts_lcens;
        vector[qrcens] lhaz_qpts_rcens;
        vector[qicens] lhaz_qpts_icenl;
        vector[qicens] lhaz_qpts_icenu;
        vector[qdelay] lhaz_qpts_delay;

        // evaluate log hazard
        if (type == 5) { // exponential model
          if (Nevent > 0) lhaz_epts_event = exponential_log_haz(eta_epts_event);
          if (qevent > 0) lhaz_qpts_event = exponential_log_haz(eta_qpts_event);
          if (qlcens > 0) lhaz_qpts_lcens = exponential_log_haz(eta_qpts_lcens);
          if (qrcens > 0) lhaz_qpts_rcens = exponential_log_haz(eta_qpts_rcens);
          if (qicens > 0) lhaz_qpts_icenl = exponential_log_haz(eta_qpts_icenl);
          if (qicens > 0) lhaz_qpts_icenu = exponential_log_haz(eta_qpts_icenu);
          if (qdelay > 0) lhaz_qpts_delay = exponential_log_haz(eta_qpts_delay);
        }
        else if (type == 1) { // weibull model
          real shape = coefs[1];
          if (Nevent > 0) lhaz_epts_event = weibull_log_haz(eta_epts_event, epts_event, shape);
          if (qevent > 0) lhaz_qpts_event = weibull_log_haz(eta_qpts_event, qpts_event, shape);
          if (qlcens > 0) lhaz_qpts_lcens = weibull_log_haz(eta_qpts_lcens, qpts_lcens, shape);
          if (qrcens > 0) lhaz_qpts_rcens = weibull_log_haz(eta_qpts_rcens, qpts_rcens, shape);
          if (qicens > 0) lhaz_qpts_icenl = weibull_log_haz(eta_qpts_icenl, qpts_icenl, shape);
          if (qicens > 0) lhaz_qpts_icenu = weibull_log_haz(eta_qpts_icenu, qpts_icenu, shape);
          if (qdelay > 0) lhaz_qpts_delay = weibull_log_haz(eta_qpts_delay, qpts_delay, shape);
        }
        else if (type == 6) { // gompertz model
          real scale = coefs[1];
          if (Nevent > 0) lhaz_epts_event = gompertz_log_haz(eta_epts_event, epts_event, scale);
          if (qevent > 0) lhaz_qpts_event = gompertz_log_haz(eta_qpts_event, qpts_event, scale);
          if (qlcens > 0) lhaz_qpts_lcens = gompertz_log_haz(eta_qpts_lcens, qpts_lcens, scale);
          if (qrcens > 0) lhaz_qpts_rcens = gompertz_log_haz(eta_qpts_rcens, qpts_rcens, scale);
          if (qicens > 0) lhaz_qpts_icenl = gompertz_log_haz(eta_qpts_icenl, qpts_icenl, scale);
          if (qicens > 0) lhaz_qpts_icenu = gompertz_log_haz(eta_qpts_icenu, qpts_icenu, scale);
          if (qdelay > 0) lhaz_qpts_delay = gompertz_log_haz(eta_qpts_delay, qpts_delay, scale);
        }
        else if (type == 4) { // M-splines, on haz scale
          if (Nevent > 0) lhaz_epts_event = mspline_log_haz(eta_epts_event, basis_epts_event, ms_coefs);
          if (qevent > 0) lhaz_qpts_event = mspline_log_haz(eta_qpts_event, basis_qpts_event, ms_coefs);
          if (qlcens > 0) lhaz_qpts_lcens = mspline_log_haz(eta_qpts_lcens, basis_qpts_lcens, ms_coefs);
          if (qrcens > 0) lhaz_qpts_rcens = mspline_log_haz(eta_qpts_rcens, basis_qpts_rcens, ms_coefs);
          if (qicens > 0) lhaz_qpts_icenl = mspline_log_haz(eta_qpts_icenl, basis_qpts_icenl, ms_coefs);
          if (qicens > 0) lhaz_qpts_icenu = mspline_log_haz(eta_qpts_icenu, basis_qpts_icenu, ms_coefs);
          if (qdelay > 0) lhaz_qpts_delay = mspline_log_haz(eta_qpts_delay, basis_qpts_delay, ms_coefs);
        }
        else if (type == 2) { // B-splines, on log haz scale
          if (Nevent > 0) lhaz_epts_event = bspline_log_haz(eta_epts_event, basis_epts_event, coefs);
          if (qevent > 0) lhaz_qpts_event = bspline_log_haz(eta_qpts_event, basis_qpts_event, coefs);
          if (qlcens > 0) lhaz_qpts_lcens = bspline_log_haz(eta_qpts_lcens, basis_qpts_lcens, coefs);
          if (qrcens > 0) lhaz_qpts_rcens = bspline_log_haz(eta_qpts_rcens, basis_qpts_rcens, coefs);
          if (qicens > 0) lhaz_qpts_icenl = bspline_log_haz(eta_qpts_icenl, basis_qpts_icenl, coefs);
          if (qicens > 0) lhaz_qpts_icenu = bspline_log_haz(eta_qpts_icenu, basis_qpts_icenu, coefs);
          if (qdelay > 0) lhaz_qpts_delay = bspline_log_haz(eta_qpts_delay, basis_qpts_delay, coefs);
        }
        else {
          reject("Bug found: invalid baseline hazard (with quadrature).");
        }

        // increment target with log-lik contributions for event submodel
        if (Nevent > 0) target +=  lhaz_epts_event;
        if (qevent > 0) target +=  quadrature_log_surv(qwts_event, lhaz_qpts_event);
        if (qlcens > 0) target +=  quadrature_log_cdf1(qwts_lcens, lhaz_qpts_lcens, qnodes, Nlcens);
        if (qrcens > 0) target +=  quadrature_log_surv(qwts_rcens, lhaz_qpts_rcens);
        if (qicens > 0) target +=  quadrature_log_cdf2(qwts_icenl, lhaz_qpts_icenl,
                                                       qwts_icenu, lhaz_qpts_icenu, qnodes, Nicens);
        if (qdelay > 0) target += -quadrature_log_surv(qwts_delay, lhaz_qpts_delay);

      }

    }

  }

  //-------- log priors

  // log priors for coefficients
  if (K > 0) {
    target += beta_custom_lpdf(z_beta | prior_dist, prior_scale, prior_df,
                               global_prior_df, local, global, mix, ool,
                               slab_df, caux);
  }

  // log prior for intercept
  if (has_intercept == 1) {
    target += gamma_custom_lpdf(gamma[1] | prior_dist_for_intercept,
                                prior_mean_for_intercept, prior_scale_for_intercept,
                                prior_df_for_intercept);
  }

  // log priors for baseline hazard parameters
  if (type == 4) {
    target += basehaz_lpdf(ms_coefs | prior_dist_for_aux, prior_conc_for_aux);
  }
  else if (nvars > 0) {
    target += basehaz_lpdf(z_coefs | prior_dist_for_aux, prior_df_for_aux);
  }

  // log priors for tve spline coefficients and their smoothing parameters
  if (S > 0) {
    target += smooth_lpdf(z_beta_tve | smooth_sd_raw,
                          prior_dist_for_smooth, prior_df_for_smooth);
  }

  // log prior for random effects
  if (t > 0) {
    target += decov_lpdf(z_b | z_T, rho, zeta, tau,
                         regularization, delta, b_prior_shape, t, p);
  }

}

generated quantities {
  // baseline hazard parameters to return
  vector[nvars] aux = (type == 4) ? ms_coefs : coefs;

  // transformed intercept
  real alpha;
  if (has_intercept == 1) {
    alpha = log_crude_event_rate - dot_product(x_bar, beta) + gamma[1];
  } else {
    alpha = log_crude_event_rate - dot_product(x_bar, beta);
  }
}
