#include "Columbia_copyright.stan"
#include "license.stan" // GPL3+

// GLM for a Gaussian, Gamma, or inverse Gaussian outcome
functions {
  #include "common_functions.stan"
  #include "continuous_likelihoods.stan"
  
  /** 
  * test function for csr_matrix_times_vector
  *
  * @param m Integer number of rows
  * @param n Integer number of columns
  * @param w Vector (see reference manual)
  * @param v Integer array (see reference manual)
  * @param u Integer array (see reference manual)
  * @param b Vector that is multiplied from the left by the CSR matrix
  * @return A vector that is the product of the CSR matrix and b
  */
  vector test_csr_matrix_times_vector(int m, int n, vector w, 
                                      int[] v, int[] u, vector b) {
    return csr_matrix_times_vector(m, n, w, v, u, b); 
  }
  
}
data {
  #include "NKX.stan"      // declares N, K, X, xbar, dense_X, nnz_x, w_x, v_x, u_x
  real lb_y; // lower bound on y
  real<lower=lb_y> ub_y; // upper bound on y
  vector<lower=lb_y, upper=ub_y>[N] y; // continuous outcome
  #include "data_glm.stan" // declares prior_PD, has_intercept, family, link, prior_dist, prior_dist_for_intercept
  #include "weights_offset.stan"  // declares has_weights, weights, has_offset, offset
  // declares prior_{mean, scale, df}, prior_{mean, scale, df}_for_intercept, prior_scale_for dispersion
  #include "hyperparameters.stan"
  // declares t, p[t], l[t], q, len_theta_L, shape, scale, {len_}concentration, {len_}regularization
  #include "glmer_stuff.stan"  
  #include "glmer_stuff2.stan" // declares num_not_zero, w, v, u
}
transformed data {
  vector[family == 3 ? N : 0] sqrt_y;
  vector[family == 3 ? N : 0] log_y;
  real sum_log_y = family == 1 ? not_a_number() : sum(log(y));
  #include "tdata_glm.stan"// defines hs, len_z_T, len_var_group, delta, pos, t_{any, all}_124
  if (family == 3) {
    sqrt_y = sqrt(y);
    log_y = log(y);
  }
}
parameters {
  real<lower=(family == 1 || link == 2 ? negative_infinity() : 0.0)> gamma[has_intercept];
  #include "parameters_glm.stan" // declares z_beta, global, local, z_b, z_T, rho, zeta, tau
  real<lower=0> dispersion_unscaled; # interpretation depends on family!
}
transformed parameters {
  real dispersion;
  #include "tparameters_glm.stan" // defines beta, b, theta_L
  
  if (prior_dist_for_dispersion == 0)
    dispersion = dispersion_unscaled;
  else {
    dispersion = prior_scale_for_dispersion * dispersion_unscaled;
    if (prior_dist_for_dispersion <= 2) // normal or student_t
      dispersion = dispersion + prior_mean_for_dispersion;
  }
    
  if (t > 0) {
    theta_L = make_theta_L(len_theta_L, p, 
                            dispersion, tau, scale, zeta, rho, z_T);
    b = make_b(z_b, theta_L, p, l);
  }
}
model {
  #include "make_eta.stan" // defines eta
  if (t > 0) eta = eta + csr_matrix_times_vector(N, q, w, v, u, b);
  if (has_intercept == 1) {
    if (family == 1 || link == 2) eta = eta + gamma[1];
    else eta = eta - min(eta) + gamma[1];
  }
  else {
    #include "eta_no_intercept.stan" // shifts eta
  }
  
  // Log-likelihood 
  if (has_weights == 0 && prior_PD == 0) { # unweighted log-likelihoods
    if (family == 1) {
      if (link == 1) 
        target += normal_lpdf(y | eta, dispersion);
      else if (link == 2) 
        target += normal_lpdf(y | exp(eta), dispersion);
      else 
        target += normal_lpdf(y | divide_real_by_vector(1, eta), dispersion);
      // divide_real_by_vector() is defined in common_functions.stan
    }
    else if (family == 2) {
      target += GammaReg(y, eta, dispersion, link, sum_log_y);
    }
    else {
      target += inv_gaussian(y, linkinv_inv_gaussian(eta, link), 
                             dispersion, sum_log_y, sqrt_y);
    }
  }
  else if (prior_PD == 0) { # weighted log-likelihoods
    vector[N] summands;
    if (family == 1) summands = pw_gauss(y, eta, dispersion, link);
    else if (family == 2) summands = pw_gamma(y, eta, dispersion, link);
    else summands = pw_inv_gaussian(y, eta, dispersion, link, log_y, sqrt_y);
    target += dot_product(weights, summands);
  }

  // Log-priors
  if (prior_dist_for_dispersion > 0 && prior_scale_for_dispersion > 0) {
    if (prior_dist_for_dispersion == 1)
      target += normal_lpdf(dispersion_unscaled | 0, 1);
    else if (prior_dist_for_dispersion == 2)
      target += student_t_lpdf(dispersion_unscaled | 
                               prior_df_for_dispersion, 0, 1);
    else 
     target += exponential_lpdf(dispersion_unscaled | 1);
  }
    
  #include "priors_glm.stan" // increments target()
  if (t > 0) decov_lp(z_b, z_T, rho, zeta, tau, 
                      regularization, delta, shape, t, p);
}
generated quantities {
  real alpha[has_intercept];
  real mean_PPD = 0;
  if (has_intercept == 1) {
    if (dense_X) alpha[1] = gamma[1] - dot_product(xbar, beta);
    else alpha[1] = gamma[1];
  }
  {
    #include "make_eta.stan" // defines eta
    if (t > 0) eta = eta + csr_matrix_times_vector(N, q, w, v, u, b);
    if (has_intercept == 1) {
      if (family == 1 || link == 2) eta = eta + gamma[1];
      else {
        real min_eta = min(eta);
        alpha[1] = alpha[1] - min_eta;
        eta = eta - min_eta + gamma[1];
      }
    }
    else {
      #include "eta_no_intercept.stan" // shifts eta
    }
    
    if (family == 1) {
      if (link > 1) eta = linkinv_gauss(eta, link);
      for (n in 1:N) mean_PPD = mean_PPD + normal_rng(eta[n], dispersion);
    }
    else if (family == 2) {
      if (link > 1) eta = linkinv_gamma(eta, link);
      for (n in 1:N) mean_PPD = mean_PPD + gamma_rng(dispersion, dispersion / eta[n]);
    }
    else if (family == 3) {
      if (link > 1) eta = linkinv_inv_gaussian(eta, link);
      for (n in 1:N) mean_PPD = mean_PPD + inv_gaussian_rng(eta[n], dispersion);
    }
    mean_PPD = mean_PPD / N;
  }
}
