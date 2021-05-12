#include /pre/Columbia_copyright.stan
#include /pre/license.stan

// GLM for a Gaussian, Gamma, inverse Gaussian, or Beta outcome
functions {
#include /functions/common_functions.stan
#include /functions/continuous_likelihoods.stan
#include /functions/SSfunctions.stan

  /**
   * Increments the log-posterior with the logarithm of a multivariate normal 
   * likelihood with a scalar standard deviation for all errors
   * Equivalent to normal_lpdf(y | intercept + X * beta + Z * b, sigma) but faster
   * @param coeff vector of coefficients (including intercept)
   * @param OLS precomputed vector of OLS coefficients (including intercept)
   * @param XtX precomputed matrix equal to crossprod(X) (including intercept)
   * @param SSR positive precomputed value of the sum of squared OLS residuals
   * @param sigma positive scalar for the standard deviation of the errors
   * @param N integer equal to the number of observations
   */
  real ll_mvn_ols(vector coeff, vector OLS, matrix XtX,
                  real SSR, real sigma, int N) {
    return -0.5 * (quad_form(XtX, coeff - OLS) + SSR) / square(sigma)
            - N * (log(sigma) + log(sqrt(2 * pi())));
  }

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
  // declares N, K, X, xbar, dense_X, nnz_x, w_x, v_x, u_x
#include /data/NKX.stan
  int<lower=0> len_y;      // length of y
  real lb_y; // lower bound on y
  real<lower=lb_y> ub_y; // upper bound on y
  vector<lower=lb_y, upper=ub_y>[len_y] y; // continuous outcome
  int<lower=1,upper=4> family; // 1 gaussian, 2 gamma, 3 inv-gaussian, 4 beta
  // declares prior_PD, has_intercept, link, prior_dist, prior_dist_for_intercept
#include /data/data_glm.stan
  // declares has_weights, weights, has_offset, offset
#include /data/weights_offset.stan
  // declares prior_{mean, scale, df}, prior_{mean, scale, df}_for_intercept, prior_{mean, scale, df}_for_aux
#include /data/hyperparameters.stan
  // declares t, p[t], l[t], q, len_theta_L, shape, scale, {len_}concentration, {len_}regularization
#include /data/glmer_stuff.stan
  // declares num_not_zero, w, v, u
#include /data/glmer_stuff2.stan
#include /data/data_betareg.stan

  int<lower=0,upper=10> SSfun; // nonlinear function indicator, 0 for identity
  vector[SSfun > 0  ? len_y : 0] input;
  vector[SSfun == 5 ? len_y : 0] Dose;
}
transformed data {
  vector[family == 3 ? len_y : 0] sqrt_y;
  vector[family == 3 ? len_y : 0] log_y;
  real sum_log_y = family == 1 ? not_a_number() : sum(log(y));
  int<lower=1> V[special_case ? t : 0, len_y] = make_V(len_y, special_case ? t : 0, v);
  int<lower=0> hs_z;                  // for tdata_betareg.stan
  int can_do_OLS = family == 1 && link == 1 && SSfun == 0 && has_offset == 0 && t == 0 && 
                   prior_PD == 0 && dense_X && N > 2 && len_y >= (has_intercept + K + K_smooth);
  vector[can_do_OLS ? has_intercept + K + K_smooth : 0] OLS;
  matrix[can_do_OLS ? has_intercept + K + K_smooth : 0, can_do_OLS ? has_intercept + K + K_smooth : 0] XtX;
  int can_do_normalidglm = K != 0 &&  // remove K!=0 after rstan includes this Stan bugfix: https://github.com/stan-dev/math/issues/1398 
                           can_do_OLS == 0 && family == 1 && link == 1 && 
                           SSfun == 0 && has_offset == 0 && dense_X && prior_PD == 0 && 
                           t == 0 && len_y < (has_intercept + K + K_smooth);
  matrix[can_do_normalidglm ? N : 0, can_do_normalidglm ? K + K_smooth : 0] XS;
  real SSR = not_a_number();
  // defines hs, len_z_T, len_var_group, delta, is_continuous, pos
#include /tdata/tdata_glm.stan
  // defines hs_z
#include /tdata/tdata_betareg.stan
  is_continuous = 1;

  if (family == 3) {
    sqrt_y = sqrt(y);
    log_y = log(y);
  }
  if (can_do_OLS) {
    matrix[N, has_intercept + K + K_smooth ] X_ = has_intercept ? append_col(rep_vector(1.0, N), 
                                                  (K_smooth > 0 ? append_col(X[1], S) : X[1])) : 
                                                  (K_smooth > 0 ? append_col(X[1], S) : X[1]);
    matrix[cols(X_), cols(X_)] R = qr_thin_R(X_);
    if (tail(diagonal(R), 1)[1] > 1e-16) {
      matrix[N, cols(R)] Q = qr_thin_Q(X_);
      XtX = crossprod(X_);
      OLS = mdivide_right_tri_low(y' * Q, R')';
      SSR = dot_self(y - X_ * OLS);
    } else can_do_OLS = 0;
  }
  if (can_do_normalidglm) {
    XS = K_smooth > 0 ? append_col(X[1], S) : X[1];
  }
}
parameters {
  real<lower=make_lower(family, link),upper=make_upper(family,link)> gamma[has_intercept];
  // declares z_beta, global, local, z_b, z_T, rho, zeta, tau
#include /parameters/parameters_glm.stan
  real<lower=0> aux_unscaled; // interpretation depends on family!
#include /parameters/parameters_betareg.stan
}
transformed parameters {
  // aux has to be defined first in the hs case
  real aux = prior_dist_for_aux == 0 ? aux_unscaled : (prior_dist_for_aux <= 2 ? 
             prior_scale_for_aux * aux_unscaled + prior_mean_for_aux :
             prior_scale_for_aux * aux_unscaled);

  vector[z_dim] omega; // used in tparameters_betareg.stan
  // defines beta, b, theta_L
#include /tparameters/tparameters_glm.stan
#include /tparameters/tparameters_betareg.stan
  
  if (prior_dist_for_aux == 0) // none
    aux = aux_unscaled;
  else {
    aux = prior_scale_for_aux * aux_unscaled;
    if (prior_dist_for_aux <= 2) // normal or student_t
      aux += prior_mean_for_aux;
  }

  if (t > 0) {
    if (special_case == 1) {
      int start = 1;
      theta_L = scale .* tau * aux;
      if (t == 1) b = theta_L[1] * z_b;
      else for (i in 1:t) {
        int end = start + l[i] - 1;
        b[start:end] = theta_L[i] * z_b[start:end];
        start = end + 1;
      }
    }
    else {
      theta_L = make_theta_L(len_theta_L, p, 
                             aux, tau, scale, zeta, rho, z_T);
      b = make_b(z_b, theta_L, p, l);
    }
  }
}
model {
  if (can_do_OLS) {
    vector[cols(XtX)] coeff = has_intercept ? append_row(to_vector(gamma), 
                                              (K_smooth > 0 ? append_row(beta, beta_smooth) : beta)) : 
                                              (K_smooth > 0 ? append_row(beta, beta_smooth) : beta);
    target += ll_mvn_ols(coeff, OLS, XtX, SSR, aux, N);
  } else if (can_do_normalidglm) {
    vector[K + K_smooth] coeff = K_smooth > 0 ? append_row(beta, beta_smooth) : beta;
    target += normal_id_glm_lpdf(y | XS, has_intercept ? gamma[1] : 0.0, coeff, aux);
  } else if (prior_PD == 0) {
    vector[link_phi > 0 ? N : 0] eta_z; // beta regression linear predictor for phi
#include /model/make_eta.stan
    if (t > 0) {
#include /model/eta_add_Zb.stan
    }
    if (has_intercept == 1) {
      if ((family == 1 || link == 2) || (family == 4 && link != 5)) eta += gamma[1];
      else if (family == 4 && link == 5) eta += gamma[1] - max(eta);
      else eta += gamma[1] - min(eta);
    }
    else {
#include /model/eta_no_intercept.stan
    }

    if (SSfun > 0) { // nlmer
      matrix[len_y, K] P = reshape_vec(eta, len_y, K);
      if (SSfun < 5) {
        if (SSfun <= 2) {
          if (SSfun == 1) target += normal_lpdf(y | SS_asymp(input, P), aux);
          else target += normal_lpdf(y | SS_asympOff(input, P), aux);
        }
        else if (SSfun == 3) target += normal_lpdf(y | SS_asympOrig(input, P), aux);
        else {
          for (i in 1:len_y) P[i,1] += exp(P[i,3]); // ordering constraint
          target += normal_lpdf(y | SS_biexp(input, P), aux);
        }
      }
      else {
        if (SSfun <= 7) {
          if (SSfun == 5) target += normal_lpdf(y | SS_fol(Dose, input, P), aux);
          else if (SSfun == 6) target += normal_lpdf(y | SS_fpl(input, P), aux);
          else target += normal_lpdf(y | SS_gompertz(input, P), aux);
        }
        else {
          if (SSfun == 8) target += normal_lpdf(y | SS_logis(input, P), aux);
          else if (SSfun == 9) target += normal_lpdf(y | SS_micmen(input, P), aux);
          else target += normal_lpdf(y | SS_weibull(input, P), aux);
        }
      }
    }
    else if (has_weights == 0) { // unweighted log-likelihoods
#include /model/make_eta_z.stan
      // adjust eta_z according to links
      if (has_intercept_z == 1) {
        if (link_phi > 1) {
          eta_z += gamma_z[1] - min(eta_z);
        }
        else {
          eta_z += gamma_z[1];
        }
      }
      else { // has_intercept_z == 0
#include /model/eta_z_no_intercept.stan
      }
      if (family == 1) {
        if (link == 1) 
          target += normal_lpdf(y | eta, aux);
        else if (link == 2) 
          target += normal_lpdf(y | exp(eta), aux);
        else 
          target += normal_lpdf(y | inv(eta), aux);
      }
      else if (family == 2) {
        target += GammaReg(y, eta, aux, link, sum_log_y);
      }
      else if (family == 3) {
        target += inv_gaussian(y, linkinv_inv_gaussian(eta, link), 
                               aux, sum_log_y, sqrt_y);
      }
      else if (family == 4 && link_phi == 0) {
        vector[N] mu;
        mu = linkinv_beta(eta, link);
        target += beta_lpdf(y | mu * aux, (1 - mu) * aux);
      }
      else if (family == 4 && link_phi > 0) {
        vector[N] mu;
        vector[N] mu_z;
        mu = linkinv_beta(eta, link);
        mu_z = linkinv_beta_z(eta_z, link_phi);
        target += beta_lpdf(y | rows_dot_product(mu, mu_z), 
                            rows_dot_product((1 - mu) , mu_z));
      }
    }
    else { // weighted log-likelihoods
      vector[N] summands;
      if (family == 1) summands = pw_gauss(y, eta, aux, link);
      else if (family == 2) summands = pw_gamma(y, eta, aux, link);
      else if (family == 3) summands = pw_inv_gaussian(y, eta, aux, link, log_y, sqrt_y);
      else if (family == 4 && link_phi == 0) summands = pw_beta(y, eta, aux, link);
      else if (family == 4 && link_phi > 0) summands = pw_beta_z(y, eta, eta_z, link, link_phi);
      target += dot_product(weights, summands);
    }
  }

  // Log-priors
  if (prior_dist_for_aux > 0 && prior_scale_for_aux > 0) {
    real log_half = -0.693147180559945286;
    if (prior_dist_for_aux == 1)
      target += normal_lpdf(aux_unscaled | 0, 1) - log_half;
    else if (prior_dist_for_aux == 2)
      target += student_t_lpdf(aux_unscaled | prior_df_for_aux, 0, 1) - log_half;
    else 
     target += exponential_lpdf(aux_unscaled | 1);
  }
    
#include /model/priors_glm.stan
#include /model/priors_betareg.stan
  if (t > 0) {
    real dummy = decov_lp(z_b, z_T, rho, zeta, tau, 
                          regularization, delta, shape, t, p);
  }
}
generated quantities {
  real mean_PPD = compute_mean_PPD ? 0 : negative_infinity();
  real alpha[has_intercept];
  real omega_int[has_intercept_z];
  
  if (has_intercept == 1) {
    if (dense_X) alpha[1] = gamma[1] - dot_product(xbar, beta);
    else alpha[1] = gamma[1];
  }
  if (has_intercept_z == 1) { 
    omega_int[1] = gamma_z[1] - dot_product(zbar, omega);  // adjust betareg intercept 
  }
  
  if (compute_mean_PPD) {
    vector[N] eta_z;
#include /model/make_eta.stan
    if (t > 0) {
#include /model/eta_add_Zb.stan
    }
    if (has_intercept == 1) {
      if (make_lower(family,link) == negative_infinity() &&
          make_upper(family,link) == positive_infinity()) eta += gamma[1];
      else if (family == 4 && link == 5) {
        real max_eta = max(eta);
        alpha[1] -= max_eta;
        eta += gamma[1] - max_eta;
      }
      else {
        real min_eta = min(eta);
        alpha[1] -= min_eta;
        eta += gamma[1] - min_eta;
      }
    }
    else {
#include /model/eta_no_intercept.stan
    }
#include /model/make_eta_z.stan
    // adjust eta_z according to links
    if (has_intercept_z == 1) {
      if (link_phi > 1) {
        omega_int[1] -= min(eta_z);
        eta_z += gamma_z[1] - min(eta_z);
      }
      else {
        eta_z += gamma_z[1];
      }
    }
    else { // has_intercept_z == 0
#include /model/eta_z_no_intercept.stan
    }
    
    if (SSfun > 0) { // nlmer
      vector[len_y] eta_nlmer;
      matrix[len_y, K] P;      
      P = reshape_vec(eta, len_y, K);
      if (SSfun < 5) {
        if (SSfun <= 2) {
          if (SSfun == 1) eta_nlmer = SS_asymp(input, P);
          else eta_nlmer = SS_asympOff(input, P);
        }
        else if (SSfun == 3) eta_nlmer = SS_asympOrig(input, P);
        else eta_nlmer = SS_biexp(input, P);
      }
      else {
        if (SSfun <= 7) {
          if (SSfun == 5) eta_nlmer = SS_fol(Dose, input, P);
          else if (SSfun == 6) eta_nlmer = SS_fpl(input, P);
          else eta_nlmer = SS_gompertz(input, P);
        }
        else {
          if (SSfun == 8) eta_nlmer = SS_logis(input, P);
          else if (SSfun == 9) eta_nlmer = SS_micmen(input, P);
          else eta_nlmer = SS_weibull(input, P);
        }
      }
      for (n in 1:len_y) mean_PPD += normal_rng(eta_nlmer[n], aux);
    }
    else if (family == 1) {
      vector[N] mu = link > 1 ? linkinv_gauss(eta, link) : eta;
      for (n in 1:len_y) mean_PPD += normal_rng(mu[n], aux);
    }
    else if (family == 2) {
      vector[N] mu = link > 1 ? linkinv_gamma(eta, link) : eta;
      for (n in 1:len_y) mean_PPD += gamma_rng(aux, aux / mu[n]);
    }
    else if (family == 3) {
      vector[N] mu = link > 1 ? linkinv_inv_gaussian(eta, link) : eta;
      for (n in 1:len_y) mean_PPD += inv_gaussian_rng(mu[n], aux);
    }
    else if (family == 4 && link_phi == 0) { 
      vector[N] mu = linkinv_beta(eta, link);
      for (n in 1:N) {
        real mu_n = mu[n];
        if (aux <= 0) mean_PPD += bernoulli_rng(0.5);
        else if (mu_n >= 1) mean_PPD += 1;
        else if (mu_n > 0)
          mean_PPD += beta_rng(mu_n * aux, (1 - mu_n) * aux);
      }
    }
    else if (family == 4 && link_phi > 0) {
      vector[N] mu = linkinv_beta(eta, link);
      vector[N] phi = linkinv_beta_z(eta_z, link_phi);
      for (n in 1:N) {
        real mu_n = mu[n];
        real aux_n = phi[n];
        if (aux_n <= 0) mean_PPD += bernoulli_rng(0.5);
        else if (mu_n >= 1) mean_PPD += 1;
        else if (mu_n > 0)
          mean_PPD += beta_rng(mu_n * aux_n, (1 - mu_n) * aux_n);
      }
    }
    mean_PPD /= len_y;
  }
}
