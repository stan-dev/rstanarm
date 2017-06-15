#include "Columbia_copyright.stan"
#include "license.stan" // GPL3+

// GLM for a binomial outcome
functions {
  #include "common_functions.stan"
  #include "binomial_likelihoods.stan"  
}
data {
  #include "NKX.stan"      // declares N, K, X, xbar, dense_X, nnz_x, w_x, v_x, u_x
  int<lower=0> y[N];       // outcome: number of successes
  int<lower=0> trials[N];  // number of trials
  #include "data_glm.stan" // declares prior_PD, has_intercept, family, link, prior_dist, prior_dist_for_intercept
  #include "weights_offset.stan"  // declares has_weights, weights, has_offset, offset
  // declares prior_{mean, scale, df}, prior_{mean, scale, df}_for_intercept, prior_scale_{mean, scale, df}_for_aux
  #include "hyperparameters.stan"
  // declares t, p[t], l[t], q, len_theta_L, shape, scale, {len_}concentration, {len_}regularization
  #include "glmer_stuff.stan"  
  #include "glmer_stuff2.stan" // declares num_not_zero, w, v, u
}
transformed data {
  real aux = not_a_number();
  int<lower=1> V[special_case ? t : 0, N] = make_V(N, special_case ? t : 0, v);
  #include "tdata_glm.stan"// defines hs, len_z_T, len_var_group, delta, pos
}
parameters {
  real<upper=(link == 4 ? 0.0 : positive_infinity())> gamma[has_intercept];
  #include "parameters_glm.stan" // declares z_beta, global, local, z_b, z_T, rho, zeta, tau
}
transformed parameters {
  #include "tparameters_glm.stan" // defines beta, b, theta_L
  if (t > 0) {
    if (special_case == 1) {
      int start = 1;
      theta_L = scale .* tau;
      if (t == 1) b = theta_L[1] * z_b;
      else for (i in 1:t) {
        int end = start + l[i] - 1;
        b[start:end] = theta_L[i] * z_b[start:end];
        start = end + 1;
      }
    }
    else {
      theta_L = make_theta_L(len_theta_L, p, 
                             1.0, tau, scale, zeta, rho, z_T);
      b = make_b(z_b, theta_L, p, l);
    }
  }
}
model {
  #include "make_eta.stan" // defines eta
  if (t > 0) {
    #include "eta_add_Zb.stan"
  }
  if (has_intercept == 1) {
    if (link != 4) eta = eta + gamma[1];
    else eta = gamma[1] + eta - max(eta);
  }
  else {
    #include "eta_no_intercept.stan" // shifts eta
  }
  
  // Log-likelihood 
  if (has_weights == 0 && prior_PD == 0) {  // unweighted log-likelihoods
    real dummy;  // irrelevant but useful for testing
    dummy = ll_binom_lp(y, trials, eta, link);
  }
  else if (prior_PD == 0) 
    target += dot_product(weights, pw_binom(y, trials, eta, link));
  
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
    vector[N] pi;
    #include "make_eta.stan" // defines eta
    if (t > 0) {
      #include "eta_add_Zb.stan"
    }
    if (has_intercept == 1) {
      if (link != 4) eta = eta + gamma[1];
      else {
        real shift;
        shift = max(eta);
        eta = gamma[1] + eta - shift;
        alpha[1] = alpha[1] - shift;
      }
    }
    else {
      #include "eta_no_intercept.stan" // shifts eta
    }
    
    pi = linkinv_binom(eta, link);
    for (n in 1:N) mean_PPD = mean_PPD + binomial_rng(trials[n], pi[n]);
    mean_PPD = mean_PPD / N;
  }
}
