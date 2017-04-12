#include "Columbia_copyright.stan"
#include "license.stan" // GPL3+

// GLM for a Bernoulli outcome
functions {
  #include "common_functions.stan"
  #include "bernoulli_likelihoods.stan"  
}
data {
  // dimensions
  int<lower=0> K;        // number of predictors
  int<lower=1> N[2];     // number of observations where y = 0 and y = 1 respectively
  vector[K] xbar;        // vector of column-means of rbind(X0, X1)
  int<lower=0,upper=1> dense_X; // flag for dense vs. sparse
  matrix[N[1],K] X0[dense_X];   // centered (by xbar) predictor matrix | y = 0
  matrix[N[2],K] X1[dense_X];   // centered (by xbar) predictor matrix | y = 1
  
  // stuff for the sparse case
  int<lower=0> nnz_X0;                       // number of non-zero elements in the implicit X0 matrix
  vector[nnz_X0] w_X0;                       // non-zero elements in the implicit X0 matrix
  int<lower=0> v_X0[nnz_X0];                 // column indices for w_X0
  int<lower=0> u_X0[dense_X ? 0 : N[1] + 1]; // where the non-zeros start in each row of X0
  int<lower=0> nnz_X1;                       // number of non-zero elements in the implicit X1 matrix
  vector[nnz_X1] w_X1;                       // non-zero elements in the implicit X1 matrix
  int<lower=0> v_X1[nnz_X1];                 // column indices for w_X1
  int<lower=0> u_X1[dense_X ? 0 : N[2] + 1]; // where the non-zeros start in each row of X1
  
  #include "data_glm.stan" // declares prior_PD, has_intercept, family, link, prior_dist, prior_dist_for_intercept

  // weights
  int<lower=0,upper=1> has_weights;  // 0 = No, 1 = Yes
  vector[has_weights ? N[1] : 0] weights0;
  vector[has_weights ? N[2] : 0] weights1;
  
  // offset
  int<lower=0,upper=1> has_offset;  // 0 = No, 1 = Yes
  vector[has_offset ? N[1] : 0] offset0;
  vector[has_offset ? N[2] : 0] offset1;
  
  // declares prior_{mean, scale, df}, prior_{mean, scale, df}_for_intercept, prior_{mean, scale, df}_for_aux
  #include "hyperparameters.stan"
  // declares t, p[t], l[t], q, len_theta_L, shape, scale, {len_}concentration, {len_}regularization
  #include "glmer_stuff.stan"  

  // more glmer stuff
  int<lower=0> num_non_zero[2];     // number of non-zero elements in the Z matrices
  vector[num_non_zero[1]] w0;       // non-zero elements in the implicit Z0 matrix
  vector[num_non_zero[2]] w1;       // non-zero elements in the implicit Z1 matrix
  int<lower=0> v0[num_non_zero[1]]; // column indices for w0
  int<lower=0> v1[num_non_zero[2]]; // column indices for w1
  int<lower=0> u0[t > 0 ? N[1] + 1 : 0];  // where the non-zeros start in each row of Z0
  int<lower=0> u1[t > 0 ? N[2] + 1 : 0];  // where the non-zeros start in each row of Z1
  int<lower=0, upper=1> special_case;     // whether we only have to deal with (1|group)
}
transformed data {
  int NN = N[1] + N[2];
  real aux = not_a_number();
  int<lower=1> V0[special_case ? t : 0,N[1]] = make_V(N[1], special_case ? t : 0, v0);
  int<lower=1> V1[special_case ? t : 0,N[2]] = make_V(N[2], special_case ? t : 0, v1);
  #include "tdata_glm.stan"// defines hs, len_z_T, len_var_group, delta, pos
}
parameters {
  real<upper=(link == 4 ? 0.0 : positive_infinity())> gamma[has_intercept];
  #include "parameters_glm.stan" // declares z_beta, global, local, z_b, z_T, rho, zeta, tau
}
transformed parameters {
  #include "tparameters_glm.stan" // defines beta, b, theta_L
  if (t > 0) {
    if (special_case) {
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
  #include "make_eta_bern.stan" // defines eta0, eta1
  if (has_intercept == 1) {
    if (link != 4) {
      eta0 = gamma[1] + eta0;
      eta1 = gamma[1] + eta1;
    }
    else {
      real shift;
      shift = fmax(max(eta0), max(eta1));
      eta0 = gamma[1] + eta0 - shift;
      eta1 = gamma[1] + eta1 - shift;
    }
  }
  // Log-likelihood 
  if (has_weights == 0 && prior_PD == 0) {  // unweighted log-likelihoods
    real dummy;  // irrelevant but useful for testing
    dummy = ll_bern_lp(eta0, eta1, link, N);
  }
  else if (prior_PD == 0) {  // weighted log-likelihoods
    target += dot_product(weights0, pw_bern(0, eta0, link));
    target += dot_product(weights1, pw_bern(1, eta1, link));
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
    vector[N[1]] pi0;
    vector[N[2]] pi1;
    #include "make_eta_bern.stan" // defines eta0, eta1
    if (has_intercept == 1) {
      if (link != 4) {
        eta0 = gamma[1] + eta0;
        eta1 = gamma[1] + eta1;
      }      
      else {
        real shift;
        shift = fmax(max(eta0), max(eta1));
        eta0 = gamma[1] + eta0 - shift;
        eta1 = gamma[1] + eta1 - shift;
        alpha[1] = alpha[1] - shift;
      }
    }
    pi0 = linkinv_bern(eta0, link);
    pi1 = linkinv_bern(eta1, link);
    for (n in 1:N[1]) mean_PPD = mean_PPD + bernoulli_rng(pi0[n]);
    for (n in 1:N[2]) mean_PPD = mean_PPD + bernoulli_rng(pi1[n]);
    mean_PPD = mean_PPD / NN;
  }
}
