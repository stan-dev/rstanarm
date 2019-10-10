#include /pre/Columbia_copyright.stan
#include /pre/license.stan

// GLM for a Bernoulli outcome
functions {
#include /functions/common_functions.stan
#include /functions/bernoulli_likelihoods.stan
}
data {
  // dimensions
  int<lower=0> K;        // number of predictors
  int<lower=0> N[2];     // number of observations where y = 0 and y = 1 respectively
  vector[K] xbar;        // vector of column-means of rbind(X0, X1)
  int<lower=0,upper=1> dense_X; // flag for dense vs. sparse
  matrix[N[1],K] X0[dense_X];   // centered (by xbar) predictor matrix | y = 0
  matrix[N[2],K] X1[dense_X];   // centered (by xbar) predictor matrix | y = 1
  
  int<lower=0, upper=1> clogit; // 1 iff the number of successes is fixed in each stratum
  int<lower=0> J; // number of strata (possibly zero)
  int<lower=1,upper=J> strata[clogit == 1 ? N[1] + N[2] : 0];

  // stuff for the sparse case
  int<lower=0> nnz_X0;                       // number of non-zero elements in the implicit X0 matrix
  vector[nnz_X0] w_X0;                       // non-zero elements in the implicit X0 matrix
  int<lower=0, upper = K - 1> v_X0[nnz_X0];  // column indices for w_X0
  // where the non-zeros start in each row of X0
  int<lower=0, upper = rows(w_X0) + 1> u_X0[dense_X ? 0 : N[1] + 1]; 
  int<lower=0> nnz_X1;                       // number of non-zero elements in the implicit X1 matrix
  vector[nnz_X1] w_X1;                       // non-zero elements in the implicit X1 matrix
  int<lower=0, upper = K - 1> v_X1[nnz_X1];  // column indices for w_X1
  // where the non-zeros start in each row of X1
  int<lower=0, upper = rows(w_X1) + 1> u_X1[dense_X ? 0 : N[2] + 1]; 
  // declares prior_PD, has_intercept, link, prior_dist, prior_dist_for_intercept
#include /data/data_glm.stan

  int<lower=0> K_smooth;
  matrix[N[1], K_smooth] S0;
  matrix[N[2], K_smooth] S1;
  int<lower=1> smooth_map[K_smooth];
  
  int<lower=5,upper=5> family;

  // weights
  int<lower=0,upper=1> has_weights;  // 0 = No, 1 = Yes
  vector[has_weights ? N[1] : 0] weights0;
  vector[has_weights ? N[2] : 0] weights1;
  
  // offset
  int<lower=0,upper=1> has_offset;  // 0 = No, 1 = Yes
  vector[has_offset ? N[1] : 0] offset0;
  vector[has_offset ? N[2] : 0] offset1;
  
  // declares prior_{mean, scale, df}, prior_{mean, scale, df}_for_intercept, prior_{mean, scale, df}_for_aux
#include /data/hyperparameters.stan
  // declares t, p[t], l[t], q, len_theta_L, shape, scale, {len_}concentration, {len_}regularization
#include /data/glmer_stuff.stan

  // more glmer stuff
  int<lower=0> num_non_zero[2];     // number of non-zero elements in the Z matrices
  vector[num_non_zero[1]] w0;       // non-zero elements in the implicit Z0 matrix
  vector[num_non_zero[2]] w1;       // non-zero elements in the implicit Z1 matrix
  int<lower=0, upper = q - 1> v0[num_non_zero[1]]; // column indices for w0
  int<lower=0, upper = q - 1> v1[num_non_zero[2]]; // column indices for w1
  // where the non-zeros start in each row of Z0
  int<lower=0, upper = rows(w0) + 1> u0[t > 0 ? N[1] + 1 : 0];  
  // where the non-zeros start in each row of Z1
  int<lower=0, upper = rows(w1) + 1> u1[t > 0 ? N[2] + 1 : 0];  
  int<lower=0, upper=1> special_case;     // whether we only have to deal with (1|group)
}
transformed data {
  int NN = N[1] + N[2];
  real aux = not_a_number();
  int<lower=1> V0[special_case ? t : 0,N[1]] = make_V(N[1], special_case ? t : 0, v0);
  int<lower=1> V1[special_case ? t : 0,N[2]] = make_V(N[2], special_case ? t : 0, v1);
  int<lower=0> successes[clogit ? J : 0];
  int<lower=0> failures[clogit ? J : 0];
  int<lower=0> observations[clogit ? J : 0];

  int can_do_bernoullilogitglm = K != 0 &&  // remove K!=0 after rstan includes this Stan bugfix: https://github.com/stan-dev/math/issues/1398
                                 link == 1 && clogit == 0 && has_offset == 0 && 
                                 prior_PD == 0 && dense_X == 1 && has_weights == 0 && t == 0;
  matrix[can_do_bernoullilogitglm ? NN : 0, can_do_bernoullilogitglm ? K + K_smooth : 0] XS;
  int y[can_do_bernoullilogitglm ? NN : 0];

  // defines hs, len_z_T, len_var_group, delta, pos
#include /tdata/tdata_glm.stan
  for (j in 1:J) {
    successes[j] = 0;
    failures[j] = 0;
  }
  if (J > 0) for (i in 1:N[2]) successes[strata[i]] += 1;
  if (J > 0) for (i in (N[2] + 1):NN) failures[strata[i]] +=  1;
  for (j in 1:J) observations[j] = failures[j] + successes[j];

  if (can_do_bernoullilogitglm) {
    XS = K_smooth > 0 ? append_col(append_row(X0[1], X1[1]), append_row(S0, S1)) : append_row(X0[1], X1[1]);
    y = append_array(rep_array(0, N[1]), rep_array(1, N[2]));
  }
}
parameters {
  real<upper=(link == 4 ? 0.0 : positive_infinity())> gamma[has_intercept];
  // declares z_beta, global, local, z_b, z_T, rho, zeta, tau
#include /parameters/parameters_glm.stan
}
transformed parameters {
  // defines beta, b, theta_L
#include /tparameters/tparameters_glm.stan
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
  if (can_do_bernoullilogitglm) {
    vector[K + K_smooth] coeff = K_smooth > 0 ? append_row(beta, beta_smooth) : beta;
    target += bernoulli_logit_glm_lpmf(y | XS, has_intercept ? gamma[1] : 0.0, coeff);
  } else if (prior_PD == 0) {
    // defines eta0, eta1
#include /model/make_eta_bern.stan
    if (has_intercept == 1) {
      if (link != 4) {
        eta0 += gamma[1];
        eta1 += gamma[1];
      }
      else {
        real shift = fmax(max(eta0), max(eta1));
        eta0 += gamma[1] - shift;
        eta1 += gamma[1] - shift;
      }
    }
    // Log-likelihood
    if (clogit) { 
      real dummy = ll_clogit_lp(eta0, eta1, successes, failures, observations);
    }
    else if (has_weights == 0) {
      real dummy = ll_bern_lp(eta0, eta1, link, N);
    }
    else {  // weighted log-likelihoods
      target += dot_product(weights0, pw_bern(0, eta0, link));
      target += dot_product(weights1, pw_bern(1, eta1, link));
    }
  }
  
#include /model/priors_glm.stan
  if (t > 0) {
    real dummy = decov_lp(z_b, z_T, rho, zeta, tau, 
                          regularization, delta, shape, t, p);
  }
}
generated quantities {
  real mean_PPD = compute_mean_PPD ? 0 : negative_infinity();
  real alpha[has_intercept];
  
  if (has_intercept == 1) {
    if (dense_X) alpha[1] = gamma[1] - dot_product(xbar, beta);
    else alpha[1] = gamma[1];
  }
  
  if (compute_mean_PPD) {
    vector[N[1]] pi0;
    vector[N[2]] pi1;
    // defines eta0, eta1
#include /model/make_eta_bern.stan
    if (has_intercept == 1) {
      if (link != 4) {
        eta0 += gamma[1];
        eta1 += gamma[1];
      }      
      else {
        real shift;
        shift = fmax(max(eta0), max(eta1));
        eta0 += gamma[1] - shift;
        eta1 += gamma[1] - shift;
        alpha[1] -= shift;
      }
    }
    if (clogit) for (j in 1:J) mean_PPD += successes[j]; // fixed by design
    else {
      pi0 = linkinv_bern(eta0, link);
      pi1 = linkinv_bern(eta1, link);
      for (n in 1:N[1]) mean_PPD += bernoulli_rng(pi0[n]);
      for (n in 1:N[2]) mean_PPD += bernoulli_rng(pi1[n]);
    }
    mean_PPD /= NN;
  }
}
