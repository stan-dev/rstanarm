#include "license.stan" // GPL3+

// GLM for a Bernoulli outcome
functions {
  #include "common_functions.stan"
  
  /** 
   * Apply inverse link function to linear predictor
   * see help(binom) in R
   *
   * @param eta Linear predictor vector
   * @param link An integer indicating the link function
   * @return A vector, i.e. inverse-link(eta)
   */
  vector linkinv_bern(vector eta, int link) {
    if (link == 1)      return(inv_logit(eta)); // logit
    else if (link == 2) return(Phi(eta)); // probit
    else if (link == 3) return(atan(eta) / pi() + 0.5);  // cauchit
    else if (link == 4) return(exp(eta)); // log
    else if (link == 5) return(inv_cloglog(eta)); // cloglog
    else reject("Invalid link");
    return eta; // never reached
  }

  /**
   * Increment with the unweighted log-likelihood
   * @param link An integer indicating the link function
   * @param eta0 A vector of linear predictors | y = 0
   * @param eta1 A vector of linear predictors | y = 1
   * @param N An integer array of length 2 giving the number of 
   *   observations where y = 0 and y = 1 respectively
   * @return lp__
   */
  real ll_bern_lp(vector eta0, vector eta1, int link, int[] N) {
    if (link < 1 || link > 5) 
      reject("Invalid link");
      
    if (link == 1) { // logit
      target += logistic_lccdf(eta0 | 0, 1);
      target += logistic_lcdf( eta1 | 0, 1);
    }
    else if (link == 2) {  // probit
      target += normal_lccdf(eta0 | 0, 1);
      target += normal_lcdf( eta1 | 0, 1);
    }
    else if (link == 3) {  // cauchit
      target += cauchy_lccdf(eta0 | 0, 1);
      target += cauchy_lcdf( eta1 | 0, 1);
    }
    else if(link == 4) {  // log
      target += log1m_exp(eta0);
      target += eta1;  // already in log form
    }
    else if(link == 5) {  // cloglog
      target += log1m_exp(-exp(eta1));
      target += -exp(eta0);
    }
    return target();
  }

  /** 
   * Pointwise (pw) log-likelihood vector
   *
   * @param y The integer outcome variable. Note that function is
   *  called separately with y = 0 and y = 1
   * @param eta Vector of linear predictions
   * @param link An integer indicating the link function
   * @return A vector
   */
  vector pw_bern(int y, vector eta, int link) {
    int N = rows(eta);
    vector[rows(eta)] ll;
    if (link < 1 || link > 5) 
      reject("Invalid link");
      
    if (link == 1) {  // logit
      for (n in 1:N) ll[n] = bernoulli_logit_lpmf(y | eta[n]);
    }
    else {  // link = probit, cauchit, log, or cloglog 
            // Note: this may not be numerically stable
      vector[N] pi;
      pi = linkinv_bern(eta, link);
      for (n in 1:N) ll[n] = bernoulli_lpmf(y | pi[n]);
    }
    return ll;
  }
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
  
  // declares prior_{mean, scale, df}, prior_{mean, scale, df}_for_intercept, prior_scale_for dispersion
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
}
transformed data {
  int NN = N[1] + N[2];
  #include "tdata_glm.stan"// defines hs, len_z_T, len_var_group, delta, pos, t_{any, all}_124
}
parameters {
  real<upper=(link == 4 ? 0.0 : positive_infinity())> gamma[has_intercept];
  #include "parameters_glm.stan" // declares z_beta, global, local, z_b, z_T, rho, zeta, tau
}
transformed parameters {
  #include "tparameters_glm.stan" // defines beta, b, theta_L
  if (t > 0) {
    theta_L = make_theta_L(len_theta_L, p, 
                            1.0, tau, scale, zeta, rho, z_T);
    b = make_b(z_b, theta_L, p, l);
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
