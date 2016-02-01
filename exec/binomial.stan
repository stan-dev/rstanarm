#include "license.stan"

// GLM for a binomial outcome
functions {
  #include "common_functions.stan"
  
  /** 
   * Apply inverse link function to linear predictor
   *
   * @param eta Linear predictor vector
   * @param link An integer indicating the link function
   * @return A vector, i.e. inverse-link(eta)
   */
  vector linkinv_binom(vector eta, int link) {
    vector[rows(eta)] pi;
    if (link < 1 || link > 5) 
      reject("Invalid link");
      
    if (link == 1)  // logit
      for(n in 1:rows(eta)) pi[n] <- inv_logit(eta[n]);
    else if (link == 2)  // probit
      for(n in 1:rows(eta)) pi[n] <- Phi(eta[n]);
    else if (link == 3)  // cauchit
      for(n in 1:rows(eta)) pi[n] <- cauchy_cdf(eta[n], 0.0, 1.0);
    else if (link == 4)  // log 
      for(n in 1:rows(eta)) pi[n] <- exp(eta[n]);
    else if (link == 5)  // cloglog
      for(n in 1:rows(eta)) pi[n] <- inv_cloglog(eta[n]);
    return pi;
  }
  
  /**
  * Increment with the unweighted log-likelihood
  * @param y An integer array indicating the number of successes
  * @param trials An integer array indicating the number of trials
  * @param eta A vector of linear predictors
  * @param link An integer indicating the link function
  * @return lp__
  */
  real ll_binom_lp(int[] y, int[] trials, vector eta, int link) {
    if (link < 1 || link > 5) 
      reject("Invalid link");
      
    if (link == 1) y ~ binomial_logit(trials, eta);
    else if (link <  4) y ~ binomial(trials, linkinv_binom(eta, link));
    else if (link == 4) {  // log
      for (n in 1:num_elements(y)) {
        increment_log_prob(y[n] * eta[n]);
        increment_log_prob( (trials[n] - y[n]) * log1m_exp(eta[n]) );
      }
    }
    else if (link == 5) {  // cloglog
      real neg_exp_eta;
      for (n in 1:num_elements(y)) {
        neg_exp_eta <- -exp(eta[n]);
        increment_log_prob(y[n] * log1m_exp(neg_exp_eta));
        increment_log_prob( (trials[n] - y[n]) * neg_exp_eta );
      }
    }
    return get_lp();
  }
  
  /** 
  * Pointwise (pw) log-likelihood vector
  *
  * @param y The integer array corresponding to the outcome variable.
  * @param link An integer indicating the link function
  * @return A vector
  */
  vector pw_binom(int[] y, int[] trials, vector eta, int link) {
    vector[rows(eta)] ll;
    if (link < 1 || link > 5) reject("Invalid link");
    if (link == 1) {  // logit
      for (n in 1:rows(eta)) 
        ll[n] <- binomial_logit_log(y[n], trials[n], eta[n]);
    }
    else {  // link = probit, cauchit, log, or cloglog (unstable)
      vector[rows(eta)] pi;
      pi <- linkinv_binom(eta, link);
      for (n in 1:rows(eta)) ll[n] <- binomial_log(y[n], trials[n], pi[n]) ;
    }
    return ll;
  }
}
data {
  #include "NKX.stan"
  int<lower=0> y[N];       // outcome: number of successes
  int<lower=0> trials[N];  // number of trials
  #include "data_glm.stan"
  #include "weights_offset.stan"
  #include "hyperparameters.stan"
  #include "glmer_stuff.stan"
  #include "glmer_stuff2.stan"
}
transformed data {
  #include "tdata_glm.stan"
}
parameters {
  real<upper=if_else(link == 4, 0, positive_infinity())> gamma[has_intercept];
  #include "parameters_glm.stan"
}
transformed parameters {
  #include "tparameters_glm.stan"
  if (t > 0) {
    theta_L <- make_theta_L(len_theta_L, p, 
                            1.0, tau, scale, zeta, rho, z_T);
    b <- make_b(z_b, theta_L, p, l);
  }
}
model {
  #include "make_eta.stan"
  if (t > 0) eta <- eta + csr_matrix_times_vector(N, q, w, v, u, b);
  if (has_intercept == 1) {
    if (link != 4) eta <- eta + gamma[1];
    else eta <- gamma[1] + eta - max(eta);
  }
  else {
    #include "eta_no_intercept.stan"
  }
  
  // Log-likelihood 
  if (has_weights == 0 && prior_PD == 0) {  // unweighted log-likelihoods
    real dummy;  // irrelevant but useful for testing
    dummy <- ll_binom_lp(y, trials, eta, link);
  }
  else if (prior_PD == 0) 
    increment_log_prob(dot_product(weights, pw_binom(y, trials, eta, link)));
  
  #include "priors_glm.stan"
  
  if (t > 0) decov_lp(z_b, z_T, rho, zeta, tau, 
                      regularization, delta, shape, t, p);
}
generated quantities {
  real alpha[has_intercept];
  real mean_PPD;
  if (has_intercept == 1) alpha[1] <- gamma[1] - dot_product(xbar, beta);
  mean_PPD <- 0;
  {
    vector[N] pi;
    #include "make_eta.stan"
    if (t > 0) eta <- eta + csr_matrix_times_vector(N, q, w, v, u, b);
    if (has_intercept == 1) {
      if (link != 4) eta <- eta + gamma[1];
      else {
        real shift;
        shift <- max(eta);
        eta <- gamma[1] + eta - shift;
        alpha[1] <- alpha[1] - shift;
      }
    }
    else {
      #include "eta_no_intercept.stan"
    }
    
    pi <- linkinv_binom(eta, link);
    for (n in 1:N) mean_PPD <- mean_PPD + binomial_rng(trials[n], pi[n]);
    mean_PPD <- mean_PPD / N;
  }
}
