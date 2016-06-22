#include "license.stan"

// GLM for a count outcome
functions {
  #include "common_functions.stan"

  vector linkinv_count(vector eta, int link) {
    vector[rows(eta)] phi;
    if (link < 1 || link > 3) 
      reject("Invalid link");
      
    if (link == 1) return exp(eta);  // log
    else if (link == 2) return eta;  // identity
    else  // link = sqrt
      for (n in 1:rows(eta)) phi[n] = square(eta[n]); 
    return phi;
  }
  
  /** 
  * Pointwise (pw) log-likelihood vector for the Poisson distribution
  *
  * @param y The integer array corresponding to the outcome variable.
  * @param link An integer indicating the link function
  * @return A vector
  */
  vector pw_pois(int[] y, vector eta, int link) {
    vector[rows(eta)] ll;
    if (link < 1 || link > 3) 
      reject("Invalid link");
      
    if (link == 1)  // log
      for (n in 1:rows(eta)) ll[n] = poisson_log_lpmf(y[n] | eta[n]);
    else {  // link = identity or sqrt
      vector[rows(eta)] phi;
      phi = linkinv_count(eta, link);
      for (n in 1:rows(eta)) ll[n] = poisson_lpmf(y[n] | phi[n]) ;
    }
    return ll;
  }
  
  /** 
  * Pointwise (pw) log-likelihood vector for the negative binomial  distribution
  *
  * @param y The integer array corresponding to the outcome variable.
  * @param link An integer indicating the link function
  * @return A vector
  */
  vector pw_nb(int[] y, vector eta, real theta, int link) {
    vector[rows(eta)] ll;
    vector[rows(eta)] rho;
    if (link < 1 || link > 3) 
      reject("Invalid link");
      
    rho = linkinv_count(eta, link);
    for (n in 1:rows(eta)) ll[n] = neg_binomial_2_lpmf(y[n] | rho[n], theta);
    return ll;
  }
}
data {
  #include "NKX.stan"
  int<lower=0> y[N];  // count outcome
  #include "data_glm.stan"
  #include "weights_offset.stan"
  #include "hyperparameters.stan"
  #include "glmer_stuff.stan"
  #include "glmer_stuff2.stan"
}
transformed data{
  real poisson_max;
  #include "tdata_glm.stan"
  poisson_max = pow(2.0, 30.0);
}
parameters {
  real<lower=if_else(link == 1, negative_infinity(), 0)> gamma[has_intercept];
  #include "parameters_glm.stan"
  real<lower=0> dispersion_unscaled[family > 1];
  vector<lower=0>[N] noise[family == 3]; // do not store this
}
transformed parameters {
  real dispersion[family > 1];
  #include "tparameters_glm.stan"
  if (family > 1 && prior_scale_for_dispersion > 0) 
    dispersion[1] = prior_scale_for_dispersion * dispersion_unscaled[1];
  else if (family > 1) dispersion[1] = dispersion_unscaled[1];
  if (t > 0) {
    if (family == 1)
      theta_L = make_theta_L(len_theta_L, p, 1.0,
                              tau, scale, zeta, rho, z_T);
    else
      theta_L = make_theta_L(len_theta_L, p, dispersion[1],
                              tau, scale, zeta, rho, z_T);
    b = make_b(z_b, theta_L, p, l);
  }
}
model {
  #include "make_eta.stan"
  if (t > 0) eta = eta + csr_matrix_times_vector(N, q, w, v, u, b);  
  if (has_intercept == 1) {
    if (link == 1) eta = eta + gamma[1];
    else eta = eta - min(eta) + gamma[1];
  }
  else {
    #include "eta_no_intercept.stan"
  }
  
  if (family == 3) {
    if      (link == 1) eta = eta + log(dispersion[1]) + log(noise[1]);
    else if (link == 2) eta = eta * dispersion[1] .* noise[1];
    else                eta = eta + sqrt(dispersion[1]) + sqrt_vec(noise[1]);
  }
  
  // Log-likelihood 
  if (has_weights == 0 && prior_PD == 0) {  // unweighted log-likelihoods
    if(family != 2) {
      if (link == 1) target += poisson_log_lpmf(y | eta);
      else target += poisson_lpmf(y | linkinv_count(eta, link));
    }
    else {
      if (link == 1) target += neg_binomial_2_log_lpmf(y | eta, dispersion[1]);
      else target += neg_binomial_2_lpmf(y | linkinv_count(eta, link), dispersion[1]);
    }
  }
  else if (family != 1 && prior_PD == 0)
    target += dot_product(weights, pw_pois(y, eta, link));
  else if (prior_PD == 0)
    target += dot_product(weights, pw_nb(y, eta, dispersion[1], link));
  
  // Log-prior for dispersion
  if (family > 1 && prior_scale_for_dispersion > 0) 
    target += cauchy_lpdf(dispersion_unscaled | 0, 1);
  
  #include "priors_glm.stan"
  
  // Log-prior for noise
  if (family == 3) target += gamma_lpdf(noise[1] | dispersion[1], 1);
  
  if (t > 0) decov_lp(z_b, z_T, rho, zeta, tau, 
                      regularization, delta, shape, t, p);
}
generated quantities {
  real alpha[has_intercept];
  real mean_PPD;
  if (has_intercept == 1) {
    if (dense_X) alpha[1] = gamma[1] - dot_product(xbar, beta);
    else alpha[1] = gamma[1];
  }
  mean_PPD = 0;
  {
    vector[N] nu;
    #include "make_eta.stan"
    if (t > 0) eta = eta + csr_matrix_times_vector(N, q, w, v, u, b);
    if (has_intercept == 1) {
      if (link == 1) eta = eta + gamma[1];
      else {
        real shift;
        shift = min(eta);
        eta = eta - shift + gamma[1];
        alpha[1] = alpha[1] - shift;
      }
    }
    else {
      #include "eta_no_intercept.stan"
    }
    
    if (family == 3) {
      if      (link == 1) eta = eta + log(dispersion[1]) + log(noise[1]);
      else if (link == 2) eta = eta * dispersion[1] .* noise[1];
      else                eta = eta + sqrt(dispersion[1]) + sqrt_vec(noise[1]);
    }
    nu = linkinv_count(eta, link);
    if (family != 2) for (n in 1:N) {
        if (nu[n] < poisson_max) mean_PPD = mean_PPD + poisson_rng(nu[n]);
        else mean_PPD = mean_PPD + normal_rng(nu[n], sqrt(nu[n]));
    }
    else for (n in 1:N) {
        real gamma_temp;
        if (is_inf(dispersion[1])) gamma_temp = nu[n];
        else gamma_temp = gamma_rng(dispersion[1], dispersion[1] / nu[n]);
        if (gamma_temp < poisson_max)
          mean_PPD = mean_PPD + poisson_rng(gamma_temp);
        else mean_PPD = mean_PPD + normal_rng(gamma_temp, sqrt(gamma_temp));
    }
    mean_PPD = mean_PPD / N;
  }
}
