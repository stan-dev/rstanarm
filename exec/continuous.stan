#include "license.stan"

# GLM for a Gaussian, Gamma, or inverse Gaussian outcome
functions {
  #include "common_functions.stan"

  /** 
   * Apply inverse link function to linear predictor
   *
   * @param eta Linear predictor vector
   * @param link An integer indicating the link function
   * @return A vector, i.e. inverse-link(eta)
   */
  vector linkinv_gauss(vector eta, int link) {
    if (link < 1 || link > 3) reject("Invalid link");
    if (link < 3)  # link = identity or log 
      return eta; # return eta for log link too bc will use lognormal
    else {# link = inverse
      vector[rows(eta)] mu;
      for(n in 1:rows(eta)) mu[n] = inv(eta[n]); 
      return mu;
    }
  }
  
  /** 
  * Apply inverse link function to linear predictor
  *
  * @param eta Linear predictor vector
  * @param link An integer indicating the link function
  * @return A vector, i.e. inverse-link(eta)
  */
  vector linkinv_gamma(vector eta, int link) {
    if (link < 1 || link > 3) reject("Invalid link");
    if (link == 1)  return eta;
    else if (link == 2) return exp(eta);
    else {
      vector[rows(eta)] mu;
      for(n in 1:rows(eta)) mu[n] = inv(eta[n]); 
      return mu;
    }
  }
  
  /** 
  * Apply inverse link function to linear predictor
  *
  * @param eta Linear predictor vector
  * @param link An integer indicating the link function
  * @return A vector, i.e. inverse-link(eta)
  */
  vector linkinv_inv_gaussian(vector eta, int link) {
    if (link < 1 || link > 4) reject("Invalid link");
    if (link == 1)  return eta;
    else if (link == 2) return exp(eta);
    else {
      vector[rows(eta)] mu;
      if (link == 3) for( n in 1:rows(eta)) mu[n] = inv(eta[n]);
      else for (n in 1:rows(eta)) mu[n] = inv_sqrt(eta[n]);      
      return mu;
    }
  }
  
  /** 
  * Pointwise (pw) log-likelihood vector
  *
  * @param y The integer array corresponding to the outcome variable.
  * @param link An integer indicating the link function
  * @return A vector
  */
  vector pw_gauss(vector y, vector eta, real sigma, int link) {
    vector[rows(eta)] ll;
    if (link < 1 || link > 3) reject("Invalid link");
    if (link == 2) # link = log
      for (n in 1:rows(eta)) ll[n] = lognormal_log(y[n], eta[n], sigma);
    else { # link = idenity or inverse
      vector[rows(eta)] mu;
      mu = linkinv_gauss(eta, link);
      for (n in 1:rows(eta)) ll[n] = normal_log(y[n], mu[n], sigma);
    }
    return ll;
  }
  
  real GammaReg_log(vector y, vector eta, real shape, 
                    int link, real sum_log_y) {
    real ret;
    if (link < 1 || link > 3) reject("Invalid link");
    ret = rows(y) * (shape * log(shape) - lgamma(shape)) +
      (shape - 1) * sum_log_y;
    if (link == 2)      # link is log
      ret = ret - shape * sum(eta) - shape * sum(y ./ exp(eta));
    else if (link == 1) # link is identity
      ret = ret - shape * sum(log(eta)) - shape * sum(y ./ eta);
    else                # link is inverse
      ret = ret + shape * sum(log(eta)) - shape * dot_product(eta, y);
    return ret;
  }
  
  /** 
  * Pointwise (pw) log-likelihood vector
  *
  * @param y The integer array corresponding to the outcome variable.
  * @param link An integer indicating the link function
  * @return A vector
  */
  vector pw_gamma(vector y, vector eta, real shape, int link) {
    vector[rows(eta)] ll;
    if (link < 1 || link > 3) reject("Invalid link");
    if (link == 3) { # link = inverse
      for (n in 1:rows(eta)) {
        ll[n] = gamma_log(y[n], shape, shape * eta[n]);
      }
    }
    else if (link == 2) { # link = log
      for (n in 1:rows(eta)) {
        ll[n] = gamma_log(y[n], shape, shape / exp(eta[n]));
      }
    }
    else { # link = identity
      for (n in 1:rows(eta)) {
        ll[n] = gamma_log(y[n], shape, shape / eta[n]);
      }
    }
    return ll;
  }
  
  /** 
  * inverse Gaussian log-PDF (for data only, excludes constants)
  *
  * @param y The vector of outcomes
  * @param eta The vector of linear predictors
  * @param lambda A positive scalar nuisance parameter
  * @param link An integer indicating the link function
  * @return A scalar
  */
  real inv_gaussian_log(vector y, vector mu, real lambda, 
                        real sum_log_y, vector sqrt_y) {
    return 0.5 * rows(y) * log(lambda / (2 * pi())) - 
      1.5 * sum_log_y - 
      0.5 * lambda * dot_self( (y - mu) ./ (mu .* sqrt_y) );
  }
  
  /** 
  * Pointwise (pw) log-likelihood vector
  *
  * @param y The integer array corresponding to the outcome variable.
  * @param eta The linear predictors
  * @param lamba A positive scalar nuisance parameter
  * @param link An integer indicating the link function
  * @param log_y A precalculated vector of the log of y
  * @param sqrt_y A precalculated vector of the square root of y
  * @return A vector of log-likelihoods
  */
  vector pw_inv_gaussian(vector y, vector eta, real lambda, 
                         int link, vector log_y, vector sqrt_y) {
    vector[rows(y)] ll;
    vector[rows(y)] mu;
    if (link < 1 || link > 4) reject("Invalid link");
    mu = linkinv_inv_gaussian(eta, link);
    for (n in 1:rows(y))
      ll[n] = -0.5 * lambda * square( (y[n] - mu[n]) / (mu[n] * sqrt_y[n]) );
    ll = ll + 0.5 * log(lambda / (2 * pi())) - 1.5 * log_y;
    return ll;
  }
  
  /** 
  * PRNG for the inverse Gaussian distribution
  *
  * Algorithm from wikipedia 
  *
  * @param mu The expectation
  * @param lambda The dispersion
  * @return A draw from the inverse Gaussian distribution
  */
  real inv_gaussian_rng(real mu, real lambda) {
    real z;
    real y;
    real x;
    real mu2;
    mu2 = square(mu);
    y = square(normal_rng(0,1));
    z = uniform_rng(0,1);
    x = mu + ( mu2 * y - mu * sqrt(4 * mu * lambda * y + mu2 * square(y)) )
      / (2 * lambda);
    if (z <= (mu / (mu + x))) return x;
    else return mu2 / x;
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
  #include "NKX.stan"
  vector[N] y; // continuous outcome
  #include "data_glm.stan"
  #include "weights_offset.stan"
  #include "hyperparameters.stan"
  #include "glmer_stuff.stan"
  #include "glmer_stuff2.stan"
}
transformed data {
  vector[N * (family == 3)] sqrt_y;
  vector[N * (family == 3)] log_y;
  real sum_log_y;
  #include "tdata_glm.stan"
  if      (family == 1) sum_log_y = not_a_number();
  else if (family == 2) sum_log_y = sum(log(y));
  else {
    for (n in 1:N) sqrt_y[n] = sqrt(y[n]);
    log_y = log(y);
    sum_log_y = sum(log_y);
  }
}
parameters {
  real<lower=if_else(family == 1 || link == 2, 
                     negative_infinity(), 0)> gamma[has_intercept];
  #include "parameters_glm.stan"
  real<lower=0> dispersion_unscaled; # interpretation depends on family!
}
transformed parameters {
  real dispersion;
  #include "tparameters_glm.stan"
  if (prior_scale_for_dispersion > 0)
    dispersion =  prior_scale_for_dispersion * dispersion_unscaled;
  else dispersion = dispersion_unscaled;
  if (t > 0) {
    theta_L = make_theta_L(len_theta_L, p, 
                            dispersion, tau, scale, zeta, rho, z_T);
    b = make_b(z_b, theta_L, p, l);
  }
}
model {
  #include "make_eta.stan"
  if (t > 0) eta = eta + csr_matrix_times_vector(N, q, w, v, u, b);
  if (has_intercept == 1) {
    if (family == 1 || link == 2) eta = eta + gamma[1];
    else eta = eta - min(eta) + gamma[1];
  }
  else {
    #include "eta_no_intercept.stan"
  }
  
  // Log-likelihood 
  if (has_weights == 0 && prior_PD == 0) { # unweighted log-likelihoods
    if (family == 1) {
      if (link == 1)      y ~ normal(eta, dispersion);
      else if (link == 2) y ~ lognormal(eta, dispersion);
      else y ~ normal(divide_real_by_vector(1, eta), dispersion);
      // divide_real_by_vector() is defined in common_functions.stan
    }
    else if (family == 2) {
      y ~ GammaReg(eta, dispersion, link, sum_log_y);
    }
    else {
      y ~ inv_gaussian(linkinv_inv_gaussian(eta, link), 
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

  // Log-prior for scale
  if (prior_scale_for_dispersion > 0) dispersion_unscaled ~ cauchy(0, 1);
  #include "priors_glm.stan"
  if (t > 0) decov_lp(z_b, z_T, rho, zeta, tau, 
                      regularization, delta, shape, t, p);
}
generated quantities {
  real alpha[has_intercept];
  real mean_PPD;
  mean_PPD = 0;
  if (has_intercept == 1)
    alpha[1] = gamma[1] - dot_product(xbar, beta);
  {
    #include "make_eta.stan"
    if (t > 0) eta = eta + csr_matrix_times_vector(N, q, w, v, u, b);
    if (has_intercept == 1) {
      if (family == 1 || link == 2) eta = eta + gamma[1];
      else {
        real min_eta;
        min_eta = min(eta);
        alpha[1] = alpha[1] - min_eta;
        eta = eta - min_eta + gamma[1];
      }
    }
    else {
      #include "eta_no_intercept.stan"
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
