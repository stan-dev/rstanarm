#include "Columbia_copyright.stan"
#include "license.stan" // GPL3+

// GLM for a count outcome
functions {
  #include "common_functions.stan"
  #include "count_likelihoods.stan"
}
data {
  #include "NKX.stan"      // declares N, K, X, xbar, dense_X, nnz_x, w_x, v_x, u_x
  int<lower=0> y[N];  // count outcome
  #include "data_glm.stan" // declares prior_PD, has_intercept, family, link, prior_dist, prior_dist_for_intercept
  #include "weights_offset.stan"  // declares has_weights, weights, has_offset, offset
  // declares prior_{mean, scale, df}, prior_{mean, scale, df}_for_intercept, prior_{mean, scale, df}_for_aux
  #include "hyperparameters.stan"
  // declares t, p[t], l[t], q, len_theta_L, shape, scale, {len_}concentration, {len_}regularization
  #include "glmer_stuff.stan"  
  #include "glmer_stuff2.stan" // declares num_not_zero, w, v, u
}
transformed data{
  real poisson_max = pow(2.0, 30.0);
  int<lower=1> V[special_case ? t : 0, N] = make_V(N, special_case ? t : 0, v);
  #include "tdata_glm.stan"// defines hs, len_z_T, len_var_group, delta, pos
}
parameters {
  real<lower=(link == 1 ? negative_infinity() : 0.0)> gamma[has_intercept];
  #include "parameters_glm.stan" // declares z_beta, global, local, z_b, z_T, rho, zeta, tau
  real<lower=0> aux_unscaled[family > 1];
  vector<lower=0>[N] noise[family == 3]; // do not store this
}
transformed parameters {
  real aux = negative_infinity(); // be careful with this in the family = 1 case
  #include "tparameters_glm.stan" // defines beta, b, theta_L
 
  if (family > 1 && (prior_dist_for_aux == 0 || prior_scale_for_aux <= 0))
    aux = aux_unscaled[1];
  else if (family > 1) {
    aux = prior_scale_for_aux * aux_unscaled[1];
    if (prior_dist_for_aux <= 2) // normal or student_t
      aux = aux + prior_mean_for_aux;
  }
 
  if (t > 0) {
    if (special_case == 1) {
      int start = 1;
      theta_L = scale .* (family == 1 ? tau : tau * aux);
      if (t == 1) b = theta_L[1] * z_b;
      else for (i in 1:t) {
        int end = start + l[i] - 1;
        b[start:end] = theta_L[i] * z_b[start:end];
        start = end + 1;
      }
    }
    else {
      if (family == 1)
        theta_L = make_theta_L(len_theta_L, p, 1.0,
                               tau, scale, zeta, rho, z_T);
      else
        theta_L = make_theta_L(len_theta_L, p, aux,
                               tau, scale, zeta, rho, z_T);
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
    if (link == 1) eta = eta + gamma[1];
    else eta = eta - min(eta) + gamma[1];
  }
  else {
    #include "eta_no_intercept.stan" // shifts eta
  }
  
  if (family == 3) {
    if      (link == 1) eta = eta + log(aux) + log(noise[1]);
    else if (link == 2) eta = eta * aux .* noise[1];
    else                eta = eta + sqrt(aux) + sqrt(noise[1]);
  }
  
  // Log-likelihood 
  if (has_weights == 0 && prior_PD == 0) {  // unweighted log-likelihoods
    if(family != 2) {
      if (link == 1) target += poisson_log_lpmf(y | eta);
      else target += poisson_lpmf(y | linkinv_count(eta, link));
    }
    else {
      if (link == 1) target += neg_binomial_2_log_lpmf(y | eta, aux);
      else target += neg_binomial_2_lpmf(y | linkinv_count(eta, link), aux);
    }
  }
  else if (family != 2 && prior_PD == 0)
    target += dot_product(weights, pw_pois(y, eta, link));
  else if (prior_PD == 0)
    target += dot_product(weights, pw_nb(y, eta, aux, link));
  
  // Log-prior for aux
  if (family > 1 && 
      prior_dist_for_aux > 0 && prior_scale_for_aux > 0) {
    if (prior_dist_for_aux == 1)
      target += normal_lpdf(aux_unscaled | 0, 1);
    else if (prior_dist_for_aux == 2)
      target += student_t_lpdf(aux_unscaled | prior_df_for_aux, 0, 1);
    else 
      target += exponential_lpdf(aux_unscaled | 1);
  }
  
  #include "priors_glm.stan" // increments target()
  
  // Log-prior for noise
  if (family == 3) target += gamma_lpdf(noise[1] | aux, 1);
  
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
    vector[N] nu;
    #include "make_eta.stan" // defines eta
    if (t > 0) {
      #include "eta_add_Zb.stan"
    }
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
      #include "eta_no_intercept.stan" // shifts eta
    }
    
    if (family == 3) {
      if      (link == 1) eta = eta + log(aux) + log(noise[1]);
      else if (link == 2) eta = eta * aux .* noise[1];
      else                eta = eta + sqrt(aux) + sqrt(noise[1]);
    }
    nu = linkinv_count(eta, link);
    if (family != 2) for (n in 1:N) {
        if (nu[n] < poisson_max) mean_PPD = mean_PPD + poisson_rng(nu[n]);
        else mean_PPD = mean_PPD + normal_rng(nu[n], sqrt(nu[n]));
    }
    else for (n in 1:N) {
        real gamma_temp;
        if (is_inf(aux)) gamma_temp = nu[n];
        else gamma_temp = gamma_rng(aux, aux / nu[n]);
        if (gamma_temp < poisson_max)
          mean_PPD = mean_PPD + poisson_rng(gamma_temp);
        else mean_PPD = mean_PPD + normal_rng(gamma_temp, sqrt(gamma_temp));
    }
    mean_PPD = mean_PPD / N;
  }
}
