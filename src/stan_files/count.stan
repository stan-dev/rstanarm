#include /pre/Columbia_copyright.stan
#include /pre/license.stan

// GLM for a count outcome
functions {
#include /functions/common_functions.stan
#include /functions/count_likelihoods.stan
}
data {
  // declares N, K, X, xbar, dense_X, nnz_x, w_x, v_x, u_x
#include /data/NKX.stan
  int<lower=0> y[N];  // count outcome
  // declares prior_PD, has_intercept, link, prior_dist, prior_dist_for_intercept
#include /data/data_glm.stan
  // declares has_weights, weights, has_offset, offset_
#include /data/weights_offset.stan
  int<lower=6,upper=7> family; // 6 poisson, 7 neg-binom, (8 poisson with gamma noise at some point?)
  // declares prior_{mean, scale, df}, prior_{mean, scale, df}_for_intercept, prior_{mean, scale, df}_for_aux
#include /data/hyperparameters.stan
  // declares t, p[t], l[t], q, len_theta_L, shape, scale, {len_}concentration, {len_}regularization
#include /data/glmer_stuff.stan
  // declares num_not_zero, w, v, u
#include /data/glmer_stuff2.stan
}
transformed data{
  real poisson_max = pow(2.0, 30.0);
  int<lower=1> V[special_case ? t : 0, N] = make_V(N, special_case ? t : 0, v);
  
  int can_do_countlogglm = K != 0 &&  // remove K!=0 after rstan includes this Stan bugfix: https://github.com/stan-dev/math/issues/1398
                           link == 1 && prior_PD == 0 && 
                           dense_X == 1 && has_weights == 0 && t == 0; 
  matrix[can_do_countlogglm ? N : 0, can_do_countlogglm ? K + K_smooth : 0] XS;
  
  // defines hs, len_z_T, len_var_group, delta, pos
#include /tdata/tdata_glm.stan

  if (can_do_countlogglm) {
    XS = K_smooth > 0 ? append_col(X[1], S) : X[1];
  }
}
parameters {
  real<lower=(link == 1 ? negative_infinity() : 0.0)> gamma[has_intercept];
  // declares z_beta, global, local, z_b, z_T, rho, zeta, tau
#include /parameters/parameters_glm.stan
  real<lower=0> aux_unscaled[family > 6];
  vector<lower=0>[N] noise[family == 8]; // do not store this
}
transformed parameters {
  real aux = negative_infinity(); // be careful with this in the family = 6 case
  // defines beta, b, theta_L
#include /tparameters/tparameters_glm.stan

  if (family > 6 && (prior_dist_for_aux == 0 || prior_scale_for_aux <= 0))
    aux = aux_unscaled[1];
  else if (family > 6) {
    aux = prior_scale_for_aux * aux_unscaled[1];
    if (prior_dist_for_aux <= 2) // normal or student_t
      aux += prior_mean_for_aux;
  }
 
  if (t > 0) {
    if (special_case == 1) {
      int start = 1;
      theta_L = scale .* (family == 6 ? tau : tau * aux);
      if (t == 1) b = theta_L[1] * z_b;
      else for (i in 1:t) {
        int end = start + l[i] - 1;
        b[start:end] = theta_L[i] * z_b[start:end];
        start = end + 1;
      }
    }
    else {
      if (family == 6)
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
  if (can_do_countlogglm) {
    vector[K + K_smooth] coeff = K_smooth > 0 ? append_row(beta, beta_smooth) : beta;
    if (family != 7) {
      if (has_offset) {
        target += poisson_log_glm_lpmf(y | XS, has_intercept ? offset_ + gamma[1] : offset_, coeff);
      } else {
        target += poisson_log_glm_lpmf(y | XS, has_intercept ? gamma[1] : 0.0, coeff);
      }
    } else {
      if (has_offset) {
        target += neg_binomial_2_log_glm_lpmf(y | XS, has_intercept ? offset_ + gamma[1] : offset_, coeff, aux);
      } else {
        target += neg_binomial_2_log_glm_lpmf(y | XS, has_intercept ? gamma[1] : 0.0, coeff, aux);
      }
    }
  } else if (prior_PD == 0) {
#include /model/make_eta.stan
    if (t > 0) {
#include /model/eta_add_Zb.stan
    }
    if (has_intercept == 1) {
      if (link == 1) eta += gamma[1];
      else eta += gamma[1] - min(eta);
    }
    else {
#include /model/eta_no_intercept.stan
    }
  
    if (family == 8) {
      if      (link == 1) eta += log(aux) + log(noise[1]);
      else if (link == 2) {
        eta *= aux;
        eta .*= noise[1];
      }
      else                eta += sqrt(aux) + sqrt(noise[1]);
    }
  
    // Log-likelihood
    if (has_weights == 0) {  // unweighted log-likelihoods
      if (family != 7) {
        if (link == 1) target += poisson_log_lpmf(y | eta);
        else target += poisson_lpmf(y | linkinv_count(eta, link));
      }
      else {
        if (link == 1) target += neg_binomial_2_log_lpmf(y | eta, aux);
        else target += neg_binomial_2_lpmf(y | linkinv_count(eta, link), aux);
      }
    }
    else if (family != 7)
      target += dot_product(weights, pw_pois(y, eta, link));
    else
      target += dot_product(weights, pw_nb(y, eta, aux, link));
  }
  
  // Log-prior for aux
  if (family > 6 && 
      prior_dist_for_aux > 0 && prior_scale_for_aux > 0) {
    real log_half = -0.693147180559945286;    
    if (prior_dist_for_aux == 1)
      target += normal_lpdf(aux_unscaled | 0, 1) - log_half;
    else if (prior_dist_for_aux == 2)
      target += student_t_lpdf(aux_unscaled | prior_df_for_aux, 0, 1) - log_half;
    else 
      target += exponential_lpdf(aux_unscaled | 1);
  }
  
#include /model/priors_glm.stan
  
  // Log-prior for noise
  if (family == 8) target += gamma_lpdf(noise[1] | aux, 1);
  
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
    vector[N] nu;
#include /model/make_eta.stan
    if (t > 0) {
#include /model/eta_add_Zb.stan
    }
    if (has_intercept == 1) {
      if (link == 1) eta += gamma[1];
      else {
        real shift = min(eta);
        eta += gamma[1] - shift;
        alpha[1] -= shift;
      }
    }
    else {
#include /model/eta_no_intercept.stan
    }

    if (family == 8) {
      if      (link == 1) eta += log(aux) + log(noise[1]);
      else if (link == 2) {
        eta *= aux;
        eta .*= noise[1];
      }
      else                eta += sqrt(aux) + sqrt(noise[1]);
    }
    nu = linkinv_count(eta, link);
    if (family != 7) for (n in 1:N) {
        if (nu[n] < poisson_max) mean_PPD += poisson_rng(nu[n]);
        else mean_PPD += normal_rng(nu[n], sqrt(nu[n]));
    }
    else for (n in 1:N) {
        real gamma_temp;
        if (is_inf(aux)) gamma_temp = nu[n];
        else gamma_temp = gamma_rng(aux, aux / nu[n]);
        if (gamma_temp < poisson_max)
          mean_PPD += poisson_rng(gamma_temp);
        else mean_PPD += normal_rng(gamma_temp, sqrt(gamma_temp));
    }
    mean_PPD /= N;
  }
}
