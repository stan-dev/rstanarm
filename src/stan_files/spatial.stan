#include /pre/Columbia_copyright.stan
#include /pre/license.stan
// CAR SPATIAL MODELS
functions {
#include /functions/continuous_likelihoods.stan
#include /functions/binomial_likelihoods.stan
#include /functions/count_likelihoods.stan
#include /functions/common_functions.stan
}
data {
  int<lower=0> N;                 // number of regions
  int<lower=0> K;                 // number of predictors (inc intercept)
  matrix[N,K] X;                  // model matrix
  int<lower=1,upper=7> family;    // 1 gaussian; 2 gamma; 5 binomial; 6 poisson; 7 neg_binomial_2
  int<lower=1,upper=5> link;
  int<lower=0,upper=1> is_continuous;
  int<lower=0,upper=1> has_aux;
  int<lower=1,upper=3> model_type;       // Besag = 1; BYM = 2; BYM2 = 3
  int<lower=0,upper=1> has_intercept;
  vector[K] xbar;
  int<lower=0> trials[family == 5 ? N : 0];
  int y_int[is_continuous == 1 ? 0 : N];
  vector[is_continuous == 1 ? N : 0] y_real;
  // pairwise difference version of CAR
  int E_n;                        // number of adjacency pairs
  int edges[E_n, 2];              // adjacency pairs
  // matrix version of CAR
  int<lower=1,upper=2> order;
  int<lower=0> Q_n[order == 2];
  vector[order == 2 ? Q_n[1] : 0] w;
  int v[order == 2 ? Q_n[1] : 0];
  int u[order == 2 ? N+1 : 0];
  // prior stuff
  int<lower=0,upper=1> prior_dist_rho;
  real prior_mean_rho;
  real<lower=0> prior_scale_rho;
  real<lower=0> prior_df_rho;
  real<lower=0> shape1_rho;
  real<lower=0> shape2_rho;
  real scaling_factor;
  int<lower=0> prior_dist_for_intercept;
  int<lower=0> prior_dist;
  int<lower=0> prior_dist_tau;
  real prior_mean_for_intercept;
  real<lower=0> prior_scale_for_intercept;
  real<lower=0> prior_df_for_intercept;
  vector[K] prior_mean;
  vector<lower=0>[K] prior_scale;
  vector<lower=0>[K] prior_df;
  real prior_mean_tau;
  real<lower=0> prior_scale_tau;
  real<lower=0> prior_df_tau;
  real<lower=0> global_prior_df;
  real<lower=0> global_prior_scale;
  int<lower=2> num_normals[prior_dist == 7 ? K : 0];
  int<lower=0,upper=3> prior_dist_for_aux;
  real<lower=0> prior_mean_for_aux;
  real<lower=0> prior_scale_for_aux;
  real<lower=0> prior_df_for_aux;
  real<lower=0> slab_scale;  // for hs prior only
}
transformed data {
  real poisson_max = pow(2.0, 30.0);
  int<lower=0> hs;
  real sum_log_y = family == 1 ? not_a_number() : sum(log(y_real));
  if (prior_dist <= 2) hs = 0;
  else if (prior_dist == 3) hs = 2;
  else if (prior_dist == 4) hs = 4;
  else hs = 0;
}
parameters {
  real<lower=make_lower(family, link), upper=make_upper(family,link)> gamma[has_intercept];
  // real<lower=negative_infinity(), upper=0> gamma[has_intercept];
  vector[K] z_beta;
  vector[model_type == 1? 0 : N] theta_raw;
  vector[N-1] phi_raw;
  // interpretation of rho and tau depends on model_type!
  real<lower=0,upper=(model_type == 3 ? 1 : positive_infinity())> rho[model_type != 1];
  real<lower=0> tau;
  real<lower=0> global[hs];
  vector<lower=0>[K] local[hs];
  vector<lower=0>[K] mix[prior_dist == 5 || prior_dist == 6];
  real<lower=0> aux_unscaled[has_aux]; // interpretation depends on family!
  real<lower=0> one_over_lambda[prior_dist == 6];
  real<lower=0> caux[hs > 0];
}
transformed parameters {
  vector[K] beta;             // predictors on covariates (including intercept)
  // aux has to be defined first in the hs case
  real aux = has_aux == 0 ? 0 : (prior_dist_for_aux == 0 ? aux_unscaled[1] : (prior_dist_for_aux <= 2 ?
             prior_scale_for_aux * aux_unscaled[1] + prior_mean_for_aux :
             prior_scale_for_aux * aux_unscaled[1]));
  vector[N] phi;          // non-centered random effect (spatial)
  vector[N] psi;
  phi[1:(N - 1)] = phi_raw;
  phi[N] = -sum(phi_raw);
  /* precision form
  if (model_type == 1)
    psi = phi * sqrt(inv(tau));
  else if (model_type == 2)
    psi = phi * sqrt(inv(rho[1])) + theta_raw * sqrt(inv(tau));
  else if (model_type == 3)
    psi = tau*(sqrt(1-rho[1])*theta_raw + sqrt(rho[1]/scaling_factor)*phi);
  */
  if (model_type == 1)
     psi = phi * tau;
  else if (model_type == 2)
     psi = phi * rho[1] + theta_raw * tau;
  else if (model_type == 3)
     psi = tau*(sqrt(1-rho[1])*theta_raw + sqrt(rho[1]/scaling_factor)*phi);
  // for regression coefficients
  // "tparameters.stan"
  if      (prior_dist == 0) beta = z_beta;
  else if (prior_dist == 1) beta = z_beta .* prior_scale + prior_mean;
  else if (prior_dist == 2) for (k in 1:K) {
    beta[k] = CFt(z_beta[k], prior_df[k]) * prior_scale[k] + prior_mean[k];
  }
  else if (prior_dist == 3) {
    real c2 = square(slab_scale) * caux[1];
    if (is_continuous == 1 && family == 1)
      beta = hs_prior(z_beta, global, local, global_prior_scale, aux, c2);
    else beta = hs_prior(z_beta, global, local, global_prior_scale, 1, c2);
  }
  else if (prior_dist == 4) {
    real c2 = square(slab_scale) * caux[1];
    if (is_continuous == 1 && family == 1)
      beta = hsplus_prior(z_beta, global, local, global_prior_scale, aux, c2);
    else beta = hsplus_prior(z_beta, global, local, global_prior_scale, 1, c2);
  }
  else if (prior_dist == 5) // laplace
    beta = prior_mean + prior_scale .* sqrt(2 * mix[1]) .* z_beta;
  else if (prior_dist == 6) // lasso
    beta = prior_mean + one_over_lambda[1] * prior_scale .* sqrt(2 * mix[1]) .* z_beta;
  else if (prior_dist == 7) { // product_normal
    int z_pos = 1;
    for (k in 1:K) {
      beta[k] = z_beta[z_pos];
      z_pos = z_pos + 1;
      for (n in 2:num_normals[k]) {
        beta[k] = beta[k] * z_beta[z_pos];
        z_pos = z_pos + 1;
      }
      beta[k] = beta[k] * prior_scale[k] ^ num_normals[k] + prior_mean[k];
    }
  }
}
model {
  vector[N] eta;   // linear predictor + spatial random effects
  // deal with intercept
  if (has_intercept == 1) {
    eta = X * beta + psi;
    if ((family == 5 && link == 4))  // binomial
      eta -= max(eta);
    else if ((family == 6 && link != 1) ||
      (family == 7 && link != 1) ||
      (family == 2 && (link == 1 || link == 3))) // poisson, neg_binomial_2, and gamma
      eta -= min(eta);
    eta += gamma[1];
  }
  else
    eta = X * beta + psi;
  // likelihoods
  if (family == 1) {
    eta = linkinv_gauss(eta, link);
    target+= normal_lpdf(y_real | eta, aux);
  }
  else if (family == 6) {
    eta = linkinv_count(eta, link);
    target+= poisson_lpmf(y_int | eta);
  }
  else if (family == 7) {
    eta = linkinv_count(eta, link);
    target += neg_binomial_2_lpmf(y_int | eta, aux);
  }
  else if (family == 5) {
    eta = linkinv_binom(eta, link);
    target+= binomial_lpmf(y_int | trials, eta);
  }
  else if (family == 2) {
    target += GammaReg(y_real, eta, aux, link, sum_log_y);
  }
  // prior on spatial parameter vector (GMRF)
  if (order == 1)
    target += -0.5 * dot_self(phi[edges[,1]] - phi[edges[,2]]);
  else if (order == 2)
    target+= -0.5 * dot_product(phi, csr_matrix_times_vector(N, N, w, v, u, phi));
  // priors on coefficients
  // "priors.stan"
  // Log-priors for coefficients
       if (prior_dist == 1) target += normal_lpdf(z_beta | 0, 1);
  else if (prior_dist == 2) target += normal_lpdf(z_beta | 0, 1); // Student t
  else if (prior_dist == 3) { // hs
    real log_half = -0.693147180559945286;
    target += normal_lpdf(z_beta | 0, 1);
    target += normal_lpdf(local[1] | 0, 1) - log_half;
    target += inv_gamma_lpdf(local[2] | 0.5 * prior_df, 0.5 * prior_df);
    target += normal_lpdf(global[1] | 0, 1) - log_half;
    target += inv_gamma_lpdf(global[2] | 0.5 * global_prior_df, 0.5 * global_prior_df);
  }
  else if (prior_dist == 4) { // hs+
    real log_half = -0.693147180559945286;
    target += normal_lpdf(z_beta | 0, 1);
    target += normal_lpdf(local[1] | 0, 1) - log_half;
    target += inv_gamma_lpdf(local[2] | 0.5 * prior_df, 0.5 * prior_df);
    target += normal_lpdf(local[3] | 0, 1) - log_half;
    // unorthodox useage of prior_scale as another df hyperparameter
    target += inv_gamma_lpdf(local[4] | 0.5 * prior_scale, 0.5 * prior_scale);
    target += normal_lpdf(global[1] | 0, 1) - log_half;
    target += inv_gamma_lpdf(global[2] | 0.5 * global_prior_df, 0.5 * global_prior_df);
  }
  else if (prior_dist == 5) { // laplace
    target += normal_lpdf(z_beta | 0, 1);
    target += exponential_lpdf(mix[1] | 1);
  }
  else if (prior_dist == 6) { // lasso
    target += normal_lpdf(z_beta | 0, 1);
    target += exponential_lpdf(mix[1] | 1);
    target += chi_square_lpdf(one_over_lambda[1] | prior_df[1]);
  }
  else if (prior_dist == 7) { // product_normal
    target += normal_lpdf(z_beta | 0, 1);
  }
  /* else prior_dist is 0 and nothing is added */

  // Log-prior for intercept
  if (has_intercept == 1) {
    if (prior_dist_for_intercept == 1)  // normal
      target += normal_lpdf(gamma | prior_mean_for_intercept, prior_scale_for_intercept);
    else if (prior_dist_for_intercept == 2)  // student_t
      target += student_t_lpdf(gamma | prior_df_for_intercept, prior_mean_for_intercept,
                               prior_scale_for_intercept);
    /* else prior_dist is 0 and nothing is added */
  }
  // model specific priors
  if (model_type == 2) {
    target+= normal_lpdf(theta_raw | 0, 1);  // unstructured (random) effect
    // prior on overall spatial variation
    if (prior_dist_rho == 1)
      target+= normal_lpdf(rho[1] | prior_mean_rho, prior_scale_rho);
    else if (prior_dist_rho == 2)
      target+= student_t_lpdf(rho[1] | prior_df_rho, prior_mean_rho, prior_scale_rho);
    else if (prior_dist_rho == 3)
      target+= exponential_lpdf(rho[1] | prior_scale_rho);
    /* else prior_dist_tau is 0 and nothing is added */
  }
  else if (model_type == 3) {  // BYM
    target+= normal_lpdf(theta_raw | 0, 1);  // unstructured (random) effect
    if (prior_dist_rho == 1)
      target+= beta_lpdf(rho[1] | shape1_rho, shape2_rho);
    /* else prior_dist_tau is 0 and nothing is added */
  }
  // prior on overall spatial variation
  if (prior_dist_tau == 1)
    target+= normal_lpdf(tau | prior_mean_tau, prior_scale_tau);
  else if (prior_dist_tau == 2)
    target+= student_t_lpdf(tau | prior_df_tau, prior_mean_tau, prior_scale_tau);
  else if (prior_dist_tau == 3)
    target+= exponential_lpdf(tau | prior_scale_tau);
  /* else prior_dist_tau is 0 and nothing is added */
  // priors on auxilliary parameters (Log-priors)
  if (has_aux == 1) {
    // "priors_aux.stan"
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
  }
}
generated quantities {
  real mean_PPD = 0;
  real alpha[has_intercept];
  {
    vector[N] lp_tf;
    vector[N] eta;
    eta = X * beta + psi;
    if (has_intercept == 1) {
      alpha[1] = gamma[1] - dot_product(beta, xbar);
      if ((family == 5 && link == 4)) { // binomial
        eta -= max(eta);
        alpha[1] -= max(eta);
      }
      else if ((family == 6 && link != 1) ||
        (family == 7 && link != 1) ||
        (family == 2 && (link == 1 || link == 3))) {// poisson, neg_binomial_2, and gamma
        eta -= min(eta);
        alpha[1] -= min(eta);
      }
      eta += gamma[1];
    }
    if (family == 1) {
      eta = linkinv_gauss(eta, link);
      for (n in 1:N) mean_PPD += normal_rng(eta[n], aux);
    }
    else if (family == 6) {
      eta = linkinv_count(eta, link);
      if (link == 2)
      for (n in 1:N) {
        if (eta[n] < poisson_max)
          mean_PPD += poisson_rng(eta[n]);
        else
          mean_PPD += normal_rng(eta[n], sqrt(eta[n]));
      }
    }
    else if (family == 7) {
      eta = linkinv_count(eta, link);
      for (n in 1:N) {
          real gamma_temp;
          if (is_inf(aux)) gamma_temp = eta[n];
          else gamma_temp = gamma_rng(aux, aux / eta[n]);
          if (gamma_temp < poisson_max)
            mean_PPD += poisson_rng(gamma_temp);
          else mean_PPD += normal_rng(gamma_temp, sqrt(gamma_temp));
      }
    }
    else if (family == 5) {
      lp_tf = linkinv_binom(eta, link);
      for (n in 1:N) mean_PPD += binomial_rng(trials[n], lp_tf[n]);
    }
    else if (family == 2) {
      if (link > 1) eta = linkinv_gamma(eta, link);
      for (n in 1:N) mean_PPD += gamma_rng(aux, aux / eta[n]);
    }
  }
  mean_PPD = mean_PPD / N;
}
