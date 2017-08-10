#include "Columbia_copyright.stan"
#include "license.stan" // GPL3+
// CAR SPATIAL MODELS
functions {
  #include "continuous_likelihoods.stan"
  #include "binomial_likelihoods.stan"
  #include "count_likelihoods.stan"
  #include "common_functions.stan"
}
data {
  int<lower=0> N;                 // number of regions
  int<lower=0> K;                 // number of predictors (inc intercept)
  matrix[N,K] X;                  // model matrix
  vector[K] xbar;
  int<lower=0> trials[N];         // binomial trials (0 1d array if not applicable)
  int y_int[N];                   // outcome
  vector[N] y_real;                 // outcome
  int<lower=1,upper=5> family;    // family (1 = Gaussian, 2 = Poisson, 3 = Binomial)
  int link;
  int E_n;                        // number of adjacency pairs
  int edges[E_n, 2];              // adjacency pairs
  real<lower=0> shape1_rho;        // priors
  real<lower=0> shape2_rho;        // priors
  int<lower=0,upper=1> has_intercept;
  int<lower=1,upper=2> mod;       // 1 = besag (icar); 2 = bym
  real scaling_factor;
  int<lower=0> prior_dist_for_intercept;
  int<lower=0> prior_dist;
  int<lower=0> prior_dist_tau;
  int<lower=0> prior_dist_rho;
  real prior_mean_for_intercept;
  real<lower=0> prior_scale_for_intercept;
  real<lower=0> prior_df_for_intercept;
  vector[K] prior_mean;
  vector<lower=0>[K] prior_scale;
  vector<lower=0>[K] prior_df;
  real prior_mean_tau;
  real<lower=0> prior_scale_tau;
  real<lower=0> prior_df_tau;
  real prior_mean_rho;
  real<lower=0> prior_scale_rho;
  real<lower=0> prior_df_rho;
  real<lower=0> global_prior_df;    // for hs priors only
  real<lower=0> global_prior_scale; // for hs priors only
  int<lower=2> num_normals[prior_dist == 7 ? K : 0];
  int<lower=0,upper=3> prior_dist_for_aux;
  real<lower=0> prior_mean_for_aux;
  real<lower=0> prior_scale_for_aux;
  real<lower=0> prior_df_for_aux;
}
transformed data {
  real poisson_max = 30 * log(2);
  int<lower=0> hs;
  int<lower=0,upper=1> is_continuous;
  real sum_log_y = family == 1 ? not_a_number() : sum(log(y_real));
  if (prior_dist <= 2) hs = 0;
  else if (prior_dist == 3) hs = 2;
  else if (prior_dist == 4) hs = 4;
  else hs = 0;
  if (family == 1) is_continuous = 1;
  else is_continuous = 0;
}
parameters {
  real gamma[has_intercept];  // raw intercept
  vector[K] z_beta;
  vector[mod == 2? N : 0] theta_raw;        // used for random effect (non-spatial)
  vector[N-1] phi_raw;        // used for random effect (spatial)
  real<lower=0,upper=(mod == 2? 1: positive_infinity())> rho;          // variance i.e. rho^2
  real<lower=0> tau[mod == 2? 1 : 0];        // variance i.e. tau^2
  real<lower=0> global[hs];
  vector<lower=0>[K] local[hs];
  vector<lower=0>[K] mix[prior_dist == 5 || prior_dist == 6];
  real<lower=0> aux_unscaled; # interpretation depends on family!
  real<lower=0> one_over_lambda[prior_dist == 6];
}
transformed parameters {
  vector[K] beta;             // predictors on covariates (including intercept)
  // aux has to be defined first in the hs case
  real aux = prior_dist_for_aux == 0 ? aux_unscaled : (prior_dist_for_aux <= 2 ?
             prior_scale_for_aux * aux_unscaled + prior_mean_for_aux :
             prior_scale_for_aux * aux_unscaled);
  vector[N] phi;          // non-centered random effect (spatial)
  vector[N] psi;
  phi[1:(N - 1)] = phi_raw;
  phi[N] = -sum(phi_raw);
  if (mod == 1)
    psi = phi * sqrt(inv(rho));
  else if (mod == 2)
    psi = tau[1]*(sqrt(1-rho)*theta_raw[1] + sqrt(rho/scaling_factor)*phi);
  // for regression coefficients
  #include "tparameters.stan"
}
model {
  vector[N] eta;   // linear predictor + spatial random effects
  // deal with intercept
  if (has_intercept == 1)
    eta = gamma[1] + X * beta + psi;
  else
    eta = X * beta + psi;
  // likelihoods
  if (family == 1) {
    eta = linkinv_gauss(eta, link);
    target+= normal_lpdf(y_real | eta, aux);
  }
  else if (family == 2) {
    eta = linkinv_count(eta, link);
    target+= poisson_lpmf(y_int | eta);
  }
  else if (family == 3) {
    eta = linkinv_count(eta, link);
    target += neg_binomial_2_lpmf(y_int | eta, aux);
  }
  else if (family == 4) {
    eta = linkinv_binom(eta, link);
    target+= binomial_lpmf(y_int | trials, inv_logit(eta));
  }
  else if (family == 5) {
    target += GammaReg(y_real, eta, aux, link, sum_log_y);
  }
  // prior on spatial parameter vector (GMRF)
  target += -0.5 * dot_self(phi[edges[,1]] - phi[edges[,2]]);
  // non-centered parameterization on beta
  target+= normal_lpdf(z_beta | 0, 1);
  // priors on coefficients
  #include "priors.stan"
  // model specific priors
  if (mod == 2) { // BYM
    target+= normal_lpdf(theta_raw[1] | 0, 1);  // unstructured (random) effect
    if (prior_dist_rho == 1)
      target+= beta_lpdf(rho | shape1_rho, shape2_rho);
    /* else prior_dist_rho is 0 and nothing is added */
    if (prior_dist_tau == 1)
      target+= normal_lpdf(tau | prior_mean_tau, prior_scale_tau);
    else if (prior_dist_tau == 2)
      target+= student_t_lpdf(tau | prior_df_tau, prior_mean_tau, prior_scale_tau);
    /* else prior_dist_tau is 0 and nothing is added */
  }
  else { // besag
    if (prior_dist_rho == 1)
      target+= normal_lpdf(rho | prior_mean_rho, prior_scale_rho);
    else if (prior_dist_rho == 2)
      target+= student_t_lpdf(rho | prior_df_rho, prior_mean_rho, prior_scale_rho);
    /* else prior_dist_rho is 0 and nothing is added */
  }
  // priors on aux
  #include "priors_aux.stan"
}
generated quantities {
  real mean_PPD = 0;
  real alpha[has_intercept];
  if (has_intercept)
    alpha[1] = gamma[1] - dot_product(beta, xbar);
  {
    vector[N] eta;
    if (has_intercept == 1) {
      eta = alpha[1] + X * beta + psi;
    }
    else {
      eta = X * beta + psi;
    }
    if (family == 1) {
      eta = linkinv_gauss(eta, link);
      for (n in 1:N) mean_PPD = mean_PPD + normal_rng(eta[n], aux);
    }
    else if (family == 2) {
      eta = linkinv_count(eta, link);
      for (n in 1:N) {
        if (eta[n] < poisson_max)
          mean_PPD = mean_PPD + poisson_rng(eta[n]);
        else
          mean_PPD = mean_PPD + normal_rng(eta[n], sqrt(eta[n]));
      }
    }
    else if (family == 3) {
      eta = linkinv_count(eta, link);
      for (n in 1:N) {
          real gamma_temp;
          if (is_inf(aux)) gamma_temp = eta[n];
          else gamma_temp = gamma_rng(aux, aux / eta[n]);
          if (gamma_temp < poisson_max)
            mean_PPD = mean_PPD + poisson_rng(gamma_temp);
          else mean_PPD = mean_PPD + normal_rng(gamma_temp, sqrt(gamma_temp));
      }
    }
    else if (family == 4) {
      eta = linkinv_binom(eta, link);
      for (n in 1:N) mean_PPD = mean_PPD + binomial_rng(trials[n], inv_logit(eta[n]));
    }
    else if (family == 5) {
      if (link > 1) eta = linkinv_gamma(eta, link);
      for (n in 1:N) mean_PPD = mean_PPD + gamma_rng(aux, aux / eta[n]);
    }
  }
  mean_PPD = mean_PPD / N;
}
