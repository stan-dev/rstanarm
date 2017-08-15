#include "Columbia_copyright.stan"
#include "license.stan" // GPL3+
// CAR SPATIAL MODELS
functions {
  #include "continuous_likelihoods.stan"
  #include "binomial_likelihoods.stan"
  #include "count_likelihoods.stan"
  #include "common_functions.stan"
  /*
   * Calculate lower bound on intercept
   *
   * @param family Integer family code
   * @param link Integer link code
   * @return real lower bound
   */
  real make_lower(int family, int link) {
    if (family == 1) return negative_infinity(); // Gaussian
    if (family == 5) { // Gamma
      if (link == 2) return negative_infinity(); // log
      return 0; // identity or inverse
    }
    if (family == 2 || family == 3) { // Poisson or nb2
      if (link == 1) return negative_infinity(); // log
      return 0.0; // identity or sqrt
    }
    return negative_infinity();
  }

  /*
   * Calculate upper bound on intercept
   *
   * @param family Integer family code
   * @param link Integer link code
   * @return real upper bound
   */
  real make_upper(int family, int link) {
    if (family == 4 && link == 4) return 0.0;  // binomial; log
    return positive_infinity();
  }
}
data {
  int<lower=0> N;                 // number of regions
  int<lower=0> K;                 // number of predictors (inc intercept)
  matrix[N,K] X;                  // model matrix
  int<lower=1,upper=5> family;
  int link;
  int<lower=0,upper=1> is_continuous;
  int<lower=0,upper=1> has_aux;
  int<lower=1,upper=3> model_type;       // Besag = 1; BYM = 2; BYM2 = 3
  int<lower=0,upper=1> has_intercept;
  vector[K] xbar;
  int<lower=0> trials[family == 4 ? N : 0];
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
  vector[K] z_beta;
  vector[model_type == 1? 0 : N] theta_raw;
  vector[N-1] phi_raw;
  // interpretation of rho and tau depends on model_type!
  real<lower=0,upper=(model_type == 3 ? 1 : positive_infinity())> rho[model_type != 1];
  real<lower=0> tau;
  real<lower=0> global[hs];
  vector<lower=0>[K] local[hs];
  vector<lower=0>[K] mix[prior_dist == 5 || prior_dist == 6];
  real<lower=0> aux_unscaled[has_aux]; # interpretation depends on family!
  real<lower=0> one_over_lambda[prior_dist == 6];
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
    target+= binomial_lpmf(y_int | trials, eta);
  }
  else if (family == 5) {
    target += GammaReg(y_real, eta, aux, link, sum_log_y);
  }
  // prior on spatial parameter vector (GMRF)
  if (order == 1)
    target += -0.5 * dot_self(phi[edges[,1]] - phi[edges[,2]]);
  else if (order == 2)
    target+= -0.5 * dot_product(phi, csr_matrix_times_vector(N, N, w, v, u, phi));
  // priors on coefficients
  #include "priors.stan"
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
    #include "priors_aux.stan"
  }
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
