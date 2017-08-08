// CAR SPATIAL MODELS
functions {
  #include "continuous_likelihoods.stan"
  #include "binomial_likelihoods.stan"
  #include "count_likelihoods.stan"
  // add negative_binomial likelihood
  // add gamma likelihood
}
data {
  int<lower=0> N;                 // number of regions
  int<lower=0> K;                 // number of predictors (inc intercept)
  matrix[N,K] X;                  // model matrix
  vector[K] xbar;
  int<lower=0> trials[N];         // binomial trials (0 1d array if not applicable)
  int y_int[N];                   // outcome
  real y_real[N];                 // outcome
  int<lower=1,upper=3> family;    // family (1 = Gaussian, 2 = Poisson, 3 = Binomial)
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
  int<lower=0> prior_dist_sigma;
  real prior_mean_for_intercept;
  real<lower=0> prior_scale_for_intercept;
  real<lower=0> prior_df_for_intercept;
  vector[K] prior_mean;
  vector<lower=0>[K] prior_scale;
  real<lower=0> prior_df[K];
  real prior_mean_tau;
  real<lower=0> prior_scale_tau;
  real<lower=0> prior_df_tau;
  real prior_mean_rho;
  real<lower=0> prior_scale_rho;
  real<lower=0> prior_df_rho;
  real prior_mean_sigma;
  real<lower=0> prior_scale_sigma;
  real<lower=0> prior_df_sigma;
}
transformed data {
  real poisson_max = 30 * log(2);
}
parameters {
  real gamma[has_intercept];  // raw intercept
  vector[K] beta;             // predictors on covariates (including intercept)
  vector[N] theta_raw[mod == 2? 1 : 0];        // used for random effect (non-spatial)
  vector[N-1] phi_raw;        // used for random effect (spatial)
  real<lower=0,upper=(mod == 2? 1: positive_infinity())> rho;          // variance i.e. rho^2
  real<lower=0> tau[mod == 2? 1 : 0];        // variance i.e. tau^2
  real<lower=0> sigma[family == 1? 1 : 0];  // applies only if family is gaussian
}
transformed parameters {
  vector[N] phi;          // non-centered random effect (spatial)
  vector[N] psi;
  phi[1:(N - 1)] = phi_raw;
  phi[N] = -sum(phi_raw);
  if (mod == 1)
    psi = phi * sqrt(inv(rho));
  else if (mod == 2)
    psi = tau[1]*(sqrt(1-rho)*theta_raw[1] + sqrt(rho/scaling_factor)*phi);
    // psi = tau[1]*(sqrt(rho)*theta_raw[1] + sqrt(1-rho)*scaling_factor*phi);
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
    target+= normal_lpdf(y_real | eta, sigma[1]);
  }
  else if (family == 2) {
    eta = linkinv_count(eta, link);
    target+= poisson_log_lpmf(y_int | eta);
  }
  else if (family == 3) {
    eta = linkinv_binom(eta, link);
    target+= binomial_lpmf(y_int | trials, inv_logit(eta));
  }
  // prior on spatial parameter vector
  target += -0.5 * dot_self(phi[edges[,1]] - phi[edges[,2]]);
  // priors on coefficients
  if (has_intercept == 1) {
    if (prior_dist_for_intercept == 1)
      target+= normal_lpdf(gamma | prior_mean_for_intercept, prior_scale_for_intercept);
    else if (prior_dist_for_intercept == 2)
      target+= student_t_lpdf(gamma | prior_df_for_intercept, prior_mean_for_intercept, prior_scale_for_intercept);
    /* else prior_dist_intercept is 0 and nothing is added */
  }
  if (K > 0) {
    if (prior_dist == 1)
      target+= normal_lpdf(beta | prior_mean, prior_scale);
    else if (prior_dist == 2)
      target+= student_t_lpdf(beta | prior_df, prior_mean, prior_scale);
    /* else prior_dist is 0 and nothing is added */
  }
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
  if (family == 1) { // prior on sd if outcome is gaussian
    if (prior_dist_sigma == 1)
      target+= normal_lpdf(sigma[1] | prior_mean_sigma, prior_scale_sigma);
    else if (prior_dist_sigma == 2)
      target+= student_t_lpdf(sigma[1] | prior_df_sigma, prior_mean_sigma, prior_scale_sigma);
    /* else prior_dist_sigma is 0 and nothing is added */
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
      for (n in 1:N) mean_PPD = mean_PPD + normal_rng(eta[n], sigma[1]);
    }
    else if (family == 2) {
      eta = linkinv_count(eta, link);
      for (n in 1:N) {
        if (eta[n] < poisson_max)
          mean_PPD = mean_PPD + poisson_log_rng(eta[n]);
        else
          mean_PPD = mean_PPD + normal_rng(eta[n], sqrt(eta[n]));
      }
    }
    else if (family == 3) {
      eta = linkinv_binom(eta, link);
      for (n in 1:N) mean_PPD = mean_PPD + binomial_rng(trials[n], inv_logit(eta[n]));
    }
  }
  mean_PPD = mean_PPD / N;
}
