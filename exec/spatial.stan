// CAR (BYM version i.e. IAR)
data {
  int<lower=0> N;                 // number of regions
  int<lower=0> K;                 // number of predictors (inc intercept)
  matrix[N,K] X;                  // model matrix
  int<lower=0> trials[N];         // binomial trials (0 1d array if not applicable)
  int y_int[N];                   // outcome
  real y_real[N];                 // outcome
  int<lower=1,upper=3> family;    // family (1 = Gaussian, 2 = Poisson, 3 = Binomial)
  int E_n;                        // number of adjacency pairs
  int edges[E_n, 2];              // adjacency pairs
  real loc_beta[K];               // priors
  real<lower=0> scale_beta[K];    // priors
  real<lower=0> shape_tau;        // priors
  real<lower=0> scale_tau;        // priors
  real<lower=0> shape_sigma;      // priors
  real<lower=0> scale_sigma;      // priors
  real<lower=0> shape_nu;         // priors
  real<lower=0> scale_nu;         // priors
  real<lower=0> loc_alpha;        // priors
  real<lower=0> scale_alpha;      // priors
  int<lower=0,upper=1> has_intercept;
  int<lower=1,upper=2> mod;
}
transformed data {
  real poisson_max = pow(2.0, 30.0);
}
parameters {
  real alpha[has_intercept];  // intercept
  vector[K] beta;             // predictors on covariates (including intercept)
  vector[N] theta_raw[mod == 1? 1 : 0];        // used for random effect (non-spatial)
  vector[N-1] phi_raw;        // used for random effect (spatial)
  real<lower=0> tau;          // variance i.e. tau^2
  real<lower=0> sigma;        // variance i.e. sigma^2
  real<lower=0> nu[family == 1? 1 : 0];  // applies only if family is gaussian
}
transformed parameters {
  vector[N] theta[mod == 1? 1 : 0];        // non-centered random effect (non-spatial)
  vector[N] phi;          // non-centered random effect (spatial)
  theta[1] = sqrt(sigma) * theta_raw[1];
  phi[1:(N - 1)] = phi_raw;
  phi[N] = -sum(phi_raw);
  phi = phi * sqrt(tau);
}
model {
  vector[N] eta;   // linear predictor + spatial random effects
  // model
  if (has_intercept == 1) {
    if (mod == 1)
      eta = alpha[1] + X * beta + phi + theta[1];
    else if (mod == 2)
      eta = alpha[1] + X * beta + phi;
  }
  else {
    if (mod == 1)
      eta = X * beta + phi + theta[1];
    else if (mod == 2)
      eta = X * beta + phi;
  }
  if (family == 1) {
    target+= normal_lpdf(y_real | eta, nu[1]);
    target+= inv_gamma_lpdf(nu[1] | shape_nu, scale_nu);
  }
  else if (family == 2) {
    target+= poisson_log_lpmf(y_int | eta);
  }
  else if (family == 3) {
    target+= binomial_lpmf(y_int | trials, inv_logit(eta));
  }
  // priors
  if (mod == 1)
    target+= - N * log(sqrt(tau)) - 0.5 * inv(tau) * dot_self(phi[edges[,1]] - phi[edges[,2]]);
  else if (mod == 2)
    reject("leroux model not supported")
  if (has_intercept == 1)
    target += normal_lpdf(alpha | loc_alpha, scale_alpha);
  target+= normal_lpdf(beta | loc_beta , scale_beta);
  target+= normal_lpdf(theta_raw[1] | 0, 1);
  target+= normal_lpdf(phi_raw | 0, 1);
  target+= inv_gamma_lpdf(tau | shape_tau, scale_tau);
  target+= inv_gamma_lpdf(sigma | shape_sigma, scale_sigma);
}
generated quantities {
  real mean_PPD = 0;
  vector[N] psi;
  if (mod == 1)
    psi = theta[1] + phi;
  else if (mod == 2)
    psi = phi;
  {
    vector[N] eta;
    if (has_intercept == 1) {
      if (mod == 1)
        eta = alpha[1] + X * beta + phi + theta[1];
      else if (mod == 2)
        eta = alpha[1] + X * beta + phi;
    }
    else {
      if (mod == 1)
        eta = X * beta + phi + theta[1];
      else if (mod == 2)
        eta = X * beta + phi;
    }
    for (n in 1:N) {
      if (family == 1)
        mean_PPD = mean_PPD + normal_rng(eta[n], nu[1]);
      else if (family == 2) {
        if (eta[n] < poisson_max)
          mean_PPD = mean_PPD + poisson_log_rng(eta[n]);
        else
          mean_PPD = mean_PPD + normal_rng(eta[n], sqrt(eta[n]));
      }
      else if (family == 3)
        mean_PPD = mean_PPD + binomial_rng(trials[n], inv_logit(eta[n]));
    }
  }
  mean_PPD = mean_PPD / N;
}
