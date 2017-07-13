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
  real<lower=0> shape1_tau;        // priors
  real<lower=0> shape2_tau;        // priors
  real<lower=0> loc_sigma;      // priors
  real<lower=0> scale_sigma;      // priors
  real<lower=0> loc_nu;         // priors
  real<lower=0> scale_nu;         // priors
  real<lower=0> loc_alpha;        // priors
  real<lower=0> scale_alpha;      // priors
  int<lower=0,upper=1> has_intercept;
  int<lower=1,upper=2> mod;       // 1 = icar; 2 = bym
  real scaling_factor;
}
transformed data {
  real poisson_max = pow(2.0, 30.0);
}
parameters {
  real alpha[has_intercept];  // intercept
  vector[K] beta;             // predictors on covariates (including intercept)
  vector[N] theta_raw[mod == 2? 1 : 0];        // used for random effect (non-spatial)
  vector[N-1] phi_raw;        // used for random effect (spatial)
  real<lower=0,upper=1> tau[mod == 2? 1 : 0];          // variance i.e. tau^2
  real<lower=0> sigma[mod == 2? 1 : 0];        // variance i.e. sigma^2
  real<lower=0> nu[family == 1? 1 : 0];  // applies only if family is gaussian
}
transformed parameters {
  vector[N] phi;          // non-centered random effect (spatial)
  vector[N] psi;
  phi[1:(N - 1)] = phi_raw;
  phi[N] = -sum(phi_raw);
  if (mod == 1)
    psi = phi;
  else if (mod == 2)
    psi = sigma[1] *(sqrt(tau[1])*theta_raw[1] + sqrt(1-tau[1])*scaling_factor*phi);
}
model {
  vector[N] eta;   // linear predictor + spatial random effects
  // model
  if (has_intercept == 1) {
    eta = alpha[1] + X * beta + psi;
  }
  else {
    eta = X * beta + psi;
  }
  if (family == 1) {
    target+= normal_lpdf(y_real | eta, nu[1]);
    target+= normal_lpdf(nu[1] | loc_nu, scale_nu);
  }
  else if (family == 2) {
    target+= poisson_log_lpmf(y_int | eta);
  }
  else if (family == 3) {
    target+= binomial_lpmf(y_int | trials, inv_logit(eta));
  }
  // priors
  target += -0.5 * dot_self(phi[edges[,1]] - phi[edges[,2]]);
  if (has_intercept == 1)
    target += normal_lpdf(alpha | loc_alpha, scale_alpha);
  target+= normal_lpdf(beta | loc_beta , scale_beta);
  if (mod == 2) {
    target+= normal_lpdf(theta_raw[1] | 0, 1);
    target+= normal_lpdf(sigma | loc_sigma, scale_sigma);
    target+= beta_lpdf(tau | shape1_tau, shape2_tau);
  }
}
generated quantities {
  real mean_PPD = 0;
  {
    vector[N] eta;
    if (has_intercept == 1) {
      eta = alpha[1] + X * beta + psi;
    }
    else {
      eta = X * beta + psi;
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
