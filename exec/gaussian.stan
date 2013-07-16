# GLM for a Gaussian outcome via Stan
data {
  int<lower=1> N; # number of observations
  int<lower=1> K; # number of predictors
  matrix[N,K]  X; # predictor matrix
  vector[N]    y; # continuous outcome

  # link function from location to linear predictor
  int<lower=1,upper=3> link; # 1 = identity, 2 = log, 3 = inverse

  int<lower=0,upper=1> has_weights; # 0 = No (weights is a ones vector), 1 = Yes
  vector[N] weights;

  int<lower=0,upper=1> has_offset;  # 0 = No (offset is a zero vector), 1 = Yes
  vector[N] offset;

  # values for hyperparameters
  vector<lower=0>[K-1] prior_scale;
  real<lower=0> prior_scale_for_intercept;
  vector[K-1] prior_mean;
  real prior_mean_for_intercept;
  vector<lower=0>[K-1] prior_df;
  real<lower=0> prior_df_for_intercept;
  real<lower=0> prior_scale_for_dispersion;
}
parameters {
  vector[K] beta;
  real<lower=0> sigma;
}
model {
  vector[N] eta;
  eta <- X * beta;
  if(has_offset == 1) eta <- eta + offset;
  
  if(has_weights == 0) { # unweighted log-likelihoods
    if     (link == 1) y ~ normal(eta, sigma);
    else if(link == 2) y ~ lognormal(eta, sigma);
    else { #link == 3
      for(n in 1:N) eta[n] <- 1.0 / eta[n];
      y ~ normal(eta, sigma);
    }
  }
  else { # weighted log-likelihoods
    vector[N] summands;
    if(link == 1) {
      for(n in 1:N) summands[n] <- normal_log(y[n], eta[n], sigma);
    }
    else if(link == 2) {
      for(n in 1:N) summands[n] <- lognormal_log(y[n], eta[n], sigma);
    }
    else { # link == 3
      for(n in 1:N) summands[n] <- normal_log(y[n], 1.0 / eta[n], sigma);
    }
    lp__ <- lp__ + dot_product(weights, summands);
  }

  # log-priors
  sigma ~ cauchy(0,prior_scale_for_dispersion);
  beta[1] ~ student_t(prior_df_for_intercept, 
                      prior_mean_for_intercept, 
                      prior_scale_for_intercept);
  for(k in 2:K) beta[k] ~ student_t(prior_df[k - 1],
                                    prior_mean[k - 1], 
                                    prior_scale[k - 1]);
}

