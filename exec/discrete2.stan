# GLM for a discrete outcome with Gaussian or t priors
data {
  # dimensions
  int<lower=1> N; # number of observations
  int<lower=1> K; # number of predictors
  
  # data
  matrix[N,K]  X; # predictor matrix
  int<lower=0> y[N]; # discrete outcome
  
  # likelihood
  int<lower=1> family; # 1 = binomial, other = poisson
  
  # link function from location to linear predictor
  int<lower=1,upper=5> link; # values depend on family
  
  # weights
  int<lower=0,upper=1> has_weights; # 0 = No (weights is a ones vector), 1 = Yes
  vector[N] weights;
  
  # offset
  int<lower=0,upper=1> has_offset;  # 0 = No (offset is a zero vector), 1 = Yes
  vector[N] offset;
  
  # prior family
  int<lower=1,upper=2> prior_dist; # 1 = normal, 2 = student_t
  int<lower=1,upper=2> prior_dist_for_intercept; # 1 = normal, 2 = student_t
  
  # hyperparameter values
  vector<lower=0>[K-1] prior_scale;
  real<lower=0> prior_scale_for_intercept;
  vector[K-1] prior_mean;
  real prior_mean_for_intercept;
  vector<lower=0>[K-1] prior_df;
  real<lower=0> prior_df_for_intercept;
}
parameters {
  vector[K] beta;
}
model {
  vector[N] eta;
  eta <- X * beta;
  if (has_offset == 1) 
    eta <- eta + offset;
  
  if (has_weights == 0) { # unweighted log-likelihoods
    if (family == 1) {    # bernoulli case
      if (link == 1) y ~ bernoulli_logit(eta);
      else {
        vector[N] pi;
        if (link == 2)      for(n in 1:N) pi[n] <- Phi_approx(eta[n]);
        else if (link == 3) for(n in 1:N) pi[n] <- cauchy_cdf(eta[n], 0.0, 1.0);
        else if (link == 4) for(n in 1:N) pi[n] <- exp(eta[n]); # may be > 1
        else if (link == 5) for(n in 1:N) pi[n] <- inv_cloglog(eta[n]);
        y ~ bernoulli(pi);
      }
    }
    else { # poisson case
      y ~ poisson_log(eta);
    }
  }
  else { # weighted log-likelihoods
    vector[N] summands;
    if (family == 1) { # bernoulli case
      if (link == 1) {
        for (n in 1:N) 
          summands[n] <- bernoulli_logit_log(y[n], eta[n]);
      }
      else if (link == 2) {
        for (n in 1:N) 
          summands[n] <- bernoulli_log(y[n], Phi(eta[n])) ;
      }
      else if (link == 3) {
        for (n in 1:N) 
          summands[n] <- bernoulli_log(y[n], cauchy_cdf(eta[n], 0.0, 1.0));
      }
      else if (link == 4) {
        for (n in 1:N) 
          summands[n] <- bernoulli_log(y[n], exp(eta[n]));
      }
      else if (link == 5) {
        for (n in 1:N) 
          summands[n] <- bernoulli_log(y[n], inv_cloglog(eta[n]));
      }
    }
    else { # poisson case
      for (n in 1:N) 
        summands[n] <- poisson_log_log(y[n], eta[n]);
    }
    increment_log_prob(dot_product(weights, summands));
  }
  
  # log-priors
  if (prior_dist_for_intercept == 1) { # normal
    beta[1] ~ normal(prior_mean_for_intercept, prior_scale_for_intercept);
  }
  else { # student_t
    beta[1] ~ student_t(prior_df_for_intercept, prior_mean_for_intercept, prior_scale_for_intercept);
  }
  
  if (prior_dist == 1) { # normal
    tail(beta, K-1) ~ normal(prior_mean, prior_scale);  
  } 
  else { # student_t
    tail(beta, K-1) ~ student_t(prior_df, prior_mean, prior_scale);
  }
}
/*
generated quantities {
  vector[N] log_lik;
  vector[N] eta;
  
  eta <- X * beta;
  if (has_offset == 1) 
    eta <- eta + offset;
  
  if (family == 1) { # bernoulli
      vector[N] pi;
      if (link == 1)      for (n in 1:N) pi[n] <- inv_logit(eta[n]);
      else if (link == 2) for (n in 1:N) pi[n] <- Phi_approx(eta[n]);
      else if (link == 3) for (n in 1:N) pi[n] <- cauchy_cdf(eta[n], 0.0, 1.0);
      else if (link == 4) for (n in 1:N) pi[n] <- exp(eta[n]); # may be > 1
      else if (link == 5) for (n in 1:N) pi[n] <- inv_cloglog(eta[n]);
      
      for (n in 1:N) 
        log_lik[n] <- bernoulli_log(y[n], pi[n]);
  }
  else { # poisson
    for (n in 1:N) 
      log_lik[n] <- poisson_log_log(y[n], eta[n]);
  }
  
  if (has_weights == 1)
    log_lik <- log_lik .* weights;
}
*/