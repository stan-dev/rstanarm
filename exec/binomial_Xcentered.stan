# GLM for a binomial outcome with Gaussian or t priors
functions {
  /** 
   * Apply inverse link function to linear predictor
   *
   * @param eta Linear predictor vector
   * @param link An integer indicating the link function
   * @return A vector, i.e. inverse-link(eta)
   */
  vector linkinv_binom(vector eta, int link) {
    vector[rows(eta)] pi;
    if (link < 1 || link > 5) reject("Invalid link");
    if      (link == 1)
      for(n in 1:rows(eta)) pi[n] <- inv_logit(eta[n]);
    else if (link == 2) 
      for(n in 1:rows(eta)) pi[n] <- Phi(eta[n]);
    else if (link == 3) 
      for(n in 1:rows(eta)) pi[n] <- cauchy_cdf(eta[n], 0.0, 1.0);
    else if (link == 4) 
      for(n in 1:rows(eta)) pi[n] <- exp(eta[n]);
    else if (link == 5) 
      for(n in 1:rows(eta)) pi[n] <- inv_cloglog(eta[n]);
    return pi;
  }

  /**
   * Increment with the unweighted log-likelihood
   * @param y An integer array indicating the number of successes
   * @param trials An integer array indicating the number of trials
   * @param eta A vector of linear predictors
   * @param link An integer indicating the link function
   * @return lp__
   */
  real ll_binom_lp(int[] y, int[] trials, vector eta, int link) {
    if (link < 1 || link > 5) reject("Invalid link");
    if      (link == 1) y ~ binomial_logit(trials, eta);
    else if (link <  4) y ~ binomial(trials, linkinv_binom(eta, link));
    else if (link == 4) { // log link
      for (n in 1:num_elements(y)) {
        increment_log_prob(y[n] * eta[n]);
        increment_log_prob( (trials[n] - y[n]) * log1m_exp(eta[n]) );
      }
    }
    else if(link == 5) { // cloglog link
      real neg_exp_eta;
      for (n in 1:num_elements(y)) {
        neg_exp_eta <- -exp(eta[n]);
        increment_log_prob(y[n] * log1m_exp(neg_exp_eta));
        increment_log_prob( (trials[n] - y[n]) * neg_exp_eta );
      }
    }
    return get_lp();
  }
  
  /** 
   * Pointwise (pw) log-likelihood vector
   *
   * @param y The integer array corresponding to the outcome variable.
   * @param link An integer indicating the link function
   * @return A vector
   */
  vector pw_binom(int[] y, int[] trials, vector eta, int link) {
    vector[rows(eta)] ll;
    if (link < 1 || link > 5) reject("Invalid link");
    if (link == 1) { # link = logit
      for (n in 1:rows(eta)) 
        ll[n] <- binomial_logit_log(y[n], trials[n], eta[n]);
    }
    else { # link = probit, cauchit, log, or cloglog (unstable)
      vector[rows(eta)] pi;
      pi <- linkinv_binom(eta, link);
      for (n in 1:rows(eta)) ll[n] <- binomial_log(y[n], trials[n], pi[n]) ;
    }
    return ll;
  }

  /** 
   * Upper bound on the intercept, which is infinity except for log link
   *
   * @param link An integer indicating the link function
   * @param X A matrix of predictors | y = 0
   * @param beta A vector of coefficients
   * @param has_offset An integer indicating an offset
   * @param offset A vector of offsets
   * @return A scalar upper bound on the intercept
   */
  real make_upper_binomial(int link, matrix X, vector beta, 
                           int has_offset, vector offset) {
    real maximum;
    if (link != 4) return positive_infinity();
    if (has_offset == 0) maximum <- max(X * beta);
    else maximum <- max(X * beta + offset);
    return -maximum;
  }
}
data {
  # dimensions
  int<lower=1> N; # number of observations
  int<lower=1> K; # number of predictors
  
  # data
  vector[K] xbar;         # predictor means
  matrix[N,K]  X;         # centered predictor matrix
  int<lower=0> y[N];      # outcome: number of successes
  int<lower=0> trials[N]; # number of trials
  
  # intercept
  int<lower=0,upper=1> has_intercept; # 0 = no, 1 = yes
  
  # link function from location to linear predictor
  int<lower=1,upper=5> link;
  
  # weights
  int<lower=0,upper=1> has_weights; # 0 = No, 1 = Yes
  vector[N * has_weights] weights;
  
  # offset
  int<lower=0,upper=1> has_offset;  # 0 = No, 1 = Yes
  vector[N * has_offset] offset;
  
  # prior family (zero indicates no prior!!!)
  int<lower=0,upper=2> prior_dist;               # 1 = normal, 2 = student_t
  int<lower=0,upper=2> prior_dist_for_intercept; # 1 = normal, 2 = student_t
  
  # hyperparameter values
  vector<lower=0>[K] prior_scale;
  real<lower=0> prior_scale_for_intercept;
  vector[K] prior_mean;
  real prior_mean_for_intercept;
  vector<lower=0>[K] prior_df;
  real<lower=0> prior_df_for_intercept;
}
parameters {
  vector[K] beta; # coefficients
  real<lower=negative_infinity(),upper=make_upper_binomial(link,
       X, beta, has_offset, offset)> gamma[has_intercept];
}
model {
  vector[N] eta; # linear predictor
  eta <- X * beta;
  if (has_intercept == 1) eta <- eta + gamma[1];
  if (has_offset == 1)    eta <- eta + offset;
  
  // Log-likelihood 
  if (has_weights == 0) { # unweighted log-likelihoods
    real dummy; # irrelevant but useful for testing
    dummy <- ll_binom_lp(y, trials, eta, link);
  }
  else increment_log_prob(dot_product(weights, 
                          pw_binom(y, trials, eta, link)));
  
  // Log-priors for coefficients 
  if (prior_dist == 1) # normal
    beta ~ normal(prior_mean, prior_scale);  
  else if (prior_dist == 2) # student_t
    beta ~ student_t(prior_df, prior_mean, prior_scale);
  /* else prior_dist = 0 and nothing is added */
  
  // Log-prior for intercept  
  if (has_intercept == 1) {
    if (link != 4) {
      if (prior_dist_for_intercept == 1) # normal
        gamma ~ normal(prior_mean_for_intercept, prior_scale_for_intercept);
      else if (prior_dist_for_intercept == 2) # student_t
        gamma ~ student_t(prior_df_for_intercept, prior_mean_for_intercept, 
                          prior_scale_for_intercept);
      /* else prior_dist = 0 and nothing is added */
    }
    else {
      real maximum;
      maximum <- -max(eta);
      if (prior_dist_for_intercept == 1) # normal
        gamma[1] ~ normal(prior_mean_for_intercept, 
                          prior_scale_for_intercept) T[,maximum];
      else if (prior_dist_for_intercept == 2) # student_t
        gamma[1] ~ student_t(prior_df_for_intercept, prior_mean_for_intercept, 
                             prior_scale_for_intercept) T[,maximum];
      /* else prior_dist = 0 and nothing is added */
    }
  }
}
generated quantities {
  real alpha[has_intercept];
  real mean_PPD;
  if (has_intercept == 1) alpha[1] <- gamma[1] - dot_product(xbar, beta);
  mean_PPD <- 0;
  {
    vector[N] eta; 
    vector[N] pi;
    eta <- X * beta;
    if (has_intercept == 1) eta <- eta + gamma[1];
    if (has_offset == 1)    eta <- eta + offset;
    pi <- linkinv_binom(eta, link);
    for (n in 1:N) mean_PPD <- mean_PPD + binomial_rng(trials[n], pi[n]);
    mean_PPD <- mean_PPD / N;
  }
}
