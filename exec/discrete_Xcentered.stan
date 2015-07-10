# GLM for a discrete outcome with Gaussian or t priors
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
    if (link < 2 || link > 5) reject("Invalid link");
    if (link == 2) 
      for(n in 1:rows(eta)) pi[n] <- Phi_approx(eta[n]);
    else if (link == 3) 
      for(n in 1:rows(eta)) pi[n] <- cauchy_cdf(eta[n], 0.0, 1.0);
    else if (link == 4) 
      for(n in 1:rows(eta)) pi[n] <- exp(eta[n]); # may be > 1
    else if (link == 5) 
      for(n in 1:rows(eta)) pi[n] <- inv_cloglog(eta[n]);
    return pi;
  }
  vector linkinv_pois(vector eta, int link) {
    vector[rows(eta)] phi;
    if (link < 2 || link > 3) reject("Invalid link");
    if (link == 2) return(eta); # link = identity
    else  # link = sqrt
      for(n in 1:rows(eta)) phi[n] <- square(eta[n]); 
    return phi;
  }
  
  /** 
   * Pointwise (pw) log-likelihood vector
   *
   * @param y The integer array corresponding to the outcome variable.
   * @param link An integer indicating the link function
   * @return A vector
   */
  vector pw_bern(int[] y, vector eta, int link) {
    vector[rows(eta)] ll;
    vector[rows(eta)] pi;
    if (link < 1 || link > 5) 
      reject("Invalid link");
    if (link == 1) { # link = logit
      for (n in 1:rows(eta)) ll[n] <- bernoulli_logit_log(y[n], eta[n]);
    }
    else { # link = probit, cauchit, log, or cloglog
      pi <- linkinv_binom(eta, link);
      for (n in 1:rows(eta)) ll[n] <- bernoulli_log(y[n], pi[n]) ;
    }
    return ll;
  }
  vector pw_pois(int[] y, vector eta, int link) {
    vector[rows(eta)] ll;
    vector[rows(eta)] phi;
    if (link < 1 || link > 3) reject("Invalid link");
    if (link == 1) # link = log
      for (n in 1:rows(eta)) ll[n] <- poisson_log_log(y[n], eta[n]);
    else { # link = identity or sqrt
      phi <- linkinv_pois(eta, link);
      for (n in 1:rows(eta)) ll[n] <- poisson_log(y[n], phi[n]) ;
    }
    return ll;
  }

  /** 
   * Column means of a matrix
   *
   * @param X matrix
   * @return row_vector (of length cols(X)) of column means
   */
  row_vector column_means(matrix X) {
    row_vector[cols(X)] out;
    for (k in 1:cols(X))
      out[k] <- mean(col(X, k));
    return out;
  }
  
  /** 
   * Center predictors
   *
   * @param X matrix containing the predictors to be centered
   * @return matrix with same dimensions as X
   */
  matrix center_predictors(matrix X) {
    matrix[rows(X), cols(X)] Xcent;
    row_vector[cols(X)] colmeans;
    colmeans <- column_means(X);
    for (n in 1:rows(X))
      Xcent[n] <- X[n] - colmeans;
    return Xcent;
  }
}
data {
  # dimensions
  int<lower=1> N; # number of observations
  int<lower=1> K; # number of predictors
  
  # data
  matrix[N,K]  X; # predictor matrix
  int<lower=0> y[N]; # discrete outcome
  
  # intercept
  int<lower=0,upper=1> has_intercept; # 1 = yes
  
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
  vector<lower=0>[K] prior_scale;
  real<lower=0> prior_scale_for_intercept;
  vector[K] prior_mean;
  real prior_mean_for_intercept;
  vector<lower=0>[K] prior_df;
  real<lower=0> prior_df_for_intercept;
}
transformed data {
  matrix[N, K] Xcent;
  Xcent <- center_predictors(X);
}
parameters {
  real alpha0[has_intercept]; # intercept (if has_intercept == 1)
  vector[K] beta; # coefficients
}
model {
  vector[N] eta; # linear predictor
  eta <- Xcent * beta;
  if (has_intercept == 1)
    eta <- eta + alpha0[1];
  if (has_offset == 1) 
    eta <- eta + offset;
  
  // Log-likelihood 
  if (has_weights == 0) { # unweighted log-likelihoods
    if (family == 1) {    # bernoulli case
      if (link == 1) y ~ bernoulli_logit(eta);
      else {
        vector[N] pi;
        pi <- linkinv_binom(eta, link);
        y ~ bernoulli(pi);
      }
    }
    else { # poisson case
      if (link == 1) y ~ poisson_log(eta);
      else {
        vector[N] phi; 
        phi <- linkinv_pois(eta, link);
        y ~ poisson(phi);
      }
    }
  }
  else { # weighted log-likelihoods
    vector[N] summands;
    if (family == 1) # bernoulli case
      summands <- pw_bern(y, eta, link);
    else # poisson case
      summands <- pw_pois(y, eta, link);
    increment_log_prob(dot_product(weights, summands));
  }
  
  // Log-priors for coefficients 
  if (prior_dist == 1) # normal
    beta ~ normal(prior_mean, prior_scale);  
  else # student_t
    beta ~ student_t(prior_df, prior_mean, prior_scale);
   
  // Log-prior for intercept  
  if (has_intercept == 1) {
    if (prior_dist_for_intercept == 1) # normal
      alpha0 ~ normal(prior_mean_for_intercept, prior_scale_for_intercept);
    else # student_t
      alpha0 ~ student_t(prior_df_for_intercept, prior_mean_for_intercept, 
                        prior_scale_for_intercept);
  }
}
generated quantities {
  real alpha[has_intercept];
  if (has_intercept == 1)
    alpha[1] <- alpha0[1] - column_means(X) * beta;
}
