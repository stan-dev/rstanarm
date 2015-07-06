# GLM for a Gaussian outcome with Gaussian or t priors
functions {
  /** 
   * Apply inverse link function to linear predictor
   *
   * @param eta Linear predictor vector
   * @param link An integer indicating the link function
   * @return A vector, i.e. inverse-link(eta)
   */
  vector linkinv_gaus(vector eta, int link) {
    if (link > 3) reject("Invalid link");
    if (link == 1 || link == 2) # link = identity or log 
      return(eta); # return eta for log link too bc will use lognormal
    else {# link = inverse
      vector[rows(eta)] mu;
      for(n in 1:rows(eta)) mu[n] <- inv(eta[n]); 
      return mu;
    }
  }
  
  /** 
   * Pointwise (pw) log-likelihood vector
   *
   * @param y The integer array corresponding to the outcome variable.
   * @param link An integer indicating the link function
   * @return A vector
   */
  vector pw_gaus(vector y, vector eta, real sigma, int link) {
    vector[rows(eta)] ll;
    if (link > 3) 
      reject("Invalid link");
    if (link == 2) # link = log
      for (n in 1:rows(eta)) ll[n] <- lognormal_log(y[n], eta[n], sigma);
    else { # link = idenity or inverse
      vector[rows(eta)] mu;
      mu <- linkinv_gaus(eta, link);
      for (n in 1:rows(eta)) ll[n] <- normal_log(y[n], mu[n], sigma);
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
  vector[N]    y; # continuous outcome
  
  # intercept
  int<lower=0,upper=1> has_intercept; # 1 = yes
  
  # link function from location to linear predictor
  int<lower=1,upper=3> link; # 1 = identity, 2 = log, 3 = inverse
  
  # weights
  int<lower=0,upper=1> has_weights; # 0 = No (weights is a ones vector), 1 = Yes
  vector[N] weights;
  
  # offset
  int<lower=0,upper=1> has_offset;  # 0 = No (offset is a zero vector), 1 = Yes
  vector[N] offset;
  
  # prior distributions
  int<lower=1,upper=2> prior_dist; # 1 = normal, 2 = student_t
  int<lower=1,upper=2> prior_dist_for_intercept; # 1 = normal, 2 = student_t
  
  # hyperparameter values
  vector<lower=0>[K] prior_scale;
  real<lower=0> prior_scale_for_intercept;
  vector[K] prior_mean;
  real prior_mean_for_intercept;
  vector<lower=0>[K] prior_df;
  real<lower=0> prior_df_for_intercept;
  real<lower=0> prior_scale_for_dispersion;
}
transformed data {
  matrix[N, K] Xcent;
  Xcent <- center_predictors(X);
}
parameters {
  real alpha0[has_intercept];
  vector[K] theta;
  real<lower=0> sigma;
}
model {
  vector[N] eta; # linear predictor
  eta <- Xcent * theta;
  if (has_intercept == 1)
    eta <- eta + alpha0[1];
  if (has_offset == 1) 
    eta <- eta + offset;
  
  // Log-likelihood 
  if (has_weights == 0) { # unweighted log-likelihoods
    vector[N] mu;
    mu <- linkinv_gaus(eta, link);
    if (link == 2)
      y ~ lognormal(mu, sigma);
    else 
      y ~ normal(mu, sigma);
  }
  else { # weighted log-likelihoods
    vector[N] summands;
    summands <- pw_gaus(y, eta, sigma, link);
    increment_log_prob(dot_product(weights, summands));
  }
  
  // Log-prior for scale
  sigma ~ cauchy(0, prior_scale_for_dispersion);
  
  // Log-priors for coefficients
  if (prior_dist == 1) # normal
    theta ~ normal(prior_mean, prior_scale);  
  else # student_t
    theta ~ student_t(prior_df, prior_mean, prior_scale);
  
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
  vector[K + has_intercept] beta;
  if (has_intercept == 0)
    beta <- theta;
  else {
    beta[1] <- alpha0[1] - column_means(X) * theta;
    for (k in 1:K)
      beta[k + 1] <- theta[k];
  }
}
