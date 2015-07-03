# GLM for a discrete outcome with Gaussian or t priors
functions {
  /** 
   * Combine intercept and coefficients into a vector or return coefficients
   * if no intercept. 
   *
   * @param alpha Array of size 0 or 1 depending on if model has intercept
   * @param theta Vector of coefficients (not including intercept) of size 
   *              K-1 or K depending on if model has intercept
   * @param K Number of columns of the predictor matrix X 
   * @return If the model has an intercept then the function combines alpha
   *         and theta into a vector. If the model doesn't have 
   *         an intercept then theta itself is returned.
   */
  vector coefficient_vector(real[] alpha, vector theta, int K) {
    vector[K] beta;
    int S;
    S <- size(alpha);
    if (S == 0) {
      if (K != rows(theta)) {
        reject("Dimension mismatch");
      }
      return theta;
    }
    else 
      if (K != 1 + rows(theta)) {
        reject("Dimension mismatch");
      }
      beta[1] <- alpha[1];
      for (k in 2:K) 
        beta[k] <- theta[k-1];
      return beta;
  }
  
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
    if (link == 2) phi <- eta; # link = identity
    else { # link = sqrt
      for(n in 1:rows(eta)) 
        phi[n] <- square(eta[n]); 
    }
    return phi;
  }
  
  /** 
   * Pointwise (pw) log-likelihood vector
   *
   * @param y The integer array corresponding to the outcome variable.
   * @param link An integer indicating the link function
   * @return A vector
   */
  vector pw_binom(int[] y, vector eta, int link) {
    vector[rows(eta)] ll;
    vector[rows(eta)] pi;
    if (link < 1 || link > 5) reject("Invalid link");
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
    if (link == 1) { # link = log
      for (n in 1:rows(eta)) ll[n] <- poisson_log_log(y[n], eta[n]);
    }
    else { # link = identity or sqrt
      phi <- linkinv_pois(eta, link);
      for (n in 1:rows(eta)) ll[n] <- poisson_log(y[n], phi[n]) ;
    }
    return ll;
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
  vector<lower=0>[K - has_intercept] prior_scale;
  real<lower=0> prior_scale_for_intercept;
  vector[K - has_intercept] prior_mean;
  real prior_mean_for_intercept;
  vector<lower=0>[K - has_intercept] prior_df;
  real<lower=0> prior_df_for_intercept;
}
parameters {
  real alpha[has_intercept];
  vector[K - has_intercept] theta;
}
model {
  vector[N] eta;
  vector[K] beta;
  beta <- coefficient_vector(alpha, theta, K);  
  eta <- X * beta;
  if (has_offset == 1) eta <- eta + offset;
  
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
      summands <- pw_binom(y, eta, link);
    else # poisson case
      summands <- pw_pois(y, eta, link);
    increment_log_prob(dot_product(weights, summands));
  }
  
  // Log-priors for coefficients 
  if (prior_dist == 1) # normal
    theta ~ normal(prior_mean, prior_scale);  
  else # student_t
    theta ~ student_t(prior_df, prior_mean, prior_scale);
   
  // Log-prior for intercept  
  if (has_intercept == 1) {
    if (prior_dist_for_intercept == 1) # normal
      alpha ~ normal(prior_mean_for_intercept, prior_scale_for_intercept);
    else # student_t
      alpha ~ student_t(prior_df_for_intercept, prior_mean_for_intercept, 
                        prior_scale_for_intercept);
  }
}
generated quantities {
  vector[K] beta;
  beta <- coefficient_vector(alpha, theta, K);  
}