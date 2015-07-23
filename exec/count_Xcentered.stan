# GLM for a count outcome with optional Gaussian or t priors
functions {
  vector linkinv_count(vector eta, int link) {
    vector[rows(eta)] phi;
    if (link < 1 || link > 3) reject("Invalid link");
    if      (link == 1) return(exp(eta));
    else if (link == 2) return(eta); # link = identity
    else  # link = sqrt
      for(n in 1:rows(eta)) phi[n] <- square(eta[n]); 
    return phi;
  }
  
  /** 
   * Pointwise (pw) log-likelihood vector for the Poisson distribution
   *
   * @param y The integer array corresponding to the outcome variable.
   * @param link An integer indicating the link function
   * @return A vector
   */
  vector pw_pois(int[] y, vector eta, int link) {
    vector[rows(eta)] ll;
    if (link < 1 || link > 3) reject("Invalid link");
    if (link == 1) # link = log
      for (n in 1:rows(eta)) ll[n] <- poisson_log_log(y[n], eta[n]);
    else { # link = identity or sqrt
      vector[rows(eta)] phi;
      phi <- linkinv_count(eta, link);
      for (n in 1:rows(eta)) ll[n] <- poisson_log(y[n], phi[n]) ;
    }
    return ll;
  }

  /** 
   * Pointwise (pw) log-likelihood vector for the negative binomial  distribution
   *
   * @param y The integer array corresponding to the outcome variable.
   * @param link An integer indicating the link function
   * @return A vector
   */
  vector pw_nb(int[] y, vector eta, real theta, int link) {
    vector[rows(eta)] ll;
    vector[rows(eta)] rho;
    if (link < 1 || link > 3) reject("Invalid link");
    rho <- linkinv_count(eta, link);
    for (n in 1:rows(eta)) ll[n] <- neg_binomial_log(y[n], rho[n], theta);
    return ll;
  }
  
  /** 
   * Lower bound on the intercept, which is negative infinity 
   * except for identity link
   *
   * @param link An integer indicating the link function
   * @param X A matrix of predictors | y = 0
   * @param beta A vector of coefficients
   * @param has_offset An integer indicating an offset
   * @param offset A vector of offsets
   * @return A scalar lower bound on the intercept
   */
  real make_lower_count(int link, matrix X, vector beta, 
                        int has_offset, vector offset) {
    real minimum;
    if (link != 2) return negative_infinity();
    if (has_offset == 0) minimum <- min(X * beta);
    else minimum <- min(X * beta + offset);
    return -minimum;
  }

  /** 
   * Elementwise square root
   *
   * @param y A vector of non-negative numbers
   * @return A vector of square roots
   */
  vector sqrt_vec(vector y) {
    vector[rows(y)] out;
    for (i in 1:rows(y)) out[i] <- sqrt(out[i]);
    return(out);
  }
  
}
data {
  # dimensions
  int<lower=1> N; # number of observations
  int<lower=1> K; # number of predictors
  
  # data
  vector[K] xbar;    # predictor means
  matrix[N,K]  X;    # centered predictor matrix
  int<lower=0> y[N]; # count outcome
  
  # intercept
  int<lower=0,upper=1> has_intercept; # 0 = no, 1 = yes
  
  # likelihood
  int<lower=1> family; # 1 = poisson, 2 = negative binomial, 3 poisson mixture
  
  # link function from location to linear predictor
  int<lower=1,upper=3> link; # 1 = log, 2 = identity, 3 = sqrt
  /* NOTE: MASS::negative.binomial switches 2 and 3 but we follow poisson */
  
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
  real<lower=0> prior_scale_for_dispersion;
}
parameters {
  vector[K] beta; # coefficients
  real<lower=make_lower_count(link, X, beta, 
                              has_offset, offset)> gamma[has_intercept];
  real<lower=0> theta_unscaled[family > 1];
  vector<lower=0>[N] noise[family == 3]; // do not store this
}
transformed parameters {
  real theta[family > 1];
  if (family > 1) theta[1] <- prior_scale_for_dispersion * 
                              theta_unscaled[1];
}
model {
  vector[N] eta; # linear predictor
  eta <- X * beta;
  if (has_intercept == 1) eta <- eta + gamma[1];
  if (has_offset == 1)    eta <- eta + offset;
  if (family == 3) {
    if      (link == 1) eta <- eta + log(theta[1]) + log(noise[1]);
    else if (link == 2) eta <- eta * theta[1] .* noise[1];
    else                eta <- eta + sqrt(theta[1]) + sqrt_vec(noise[1]);
  }
  
  // Log-likelihood 
  if (has_weights == 0) { # unweighted log-likelihoods
    if(family != 2) {
      if (link == 1) y ~ poisson_log(eta);
      else y ~ poisson(linkinv_count(eta, link));
    }
    else y ~ neg_binomial(linkinv_count(eta, link), theta[1]);
  }
  else if (family != 1)
    increment_log_prob(dot_product(weights, pw_pois(y, eta, link)));
  else
    increment_log_prob(dot_product(weights, pw_nb(y, eta, theta[1], link)));
  
  // Log-priors for coefficients 
  if (prior_dist == 1) # normal
    beta ~ normal(prior_mean, prior_scale);  
  else if (prior_dist == 2) # student_t
    beta ~ student_t(prior_df, prior_mean, prior_scale);
  /* else prior_dist = 0 and nothing is added */
   
  // Log-prior for intercept  
  if (has_intercept == 1) {
    if (link != 2) {
      if (prior_dist_for_intercept == 1) # normal
        gamma ~ normal(prior_mean_for_intercept, prior_scale_for_intercept);
      else if (prior_dist_for_intercept == 2) # student_t
        gamma ~ student_t(prior_df_for_intercept, prior_mean_for_intercept, 
                          prior_scale_for_intercept);
      /* else prior_dist = 0 and nothing is added */
    }
    else {
      real minimum;
      minimum <- -min(eta);
      if (prior_dist_for_intercept == 1) # normal
        gamma[1] ~ normal(prior_mean_for_intercept, 
                          prior_scale_for_intercept) T[minimum,];
      else if (prior_dist_for_intercept == 2) # student_t
        gamma[1] ~ student_t(prior_df_for_intercept, prior_mean_for_intercept, 
                             prior_scale_for_intercept) T[minimum,];
      /* else prior_dist = 0 and nothing is added */
    }
  }
  
  // Log-prior for dispersion
  if (family > 1) theta_unscaled ~ cauchy(0, 1);

  // Log-prior for noise
  if (family == 3) noise[1] ~ gamma(theta[1], 1);
}
generated quantities {
  real alpha[has_intercept];
  real mean_PPD;
  if (has_intercept == 1) alpha[1] <- gamma[1] - dot_product(xbar, beta);
  mean_PPD <- 0;
  {
    vector[N] eta;
    vector[N] rho;
    eta <- X * beta;
    if (has_intercept == 1) eta <- eta + gamma[1];
    if (has_offset == 1)    eta <- eta + offset;
    if (family == 3) {
      if      (link == 1) eta <- eta + log(theta[1]) + log(noise[1]);
      else if (link == 2) eta <- eta * theta[1] .* noise[1];
      else                eta <- eta + sqrt(theta[1]) + sqrt_vec(noise[1]);
    }
    rho <- linkinv_count(eta, link);
    if (family != 2) 
      for (n in 1:N) mean_PPD <- mean_PPD + poisson_rng(rho[n]);
    else
      for (n in 1:N) mean_PPD <- mean_PPD + neg_binomial_rng(rho[n], theta[1]);
    mean_PPD <- mean_PPD / N;
  }
}
