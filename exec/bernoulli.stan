# GLM for a Bernoulli outcome with optional Gaussian or t priors
functions {
  /** 
   * Apply inverse link function to linear predictor
   * see help(binom) in R
   *
   * @param eta Linear predictor vector
   * @param link An integer indicating the link function
   * @return A vector, i.e. inverse-link(eta)
   */
  vector linkinv_bern(vector eta, int link) {
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
   * @param link An integer indicating the link function
   * @param eta0 A vector of linear predictors | y = 0
   * @param eta1 A vector of linear predictors | y = 1
   * @param N An integer array of length 2 giving the number of 
   *   observations where y = 0 and y = 1 respectively
   * @return lp__
   */
  real ll_bern_lp(vector eta0, vector eta1, int link, int[] N) {
    if (link < 1 || link > 5) reject("Invalid link");
    if (link == 1) { // logit link
      0 ~ bernoulli_logit(eta0);
      1 ~ bernoulli_logit(eta1);
    }
    else if (link == 2) { // probit link
      increment_log_prob(normal_ccdf_log(eta0, 0, 1));
      increment_log_prob(normal_cdf_log(eta1, 0, 1));
    }
    else if (link == 3) { // cauchit link
      increment_log_prob(cauchy_ccdf_log(eta0, 0, 1));
      increment_log_prob(cauchy_cdf_log(eta1, 0, 1));
    }
    else if(link == 4) { // log link
      vector[N[1]]       log_pi0;
      for (n in 1:N[1])  log_pi0[n] <- log1m_exp(eta0[n]);
      increment_log_prob(log_pi0);
      increment_log_prob(eta1); # already in log form
    }
    else if(link == 5) { // cloglog link
      vector[N[2]]       log_pi1;
      for (n in 1:N[2])  log_pi1[n] <- log1m_exp(-exp(eta1[n]));
      increment_log_prob(log_pi1);
      increment_log_prob(-exp(eta0));
    }
    return get_lp();
  }

  /** 
   * Pointwise (pw) log-likelihood vector
   *
   * @param y The integer outcome variable. Note that function is
   *  called separately with y = 0 and y = 1
   * @param eta Vector of linear predictions
   * @param link An integer indicating the link function
   * @return A vector
   */
  vector pw_bern(int y, vector eta, int link) {
    vector[rows(eta)] ll;
    if (link < 1 || link > 5) 
      reject("Invalid link");
    if (link == 1) { # link = logit
      for (n in 1:rows(eta)) ll[n] <- bernoulli_logit_log(y, eta[n]);
    }
    else { # link = probit, cauchit, log, or cloglog 
           # Note: this may not be numerically stable
      vector[rows(eta)] pi;
      pi <- linkinv_bern(eta, link);
      for (n in 1:rows(eta)) ll[n] <- bernoulli_log(y, pi[n]) ;
    }
    return ll;
  }

  /** 
   * Upper bound on the intercept, which is infinity except for log link
   *
   * @param link An integer indicating the link function
   * @param X0 A matrix of predictors | y = 0
   * @param X1 A matrix of predictors | y = 1
   * @param beta A vector of coefficients
   * @param has_offset An integer indicating an offset
   * @param offset0 A vector of offsets | y = 0
   * @param offset1 A vector of offsets | y = 1
   * @return A scalar upper bound on the intercept
   */
  real make_upper_bernoulli(int link, matrix X0, matrix X1, 
                            vector beta, int has_offset, 
                            vector offset0, vector offset1) {
    real maximum;
    if (link != 4) return positive_infinity();
    if (has_offset == 0) maximum <- fmax( max(X0 * beta), max(X1 * beta) );
    else
      maximum <- fmax( max(X0 * beta + offset0), max(X1 * beta + offset1) );
      
    return -maximum;
  }
}
data {
  # dimensions
  int<lower=0> K;                # number of predictors
  int<lower=1> N[2];             # number of observations where y = 0 and y = 1 respectively
  vector[K] xbar;                # vector of column-means of rbind(X0, X1)
  matrix[N[1],K] X0;             # centered (by xbar) predictor matrix | y = 0
  matrix[N[2],K] X1;             # centered (by xbar) predictor matrix | y = 1
  int<lower=0,upper=1> prior_PD; # flag indicating whether to draw from the prior predictive
  
  # intercept
  int<lower=0,upper=1> has_intercept; # 0 = no, 1 = yes
  
  # link function from location to linear predictor
  int<lower=1,upper=5> link;
  
  # weights
  int<lower=0,upper=1> has_weights; # 0 = No, 1 = Yes
  vector[N[1] * has_weights] weights0;
  vector[N[2] * has_weights] weights1;
  
  # offset
  int<lower=0,upper=1> has_offset;  # 0 = No, 1 = Yes
  vector[N[1] * has_offset] offset0;
  vector[N[2] * has_offset] offset1;

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
transformed data {
  int NN;
  NN <- N[1] + N[2];
}
parameters {
  vector[K] z_beta;
  real<upper=if_else(link == 4, 0, positive_infinity())> gamma[has_intercept];
}
transformed parameters {
  vector[K] beta;
  if (prior_dist > 0) beta <- prior_mean + prior_scale .* z_beta;
  else beta <- z_beta;
}
model {
  vector[N[1]] eta0;
  vector[N[2]] eta1;
  if (K > 0) {
    eta0 <- X0 * beta;
    eta1 <- X1 * beta;
  }
  else {
    eta0 <- rep_vector(0.0, N[1]);
    eta1 <- rep_vector(0.0, N[2]);
  }
  if (has_offset == 1) {
    eta0 <- eta0 + offset0;
    eta1 <- eta1 + offset1;
  }
  if (has_intercept == 1) {
    if (link != 4) {
      eta0 <- gamma[1] + eta0;
      eta1 <- gamma[1] + eta1;
    }
    else {
      real shift;
      shift <- fmax(max(eta0), max(eta1));
      eta0 <- gamma[1] + eta0 - shift + 1;
      eta1 <- gamma[1] + eta1 - shift + 1;
    }
  }
  
  // Log-likelihood 
  if (has_weights == 0 && prior_PD == 0) { # unweighted log-likelihoods
    real dummy; # irrelevant but useful for testing
    dummy <- ll_bern_lp(eta0, eta1, link, N);
  }
  else if (prior_PD == 0) { # weighted log-likelihoods
    increment_log_prob(dot_product(weights0, pw_bern(0, eta0, link)));
    increment_log_prob(dot_product(weights1, pw_bern(1, eta1, link)));
  }
  
  // Log-priors for coefficients 
  if (prior_dist == 1) # normal
    z_beta ~ normal(0, 1);  
  else if (prior_dist == 2) # student_t
    z_beta ~ student_t(prior_df, 0, 1);
  /* else prior_dist is 0 and nothing is added */
   
  // Log-prior for intercept  
  if (has_intercept == 1) {
    if (prior_dist_for_intercept == 1) # normal
      gamma ~ normal(prior_mean_for_intercept, prior_scale_for_intercept);
    else if (prior_dist_for_intercept == 2) # student_t
      gamma ~ student_t(prior_df_for_intercept, prior_mean_for_intercept, 
                        prior_scale_for_intercept);
    /* else prior_dist = 0 and nothing is added */
  }
}
generated quantities {
  real alpha[has_intercept];
  real mean_PPD;
  if (has_intercept == 1) {
    alpha[1] <- gamma[1] - dot_product(xbar, beta);
  }
  mean_PPD <- 0;
  {
    vector[N[1]] eta0; 
    vector[N[2]] eta1;
    vector[N[1]] pi0;
    vector[N[2]] pi1;
    if (K > 0) {
      eta0 <- X0 * beta;
      eta1 <- X1 * beta;
    }
    else {
      eta0 <- rep_vector(0.0, N[1]);
      eta1 <- rep_vector(0.0, N[2]);
    }
    if (has_offset == 1) {
      eta0 <- eta0 + offset0;
      eta1 <- eta1 + offset1;
    }
    if (has_intercept == 1) {
      if (link != 4) {
        eta0 <- gamma[1] + eta0;
        eta1 <- gamma[1] + eta1;
      }      
      else {
        real shift;
        shift <- fmax(max(eta0), max(eta1));
        eta0 <- gamma[1] + eta0 - shift + 1;
        eta1 <- gamma[1] + eta1 - shift + 1;
        alpha[1] <- alpha[1] - shift;
      }
    }
    pi0 <- linkinv_bern(eta0, link);
    pi1 <- linkinv_bern(eta1, link);
    for (n in 1:N[1]) mean_PPD <- mean_PPD + bernoulli_rng(pi0[n]);
    for (n in 1:N[2]) mean_PPD <- mean_PPD + bernoulli_rng(pi1[n]);
    mean_PPD <- mean_PPD / NN;
  }
}
