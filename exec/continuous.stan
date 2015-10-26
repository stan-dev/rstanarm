# GLM for a Gaussian, Gamma, or inverse Gaussian outcome
data {
  # dimensions
  int<lower=1> N; # number of observations
  int<lower=0> K; # number of predictors
  
  # data
  vector[K] xbar;                # predictor means
  matrix[N,K]  X;                # centered predictor matrix
  vector[N]    y;                # continuous outcome
  int<lower=0,upper=1> prior_PD; # flag indicating whether to draw from the prior
  
  # intercept
  int<lower=0,upper=1> has_intercept; # 1 = yes
  
  # glmer stuff, see table 3 of
  # https://cran.r-project.org/web/packages/lme4/vignettes/lmer.pdf
  int<lower=0> t;               # num. terms (maybe 0) with a | in the glmer formula
  int<lower=1> p[t];            # num. variables on the LHS of each |
  int<lower=1> l[t];            # num. levels for the factor(s) on the RHS of each |
  int<lower=0> q;               # conceptually equals \sum_{i=1}^t p_i \times l_i
  int<lower=0> num_non_zero;    # number of non-zero elements in the Z matrix
  vector[num_non_zero] w;       # non-zero elements in the implicit Z matrix
  int<lower=0> v[num_non_zero]; # column indices for w
  int<lower=0> u[(N+1)*(t>0)];  # where the non-zeros start in each row
  int<lower=0> len_theta_L;     # length of the theta_L vector

  # family 
  int<lower=1,upper=3> family; # 1 = gaussian, 2 = Gamma, 3 = inverse Gaussian
  
  # link function from location to linear predictor
  int<lower=1,upper=4> link; # 1 = identity, 2 = log, 3 = inverse, 4 = 1 / mu^2
  
  # weights
  int<lower=0,upper=1> has_weights; # 0 = No, 1 = Yes
  vector[N * has_weights] weights;
  
  # offset
  int<lower=0,upper=1> has_offset;  # 0 = No, 1 = Yes
  vector[N * has_offset] offset;
  
  # prior family: 0 = none, 1 = normal, 2 = student_t, 3 = hs, 4 = hs_plus
  int<lower=0,upper=4> prior_dist;
  int<lower=0,upper=2> prior_dist_for_intercept;
  
  # hyperparameter values are set to 0 if there is no prior
  vector<lower=0>[K] prior_scale;
  real<lower=0> prior_scale_for_intercept;
  vector[K] prior_mean;
  real prior_mean_for_intercept;
  vector<lower=0>[K] prior_df;
  real<lower=0> prior_df_for_intercept;
  real<lower=0> prior_scale_for_dispersion;
  
  # hyperparameters for glmer stuff; if t > 0 priors are mandatory
  vector<lower=0>[t] shape; 
  vector<lower=0>[t] scale;
  int<lower=0> len_concentration;
  real<lower=0> concentration[len_concentration];
  int<lower=0> len_regularization;
  real<lower=0> regularization[len_regularization];
}
transformed data {
  int<lower=0> hs;
  int<lower=0> len_z_T;
  int<lower=0> len_var_group;
  int<lower=0> len_rho;
  real<lower=0> delta[len_concentration];
  int<lower=1> pos;
  vector[N * (family == 3)] sqrt_y;
  vector[N * (family == 3)] log_y;
  real sum_log_y;
  if      (prior_dist <= 2) hs <- 0;
  else if (prior_dist == 3) hs <- 2;
  else if (prior_dist == 4) hs <- 4;
  len_z_T <- 0;
  len_var_group <- sum(p) * (t > 0);
  len_rho <- sum(p) - t;
  pos <- 1;
  for (i in 1:t) {
    if (p[i] > 1) {
      for (j in 1:p[i]) {
        delta[pos] <- concentration[j];
        pos <- pos + 1;
      }
    }
    for (j in 3:p[i]) len_z_T <- len_z_T + p[i] - 1;
  }
  
  if (family == 1) sum_log_y <- not_a_number();
  else if (family == 2) sum_log_y <- sum(log(y));
  else {
    for (n in 1:N) sqrt_y[n] <- sqrt(y[n]);
    log_y <- log(y);
    sum_log_y <- sum(log_y);
  }
}
parameters {
  real<lower=if_else(family == 1 || link == 2, 
                     negative_infinity(), 0)> gamma[has_intercept];
  vector[K] z_beta;
  real<lower=0> global[hs];
  vector<lower=0>[K] local[hs];
  real<lower=0> dispersion_unscaled; # interpretation depends on family!
  vector[q] z_b;
  vector[len_z_T] z_T;
  vector<lower=0,upper=1>[len_rho] rho;
  vector<lower=0>[len_concentration] zeta;
  vector<lower=0>[t] tau;
}
transformed parameters {
  vector[K] beta;
  vector[q] b;
  real dispersion;
  vector[len_theta_L] theta_L;
  if (prior_dist == 0)      beta <- z_beta;
  else if (prior_dist <= 2) beta <- z_beta .* prior_scale + prior_mean;
  else if (prior_dist == 3) {
    vector[K] lambda;
    for (k in 1:K) lambda[k] <- local[1][k] * sqrt(local[2][k]);
    beta <- z_beta .* lambda * global[1]    * sqrt(global[2]);
  }
  else if (prior_dist == 4) {
    vector[K] lambda;
    vector[K] lambda_plus;
    for (k in 1:K) {
      lambda[k] <- local[1][k] * sqrt(local[2][k]);
      lambda_plus[k] <- local[3][k] * sqrt(local[4][k]);
    }
    beta <- z_beta .* lambda .* lambda_plus * global[1] * sqrt(global[2]);
  }
  if (prior_scale_for_dispersion > 0)
    dispersion <-  prior_scale_for_dispersion * dispersion_unscaled;
  else dispersion <- dispersion_unscaled;
  if (t > 0) {
    theta_L <- make_theta_L(len_theta_L, p, 
                            dispersion, tau, scale, zeta, rho, z_T);
    b <- make_b(z_b, theta_L, p, l);
  }
}
model {
  vector[N] eta; # linear predictor
  if (K > 0) eta <- X * beta;
  else eta <- rep_vector(0.0, N);
  if (has_offset == 1)    eta <- eta + offset;
  if (t > 0) eta <- eta + csr_matrix_times_vector(N, q, w, v, u, b);
  if (has_intercept == 1) {
    if (family == 1 || link == 2) eta <- eta + gamma[1];
    else eta <- eta - min(eta) + gamma[1];
  }
  
  // Log-likelihood 
  if (has_weights == 0 && prior_PD == 0) { # unweighted log-likelihoods
    if (family == 1) {
      if (link == 1)      y ~ normal(eta, dispersion);
      else if (link == 2) y ~ lognormal(eta, dispersion);
      else y ~ normal(divide_real_by_vector(1, eta), dispersion);
    }
    else if (family == 2) {
      y ~ GammaReg(eta, dispersion, link, sum_log_y);
    }
    else {
      y ~ inv_gaussian(linkinv_inv_gaussian(eta, link), 
                       dispersion, sum_log_y, sqrt_y);
    }
  }
  else if (prior_PD == 0) { # weighted log-likelihoods
    vector[N] summands;
    if (family == 1) summands <- pw_gauss(y, eta, dispersion, link);
    else if (family == 2) summands <- pw_gamma(y, eta, dispersion, link);
    else summands <- pw_inv_gaussian(y, eta, dispersion, link, log_y, sqrt_y);
    increment_log_prob(dot_product(weights, summands));
  }
  
  // Log-prior for scale
  if (prior_scale_for_dispersion > 0) dispersion_unscaled ~ cauchy(0, 1);
  
  // Log-priors for coefficients
  if      (prior_dist == 1) z_beta ~ normal(0, 1);
  else if (prior_dist == 2) z_beta ~ student_t(prior_df, 0, 1);
  else if (prior_dist == 3) { # hs
    z_beta ~ normal(0,1);
    local[1] ~ normal(0,1);
    local[2] ~ inv_gamma(0.5 * prior_df, 0.5 * prior_df);
    global[1] ~ normal(0,1);
    global[2] ~ inv_gamma(0.5, 0.5);
  }
  else if (prior_dist == 4) { # hs+
    z_beta ~ normal(0,1);
    local[1] ~ normal(0,1);
    local[2] ~ inv_gamma(0.5 * prior_df, 0.5 * prior_df);
    local[3] ~ normal(0,1);
    // unorthodox useage of prior_scale as another df hyperparameter
    local[4] ~ inv_gamma(0.5 * prior_scale, 0.5 * prior_scale);
    global[1] ~ normal(0,1);
    global[2] ~ inv_gamma(0.5, 0.5);
  }
  /* else prior_dist is 0 and nothing is added */
  
  // Log-prior for intercept  
  if (has_intercept == 1) {
    if (prior_dist_for_intercept == 1) # normal
      gamma ~ normal(prior_mean_for_intercept, prior_scale_for_intercept);
    else if (prior_dist_for_intercept == 2) # student_t
      gamma ~ student_t(prior_df_for_intercept, prior_mean_for_intercept, 
                        prior_scale_for_intercept);
    /* else prior_dist is 0 and nothing is added */
  }
  
  if (t > 0) {
    int pos_reg;
    int pos_rho;
    z_b ~ normal(0,1);
    z_T ~ normal(0,1);
    pos_reg <- 1;
    pos_rho <- 1;
    for (i in 1:t) if (p[i] > 1) {
      vector[p[i] - 1] shape1;
      vector[p[i] - 1] shape2;
      real nu;
      nu <- regularization[pos_reg] + 0.5 * (p[i] - 2);
      pos_reg <- pos_reg + 1;
      shape1[1] <- nu;
      shape2[1] <- nu;
      for (j in 2:(p[i]-1)) {
        nu <- nu - 0.5;
        shape1[j] <- 0.5 * j;
        shape2[j] <- nu;
      }
      segment(rho, pos_rho, p[i] - 1) ~ beta(shape1,shape2);
      pos_rho <- pos_rho + p[i] - 1;
    }
    zeta ~ gamma(delta, 1);
    tau ~ gamma(shape, 1);
  }
}
generated quantities {
  real alpha[has_intercept];
  real mean_PPD;
  mean_PPD <- 0;
  if (has_intercept == 1)
    alpha[1] <- gamma[1] - dot_product(xbar, beta);
  {
    vector[N] eta;
    if (K > 0) eta <- X * beta;
    else eta <- rep_vector(0.0, N);
    if (has_offset) eta <- eta + offset;
    if (t > 0) eta <- eta + csr_matrix_times_vector(N, q, w, v, u, b);
    if (has_intercept == 1) {
      if (family == 1 || link == 2) eta <- eta + gamma[1];
      else {
        real min_eta;
        min_eta <- min(eta);
        alpha[1] <- alpha[1] - min_eta;
        eta <- eta - min_eta + gamma[1];
      }
    }
    if (family == 1) {
      if (link > 1) eta <- linkinv_gauss(eta, link);
      for (n in 1:N) mean_PPD <- mean_PPD + normal_rng(eta[n], dispersion);
    }
    else if (family == 2) {
      if (link > 1) eta <- linkinv_gamma(eta, link);
      for (n in 1:N) mean_PPD <- mean_PPD + gamma_rng(dispersion, dispersion / eta[n]);
    }
    else if (family == 3) {
      if (link > 1) eta <- linkinv_inv_gaussian(eta, link);
      for (n in 1:N) mean_PPD <- mean_PPD + inv_gaussian_rng(eta[n], dispersion);
    }
    mean_PPD <- mean_PPD / N;
  }
}
