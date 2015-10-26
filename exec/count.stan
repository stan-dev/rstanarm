# GLM for a count outcome
data {
  # dimensions
  int<lower=1> N; # number of observations
  int<lower=0> K; # number of predictors
  
  # data
  vector[K] xbar;                # predictor means
  matrix[N,K]  X;                # centered predictor matrix
  int<lower=0> y[N];             # count outcome
  int<lower=0,upper=1> prior_PD; # flag indicating whether to draw from the prior predictive
  
  # intercept
  int<lower=0,upper=1> has_intercept; # 0 = no, 1 = yes
  
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
  
  # prior family: 0 = none, 1 = normal, 2 = student_t, 3 = hs, 4 = hs_plus
  int<lower=0,upper=4> prior_dist;
  int<lower=0,upper=2> prior_dist_for_intercept;
  
  # hyperparameter values
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
transformed data{
  real poisson_max;
  int<lower=0> hs;
  int<lower=0> len_z_T;
  int<lower=0> len_var_group;
  int<lower=0> len_rho;
  real<lower=0> delta[len_concentration];
  int<lower=1> pos;
  
  poisson_max <- pow(2.0, 30.0);
  if (prior_dist <=  2)     hs <- 0;
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
}
parameters {
  real<lower=0> theta_unscaled[family > 1];
  vector<lower=0>[N] noise[family == 3]; // do not store this
  vector[K] z_beta;
  real<lower=if_else(link == 1, negative_infinity(), 0)> gamma[has_intercept];
  real<lower=0> global[hs];
  vector<lower=0>[K] local[hs];
  real<lower=0> dispersion_unscaled[family > 1];
  vector[q] z_b;
  vector[len_z_T] z_T;
  vector<lower=0,upper=1>[len_rho] rho;
  vector<lower=0>[len_concentration] zeta;
  vector<lower=0>[t] tau; 
}
transformed parameters {
  real theta[family > 1];
  vector[K] beta;
  vector[q] b;
  vector[len_theta_L] theta_L;
  real dispersion[family > 1];
  
  if (family > 1 && prior_scale_for_dispersion > 0) 
    theta[1] <- prior_scale_for_dispersion * theta_unscaled[1];
  else if (family > 1) theta[1] <- theta_unscaled[1];
  
  if (prior_dist == 0) beta <- z_beta;
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
  if (t > 0) {
    theta_L <- make_theta_L(len_theta_L, p, if_else(family > 1, theta[1], 1),
                            tau, scale, zeta, rho, z_T);
    b <- make_b(z_b, theta_L, p, l);
  }
}
model {
  vector[N] eta; # linear predictor
  if (K > 0) eta <- X * beta;
  else eta <- rep_vector(0.0, N);
  if (has_offset == 1) eta <- eta + offset;
  if (t > 0) eta <- eta + csr_matrix_times_vector(N, q, w, v, u, b);  
  if (has_intercept == 1) {
    if (link == 1) eta <- eta + gamma[1];
    else eta <- eta - min(eta) + gamma[1];
  }
  if (family == 3) {
    if      (link == 1) eta <- eta + log(theta[1]) + log(noise[1]);
    else if (link == 2) eta <- eta * theta[1] .* noise[1];
    else                eta <- eta + sqrt(theta[1]) + sqrt_vec(noise[1]);
  }
  
  // Log-likelihood 
  if (has_weights == 0 && prior_PD == 0) { # unweighted log-likelihoods
    if(family != 2) {
      if (link == 1) y ~ poisson_log(eta);
      else y ~ poisson(linkinv_count(eta, link));
    }
    else {
      if (link == 1) y ~ neg_binomial_2_log(eta, theta[1]);
      else y ~ neg_binomial_2(linkinv_count(eta, link), theta[1]);
    }
  }
  else if (family != 1 && prior_PD == 0)
    increment_log_prob(dot_product(weights, pw_pois(y, eta, link)));
  else if (prior_PD == 0)
    increment_log_prob(dot_product(weights, pw_nb(y, eta, theta[1], link)));
  
  // Log-priors for coefficients
  if (prior_dist == 1) z_beta ~ normal(0, 1);
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
    /* else prior_dist = 0 and nothing is added */
  }
  
  // Log-prior for dispersion
  if (family > 1 && prior_scale_for_dispersion > 0) theta_unscaled ~ cauchy(0, 1);

  // Log-prior for noise
  if (family == 3) noise[1] ~ gamma(theta[1], 1);
  
  if (t > 0) {
    int pos_reg;
    int pos_rho;
    if (family > 1 && prior_scale_for_dispersion > 0) 
      dispersion_unscaled ~ cauchy(0, 1);
    z_b ~ normal(0,1);
    z_T ~ normal(0,1);
    pos_reg <- 1;
    pos_rho <- 1;
    for (i in 1:t) if (p[i] > 1) {
      vector[p[i] - 1] shape1;
      vector[p[i] - 1] shape2;
      real nu;
      nu <- shape[pos_reg] + 0.5 * (p[i] - 2);
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
  if (has_intercept == 1) alpha[1] <- gamma[1] - dot_product(xbar, beta);
  mean_PPD <- 0;
  {
    vector[N] eta;
    vector[N] nu;
    if (K > 0) eta <- X * beta;
    else eta <- rep_vector(0.0, N);
    if (has_offset == 1) eta <- eta + offset;
    if (t > 0) eta <- eta + csr_matrix_times_vector(N, q, w, v, u, b);
    if (has_intercept == 1) {
      if (link == 1) eta <- eta + gamma[1];
      else {
        real shift;
        shift <- min(eta);
        eta <- eta - shift + gamma[1];
        alpha[1] <- alpha[1] - shift;
      }
    }
    if (family == 3) {
      if      (link == 1) eta <- eta + log(theta[1]) + log(noise[1]);
      else if (link == 2) eta <- eta * theta[1] .* noise[1];
      else                eta <- eta + sqrt(theta[1]) + sqrt_vec(noise[1]);
    }
    nu <- linkinv_count(eta, link);
    if (family != 2) for (n in 1:N) {
        if (nu[n] < poisson_max) mean_PPD <- mean_PPD + poisson_rng(nu[n]);
        else mean_PPD <- mean_PPD + normal_rng(nu[n], sqrt(nu[n]));
    }
    else for (n in 1:N) {
        real gamma_temp;
        if (is_inf(theta[1])) gamma_temp <- nu[n];
        else gamma_temp <- gamma_rng(theta[1], theta[1] / nu[n]);
        if (gamma_temp < poisson_max)
          mean_PPD <- mean_PPD + poisson_rng(gamma_temp);
        else mean_PPD <- mean_PPD + normal_rng(gamma_temp, sqrt(gamma_temp));
    }
    mean_PPD <- mean_PPD / N;
  }
}
