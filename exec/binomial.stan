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
  
  /** 
   * Create group-specific block-diagonal Cholesky factor, see section 2.3 of
   * https://cran.r-project.org/web/packages/lme4/vignettes/lmer.pdf
   * @param p An integer array with the number variables on the LHS of each |
   * @param tau Vector of scale parameters for the decomposed covariance matrices
   * @param scale Vector of scale hyperparameters
   * @param zeta Vector of positive parameters that are normalized into simplexes
   * @param rho Vector of radii in the onion method for creating Cholesky factors
   * @param z_T Vector used in the onion method for creating Cholesky factors
   * @return A vector that corresponds to theta in lme4
   */
  vector make_theta_L_binom(int len_theta_L, int[] p,
                            vector tau, vector scale, vector zeta,
                            vector rho, vector z_T) {
    vector[len_theta_L] theta_L;
    int zeta_mark;
    int z_T_mark;
    int rho_mark;
    int theta_L_mark;
    zeta_mark <- 1;
    z_T_mark <- 1;
    rho_mark <- 1;
    theta_L_mark <- 1;
    
    // each of these is a diagonal block of the implicit Cholesky factor
    for (i in 1:size(p)) { 
      int nc;
      nc <- p[i];
      if (nc == 1) { // "block" is just a standard deviation
        theta_L[theta_L_mark] <- tau[i] * scale[i];
        theta_L_mark <- theta_L_mark + 1;
      }
      else { // block is lower-triangular               
        matrix[nc,nc] T_i; 
        real trace_T_i;
        vector[nc] pi; // variance = proportion of trace_T_i
        real std_dev;
        real T21;
        
        trace_T_i <- square(tau[i] * scale[i]) * nc;
        pi <- segment(zeta, zeta_mark, nc); // zeta ~ gamma(shape, 1)
        pi <- pi / sum(pi);                 // thus pi ~ dirichlet(shape)
        zeta_mark <- zeta_mark + nc;
        std_dev <- sqrt(pi[1] * trace_T_i);
        T_i[1,1] <- std_dev;

        // Put a correlation into T_i[2,1] and scale by std_dev
        std_dev <- sqrt(pi[2] * trace_T_i);
        T21 <- 2.0 * rho[rho_mark] - 1.0;
        rho_mark <- rho_mark + 1;
        T_i[2,2] <- std_dev * sqrt(1.0 - square(T21));
        T_i[2,1] <- std_dev * T21;
        
        for (r in 2:(nc - 1)) { // scaled onion method
          int rp1;
          vector[r] T_row;
          real scale_factor;
          T_row <- segment(z_T, z_T_mark, r);
          z_T_mark <- z_T_mark + r;
          rp1 <- r + 1;
          std_dev <- sqrt(pi[rp1] * trace_T_i);
          scale_factor <- sqrt(rho[rho_mark] / dot_self(T_row)) * std_dev;
          for(c in 1:r) T_i[rp1,c] <- T_row[c] * scale_factor;
          T_i[rp1,rp1] <- sqrt(1.0 - rho[rho_mark]) * std_dev;
          rho_mark <- rho_mark + 1;
        }
        
        // vec T_i
        for (c in 1:nc) for (r in c:nc) {
          theta_L[theta_L_mark] <- T_i[r,c];
          theta_L_mark <- theta_L_mark + 1;
        }
      }
    }
    return theta_L;
  }
  
  /** 
   * Create group-specific coefficients, see section 2.3 of
   * https://cran.r-project.org/web/packages/lme4/vignettes/lmer.pdf
   *
   * @param z_b Vector whose elements are iid normal(0,sigma) a priori
   * @param theta Vector with covariance parameters
   * @param p An integer array with the number variables on the LHS of each |
   * @param l An integer array with the number of levels for the factor(s) on 
   *   the RHS of each |
   * @return A vector of group-specific coefficients
   */
  vector make_b_binom(vector z_b, vector theta_L, int[] p, int[] l) { 
    vector[rows(z_b)] b;
    int b_mark;
    int theta_L_mark;
    b_mark <- 1;
    theta_L_mark <- 1;
    for (i in 1:size(p)) {
      int nc;
      nc <- p[i];
      if (nc == 1) {
        real theta_L_start;
        theta_L_start <- theta_L[theta_L_mark]; // needs to be positive
        for (s in b_mark:(b_mark + l[i] - 1)) 
          b[s] <- theta_L_start * z_b[s];
        b_mark <- b_mark + l[i];
        theta_L_mark <- theta_L_mark + 1;
      }
      else {
        matrix[nc,nc] T_i;
        T_i <- rep_matrix(0, nc, nc);
        for (c in 1:nc) {
          T_i[c,c] <- theta_L[theta_L_mark];    // needs to be positive
          theta_L_mark <- theta_L_mark + 1;
          for(r in (c+1):nc) {
            T_i[r,c] <- theta_L[theta_L_mark];
            theta_L_mark <- theta_L_mark + 1;
          }
        }
        for (j in 1:l[i]) {
          vector[nc] temp;
          temp <- T_i * segment(z_b, b_mark, nc);
          b_mark <- b_mark - 1;
          for (s in 1:nc) b[b_mark + s] <- temp[s];
          b_mark <- b_mark + nc + 1;
        }
      }
    }
    return b;
  }
  
}
data {
  # dimensions
  int<lower=1> N; # number of observations
  int<lower=0> K; # number of predictors
  
  # data
  vector[K] xbar;                # predictor means
  matrix[N,K]  X;                # centered predictor matrix
  int<lower=0> y[N];             # outcome: number of successes
  int<lower=0> trials[N];        # number of trials
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

  # link function from location to linear predictor
  int<lower=1,upper=5> link;
  
  # weights
  int<lower=0,upper=1> has_weights; # 0 = No, 1 = Yes
  vector[N * has_weights] weights;
  
  # offset
  int<lower=0,upper=1> has_offset;  # 0 = No, 1 = Yes
  vector[N * has_offset] offset;
  
  # prior family: 0 = none, 1 = normal, 2 = student_t, 3 = horseshoe, 4 = horseshoe_plus
  int<lower=0,upper=4> prior_dist;
  int<lower=0,upper=2> prior_dist_for_intercept;
  
  # hyperparameter values
  vector<lower=0>[K] prior_scale;
  real<lower=0> prior_scale_for_intercept;
  vector[K] prior_mean;
  real prior_mean_for_intercept;
  vector<lower=0>[K] prior_df;
  real<lower=0> prior_df_for_intercept;

  # hyperparameters for glmer stuff; if t > 0 priors are mandatory
  vector<lower=0>[t] gamma_shape; 
  vector<lower=0>[t] scale;
  int<lower=0> len_concentration;
  real<lower=0> concentration[len_concentration];
  int<lower=0> len_shape;
  real<lower=0> shape[len_shape];  
}
transformed data {
  int<lower=0> horseshoe;
  int<lower=0> len_z_T;
  int<lower=0> len_var_group;
  int<lower=0> len_rho;
  real<lower=0> delta[len_concentration];
  int<lower=1> pos;

  if (prior_dist <  2) horseshoe <- 0;
  else if (prior_dist == 3) horseshoe <- 2;
  else if (prior_dist == 4) horseshoe <- 4;
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
  vector[K] z_beta;
  real<upper=if_else(link == 4, 0, positive_infinity())> gamma[has_intercept];
  real<lower=0> global[horseshoe];
  vector<lower=0>[K] local[horseshoe];
  vector[q] z_b;
  vector[len_z_T] z_T;
  vector<lower=0,upper=1>[len_rho] rho;
  vector<lower=0>[len_concentration] zeta;
  vector<lower=0>[t] tau;  
}
transformed parameters {
  vector[K] beta;
  vector[q] b;
  vector[len_theta_L] theta_L;
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
    theta_L <- make_theta_L_binom(len_theta_L, p, 
                                  tau, scale, zeta, rho, z_T);
    b <- make_b_binom(z_b, theta_L, p, l);
  }
}
model {
  vector[N] eta; # linear predictor
  if (K > 0) eta <- X * beta;
  else eta <- rep_vector(0.0, N);
  if (has_offset == 1) eta <- eta + offset;
  if (t > 0) eta <- eta + csr_matrix_times_vector(N, q, w, v, u, b);
  if (has_intercept == 1) {
    if (link != 4) eta <- eta + gamma[1];
    else eta <- gamma[1] + eta - max(eta);
  }
  
  // Log-likelihood 
  if (has_weights == 0 && prior_PD == 0) { # unweighted log-likelihoods
    real dummy; # irrelevant but useful for testing
    dummy <- ll_binom_lp(y, trials, eta, link);
  }
  else if (prior_PD == 0) 
    increment_log_prob(dot_product(weights, pw_binom(y, trials, eta, link)));
  
  // Log-priors for coefficients
  if      (prior_dist == 1) z_beta ~ normal(0, 1);
  else if (prior_dist == 2) z_beta ~ student_t(prior_df, 0, 1);
  else if (prior_dist == 3) { # horseshoe
    z_beta ~ normal(0,1);
    local[1] ~ normal(0,1);
    local[2] ~ inv_gamma(0.5 * prior_df, 0.5 * prior_df);
    global[1] ~ normal(0,1);
    global[2] ~ inv_gamma(0.5, 0.5);
  }
  else if (prior_dist == 4) { # horseshoe+
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
  
  if (t > 0) {
    int pos_shape;
    int pos_rho;
    z_b ~ normal(0,1);
    z_T ~ normal(0,1);
    pos_shape <- 1;
    pos_rho <- 1;
    for (i in 1:t) if (p[i] > 1) {
      vector[p[i] - 1] shape1;
      vector[p[i] - 1] shape2;
      real nu;
      nu <- shape[pos_shape] + 0.5 * (p[i] - 2);
      pos_shape <- pos_shape + 1;
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
    tau ~ gamma(gamma_shape, 1);
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
    if (K > 0) eta <- X * beta;
    else eta <- rep_vector(0.0, N);
    if (has_offset == 1) eta <- eta + offset;
    if (t > 0) eta <- eta + csr_matrix_times_vector(N, q, w, v, u, b);
    if (has_intercept == 1) {
      if (link != 4) eta <- eta + gamma[1];
      else {
        real shift;
        shift <- max(eta);
        eta <- gamma[1] + eta - shift;
        alpha[1] <- alpha[1] - shift;
      }
    }
    pi <- linkinv_binom(eta, link);
    for (n in 1:N) mean_PPD <- mean_PPD + binomial_rng(trials[n], pi[n]);
    mean_PPD <- mean_PPD / N;
  }
}
