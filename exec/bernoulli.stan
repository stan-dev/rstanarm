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
  
  /** 
   * Create group-specific coefficients, see section 2.3 of
   * https://cran.r-project.org/web/packages/lme4/vignettes/lmer.pdf
   *
   * The group-specific coefficients have a centered multivariate normal prior.
   * To apply the Matt trick, we need Cholesky factors of covariance matrices.
   * Due to the lack of ragged arrays in Stan, this is a bit tedious.
   * We decompose a covariance matrix into a correlation matrix and variances.
   * We represent the variances as a scaled simplex, where the scale component
   * consists of the square root of the number of variables and an unknown.
   * The Cholesky factor of a correlation matrix can be built up via the onion
   * method, which inputs standard normally distributed variables and 
   * beta-distributed variables.
   *
   * @param u Vector whose elements are iid normal(0,1) a priori
   * @param z_T Vector used in the onion method for creating Cholesky factors
   * @param rho Vector of radii in the onion method for creating Cholesky factors
   * @param pi Vector of simplexes concatenated together
   * @param tau Vector of scale parameters
   * @param p An integer array with the number variables on the LHS of each |
   * @param l An integer array with the number of levels for the factor(s) on 
   *   the RHS of each |
   * @return A vector of group-specific coefficients
   */
  vector make_b(vector u, vector z_T, vector rho, vector var_group, 
                int[] p, int[] l) {
    vector[rows(u)] b;
    int b_mark;
    int z_T_mark;
    int rho_mark;
    int vg_mark;
    // Due to lack of ragged arrays, everything is input as long vectors
    b_mark   <- 1;
    z_T_mark <- 1;
    rho_mark <- 1;
    vg_mark <- 1;
    for (i in 1:size(p)) {
      int nc;
      nc <- p[i];
      if (nc == 1) { // just a standard deviation times a part of a vector
        real vg_i;
        vg_i <- var_group[vg_mark];
        for (s in b_mark:(b_mark + l[i] - 1)) 
          b[s] <- vg_i * u[s];
        vg_mark <- vg_mark + 1;  
        b_mark <- b_mark + l[i];
      }
      else {         // deal with a (Cholesky factor of a) covariance matrix
        matrix[nc,nc] T_i;
        real std_dev;
        T_i <- rep_matrix(0, nc, nc);
        std_dev <- sqrt(var_group[vg_mark]);
        vg_mark <- vg_mark + 1;
        T_i[1,1] <- std_dev;
        std_dev <- sqrt(var_group[vg_mark]);
        vg_mark <- vg_mark + 1;
        T_i[2,1] <- std_dev * (2.0 * rho[rho_mark] - 1.0);
        rho_mark <- rho_mark + 1;
        T_i[2,2] <- std_dev * sqrt(1.0 - square(T_i[2,1]));
        for (r in 2:(nc - 1)) { // modified onion method
          int rp1;
          vector[r] T_row;
          real scale;
          T_row <- segment(z_T, z_T_mark, r);
          z_T_mark <- z_T_mark + r;
          scale <- sqrt(rho[rho_mark] / dot_self(T_row));
          std_dev <- sqrt(var_group[vg_mark]);
          rp1 <- r + 1;
          for(c in 1:r) T_i[rp1,c] <- T_row[c] * scale * std_dev;
          T_i[rp1,rp1] <- sqrt(1.0 - rho[rho_mark]);
          rho_mark <- rho_mark + 1;
          vg_mark <- vg_mark + 1;
        }
        for (j in 1:l[i]) { // multiply Cholesky factor by relevant parts of u
          vector[nc] temp;
          temp <- T_i * segment(u, b_mark, nc);
          for (s in 1:nc) b[b_mark + s - 1] <- temp[s];
          b_mark <- b_mark + nc;
        }
      }
    }
    return b;
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
  
  # glmer stuff, see table 3 of
  # https://cran.r-project.org/web/packages/lme4/vignettes/lmer.pdf
  int<lower=0> t;    # num. terms (maybe 0) with a | in the glmer formula
  int<lower=1> p[t]; # num. variables on the LHS of each |
  int<lower=1> l[t]; # num. levels for the factor(s) on the RHS of each |
  int<lower=0> q;    # conceptually equals \sum_{i=1}^t p_i \times l_i
  matrix[N[1],q]  Z0; # uncentered design matrix for group-specific variables
  matrix[N[2],q]  Z1; # uncentered design matrix for group-specific variables

  
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
  real<lower=0> prior_scale_for_dispersion;
  
  # hyperparameters for glmer stuff; if t > 0 priors are mandatory
  vector<lower=0>[t] gamma_shape; 
  vector<lower=0>[t] scale;
  int<lower=0> len_concentration;
  real<lower=0> concentration[len_concentration];
  int<lower=0> len_shape;
  real<lower=0> shape[len_shape];
}
transformed data {
  int NN;
  int<lower=0> horseshoe;
  int<lower=0> len_z_T;
  int<lower=0> len_var_group;
  int<lower=0> len_rho;
  real<lower=0> shape1[len_shape];
  real<lower=0> shape2[len_shape];
  real<lower=0> delta[len_concentration];
  int<lower=1> pos[2];
  
  NN <- N[1] + N[2];
  if (prior_dist <  2) horseshoe <- 0;
  else if (prior_dist == 3) horseshoe <- 2;
  else if (prior_dist == 4) horseshoe <- 4;
  len_z_T <- 0;
  len_var_group <- sum(p) * (t > 0);
  len_rho <- sum(p) - t;
  pos[1] <- 1;
  pos[2] <- 1;
  for (i in 1:t) {
    real nu;
    if (p[i] > 1) {
      for (j in 1:p[i]) {
        delta[pos[2]] <- concentration[j];
        pos[2] <- pos[2] + 1;
      }
      nu <- shape[pos[1]] + 0.5 * (p[i] - 2);
      shape1[pos[1]] <- nu;
      shape2[pos[1]] <- nu;
      pos[1] <- pos[1] + 1;
    }
    if (p[i] > 2) for (j in 2:p[i]) {
      nu <- nu - 0.5;
      shape1[pos[1]] <- 0.5 * j;
      shape2[pos[1]] <- nu;
      pos[1] <- pos[1] + 1;
    }
    if (p[i] > 2) for (j in 3:p[i]) {
      len_z_T <- len_z_T + p[i] - 1;
    }
  }
}
parameters {
  vector[K] z_beta;
  real<upper=if_else(link == 4, 0, positive_infinity())> gamma[has_intercept];
  real<lower=0> global[horseshoe];
  vector<lower=0>[K] local[horseshoe];
  real<lower=0> dispersion_unscaled[t > 0]; # interpretation depends on family!
  vector[q] u;
  vector[len_z_T] z_T;
  vector<lower=0,upper=1>[len_rho] rho;
  vector<lower=0>[len_concentration] zeta;
  vector<lower=0>[t] tau;
}
transformed parameters {
  vector[K] beta;
  vector[q] b;
  vector[len_var_group] var_group;
  real dispersion[t > 0];
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
    vector[t] scaled_tau;
    int mark;
    if (prior_scale_for_dispersion > 0)
      dispersion[1] <-  prior_scale_for_dispersion * dispersion_unscaled[1];
    else dispersion[1] <- dispersion_unscaled[1];
    mark <- 1;
    scaled_tau <- tau .* scale * square(dispersion[1]);
    for (i in 1:t) {
      real trace_mat;
      trace_mat <- square(scaled_tau[i]);
      if (p[i] == 1) {
        var_group[mark] <- trace_mat;
        mark <- mark + 1;
      }
      else {
        int nc;
        vector[p[i]] temp;
        nc <- p[i];
        trace_mat <- trace_mat * sqrt(nc);
        temp <- segment(zeta, mark, nc);
        temp <- temp / sum(temp);
        for (j in 1:nc) {
          var_group[mark] <- temp[j] * trace_mat;
          mark <- mark + 1;
        }
      }
    }
    b <- make_b(u, z_T, rho, var_group, p, l);
  }
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
  if (t > 0) {
    eta0 <- eta0 + Z0 * b;
    eta1 <- eta1 + Z1 * b;
  }
  if (has_intercept == 1) {
    if (link != 4) {
      eta0 <- gamma[1] + eta0;
      eta1 <- gamma[1] + eta1;
    }
    else {
      real shift;
      shift <- fmax(max(eta0), max(eta1));
      eta0 <- gamma[1] + eta0 - shift;
      eta1 <- gamma[1] + eta1 - shift;
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
    if (prior_scale_for_dispersion > 0) dispersion_unscaled ~ cauchy(0, 1);    
    u ~ normal(0,1);
    z_T ~ normal(0,1);
    rho ~ beta(shape1,shape2);
    zeta ~ gamma(delta, 1);
    tau ~ gamma(gamma_shape, 1);
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
    if (t > 0) {
      eta0 <- eta0 + Z0 * b;
      eta1 <- eta1 + Z1 * b;
    }
    if (has_intercept == 1) {
      if (link != 4) {
        eta0 <- gamma[1] + eta0;
        eta1 <- gamma[1] + eta1;
      }      
      else {
        real shift;
        shift <- fmax(max(eta0), max(eta1));
        eta0 <- gamma[1] + eta0 - shift;
        eta1 <- gamma[1] + eta1 - shift;
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
