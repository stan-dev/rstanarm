# GLM for a Gaussian (or Gamma) outcome
functions {
  /** 
   * Apply inverse link function to linear predictor
   *
   * @param eta Linear predictor vector
   * @param link An integer indicating the link function
   * @return A vector, i.e. inverse-link(eta)
   */
  vector linkinv_gauss(vector eta, int link) {
    if (link < 1 || link > 3) reject("Invalid link");
    if (link < 3)  # link = identity or log 
      return(eta); # return eta for log link too bc will use lognormal
    else {# link = inverse
      vector[rows(eta)] mu;
      for(n in 1:rows(eta)) mu[n] <- inv(eta[n]); 
      return mu;
    }
  }

  /** 
   * Apply inverse link function to linear predictor
   *
   * @param eta Linear predictor vector
   * @param link An integer indicating the link function
   * @return A vector, i.e. inverse-link(eta)
   */
  vector linkinv_gamma(vector eta, int link) {
    if (link < 1 || link > 3) reject("Invalid link");
    if (link == 1)  return eta;
    else if (link == 2) return exp(eta);
    else {
      vector[rows(eta)] mu;
      for(n in 1:rows(eta)) mu[n] <- inv(eta[n]); 
      return mu;
    }
  }

  /** 
   * Apply inverse link function to linear predictor
   *
   * @param eta Linear predictor vector
   * @param link An integer indicating the link function
   * @return A vector, i.e. inverse-link(eta)
   */
  vector linkinv_inv_gaussian(vector eta, int link) {
    if (link < 1 || link > 4) reject("Invalid link");
    if (link == 1)  return eta;
    else if (link == 2) return exp(eta);
    else {
      vector[rows(eta)] mu;
      if (link == 3) for( n in 1:rows(eta)) mu[n] <- inv(eta[n]);
      else for (n in 1:rows(eta)) mu[n] <- inv_sqrt(eta[n]);      
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
  vector pw_gauss(vector y, vector eta, real sigma, int link) {
    vector[rows(eta)] ll;
    if (link < 1 || link > 3) reject("Invalid link");
    if (link == 2) # link = log
      for (n in 1:rows(eta)) ll[n] <- lognormal_log(y[n], eta[n], sigma);
    else { # link = idenity or inverse
      vector[rows(eta)] mu;
      mu <- linkinv_gauss(eta, link);
      for (n in 1:rows(eta)) ll[n] <- normal_log(y[n], mu[n], sigma);
    }
    return ll;
  }

  /** 
   * Divide a scalar by a vector
   *
   * @param x The scalar in the numerator
   * @param y The vector in the denominator
   * @return An elementwise vector
   */
  vector divide_real_by_vector(real x, vector y) {
    vector[rows(y)] ret;
    for (n in 1:rows(y)) ret[n] <- x / y[n];
    return ret;
  }

  real GammaReg_log(vector y, vector eta, real shape, 
                    int link, real sum_log_y) {
    real ret;
    if (link < 1 || link > 3) reject("Invalid link");
    ret <- rows(y) * (shape * log(shape) - lgamma(shape)) +
           (shape - 1) * sum_log_y;
    if (link == 2)      # link is log
      ret <- ret - shape * sum(eta) - shape * sum(y ./ exp(eta));
    else if (link == 1) # link is identity
      ret <- ret - shape * sum(log(eta)) - shape * sum(y ./ eta);
    else                # link is inverse
      ret <- ret + shape * sum(log(eta)) - shape * dot_product(eta, y);
    return ret;
  }

  /** 
   * Pointwise (pw) log-likelihood vector
   *
   * @param y The integer array corresponding to the outcome variable.
   * @param link An integer indicating the link function
   * @return A vector
   */
  vector pw_gamma(vector y, vector eta, real shape, int link) {
    vector[rows(eta)] ll;
    if (link < 1 || link > 3) reject("Invalid link");
    if (link == 3) { # link = inverse
      for (n in 1:rows(eta)) {
        ll[n] <- gamma_log(y[n], shape, shape * eta[n]);
      }
    }
    else if (link == 2) { # link = log
      for (n in 1:rows(eta)) {
        ll[n] <- gamma_log(y[n], shape, shape / exp(eta[n]));
      }
    }
    else { # link = identity
      for (n in 1:rows(eta)) {
        ll[n] <- gamma_log(y[n], shape, shape / eta[n]);
      }
    }
    return ll;
  }

  /** 
   * inverse Gaussian log-PDF (for data only, excludes constants)
   *
   * @param y The vector of outcomes
   * @param eta The vector of linear predictors
   * @param lambda A positive scalar nuisance parameter
   * @param link An integer indicating the link function
   * @return A scalar
   */
  real inv_gaussian_log(vector y, vector mu, real lambda, 
                        real sum_log_y, vector sqrt_y) {
    return 0.5 * rows(y) * log(lambda / (2 * pi())) - 
           1.5 * sum_log_y - 
           0.5 * lambda * dot_self( (y - mu) ./ (mu .* sqrt_y) );
  }

  /** 
   * Pointwise (pw) log-likelihood vector
   *
   * @param y The integer array corresponding to the outcome variable.
   * @param eta The linear predictors
   * @param lamba A positive scalar nuisance parameter
   * @param link An integer indicating the link function
   * @param log_y A precalculated vector of the log of y
   * @param sqrt_y A precalculated vector of the square root of y
   * @return A vector of log-likelihoods
   */
  vector pw_inv_gaussian(vector y, vector eta, real lambda, 
                         int link, vector log_y, vector sqrt_y) {
    vector[rows(y)] ll;
    vector[rows(y)] mu;
    if (link < 1 || link > 4) reject("Invalid link");
    mu <- linkinv_inv_gaussian(eta, link);
    for (n in 1:rows(y))
      ll[n] <- -0.5 * lambda * square( (y[n] - mu[n]) / (mu[n] * sqrt_y[n]) );
    ll <- ll + 0.5 * log(lambda / (2 * pi())) - 1.5 * log_y;
    return ll;
  }

  /** 
   * PRNG for the inverse Gaussian distribution
   *
   * Algorithm from wikipedia 
   *
   * @param mu The expectation
   * @param lambda The dispersion
   * @return A draw from the inverse Gaussian distribution
   */
  real inv_gaussian_rng(real mu, real lambda) {
    real z;
    real y;
    real x;
    real mu2;
    mu2 <- square(mu);
    y <- square(normal_rng(0,1));
    z <- uniform_rng(0,1);
    x <- mu + ( mu2 * y - mu * sqrt(4 * mu * lambda * y + mu2 * square(y)) )
         / (2 * lambda);
    if (z <= (mu / (mu + x))) return x;
    else return mu2 / x;
  }

  /** 
   * Create group-specific block-diagonal Cholesky factor, see section 2.3 of
   * https://cran.r-project.org/web/packages/lme4/vignettes/lmer.pdf
   * @param p An integer array with the number variables on the LHS of each |
   * @param dispersion Scalar standard deviation of the errors
   * @param tau Vector of scale parameters for the decomposed covariance matrices
   * @param scale Vector of scale hyperparameters
   * @param zeta Vector of positive parameters that are normalized into simplexes
   * @param rho Vector of radii in the onion method for creating Cholesky factors
   * @param z_T Vector used in the onion method for creating Cholesky factors
   * @return A vector that corresponds to theta in lme4
   */
  vector make_theta_L(int len_theta_L, int[] p, real dispersion,
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
        theta_L[theta_L_mark] <- tau[i] * scale[i] * dispersion;
        theta_L_mark <- theta_L_mark + 1;
      }
      else { // block is lower-triangular               
        matrix[nc,nc] T_i; 
        real trace_T_i;
        vector[nc] pi; // variance = proportion of trace_T_i
        real std_dev;
        real T21;
        
        trace_T_i <- square(tau[i] * scale[i] * dispersion) * nc;
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
   * @param theta A real array with covariance parameters
   * @param p An integer array with the number variables on the LHS of each |
   * @param l An integer array with the number of levels for the factor(s) on 
   *   the RHS of each |
   * @return A vector of group-specific coefficients
   */
  vector make_b(vector z_b, vector theta_L, int[] p, int[] l) {
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
  
  vector test_csr_matrix_times_vector(int m, int n, vector w, 
                                      int[] v, int[] u, vector b) {
    return csr_matrix_times_vector(m, n, w, v, u, b);            
  }
}
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
  
  # prior family: 0 = none, 1 = normal, 2 = student_t, 3 = horseshoe, 4 = horseshoe_plus
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
  vector<lower=0>[t] gamma_shape; 
  vector<lower=0>[t] scale;
  int<lower=0> len_concentration;
  real<lower=0> concentration[len_concentration];
  int<lower=0> len_shape;
  real<lower=0> shape[len_shape];
}
transformed data {
  int<lower=0> horseshoe;
  int<lower=0> len_theta_L;
  int<lower=0> len_z_T;
  int<lower=0> len_var_group;
  int<lower=0> len_rho;
  real<lower=0> delta[len_concentration];
  int<lower=1> pos;
  vector[N * (family == 3)] sqrt_y;
  vector[N * (family == 3)] log_y;
  real sum_log_y;
  if      (prior_dist <  2) horseshoe <- 0;
  else if (prior_dist == 3) horseshoe <- 2;
  else if (prior_dist == 4) horseshoe <- 4;
  len_theta_L <- 0;
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
    len_theta_L <- len_theta_L + (p[i] * (p[i] - 1)) / 2 + p[i];
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
  real<lower=0> global[horseshoe];
  vector<lower=0>[K] local[horseshoe];
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
    /* else prior_dist is 0 and nothing is added */
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
