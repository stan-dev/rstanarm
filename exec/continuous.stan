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
   * Pointwise (pw) log-likelihood vector
   *
   * @param y The integer array corresponding to the outcome variable.
   * @param link An integer indicating the link function
   * @return A vector
   */
  vector pw_gamma(vector y, vector eta, real shape, int link) {
    vector[rows(eta)] ll;
    real rate;
    if (link < 1 || link > 3) reject("Invalid link");
    if (link == 3) { # link = inverse
      for (n in 1:rows(eta)) {
        rate <- shape * eta[n];
        ll[n] <- gamma_log(y[n], shape, rate);
      }
    }
    else if (link == 2) { # link = log
      for (n in 1:rows(eta)) {
        rate <- shape / exp(eta[n]);
        ll[n] <- gamma_log(y[n], shape, rate);
      }
    }
    else { # link = identity
      for (n in 1:rows(eta)) {
        rate <- shape / eta[n];
        ll[n] <- gamma_log(y[n], shape, rate);
      }
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
   * @param link An integer indicating the link function
   * @return A vector of log-likelihoods
   */
  vector pw_inv_gaussian(vector y, vector eta, real lambda, 
                         int link, real sum_log_y, vector sqrt_y) {
    vector[rows(y)] ll;
    vector[rows(y)] mu;
    if (link < 1 || link > 4) reject("Invalid link");
    mu <- linkinv_inv_gaussian(eta, link);
    for (n in 1:rows(y))
      ll[n] <- -0.5 * lambda * square( (y[n] - mu[n]) / (mu[n] * sqrt_y[n]) );
    ll <- ll + 0.5 * log(lambda / (2 * pi())) - 1.5 * sum_log_y;
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
   * Create group-specific coefficients, see section 2.3 of
   * https://cran.r-project.org/web/packages/lme4/vignettes/lmer.pdf
   * This function is only used for checking against lme4 behavior
   *
   * @param u Vector whose elements are iid normal(0,sigma) a priori
   * @param theta A real array with covariance parameters
   * @param p An integer array with the number variables on the LHS of each |
   * @param l An integer array with the number of levels for the factor(s) on 
   *   the RHS of each |
   * @return A vector of group-specific coefficients
   */
  vector DoTheDougie(vector u, real[] theta, int[] p, int[] l) {
    vector[rows(u)] b;
    int mark;
    int start;
    mark <- 1;
    start <- 1;
    for (i in 1:size(p)) {
      int nc;
      nc <- p[i];
      if (nc == 1) {
        real theta_start;
        theta_start <- theta[start]; // needs to be positive
        for (s in mark:(mark + l[i] - 1)) 
          b[s] <- theta_start * u[s];
        mark <- mark + l[i];
        start <- start + 1;
      }
      else {
        matrix[nc,nc] T_i;
        int pos;
        T_i <- rep_matrix(0, nc, nc);
        pos <- mark;
        for (c in 1:nc) {
          T_i[c,c] <- theta[pos];    // needs to be positive
          pos <- pos + 1;
          for(r in (c+1):nc) {
            T_i[r,c] <- theta[pos];
            pos <- pos + 1;
          }
        }
        for (j in 1:l[i]) {
          vector[nc] temp;
          temp <- T_i * segment(u, start, nc);
          for (s in 1:nc) b[start + s - 1] <- temp[s];
          start <- start + nc;
        }
        mark <- mark + l[i] * p[i];
      }
    }
    return b;
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
  int<lower=0> t;    # num. terms (maybe 0) with a | in the glmer formula
  int<lower=1> p[t]; # num. variables on the LHS of each |
  int<lower=1> l[t]; # num. levels for the factor(s) on the RHS of each |
  int<lower=0> q;    # conceptually equals \sum_{i=1}^t p_i \times l_i
  matrix[N,q]  Z;    # uncentered design matrix for group-specific variables

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
  int<lower=0> len_z_T;
  int<lower=0> len_var_group;
  int<lower=0> len_rho;
  real<lower=0> shape1[len_shape];
  real<lower=0> shape2[len_shape];
  real<lower=0> delta[len_concentration];
  int<lower=1> pos[2];
  vector[N * (family == 3)] sqrt_y;
  real sum_log_y;
  if      (prior_dist <  2) horseshoe <- 0;
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
  if (family == 3) {
    for (n in 1:N) sqrt_y[n] <- sqrt(y[n]);
    sum_log_y <- sum(log(y));
  }
}
parameters {
  real<lower=if_else(family == 1 || link == 2, 
                     negative_infinity(), 0)> gamma[has_intercept];
  vector[K] z_beta;
  real<lower=0> global[horseshoe];
  vector<lower=0>[K] local[horseshoe];
  real<lower=0> dispersion_unscaled; # interpretation depends on family!
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
  real dispersion;
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
    vector[t] scaled_tau;
    int mark;
    mark <- 1;
    scaled_tau <- tau .* scale * square(dispersion);
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
  vector[N] eta; # linear predictor
  if (K > 0) eta <- X * beta;
  else eta <- rep_vector(0.0, N);
  if (has_intercept == 1) {
    if (family == 1 || link == 2) eta <- eta + gamma[1];
    else eta <- eta - min(eta) + gamma[1];
  }
  if (has_offset == 1)    eta <- eta + offset;
  if (t > 0)              eta <- eta + Z * b;
  
  
  // Log-likelihood 
  if (has_weights == 0 && prior_PD == 0) { # unweighted log-likelihoods
    vector[N] mu;
    if (family == 1) {
      mu <- linkinv_gauss(eta, link);
      if (link == 2) y ~ lognormal(mu, dispersion);
      else y ~ normal(mu, dispersion);
    }
    else if (family == 2) {
      mu <- linkinv_gamma(eta, link);
      y ~ gamma(dispersion, divide_real_by_vector(dispersion, mu));
    }
    else {
      mu <- linkinv_inv_gaussian(eta, link);
      y ~ inv_gaussian(mu, dispersion, sum_log_y, sqrt_y);
    }
  }
  else if (prior_PD == 0) { # weighted log-likelihoods
    vector[N] summands;
    if (family == 1) summands <- pw_gauss(y, eta, dispersion, link);
    else if (family == 2) summands <- pw_gamma(y, eta, dispersion, link);
    else summands <- pw_inv_gaussian(y, eta, dispersion, link, sum_log_y, sqrt_y);
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
    if (family == 2 && link != 2) {
      # nothing because of the weird constraint
    }
    else if (prior_dist_for_intercept == 1) # normal
      gamma ~ normal(prior_mean_for_intercept, prior_scale_for_intercept);
    else if (prior_dist_for_intercept == 2) # student_t
      gamma ~ student_t(prior_df_for_intercept, prior_mean_for_intercept, 
                        prior_scale_for_intercept);
    /* else prior_dist is 0 and nothing is added */
  }
  
  if (t > 0) {
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
  mean_PPD <- 0;
  if (has_intercept == 1)
    alpha[1] <- gamma[1] - dot_product(xbar, beta);
  {
    vector[N] eta;
    if (K > 0) eta <- X * beta;
    else eta <- rep_vector(0.0, N);
    if (has_intercept == 1) {
      if (family == 1 || link == 2) eta <- eta + gamma[1];
      else {
        real min_eta;
        min_eta <- min(eta);
        alpha[1] <- alpha[1] - min_eta;
        eta <- eta - min_eta + gamma[1];
      }
    }
    if (has_offset)         eta <- eta + offset;
    if (t > 0)              eta <- eta + Z * b;
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
