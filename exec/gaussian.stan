# GLM for a Gaussian outcome with Gaussian or t priors
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
   * @param u Vector whose elements are iid normal(0,sigma) a priori
   * @param z_T Vector used in the onion method for creating Cholesky factors
   * @param rho Vector of radii in the onion method for creating Cholesky factors
   * @param pi Vector of simplexes concatenated together
   * @param tau Vector of scale parameters
   * @param p An integer array with the number variables on the LHS of each |
   * @param l An integer array with the number of levels for the factor(s) on 
   *   the RHS of each |
   * @return A vector of group-specific coefficients
   */
  vector make_b(vector u, vector z_T, vector rho, vector pi, vector tau,
                int[] p, int[] l) {
    vector[rows(u)] b;
    int b_mark;
    int z_T_mark;
    int rho_mark;
    int pi_mark;
    // Due to lack of ragged arrays, everything is input as long vectors
    b_mark   <- 1;
    z_T_mark <- 1;
    rho_mark <- 1;
    pi_mark  <- 1;
    for (i in 1:size(p)) {
      int nc;
      nc <- p[i];
      if (nc == 1) { // just a standard deviation times a part of a vector
        for (s in b_mark:(b_mark + l[i] - 1)) 
          b[s] <- tau[i] * u[s];
        b_mark <- b_mark + l[i];
      }
      else {         // deal with a (Cholesky factor of a ) covariance matrix
        matrix[nc,nc] T_i; 
        real std_dev;
        real sqrt_nc_tau;
        sqrt_nc_tau <- sqrt(nc) * tau[i];
        std_dev <- sqrt(pi[pi_mark]) * sqrt_nc_tau;
        T_i <- rep_matrix(0, nc, nc);
        T_i[1,1] <- std_dev;
        pi_mark <- pi_mark + 1;
        std_dev <- sqrt(pi[pi_mark]) * sqrt_nc_tau;
        T_i[2,1] <- std_dev * (2.0 * rho[rho_mark] - 1.0);
        rho_mark <- rho_mark + 1;
        pi_mark <- pi_mark + 1;
        T_i[2,2] <- std_dev * sqrt(1.0 - square(T_i[2,1]));
        for (r in 2:(nc - 1)) { // modified onion method
          int rp1;
          vector[r] T_row;
          real scale;
          T_row <- segment(z_T, z_T_mark, r);
          z_T_mark <- z_T_mark + r;
          scale <- sqrt(rho[rho_mark] / dot_self(T_row));
          std_dev <- sqrt(pi[pi_mark]) * sqrt_nc_tau;
          rp1 <- r + 1;
          for(c in 1:r) T_i[rp1,c] <- T_row[c] * scale * std_dev;
          T_i[rp1,rp1] <- sqrt(1.0 - rho[rho_mark]);
          rho_mark <- rho_mark + 1;
          pi_mark <- pi_mark + 1;
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
  int<lower=1> K; # number of predictors
  
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

  # link function from location to linear predictor
  int<lower=1,upper=3> link; # 1 = identity, 2 = log, 3 = inverse
  
  # weights
  int<lower=0,upper=1> has_weights; # 0 = No, 1 = Yes
  vector[N * has_weights] weights;
  
  # offset
  int<lower=0,upper=1> has_offset;  # 0 = No, 1 = Yes
  vector[N * has_offset] offset;
  
  # prior family (zero indicates no prior!!!)
  int<lower=0,upper=2> prior_dist;               # 1 = normal, 2 = student_t
  int<lower=0,upper=2> prior_dist_for_intercept; # 1 = normal, 2 = student_t
  
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
  int<lower=0> len_z_T;
  int<lower=0> len_rho;
  real<lower=0> shape1[len_shape];
  real<lower=0> shape2[len_shape];
  real<lower=0> delta[len_concentration];
  int<lower=1> pos[2];
  len_z_T <- 0;
  len_rho <- sum(p) - t;
  pos[1] <- 1;
  pos[2] <- 1;
  for (i in 1:t) {
    real nu;
    if (p[i] > 1) for (j in 1:p[i]) {
      delta[pos[2]] <- concentration[j];
      pos[2] <- pos[2] + 1;
      nu <- shape[i] + 0.5  + 0.5 * (p[i] - 2);
      shape1[pos[1]] <- nu;
      shape2[pos[1]] <- nu;
    }
    if (p[i] > 2) for (j in 2:p[i]) {
      pos[1] <- pos[1] + 1;
      nu <- nu - 0.5;
      shape1[pos[1]] <- 0.5 * j;
      shape2[pos[1]] <- nu;
    }
    if (p[i] > 2) for (j in 3:p[i]) {
      len_z_T <- len_z_T + p[i] - 1;
    }
  }
}
parameters {
  real gamma[has_intercept];
  vector[K] z_beta;
  real<lower=0> sigma_unscaled;
  vector[q] u;
  vector[len_z_T] z_T;
  vector<lower=0,upper=1>[len_rho] rho;
  vector<lower=0>[len_concentration] zeta;
  vector<lower=0>[t] tau;
}
transformed parameters {
  vector[K] beta;
  real sigma;
  vector[q] b;
  if (prior_dist > 0) beta <- prior_mean + prior_scale .* z_beta;
  else beta <- z_beta;
  if (prior_scale_for_dispersion > 0)
    sigma <-  prior_scale_for_dispersion * sigma_unscaled;
  else sigma <- sigma_unscaled;
  if (t > 0) {
    vector[len_concentration] pi;
    int mark;
    mark <- 1;
    for (i in 1:t) if (p[i] > 1) {
      int nc;
      vector[p[i]] temp;
      nc <- p[i];
      temp <- segment(zeta, mark, nc);
      temp <- temp / sum(temp);
      for (j in 1:nc) {
        pi[mark] <- temp[j];
        mark <- mark + 1;
      }
    }
    b <- make_b(u * sigma, z_T, rho, pi, tau .* scale, p, l);
  }
}
model {
  vector[N] eta; # linear predictor
  eta <- X * beta;
  if (has_intercept == 1) eta <- eta + gamma[1];
  if (has_offset == 1)    eta <- eta + offset;
  if (t > 0)              eta <- eta + Z * b;
  
  
  // Log-likelihood 
  if (has_weights == 0 && prior_PD == 0) { # unweighted log-likelihoods
    vector[N] mu;
    mu <- linkinv_gauss(eta, link);
    if (link == 2)
      y ~ lognormal(mu, sigma);
    else 
      y ~ normal(mu, sigma);
  }
  else if (prior_PD == 0) { # weighted log-likelihoods
    vector[N] summands;
    summands <- pw_gauss(y, eta, sigma, link);
    increment_log_prob(dot_product(weights, summands));
  }
  
  // Log-prior for scale
  if (prior_scale_for_dispersion > 0) sigma_unscaled ~ cauchy(0, 1);
  
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
    eta <- X * beta;
    if (has_intercept == 1) eta <- eta + gamma[1];
    if (has_offset)         eta <- eta + offset;
    if (t > 0)              eta <- eta + Z * b;
    eta <- linkinv_gauss(eta, link);
    for (n in 1:N) mean_PPD <- mean_PPD + normal_rng(eta[n], sigma);
    mean_PPD <- mean_PPD / N;
  }
}
