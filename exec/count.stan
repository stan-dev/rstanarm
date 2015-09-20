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
    for (n in 1:rows(eta)) ll[n] <- neg_binomial_2_log(y[n], rho[n], theta);
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
  vector make_b_count(vector u, vector z_T, vector rho, vector var_group, 
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
transformed data{
  real poisson_max;
  int<lower=0> horseshoe;
  int<lower=0> len_z_T;
  int<lower=0> len_var_group;
  int<lower=0> len_rho;
  real<lower=0> delta[len_concentration];
  int<lower=1> pos;
  
  poisson_max <- pow(2.0, 30.0);
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
  real<lower=0> theta_unscaled[family > 1];
  vector<lower=0>[N] noise[family == 3]; // do not store this
  vector[K] z_beta;
  real<upper=if_else(link == 4, 0, positive_infinity())> gamma[has_intercept];
  real<lower=0> global[horseshoe];
  vector<lower=0>[K] local[horseshoe];
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
  vector[len_var_group] var_group;
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
    vector[t] scaled_tau;
    int mark;
    if (family > 1 && prior_scale_for_dispersion > 0)
      dispersion[1] <- prior_scale_for_dispersion * dispersion_unscaled[1];
    else if (family > 1) dispersion[1] <- dispersion_unscaled[1];
    mark <- 1;
    if (family > 1) scaled_tau <- tau .* scale * square(dispersion[1]);
    else scaled_tau <- tau .* scale;
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
    b <- make_b_count(z_b, z_T, rho, var_group, p, l);
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
  
  // Log-prior for dispersion
  if (family > 1 && prior_scale_for_dispersion > 0) theta_unscaled ~ cauchy(0, 1);

  // Log-prior for noise
  if (family == 3) noise[1] ~ gamma(theta[1], 1);
  
  if (t > 0) {
    int pos_shape;
    int pos_rho;
    if (family > 1 && prior_scale_for_dispersion > 0) 
      dispersion_unscaled ~ cauchy(0, 1);
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
