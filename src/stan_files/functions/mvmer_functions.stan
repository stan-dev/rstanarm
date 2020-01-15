  /**
  * Return the required number of local hs parameters
  *
  * @param prior_dist An integer indicating the prior distribution
  * @return An integer
  */
  int get_nvars_for_hs(int prior_dist) {
    int hs = 0;
    if (prior_dist == 3) hs = 2;
    else if (prior_dist == 4) hs = 4;
    return hs;
  }

  /**
  * Return the lower/upper bound for the specified intercept type
  *
  * @param intercept_type An integer specifying the type of intercept;
  *   0=no intercept, 1=unbounded, 2=lower bounded, 3=upper bounded
  * @return A real, corresponding to the lower bound
  */
  real lb(int intercept_type) {
    return intercept_type == 2 ? 0 : negative_infinity();
  }
  real ub(int intercept_type) {
    return intercept_type == 3 ? 0 : positive_infinity();
  }

  /**
  * Get the indices corresponding to the lower tri of a square matrix
  *
  * @param dim The number of rows in the square matrix
  * @return A vector of indices
  */
  int[] lower_tri_indices(int dim) {
    int indices[dim + choose(dim, 2)];
    int mark = 1;
    for (r in 1:dim) {
      for (c in r:dim) {
        indices[mark] = (r - 1) * dim + c;
        mark += 1;
      }
    }
    return indices;
  }

  /**
  * Scale the auxiliary parameter based on prior information
  *
  * @param aux_unscaled A real, the unscaled auxiliary parameter
  * @param prior_dist Integer, the type of prior distribution
  * @param prior_mean,prior_scale Real scalars, the mean and scale
  *   of the prior distribution
  * @return A real, corresponding to the scaled auxiliary parameter
  */
  real make_aux(real aux_unscaled, int prior_dist,
                real prior_mean, real prior_scale) {
    real aux;
    if (prior_dist == 0) // none
      aux = aux_unscaled;
    else {
      aux = prior_scale * aux_unscaled;
      if (prior_dist <= 2) // normal or student_t
        aux += prior_mean;
    }
    return aux;
  }

  /**
  * Scale the primitive population level parameters based on prior information
  *
  * @param z_beta A vector of primitive parameters
  * @param prior_dist Integer, the type of prior distribution
  * @param prior_mean,prior_scale Vectors of mean and scale parameters
  *   for the prior distributions
  * @return A vector containing the population level parameters (coefficients)
  */
  vector make_beta(vector z_beta, int prior_dist, vector prior_mean,
                   vector prior_scale, vector prior_df, real global_prior_scale,
                   real[] global, vector[] local, real[] ool, vector[] mix,
                   real[] aux, int family, real slab_scale, real[] caux) {
    vector[rows(z_beta)] beta;
    if (prior_dist == 0) beta = z_beta;
    else if (prior_dist == 1) beta = z_beta .* prior_scale + prior_mean;
    else if (prior_dist == 2) for (k in 1:rows(prior_mean)) {
      beta[k] = CFt(z_beta[k], prior_df[k]) * prior_scale[k] + prior_mean[k];
    }
    else if (prior_dist == 3) {
      real c2 = square(slab_scale) * caux[1];
      if (family == 1) // don't need is_continuous since family == 1 is gaussian in mvmer
        beta = hs_prior(z_beta, global, local, global_prior_scale, aux[1], c2);
      else
        beta = hs_prior(z_beta, global, local, global_prior_scale, 1, c2);
    }
    else if (prior_dist == 4) {
      real c2 = square(slab_scale) * caux[1];
      if (family == 1) // don't need is_continuous since family == 1 is gaussian in mvmer
        beta = hsplus_prior(z_beta, global, local, global_prior_scale, aux[1], c2);
      else
        beta = hsplus_prior(z_beta, global, local, global_prior_scale, 1, c2);
    }
    else if (prior_dist == 5) // laplace
      beta = prior_mean + prior_scale .* sqrt(2 * mix[1]) .* z_beta;
    else if (prior_dist == 6) // lasso
      beta = prior_mean + ool[1] * prior_scale .* sqrt(2 * mix[1]) .* z_beta;
    return beta;
  }

  /**
  * Create group-specific coefficients, see section 2 of
  * https://cran.r-project.org/web/packages/lme4/vignettes/lmer.pdf
  *
  * @param z_b Vector whose elements are iid normal(0,sigma) a priori
  * @param theta Vector with covariance parameters as defined in lme4
  * @param p An integer array with the number variables on the LHS of each |
  * @param l An integer array with the number of levels for the factor(s) on
  *   the RHS of each |
  * @param i The index of the grouping factor for which you want to return
  *   the group-specific coefficients for
  * @return An array of group-specific coefficients for grouping factor i
  */
  matrix make_b_matrix(vector z_b, vector theta_L, int[] p, int[] l, int i) {
    matrix[p[i],l[i]] b_matrix;
    int nc = p[i];
    int b_mark = 1;
    int theta_L_mark = 1;
    if (i > 1) {
      for (j in 1:(i-1)) {
        theta_L_mark += p[j] + choose(p[j], 2);
        b_mark += p[j] * l[j];
      }
    }
    if (nc == 1) {
      real theta_L_start = theta_L[theta_L_mark];
      for (s in b_mark:(b_mark + l[i] - 1))
        b_matrix[nc,s] = theta_L_start * z_b[s];
    }
    else {
      matrix[nc,nc] T_i = rep_matrix(0, nc, nc);
      for (c in 1:nc) {
        T_i[c,c] = theta_L[theta_L_mark];
        theta_L_mark += 1;
        for(r in (c+1):nc) {
          T_i[r,c] = theta_L[theta_L_mark];
          theta_L_mark += 1;
        }
      }
      for (j in 1:l[i]) {
        vector[nc] temp = T_i * segment(z_b, b_mark, nc);
        b_matrix[,j] = temp;
        b_mark += nc;
      }
    }
    return b_matrix';
  }

  /**
  * Evaluate the linear predictor for the glmer submodel
  *
  * @param X Design matrix for fe
  * @param Z1 Design matrix for re, for first grouping factor
  * @param Z2 Design matrix for re, for second grouping factor
  * @param Z1_id Group indexing for Z1
  * @param Z2_id Group indexing for Z2
  * @param gamma The intercept parameter
  * @param beta Vector of population level parameters
  * @param b1Mat Matrix of group level params for first grouping factor
  * @param b2Mat Matrix of group level params for second grouping factor
  * @param b1Mat_colshift,b2Mat_colshift Number of columns in b1Mat/b2Mat
  *   that correpond to group level params from prior glmer submodels
  * @param intercept_type The type of intercept parameter (0 = none,
  *   1 = unbounded, 2 = lower bound, 3 = upper bound)
  * @return A vector containing the linear predictor for the glmer submodel
  */
  vector evaluate_eta(matrix X, vector[] Z1, vector[] Z2, int[] Z1_id, int[] Z2_id,
                      real[] gamma, vector beta, matrix b1Mat, matrix b2Mat,
                      int b1Mat_colshift, int b2Mat_colshift,
                      int intercept_type) {
    int N = rows(X);    // num rows in design matrix
    int K = rows(beta); // num predictors
    int p1 = size(Z1);  // num group level params for group factor 1
    int p2 = size(Z2);  // num group level params for group factor 2
    vector[N] eta;

    if (K > 0) eta = X * beta;
    else eta = rep_vector(0.0, N);

    if (intercept_type > 0) { // submodel has an intercept
      if (intercept_type == 1) eta += gamma[1];
      else if (intercept_type == 2) eta += gamma[1] - max(eta);
      else if (intercept_type == 3) eta += gamma[1] - min(eta);
    }

    if (p1 > 0) { // submodel includes group factor 1
      for (k in 1:p1)
        for (n in 1:N)
          eta[n] += (b1Mat[Z1_id[n], k+b1Mat_colshift]) * Z1[k,n];
    }
    if (p2 > 0) { // submodel includes group factor 2
      for (k in 1:p2)
        for (n in 1:N)
          eta[n] += (b2Mat[Z2_id[n], k+b2Mat_colshift]) * Z2[k,n];
    }

    return eta;
  }

  /**
  * Evaluate mu based on eta, family and link
  *
  * @param eta Vector of linear predictors
  * @param family An integer indicating the family
  * @param link An integer indicating the link function (differs by family)
  * @return A vector
  */
  vector evaluate_mu(vector eta, int family, int link) {
    vector[rows(eta)] mu;
    if (family == 1)
      mu = linkinv_gauss(eta, link);
    else if (family == 2)
      mu = linkinv_gamma(eta, link);
    else if (family == 3)
      mu = linkinv_inv_gaussian(eta, link);
    else if (family == 4)
      mu = linkinv_bern(eta, link);
    else if (family == 5)
      mu = linkinv_binom(eta, link);
    else if (family == 6 || family == 7 || family == 8)
      mu = linkinv_count(eta, link);
    return mu;
  }

  /**
  * Increment the target with the log-likelihood for the glmer submodel
  *
  * @param z_beta A vector of primitive parameters
  * @param prior_dist Integer, the type of prior distribution
  * @param prior_mean,prior_scale Vectors of mean and scale parameters
  *   for the prior distributions
  * @return A vector containing the population level parameters (coefficients)
  */
  void glm_lp(vector y_real, int[] y_integer, vector eta, real[] aux,
              int family, int link, real sum_log_y, vector sqrt_y, vector log_y) {
    if (family == 1) {  // gaussian
      if (link == 1) target += normal_lpdf(y_real | eta, aux[1]);
      else if (link == 2) target += lognormal_lpdf(y_real | eta, aux[1]);
      else target += normal_lpdf(y_real | inv(eta), aux[1]);
    }
    else if (family == 2) {  // gamma
      target += GammaReg(y_real, eta, aux[1], link, sum_log_y);
    }
    else if (family == 3) {  // inverse gaussian
      target += inv_gaussian(y_real, linkinv_inv_gaussian(eta, link),
                             aux[1], sum_log_y, sqrt_y);
    }
    else if (family == 4) {  // bernoulli
      if (link == 1) target += bernoulli_logit_lpmf(y_integer | eta);
      else target += bernoulli_lpmf(y_integer | linkinv_bern(eta, link));
    }
    else if (family == 5) {  // binomial
      reject("Binomial with >1 trials not allowed.");
    }
    else if (family == 6 || family == 8) {  // poisson or poisson-gamma
      if (link == 1) target += poisson_log_lpmf(y_integer | eta);
      else target += poisson_lpmf(y_integer | linkinv_count(eta, link));
    }
    else if (family == 7) {  // negative binomial
        if (link == 1) target += neg_binomial_2_log_lpmf(y_integer | eta, aux[1]);
      else target += neg_binomial_2_lpmf(y_integer | linkinv_count(eta, link), aux[1]);
    }
    else reject("Invalid family.");
  }

  /**
  * Log-prior for coefficients
  *
  * @param z_beta Vector of primative coefficients
  * @param prior_dist Integer, the type of prior distribution
  * @param prior_scale Real, scale for the prior distribution
  * @param prior_df Real, df for the prior distribution
  * @param global_prior_df Real, df for the prior for the global hs parameter
  * @param local Vector of hs local parameters
  * @param global Real, the global parameter
  * @param mix Vector of shrinkage parameters
  * @param one_over_lambda Real
  * @return nothing
  */
  void beta_lp(vector z_beta, int prior_dist, vector prior_scale,
               vector prior_df, real global_prior_df, vector[] local,
               real[] global, vector[] mix, real[] one_over_lambda,
               real slab_df, real[] caux) {
    if      (prior_dist == 1) target += normal_lpdf(z_beta | 0, 1);
    else if (prior_dist == 2) target += normal_lpdf(z_beta | 0, 1); // Student t
    else if (prior_dist == 3) { // hs
      target += normal_lpdf(z_beta | 0, 1);
      target += normal_lpdf(local[1] | 0, 1);
      target += inv_gamma_lpdf(local[2] | 0.5 * prior_df, 0.5 * prior_df);
      target += normal_lpdf(global[1] | 0, 1);
      target += inv_gamma_lpdf(global[2] | 0.5 * global_prior_df, 0.5 * global_prior_df);
      target += inv_gamma_lpdf(caux | 0.5 * slab_df, 0.5 * slab_df);
    }
    else if (prior_dist == 4) { // hs+
      target += normal_lpdf(z_beta | 0, 1);
      target += normal_lpdf(local[1] | 0, 1);
      target += inv_gamma_lpdf(local[2] | 0.5 * prior_df, 0.5 * prior_df);
      target += normal_lpdf(local[3] | 0, 1);
      // unorthodox useage of prior_scale as another df hyperparameter
      target += inv_gamma_lpdf(local[4] | 0.5 * prior_scale, 0.5 * prior_scale);
      target += normal_lpdf(global[1] | 0, 1);
      target += inv_gamma_lpdf(global[2] | 0.5 * global_prior_df, 0.5 * global_prior_df);
      target += inv_gamma_lpdf(caux | 0.5 * slab_df, 0.5 * slab_df);
    }
    else if (prior_dist == 5) { // laplace
      target += normal_lpdf(z_beta | 0, 1);
      target += exponential_lpdf(mix[1] | 1);
    }
    else if (prior_dist == 6) { // lasso
      target += normal_lpdf(z_beta | 0, 1);
      target += exponential_lpdf(mix[1] | 1);
      target += chi_square_lpdf(one_over_lambda[1] | prior_df[1]);
    }
    else if (prior_dist == 7) { // product_normal
      target += normal_lpdf(z_beta | 0, 1);
    }
    /* else prior_dist is 0 and nothing is added */
  }

  /**
  * Log-prior for intercept parameters
  *
  * @param gamma Real, the intercept parameter
  * @param dist Integer, the type of prior distribution
  * @param mean_ Real, mean of prior distribution
  * @param scale Real, scale for the prior distribution
  * @param df Real, df for the prior distribution
  * @return nothing
  */
  void gamma_lp(real gamma, int dist, real mean_, real scale, real df) {
    if (dist == 1)  // normal
      target += normal_lpdf(gamma | mean_, scale);
    else if (dist == 2)  // student_t
      target += student_t_lpdf(gamma | df, mean_, scale);
    /* else dist is 0 and nothing is added */
  }

  /**
  * Log-prior for auxiliary parameters
  *
  * @param aux_unscaled Vector (potentially of length 1) of unscaled
  *   auxiliary parameter(s)
  * @param dist Integer specifying the type of prior distribution
  * @param scale Real specifying the scale for the prior distribution
  * @param df Real specifying the df for the prior distribution
  * @return nothing
  */
  void aux_lp(real aux_unscaled, int dist, real scale, real df) {
    if (dist > 0 && scale > 0) {
      if (dist == 1)
        target += normal_lpdf(aux_unscaled | 0, 1);
      else if (dist == 2)
        target += student_t_lpdf(aux_unscaled | df, 0, 1);
      else
        target += exponential_lpdf(aux_unscaled | 1);
    }
  }

  /**
  * Evaluate the mean of the posterior predictive distribution
  *
  * @param mu Vector containing the mean of the posterior predictive
  *   distribution for each observation (ie. the linear predictor after
  *   applying the inverse link function).
  * @param real The auxiliary parameter for the glmer submodel. This will be
  *   an empty array if the submodel does not have an auxiliary parameter
  * @param family An integer specifying the family
  * @return A real, the mean of the posterior predictive distribution
  */
  real mean_PPD_rng(vector mu, real[] aux, int family) {
    int N = rows(mu);
    real mean_PPD = 0;
    if (family == 1) { // gaussian
      for (n in 1:N)
        mean_PPD += normal_rng(mu[n], aux[1]);
    }
    else if (family == 2) {  // gamma
      for (n in 1:N)
        mean_PPD += gamma_rng(aux[1], aux[1] / mu[n]);
    }
    else if (family == 3) {  // inverse gaussian
      for (n in 1:N)
        mean_PPD += inv_gaussian_rng(mu[n], aux[1]);
    }
    else if (family == 4) {  // bernoulli
      for (n in 1:N)
        mean_PPD += bernoulli_rng(mu[n]);
    }
    else if (family == 5) {  // binomial
      reject("Binomial with >1 trials not allowed.");
    }
    else if (family == 6 || family == 8) {
      real poisson_max = pow(2.0, 30.0);
      for (n in 1:N) {  // poisson or poisson-gamma
        if (mu[n] < poisson_max)
          mean_PPD += poisson_rng(mu[n]);
        else
          mean_PPD += normal_rng(mu[n], sqrt(mu[n]));
      }
    }
    else if (family == 7) {
      real poisson_max = pow(2.0, 30.0);
      for (n in 1:N) {  // negative binomial
        real gamma_temp;
        if (is_inf(aux[1]))
          gamma_temp = mu[n];
        else
          gamma_temp = gamma_rng(aux[1], aux[1] / mu[n]);
        if (gamma_temp < poisson_max)
          mean_PPD += poisson_rng(gamma_temp);
        else
          mean_PPD += normal_rng(gamma_temp, sqrt(gamma_temp));
      }
    }
    mean_PPD /= N;
    return mean_PPD;
  }
