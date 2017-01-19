  /** 
   * Apply inverse link function to linear predictor
   *
   * @param eta Linear predictor vector
   * @param link An integer indicating the link function
   * @return A vector, i.e. inverse-link(eta)
   */
  vector linkinv_gauss(vector eta, int link) {
    if (link == 1)      return eta;
    else if (link == 2) return exp(eta); 
    else if (link == 3) return inv(eta);
    else reject("Invalid link");
    return eta; // never reached
  }

  /** 
  * Pointwise (pw) log-likelihood vector
  *
  * @param y The integer array corresponding to the outcome variable.
  * @param link An integer indicating the link function
  * @return A vector
  */
  vector pw_gauss(vector y, vector eta, real sigma, int link) {
    return -0.5 * log(6.283185307179586232 * sigma) - 
            0.5 * square((y - linkinv_gauss(eta, link)) / sigma);
  }

  /** 
  * Apply inverse link function to linear predictor
  *
  * @param eta Linear predictor vector
  * @param link An integer indicating the link function
  * @return A vector, i.e. inverse-link(eta)
  */
  vector linkinv_gamma(vector eta, int link) {
    if (link == 1)      return eta;
    else if (link == 2) return exp(eta);
    else if (link == 3) return inv(eta);
    else reject("Invalid link");
    return eta; // never reached
  }

  /** 
  * Pointwise (pw) log-likelihood vector
  *
  * @param y A vector corresponding to the outcome variable.
  * @param eta A vector of linear predictors
  * @param shape A real number for the shape parameter
  * @param link An integer indicating the link function
  * @param sum_log_y A scalar equal to the sum of log(y)
  * @return A scalar log-likelihood
  */
  real GammaReg(vector y, vector eta, real shape, 
                int link, real sum_log_y) {
    real ret;
    if (link < 1 || link > 3) reject("Invalid link");
    ret = rows(y) * (shape * log(shape) - lgamma(shape)) +
      (shape - 1) * sum_log_y;
    if (link == 2)      // link is log
      ret = ret - shape * sum(eta) - shape * sum(y ./ exp(eta));
    else if (link == 1) // link is identity
      ret = ret - shape * sum(log(eta)) - shape * sum(y ./ eta);
    else                // link is inverse
      ret = ret + shape * sum(log(eta)) - shape * dot_product(eta, y);
    return ret;
  }
  
  /** 
  * Pointwise (pw) log-likelihood vector
  *
  * @param y A vector corresponding to the outcome variable.
  * @param shape A real number for the shape parameter
  * @param link An integer indicating the link function
  * @return A vector
  */
  vector pw_gamma(vector y, vector eta, real shape, int link) {
    int N = rows(eta);
    vector[N] ll;
    if (link == 3) { // link = inverse
      for (n in 1:N) {
        ll[n] = gamma_lpdf(y[n] | shape, shape * eta[n]);
      }
    }
    else if (link == 2) { // link = log
      for (n in 1:N) {
        ll[n] = gamma_lpdf(y[n] | shape, shape / exp(eta[n]));
      }
    }
    else if (link == 1) { // link = identity
      for (n in 1:N) {
        ll[n] = gamma_lpdf(y[n] | shape, shape / eta[n]);
      }
    }
    else reject("Invalid link");
    return ll;
  }

  /** 
  * Apply inverse link function to linear predictor
  *
  * @param eta Linear predictor vector
  * @param link An integer indicating the link function
  * @return A vector, i.e. inverse-link(eta)
  */
  vector linkinv_inv_gaussian(vector eta, int link) {
    if (link == 1)      return eta;
    else if (link == 2) return exp(eta);
    else if (link == 3) return inv(eta);
    else if (link == 4) return inv_sqrt(eta);
    else reject("Invalid link");
    return eta; // never reached
  }

  /** 
  * inverse Gaussian log-PDF (for data only, excludes constants)
  *
  * @param y The vector of outcomes
  * @param mu The vector of conditional means
  * @param lambda A positive scalar dispersion parameter
  * @param sum_log_y A scalar equal to the sum of log(y)
  * @param sqrt_y A vector equal to sqrt(y)
  * @return A scalar
  */
  real inv_gaussian(vector y, vector mu, real lambda, 
                    real sum_log_y, vector sqrt_y) {
    return 0.5 * rows(y) * log(lambda / 6.283185307179586232) - 
      1.5 * sum_log_y - 
      0.5 * lambda * dot_self( (y - mu) ./ (mu .* sqrt_y) );
  }
  
  /** 
  * Pointwise (pw) log-likelihood vector
  *
  * @param y The integer array corresponding to the outcome variable.
  * @param eta The linear predictors
  * @param lamba A positive scalar dispersion parameter
  * @param link An integer indicating the link function
  * @param log_y A precalculated vector of the log of y
  * @param sqrt_y A precalculated vector of the square root of y
  * @return A vector of log-likelihoods
  */
  vector pw_inv_gaussian(vector y, vector eta, real lambda, 
                         int link, vector log_y, vector sqrt_y) {
    vector[rows(y)] mu;
    mu = linkinv_inv_gaussian(eta, link); // link checked
    return -0.5 * lambda * square( (y - mu) ./ (mu .* sqrt_y) ) +
            0.5 * log(lambda / 6.283185307179586232) - 1.5 * log_y;
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
    real mu2 = square(mu);
    // compound declare & define does not work with _rng
    real z;
    real y;
    real x;
    z = uniform_rng(0,1);
    y = square(normal_rng(0,1));
    x = mu + ( mu2 * y - mu * sqrt(4 * mu * lambda * y + mu2 * square(y)) )
      / (2 * lambda);
    if (z <= (mu / (mu + x))) return x;
    else return mu2 / x;
  }
  
  /** 
  * Apply inverse link function to linear predictor for beta models
  *
  * @param eta Linear predictor vector
  * @param link An integer indicating the link function
  * @return A vector, i.e. inverse-link(eta)
  */
  vector linkinv_beta(vector eta, int link) {
    vector[rows(eta)] mu;
    if (link < 1 || link > 6) reject("Invalid link");
    if (link == 1)  // logit
      for(n in 1:rows(eta)) mu[n] = inv_logit(eta[n]);
    else if (link == 2)  // probit
      for(n in 1:rows(eta)) mu[n] = Phi(eta[n]);
    else if (link == 3)  // cloglog
      for(n in 1:rows(eta)) mu[n] = inv_cloglog(eta[n]);
    else if (link == 4) // cauchy
      for(n in 1:rows(eta)) mu[n] = cauchy_cdf(eta[n], 0.0, 1.0);
    else if (link == 5)  // log 
      for(n in 1:rows(eta)) {
          mu[n] = exp(eta[n]);
          if (mu[n] < 0 || mu[n] > 1)
            reject("mu needs to be between 0 and 1");
      }
    else if (link == 6) // loglog
      for(n in 1:rows(eta)) mu[n] = 1-inv_cloglog(-eta[n]); 
      
    return mu;
  }
  
  /** 
  * Apply inverse link function to linear predictor for dispersion for beta models
  *
  * @param eta Linear predictor vector
  * @param link An integer indicating the link function
  * @return A vector, i.e. inverse-link(eta)
  */
  vector linkinv_beta_z(vector eta, int link) {
    vector[rows(eta)] mu;
    if (link < 1 || link > 3) reject("Invalid link");
    if (link == 1)        // log
      for(n in 1:rows(eta)) mu[n] = exp(eta[n]);
    else if (link == 2)   // identity
      return eta;
    else if (link == 3)   // sqrt
      for(n in 1:rows(eta)) mu[n] = square(eta[n]);
    return mu;
  }
  
  /** 
  * Pointwise (pw) log-likelihood vector for beta models
  *
  * @param y The vector of outcomes
  * @param eta The linear predictors
  * @param dispersion Positive dispersion parameter
  * @param link An integer indicating the link function
  * @return A vector of log-likelihoods
  */
  vector pw_beta(vector y, vector eta, real dispersion, int link) {
    vector[rows(y)] ll;
    vector[rows(y)] mu;
    vector[rows(y)] shape1;
    vector[rows(y)] shape2;
    if (link < 1 || link > 6) reject("Invalid link");
    mu = linkinv_beta(eta, link);
    shape1 = mu * dispersion;
    shape2 = (1 - mu) * dispersion;
    for (n in 1:rows(y)) {
      ll[n] = beta_lpdf(y[n] | shape1[n], shape2[n]);
    }
    return ll;
  }

  /** 
  * Pointwise (pw) log-likelihood vector for beta models with z variables
  *
  * @param y The vector of outcomes
  * @param eta The linear predictors (for y)
  * @param eta_z The linear predictors (for dispersion)
  * @param link An integer indicating the link function passed to linkinv_beta
  * @param link_phi An integer indicating the link function passed to linkinv_beta_z
  * @return A vector of log-likelihoods
  */
  vector pw_beta_z(vector y, vector eta, vector eta_z, int link, int link_phi) {
    vector[rows(y)] ll;
    vector[rows(y)] mu;
    vector[rows(y)] mu_z;
    if (link < 1 || link > 6) reject("Invalid link");
    if (link_phi < 1 || link_phi > 3) reject("Invalid link");
    mu = linkinv_beta(eta, link);
    mu_z = linkinv_beta_z(eta_z, link_phi);
    for (n in 1:rows(y)) {
      ll[n] = beta_lpdf(y[n] | mu[n] * mu_z[n], (1-mu[n]) * mu_z[n]);
    }
    return ll;
  }
