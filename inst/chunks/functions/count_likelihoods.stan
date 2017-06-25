  /** 
   * Apply inverse link function to linear predictor
   * see help(poisson) in R
   *
   * @param eta Linear predictor vector
   * @param link An integer indicating the link function
   * @return A vector, i.e. inverse-link(eta)
   */
  vector linkinv_count(vector eta, int link) {
    if (link == 1)      return exp(eta);     // log
    else if (link == 2) return eta;          // identity
    else if (link == 3) return(square(eta)); // sqrt
    else reject("Invalid link");
    return eta; // never reached
  }
  
  /** 
  * Pointwise (pw) log-likelihood vector for the Poisson distribution
  *
  * @param y The integer array corresponding to the outcome variable.
  * @param eta The vector of linear predictors
  * @param link An integer indicating the link function
  * @return A vector
  */
  vector pw_pois(int[] y, vector eta, int link) {
    int N = rows(eta);
    vector[N] ll;
    if (link < 1 || link > 3) 
      reject("Invalid link");
      
    if (link == 1)  // log
      for (n in 1:N) ll[n] = poisson_log_lpmf(y[n] | eta[n]);
    else {  // link = identity or sqrt
      vector[N] phi = linkinv_count(eta, link);
      for (n in 1:N) ll[n] = poisson_lpmf(y[n] | phi[n]) ;
    }
    return ll;
  }
  
  /** 
  * Pointwise (pw) log-likelihood vector for the negative binomial  distribution
  *
  * @param y The integer array corresponding to the outcome variable.
  * @param eta The vector of linear predictors
  * @param theta The reciprocal_dispersion parameter
  * @param link An integer indicating the link function
  * @return A vector
  */
  vector pw_nb(int[] y, vector eta, real theta, int link) {
    int N = rows(eta);
    vector[N] rho = linkinv_count(eta, link);
    vector[N] ll;
    for (n in 1:N) ll[n] = neg_binomial_2_lpmf(y[n] | rho[n], theta);
    return ll;
  }
