  /**
  * Log hazard for exponential distribution
  *
  * @param eta Vector, linear predictor
  * @return A vector
  */
  vector exponential_log_haz(vector eta) {
    return eta;
  }

  /**
  * Log hazard for Weibull distribution
  *
  * @param eta Vector, linear predictor
  * @param t Vector, event or censoring times
  * @param shape Real, Weibull shape
  * @return A vector
  */
  vector weibull_log_haz(vector eta, vector t, real shape) {
    vector[rows(eta)] res;
    res = log(shape) + (shape - 1) * log(t) + eta;
    return res;
  }

  /**
  * Log hazard for Gompertz distribution
  *
  * @param eta Vector, linear predictor
  * @param t Vector, event or censoring times
  * @param scale Real, Gompertz scale
  * @return A vector
  */
  vector gompertz_log_haz(vector eta, vector t, real scale) {
    vector[rows(eta)] res;
    res = scale * t + eta;
    return res;
  }

  /**
  * Log hazard for M-spline model
  *
  * @param eta Vector, linear predictor
  * @param t Vector, event or censoring times
  * @param coefs Vector, M-spline coefficients
  * @return A vector
  */
  vector mspline_log_haz(vector eta, matrix basis, vector coefs) {
    vector[rows(eta)] res;
    res = log(basis * coefs) + eta;
    return res;
  }

  /**
  * Log hazard for B-spline model
  *
  * @param eta Vector, linear predictor
  * @param t Vector, event or censoring times
  * @param coefs Vector, B-spline coefficients
  * @return A vector
  */
  vector bspline_log_haz(vector eta, matrix basis, vector coefs) {
    vector[rows(eta)] res;
    res = basis * coefs + eta;
    return res;
  }

  /**
  * Evaluate log survival or log CDF from the log hazard evaluated at
  * quadrature points and a corresponding vector of quadrature weights
  *
  * @param qwts Vector, the quadrature weights
  * @param log_hazard Vector, log hazard at the quadrature points
  * @param qnodes Integer, the number of quadrature points for each individual
  * @return A vector
  */
  real quadrature_log_surv(vector qwts, vector log_hazard) {
    real res;
    res = - dot_product(qwts, exp(log_hazard)); // sum across all individuals
    return res;
  }

  vector quadrature_log_cdf(vector qwts, vector log_hazard, int qnodes) {
    int M = rows(log_hazard);
    int N = M / qnodes; // num of individuals
    vector[M] hazard = exp(log_hazard);
    matrix[N,qnodes] qwts_mat = to_matrix(qwts,   N, qnodes);
    matrix[N,qnodes] haz_mat  = to_matrix(hazard, N, qnodes);
    vector[N] chaz = rows_dot_product(qwts_mat, haz_mat);
    vector[N] res;
    res = log(1 - exp(- chaz));
    return res;
  }

  vector quadrature_log_cdf2(vector qwts_lower, vector log_hazard_lower,
                             vector qwts_upper, vector log_hazard_upper,
                             int qnodes) {
    int M = rows(log_hazard_lower);
    int N = M / qnodes; // num of individuals
    vector[M] hazard_lower = exp(log_hazard_lower);
    vector[M] hazard_upper = exp(log_hazard_upper);
    matrix[N,qnodes] qwts_lower_mat = to_matrix(qwts_lower,   N, qnodes);
    matrix[N,qnodes] qwts_upper_mat = to_matrix(qwts_lower,   N, qnodes);
    matrix[N,qnodes] haz_lower_mat  = to_matrix(hazard_lower, N, qnodes);
    matrix[N,qnodes] haz_upper_mat  = to_matrix(hazard_upper, N, qnodes);
    vector[N] chaz_lower = rows_dot_product(qwts_lower_mat, haz_lower_mat);
    vector[N] chaz_upper = rows_dot_product(qwts_upper_mat, haz_upper_mat);
    vector[N] surv_lower = exp(- chaz_lower);
    vector[N] surv_upper = exp(- chaz_upper);
    vector[N] res;
    res = log(surv_lower - surv_upper);
    return res;
  }
