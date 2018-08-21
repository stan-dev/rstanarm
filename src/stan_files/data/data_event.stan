  // prior family:
  //   0 = none
  //   1 = normal
  //   2 = student_t
  //   3 = hs
  //   4 = hs_plus
  //   5 = laplace
  //   6 = lasso
  int<lower=0,upper=6> e_prior_dist;
  int<lower=0,upper=2> e_prior_dist_for_intercept;

  // prior family:
  //   0 = none
  //   1 = normal
  //   2 = student_t
  //   3 = exponential
  int<lower=0,upper=3> e_prior_dist_for_aux; // prior for basehaz params

  // baseline hazard type:
  //   1 = weibull
  //   2 = B-splines
  //   3 = piecewise
  int<lower=1,upper=3> basehaz_type;

  // dimensions
  int<lower=0> e_K;                // num. predictors in event submodel
  int<lower=0> len_epts;           // num. events (ie. not censored)
  int<lower=0> len_qpts;           // num. rows used for quadrature
  int<lower=0> len_ipts;           // num. rows used for quadrature for interval cens
  int<lower=0> basehaz_nvars;      // num. aux parameters for baseline hazard
  int<lower=0> qnodes;             // num. nodes for GK quadrature

  // response and time variables
  vector[len_epts] epts;           // time of events
  vector[len_qpts] qpts;           // time at quadpoints
  vector[len_ipts] ipts;           // time at quadpoints for interval censoring

  // predictor matrices
  matrix[len_epts, e_K] e_x_epts;  // for rows with events
  matrix[len_qpts, e_K] e_x_qpts;  // for rows at quadpoints
  matrix[len_ipts, e_K] e_x_ipts;  // for rows at quadpoints for interval censoring

  // predictor means
  vector[e_K] e_xbar;

  // design matrices for baseline hazard
  matrix[len_epts, basehaz_nvars] basis_epts;  // spline basis for rows with events
  matrix[len_qpts, basehaz_nvars] basis_qpts;  // spline basis for rows at quadpoints
  matrix[len_ipts, basehaz_nvars] basis_ipts;  // spline basis for rows at quadpoints
                                               //   for interval censoring

  // GK quadrature weights, with (b-a)/2 scaling already incorporated
  vector[len_qpts] qwts;
  vector[len_ipts] iwts;

  // weights, set to zero if not used
  vector[len_epts] e_weights_epts;
  vector[len_qpts] e_weights_qpts;
  vector[len_ipts] e_weights_ipts;

  // constant shift for log baseline hazard
  real norm_const;

  // flags
  int<lower=0,upper=1> e_has_intercept; // basehaz requires intercept
