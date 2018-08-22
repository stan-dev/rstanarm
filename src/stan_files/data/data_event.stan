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
  int<lower=0> basehaz_nvars;      // num. aux parameters for baseline hazard
  int<lower=0> qnodes;             // num. nodes for GK quadrature
  int<lower=0> len_epts;           // num. epts (event times)
  int<lower=0> len_qpts;           // num. qpts (quadrature points)
  int<lower=0> len_ipts;           // num. ipts (qpts for interval cens.)
  int<lower=0> len_cpts;           // = len_epts + len_qpts + len_ipts
  int idx_cpts[3,2];               // index for breaking cpts into epts,qpts,ipts

  // response and time variables
  vector[len_epts] epts;           // time of events
  vector[len_qpts] qpts;           // time at quadpoints
  vector[len_ipts] ipts;           // time at quadpoints for interval censoring

  // predictor matrices
  matrix[len_cpts, e_K] e_x;             // predictor matrix
  vector[e_K] e_xbar;                    // predictor means

  // spline basis for baseline hazard
  matrix[len_epts, basehaz_nvars] basis_epts;
  matrix[len_qpts, basehaz_nvars] basis_qpts;
  matrix[len_ipts, basehaz_nvars] basis_ipts;

  // GK quadrature weights, with (b-a)/2 scaling already incorporated
  vector[len_qpts] qwts;
  vector[len_ipts] iwts;

  // constant shift for log baseline hazard
  real norm_const;

  // flags
  int<lower=0,upper=1> e_has_intercept; // basehaz requires intercept
