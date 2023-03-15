  // intercepts
  real<lower=lb(intercept_type[1]),upper=ub(intercept_type[1])>
    yGamma1[intercept_type[1] > 0];
  real<lower=lb(intercept_type[2]),upper=ub(intercept_type[2])>
    yGamma2[intercept_type[2] > 0];
  real<lower=lb(intercept_type[3]),upper=ub(intercept_type[3])>
    yGamma3[intercept_type[3] > 0];

  // population level primitive params
  vector[yK[1]] z_yBeta1;
  vector[yK[2]] z_yBeta2;
  vector[yK[3]] z_yBeta3;

  // group level params, decov prior
  vector[prior_dist_for_cov == 1 ? q : 0] z_b;
  vector[prior_dist_for_cov == 1 ? len_z_T : 0] z_T;
  vector<lower=0,upper=1>[prior_dist_for_cov == 1 ? len_rho : 0] rho;
  vector<lower=0>[prior_dist_for_cov == 1 ? len_concentration : 0] zeta;
  vector<lower=0>[prior_dist_for_cov == 1 ? t : 0] tau;

  // group level params for first grouping factor
    // group-level sds
    vector<lower=0>[prior_dist_for_cov == 2 ? bK1 : 0] bSd1;
    // unscaled group-level params
    matrix[prior_dist_for_cov == 2 && bK1 >  0 ? bK1 : 0, bK1 >  0 ? bN1 : 0] z_bMat1;
    // cholesky factor of corr matrix (if > 1 random effect)
    cholesky_factor_corr[prior_dist_for_cov == 2 && bK1 > 1 ? bK1 : 0] bCholesky1;

  // group level params for second grouping factor
    // group-level sds
    vector<lower=0>[prior_dist_for_cov == 2 ? bK2 : 0] bSd2;
    // unscaled group-level params
    matrix[prior_dist_for_cov == 2 && bK2 >  0 ? bK2 : 0, bK2 >  0 ? bN2 : 0] z_bMat2;
    // cholesky factor of corr matrix (if > 1 random effect)
    cholesky_factor_corr[prior_dist_for_cov == 2 && bK2 > 1 ? bK2 : 0] bCholesky2;

  // auxiliary params, interpretation depends on family
  real<lower=0> yAux1_unscaled[has_aux[1]];
  real<lower=0> yAux2_unscaled[has_aux[2]];
  real<lower=0> yAux3_unscaled[has_aux[3]];

  // params for priors
  real<lower=0> yGlobal1[yHs1];
  real<lower=0> yGlobal2[yHs2];
  real<lower=0> yGlobal3[yHs3];
  vector<lower=0>[yK[1]] yLocal1[yHs1];
  vector<lower=0>[yK[2]] yLocal2[yHs2];
  vector<lower=0>[yK[3]] yLocal3[yHs3];
  real<lower=0> y_caux1[yHs1 > 0];
  real<lower=0> y_caux2[yHs2 > 0];
  real<lower=0> y_caux3[yHs3 > 0];
  real<lower=0> yOol1[y_prior_dist[1] == 6]; // one_over_lambda
  real<lower=0> yOol2[y_prior_dist[2] == 6];
  real<lower=0> yOol3[y_prior_dist[3] == 6];
  vector<lower=0>[yK[1]] yMix1[y_prior_dist[1] == 5 || y_prior_dist[1] == 6];
  vector<lower=0>[yK[2]] yMix2[y_prior_dist[2] == 5 || y_prior_dist[2] == 6];
  vector<lower=0>[yK[3]] yMix3[y_prior_dist[3] == 5 || y_prior_dist[3] == 6];
