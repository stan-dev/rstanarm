  // Log-priors, auxiliary params
  if (has_aux[1] == 1)
    aux_lp(yAux1_unscaled[1], y_prior_dist_for_aux[1],
           y_prior_scale_for_aux[1], y_prior_df_for_aux[1]);
  if (M > 1 && has_aux[2] == 1)
    aux_lp(yAux2_unscaled[1], y_prior_dist_for_aux[2],
           y_prior_scale_for_aux[2], y_prior_df_for_aux[2]);
  if (M > 2 && has_aux[3] == 1)
    aux_lp(yAux3_unscaled[1], y_prior_dist_for_aux[3],
           y_prior_scale_for_aux[3], y_prior_df_for_aux[3]);

  // Log priors, intercepts
  if (intercept_type[1] > 0)
    gamma_lp(yGamma1[1], y_prior_dist_for_intercept[1], y_prior_mean_for_intercept[1],
             y_prior_scale_for_intercept[1], y_prior_df_for_intercept[1]);
  if (M > 1 && intercept_type[2] > 0)
    gamma_lp(yGamma2[1], y_prior_dist_for_intercept[2], y_prior_mean_for_intercept[2],
             y_prior_scale_for_intercept[2], y_prior_df_for_intercept[2]);
  if (M > 2 && intercept_type[3] > 0)
    gamma_lp(yGamma3[1], y_prior_dist_for_intercept[3], y_prior_mean_for_intercept[3],
             y_prior_scale_for_intercept[3], y_prior_df_for_intercept[3]);

  // Log priors, population level params
  if (yK[1] > 0)
    beta_lp(z_yBeta1, y_prior_dist[1], y_prior_scale1, y_prior_df1,
            y_global_prior_df[1], yLocal1, yGlobal1, yMix1, yOol1,
            y_slab_df[1], y_caux1);
  if (M > 1 && yK[2] > 0)
    beta_lp(z_yBeta2, y_prior_dist[2], y_prior_scale2, y_prior_df2,
            y_global_prior_df[2], yLocal2, yGlobal2, yMix2, yOol2,
            y_slab_df[2], y_caux2);
  if (M > 2 && yK[3] > 0)
    beta_lp(z_yBeta3, y_prior_dist[3], y_prior_scale3, y_prior_df3,
            y_global_prior_df[3], yLocal3, yGlobal3, yMix3, yOol3,
            y_slab_df[3], y_caux3);

  // Log priors, group level terms
  if (prior_dist_for_cov == 1) { // decov
    real dummy = decov_lp(z_b, z_T, rho, zeta, tau, b_prior_regularization,
                          delta, b_prior_shape, t, p);
  }
  else if (prior_dist_for_cov == 2) { // lkj
    if (bK1 > 0) {
      // sds for group factor 1
      target += student_t_lpdf(bSd1 | b1_prior_df, 0, b1_prior_scale);
      // primitive coefs for group factor 1
      target += normal_lpdf(to_vector(z_bMat1) | 0, 1);
      // corr matrix for group factor 1
      if (bK1 > 1)
        target += lkj_corr_cholesky_lpdf(bCholesky1 | b1_prior_regularization);
    }
    if (bK2 > 0) {
      // sds for group factor 2
      target += student_t_lpdf(bSd2 | b2_prior_df, 0, b2_prior_scale);
      // primitive coefs for group factor 2
      target += normal_lpdf(to_vector(z_bMat2) | 0, 1);
      // corr matrix for group factor 2
      if (bK2 > 1)
        target += lkj_corr_cholesky_lpdf(bCholesky2 | b2_prior_regularization);
    }
  }
