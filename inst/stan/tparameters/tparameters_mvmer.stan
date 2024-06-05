  vector[yK[1]] yBeta1; // population level params
  vector[yK[2]] yBeta2;
  vector[yK[3]] yBeta3;
  array[has_aux[1]] real yAux1; // auxiliary params
  array[has_aux[2]] real yAux2;
  array[has_aux[3]] real yAux3;
  vector[len_theta_L] theta_L; // cov matrix for decov prior
  real yAuxMaximum = 1.0; // used for scaling in theta_L

  // group level params
  matrix[bK1 >  0 ? bN1 : 0, bK1] bMat1; // for grouping factor 1
  matrix[bK2 >  0 ? bN2 : 0, bK2] bMat2; // for grouping factor 2

  // population level params, auxiliary params
  if (has_aux[1] == 1) {
    yAux1[1] = make_aux(yAux1_unscaled[1], y_prior_dist_for_aux[1],
                        y_prior_mean_for_aux[1], y_prior_scale_for_aux[1]);
    if (yAux1[1] > yAuxMaximum)
      yAuxMaximum = yAux1[1];
  }

  if (yK[1] > 0)
    yBeta1 = make_beta(z_yBeta1, y_prior_dist[1], y_prior_mean1,
                       y_prior_scale1, y_prior_df1, y_global_prior_scale[1],
                       yGlobal1, yLocal1, yOol1, yMix1, yAux1, family[1],
                       y_slab_scale[1], y_caux1);
  if (M > 1) {
    if (has_aux[2] == 1) {
      yAux2[1] = make_aux(yAux2_unscaled[1], y_prior_dist_for_aux[2],
                          y_prior_mean_for_aux[2], y_prior_scale_for_aux[2]);
      if (yAux2[1] > yAuxMaximum)
        yAuxMaximum = yAux2[1];
    }
    if (yK[2] > 0)
      yBeta2 = make_beta(z_yBeta2, y_prior_dist[2], y_prior_mean2,
                         y_prior_scale2, y_prior_df2, y_global_prior_scale[2],
                         yGlobal2, yLocal2, yOol2, yMix2, yAux2, family[2],
                         y_slab_scale[2], y_caux2);
  }
  if (M > 2) {
    if (has_aux[3] == 1) {
      yAux3[1] = make_aux(yAux3_unscaled[1], y_prior_dist_for_aux[3],
                          y_prior_mean_for_aux[3], y_prior_scale_for_aux[3]);
      if (yAux3[1] > yAuxMaximum)
        yAuxMaximum = yAux3[1];
    }
    if (yK[3] > 0)
      yBeta3 = make_beta(z_yBeta3, y_prior_dist[3], y_prior_mean3,
                         y_prior_scale3, y_prior_df3, y_global_prior_scale[3],
                         yGlobal3, yLocal3, yOol3, yMix3, yAux3, family[3],
                         y_slab_scale[3], y_caux3);
  }

  // group level params, under decov prior
  if (prior_dist_for_cov == 1) {
    int mark = 1;
    // cov matrix
    theta_L = make_theta_L(len_theta_L, p, yAuxMaximum, tau,
                           b_prior_scale, zeta, rho, z_T);
    // group-level params for first grouping factor
    if (bK1 > 0)
      bMat1 = make_b_matrix(z_b, theta_L, p, l, 1);
    // group level params for second grouping factor
    if (bK2 > 0)
      bMat2 = make_b_matrix(z_b, theta_L, p, l, 2);
  }

  // group-level params, under lkj prior
  else if (prior_dist_for_cov == 2) {
    // group-level params for first grouping factor
    if (bK1 == 1)
      bMat1 = (bSd1[1] * z_bMat1)';
    else if (bK1 > 1)
      bMat1 = (diag_pre_multiply(bSd1, bCholesky1) * z_bMat1)';
    // group level params for second grouping factor
    if (bK2 == 1)
      bMat2 = (bSd2[1] * z_bMat2)';
    else if (bK2 > 1)
      bMat2 = (diag_pre_multiply(bSd2, bCholesky2) * z_bMat2)';
  }
