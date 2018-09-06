  vector[yNeta[1]] yEta1; // linear predictor
  vector[yNeta[2]] yEta2;
  vector[yNeta[3]] yEta3;

  // Linear predictor for submodel 1
  if (M > 0) {
    int bMat1_colshift = 0; // column shift in bMat1
    int bMat2_colshift = 0; // column shift in bMat2
    yEta1 = evaluate_eta(yX1, y1_Z1, y1_Z2, y1_Z1_id, y1_Z2_id, yGamma1, yBeta1,
                         bMat1, bMat2, bMat1_colshift, bMat2_colshift, intercept_type[1]);
  }

  // Linear predictor for submodel 2
  if (M > 1) {
    int bMat1_colshift = bK1_len[1]; // column shift in bMat1
    int bMat2_colshift = bK2_len[1]; // column shift in bMat2
    yEta2 = evaluate_eta(yX2, y2_Z1, y2_Z2, y2_Z1_id, y2_Z2_id, yGamma2, yBeta2,
                         bMat1, bMat2, bMat1_colshift, bMat2_colshift, intercept_type[2]);
  }

  // Linear predictor for submodel 3
  if (M > 2) {
    int bMat1_colshift = sum(bK1_len[1:2]); // column shift in bMat1
    int bMat2_colshift = sum(bK2_len[1:2]); // column shift in bMat2
    yEta3 = evaluate_eta(yX3, y3_Z1, y3_Z2, y3_Z1_id, y3_Z2_id, yGamma3, yBeta3,
                         bMat1, bMat2, bMat1_colshift, bMat2_colshift, intercept_type[3]);
  }

  // Log-likelihoods
  if (prior_PD == 0) {
    glm_lp(yReal1, yInt1, yEta1, yAux1, family[1], link[1], sum_log_y1, sqrt_y1, log_y1);
    if (M > 1)
      glm_lp(yReal2, yInt2, yEta2, yAux2, family[2], link[2], sum_log_y2, sqrt_y2, log_y2);
    if (M > 2)
      glm_lp(yReal3, yInt3, yEta3, yAux3, family[3], link[3], sum_log_y3, sqrt_y3, log_y3);
  }
