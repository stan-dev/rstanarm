  // prior family: 0 = none, 1 = normal, 2 = student_t, 3 = hs, 4 = hs_plus,
  //   5 = laplace, 6 = lasso
  int<lower=0,upper=6> a_prior_dist;

  //--- dimensions for association structure

    // num. of association parameters
    int<lower=0> a_K;

    // used for centering assoc terms
    vector[a_K] a_xbar;
    
    // used for scaling assoc terms
    vector[a_K] a_scale;

    // 0 = no assoc structure, 1 = any assoc structure
    int<lower=0,upper=1> assoc;

    // which components are required to build association terms
    array[6,20] int<lower=0,upper=1> assoc_uses;

    // which association terms does each submodel use
    array[16,M] int<lower=0,upper=1> has_assoc;

    // num. of shared random effects
    int<lower=0> sum_size_which_b;

    // num. of shared random effects for each long submodel
    array[M] int<lower=0> size_which_b;

    // which random effects are shared for each long submodel
    array[sum_size_which_b] int<lower=1> which_b_zindex;

    // num. of shared random effects incl fixed component
    int<lower=0> sum_size_which_coef;

    // num. of shared random effects incl fixed component for each long submodel
    array[M] int<lower=0> size_which_coef;

    // which random effects are shared incl fixed component
    array[sum_size_which_coef] int<lower=1> which_coef_zindex;

    // which fixed effects are shared
    array[sum_size_which_coef] int<lower=1> which_coef_xindex;

    // total num pars used in assoc*assoc interactions
    int<lower=0> sum_size_which_interactions;

    // num pars used in assoc*assoc interactions, by submodel
    //   and by evev/evmv/mvev/mvmv interactions
    array[M*4] int<lower=0,upper=sum_size_which_interactions> size_which_interactions;

    // which terms to interact with
    array[sum_size_which_interactions] int<lower=1> which_interactions;

  //---- data for calculating eta in GK quadrature

    array[20] int<lower=0> nrow_y_Xq;     // num. rows in long. predictor matrix at quadpoints

    // fe design matrix at quadpoints
    matrix[assoc_uses[1,1] == 1 ? nrow_y_Xq[1] : 0, yK[1]] y1_xq_eta;
    matrix[assoc_uses[1,2] == 1 ? nrow_y_Xq[2] : 0, yK[2]] y2_xq_eta;
    matrix[assoc_uses[1,3] == 1 ? nrow_y_Xq[3] : 0, yK[3]] y3_xq_eta;
    matrix[assoc_uses[1,4] == 1 ? nrow_y_Xq[4] : 0, yK[4]] y4_xq_eta;
    matrix[assoc_uses[1,5] == 1 ? nrow_y_Xq[5] : 0, yK[5]] y5_xq_eta;
    matrix[assoc_uses[1,6] == 1 ? nrow_y_Xq[6] : 0, yK[6]] y6_xq_eta;
    matrix[assoc_uses[1,7] == 1 ? nrow_y_Xq[7] : 0, yK[7]] y7_xq_eta;
    matrix[assoc_uses[1,8] == 1 ? nrow_y_Xq[8] : 0, yK[8]] y8_xq_eta;
    matrix[assoc_uses[1,9] == 1 ? nrow_y_Xq[9] : 0, yK[9]] y9_xq_eta;
    matrix[assoc_uses[1,10] == 1 ? nrow_y_Xq[10] : 0, yK[10]] y10_xq_eta;
    matrix[assoc_uses[1,11] == 1 ? nrow_y_Xq[11] : 0, yK[11]] y11_xq_eta;
    matrix[assoc_uses[1,12] == 1 ? nrow_y_Xq[12] : 0, yK[12]] y12_xq_eta;
    matrix[assoc_uses[1,13] == 1 ? nrow_y_Xq[13] : 0, yK[13]] y13_xq_eta;
    matrix[assoc_uses[1,14] == 1 ? nrow_y_Xq[14] : 0, yK[14]] y14_xq_eta;
    matrix[assoc_uses[1,15] == 1 ? nrow_y_Xq[15] : 0, yK[15]] y15_xq_eta;
    matrix[assoc_uses[1,16] == 1 ? nrow_y_Xq[16] : 0, yK[16]] y16_xq_eta;
    matrix[assoc_uses[1,17] == 1 ? nrow_y_Xq[17] : 0, yK[17]] y17_xq_eta;
    matrix[assoc_uses[1,18] == 1 ? nrow_y_Xq[18] : 0, yK[18]] y18_xq_eta;
    matrix[assoc_uses[1,19] == 1 ? nrow_y_Xq[19] : 0, yK[19]] y19_xq_eta;
    matrix[assoc_uses[1,20] == 1 ? nrow_y_Xq[20] : 0, yK[20]] y20_xq_eta;
    
    // offset values at quadpoints
    vector[has_offset[1] && assoc_uses[1,1] == 1 ? nrow_y_Xq[1] : 0] y1_offset_eta;
    vector[has_offset[2] && assoc_uses[1,2] == 1 ? nrow_y_Xq[2] : 0] y2_offset_eta;
    vector[has_offset[3] && assoc_uses[1,3] == 1 ? nrow_y_Xq[3] : 0] y3_offset_eta;
    vector[has_offset[4] && assoc_uses[1,4] == 1 ? nrow_y_Xq[4] : 0] y4_offset_eta;
    vector[has_offset[5] && assoc_uses[1,5] == 1 ? nrow_y_Xq[5] : 0] y5_offset_eta;
    vector[has_offset[6] && assoc_uses[1,6] == 1 ? nrow_y_Xq[6] : 0] y6_offset_eta;
    vector[has_offset[7] && assoc_uses[1,7] == 1 ? nrow_y_Xq[7] : 0] y7_offset_eta;
    vector[has_offset[8] && assoc_uses[1,8] == 1 ? nrow_y_Xq[8] : 0] y8_offset_eta;
    vector[has_offset[9] && assoc_uses[1,9] == 1 ? nrow_y_Xq[9] : 0] y9_offset_eta;
    vector[has_offset[10] && assoc_uses[1,10] == 1 ? nrow_y_Xq[10] : 0] y10_offset_eta;
    vector[has_offset[11] && assoc_uses[1,11] == 1 ? nrow_y_Xq[11] : 0] y11_offset_eta;
    vector[has_offset[12] && assoc_uses[1,12] == 1 ? nrow_y_Xq[12] : 0] y12_offset_eta;
    vector[has_offset[13] && assoc_uses[1,13] == 1 ? nrow_y_Xq[13] : 0] y13_offset_eta;
    vector[has_offset[14] && assoc_uses[1,14] == 1 ? nrow_y_Xq[14] : 0] y14_offset_eta;
    vector[has_offset[15] && assoc_uses[1,15] == 1 ? nrow_y_Xq[15] : 0] y15_offset_eta;
    vector[has_offset[16] && assoc_uses[1,16] == 1 ? nrow_y_Xq[16] : 0] y16_offset_eta;
    vector[has_offset[17] && assoc_uses[1,17] == 1 ? nrow_y_Xq[17] : 0] y17_offset_eta;
    vector[has_offset[18] && assoc_uses[1,18] == 1 ? nrow_y_Xq[18] : 0] y18_offset_eta;
    vector[has_offset[19] && assoc_uses[1,19] == 1 ? nrow_y_Xq[19] : 0] y19_offset_eta;
    vector[has_offset[20] && assoc_uses[1,20] == 1 ? nrow_y_Xq[20] : 0] y20_offset_eta;

    // re design matrix at quadpoints, group factor 1
    array[bK1_len[1]] vector[assoc_uses[1,1] == 1 && bK1_len[1] > 0 ? nrow_y_Xq[1] : 0] y1_z1q_eta;
    array[bK1_len[2]] vector[assoc_uses[1,2] == 1 && bK1_len[2] > 0 ? nrow_y_Xq[2] : 0] y2_z1q_eta;
    array[bK1_len[3]] vector[assoc_uses[1,3] == 1 && bK1_len[3] > 0 ? nrow_y_Xq[3] : 0] y3_z1q_eta;
    array[bK1_len[4]] vector[assoc_uses[1,4] == 1 && bK1_len[4] > 0 ? nrow_y_Xq[4] : 0] y4_z1q_eta;
    array[bK1_len[5]] vector[assoc_uses[1,5] == 1 && bK1_len[5] > 0 ? nrow_y_Xq[5] : 0] y5_z1q_eta;
    array[bK1_len[6]] vector[assoc_uses[1,6] == 1 && bK1_len[6] > 0 ? nrow_y_Xq[6] : 0] y6_z1q_eta;
    array[bK1_len[7]] vector[assoc_uses[1,7] == 1 && bK1_len[7] > 0 ? nrow_y_Xq[7] : 0] y7_z1q_eta;
    array[bK1_len[8]] vector[assoc_uses[1,8] == 1 && bK1_len[8] > 0 ? nrow_y_Xq[8] : 0] y8_z1q_eta;
    array[bK1_len[9]] vector[assoc_uses[1,9] == 1 && bK1_len[9] > 0 ? nrow_y_Xq[9] : 0] y9_z1q_eta;
    array[bK1_len[10]] vector[assoc_uses[1,10] == 1 && bK1_len[10] > 0 ? nrow_y_Xq[10] : 0] y10_z1q_eta;
    array[bK1_len[11]] vector[assoc_uses[1,11] == 1 && bK1_len[11] > 0 ? nrow_y_Xq[11] : 0] y11_z1q_eta;
    array[bK1_len[12]] vector[assoc_uses[1,12] == 1 && bK1_len[12] > 0 ? nrow_y_Xq[12] : 0] y12_z1q_eta;
    array[bK1_len[13]] vector[assoc_uses[1,13] == 1 && bK1_len[13] > 0 ? nrow_y_Xq[13] : 0] y13_z1q_eta;
    array[bK1_len[14]] vector[assoc_uses[1,14] == 1 && bK1_len[14] > 0 ? nrow_y_Xq[14] : 0] y14_z1q_eta;
    array[bK1_len[15]] vector[assoc_uses[1,15] == 1 && bK1_len[15] > 0 ? nrow_y_Xq[15] : 0] y15_z1q_eta;
    array[bK1_len[16]] vector[assoc_uses[1,16] == 1 && bK1_len[16] > 0 ? nrow_y_Xq[16] : 0] y16_z1q_eta;
    array[bK1_len[17]] vector[assoc_uses[1,17] == 1 && bK1_len[17] > 0 ? nrow_y_Xq[17] : 0] y17_z1q_eta;
    array[bK1_len[18]] vector[assoc_uses[1,18] == 1 && bK1_len[18] > 0 ? nrow_y_Xq[18] : 0] y18_z1q_eta;
    array[bK1_len[19]] vector[assoc_uses[1,19] == 1 && bK1_len[19] > 0 ? nrow_y_Xq[19] : 0] y19_z1q_eta;
    array[bK1_len[20]] vector[assoc_uses[1,20] == 1 && bK1_len[20] > 0 ? nrow_y_Xq[20] : 0] y20_z1q_eta;
    array[assoc_uses[1,1] == 1 && bK1_len[1] > 0 ? nrow_y_Xq[1] : 0] int<lower=0> y1_z1q_id_eta;
    array[assoc_uses[1,2] == 1 && bK1_len[2] > 0 ? nrow_y_Xq[2] : 0] int<lower=0> y2_z1q_id_eta;
    array[assoc_uses[1,3] == 1 && bK1_len[3] > 0 ? nrow_y_Xq[3] : 0] int<lower=0> y3_z1q_id_eta;
    array[assoc_uses[1,4] == 1 && bK1_len[4] > 0 ? nrow_y_Xq[4] : 0] int<lower=0> y4_z1q_id_eta;
    array[assoc_uses[1,5] == 1 && bK1_len[5] > 0 ? nrow_y_Xq[5] : 0] int<lower=0> y5_z1q_id_eta;
    array[assoc_uses[1,6] == 1 && bK1_len[6] > 0 ? nrow_y_Xq[6] : 0] int<lower=0> y6_z1q_id_eta;
    array[assoc_uses[1,7] == 1 && bK1_len[7] > 0 ? nrow_y_Xq[7] : 0] int<lower=0> y7_z1q_id_eta;
    array[assoc_uses[1,8] == 1 && bK1_len[8] > 0 ? nrow_y_Xq[8] : 0] int<lower=0> y8_z1q_id_eta;
    array[assoc_uses[1,9] == 1 && bK1_len[9] > 0 ? nrow_y_Xq[9] : 0] int<lower=0> y9_z1q_id_eta;
    array[assoc_uses[1,10] == 1 && bK1_len[10] > 0 ? nrow_y_Xq[10] : 0] int<lower=0> y10_z1q_id_eta;
    array[assoc_uses[1,11] == 1 && bK1_len[11] > 0 ? nrow_y_Xq[11] : 0] int<lower=0> y11_z1q_id_eta;
    array[assoc_uses[1,12] == 1 && bK1_len[12] > 0 ? nrow_y_Xq[12] : 0] int<lower=0> y12_z1q_id_eta;
    array[assoc_uses[1,13] == 1 && bK1_len[13] > 0 ? nrow_y_Xq[13] : 0] int<lower=0> y13_z1q_id_eta;
    array[assoc_uses[1,14] == 1 && bK1_len[14] > 0 ? nrow_y_Xq[14] : 0] int<lower=0> y14_z1q_id_eta;
    array[assoc_uses[1,15] == 1 && bK1_len[15] > 0 ? nrow_y_Xq[15] : 0] int<lower=0> y15_z1q_id_eta;
    array[assoc_uses[1,16] == 1 && bK1_len[16] > 0 ? nrow_y_Xq[16] : 0] int<lower=0> y16_z1q_id_eta;
    array[assoc_uses[1,17] == 1 && bK1_len[17] > 0 ? nrow_y_Xq[17] : 0] int<lower=0> y17_z1q_id_eta;
    array[assoc_uses[1,18] == 1 && bK1_len[18] > 0 ? nrow_y_Xq[18] : 0] int<lower=0> y18_z1q_id_eta;
    array[assoc_uses[1,19] == 1 && bK1_len[19] > 0 ? nrow_y_Xq[19] : 0] int<lower=0> y19_z1q_id_eta;
    array[assoc_uses[1,20] == 1 && bK1_len[20] > 0 ? nrow_y_Xq[20] : 0] int<lower=0> y20_z1q_id_eta;

    // re design matrix at quadpoints, group factor 2
    array[bK2_len[1]] vector[assoc_uses[1,1] == 1 && bK2_len[1] > 0 ? nrow_y_Xq[1] : 0] y1_z2q_eta;
    array[bK2_len[2]] vector[assoc_uses[1,2] == 1 && bK2_len[2] > 0 ? nrow_y_Xq[2] : 0] y2_z2q_eta;
    array[bK2_len[3]] vector[assoc_uses[1,3] == 1 && bK2_len[3] > 0 ? nrow_y_Xq[3] : 0] y3_z2q_eta;
    array[bK2_len[4]] vector[assoc_uses[1,4] == 1 && bK2_len[4] > 0 ? nrow_y_Xq[4] : 0] y4_z2q_eta;
    array[bK2_len[5]] vector[assoc_uses[1,5] == 1 && bK2_len[5] > 0 ? nrow_y_Xq[5] : 0] y5_z2q_eta;
    array[bK2_len[6]] vector[assoc_uses[1,6] == 1 && bK2_len[6] > 0 ? nrow_y_Xq[6] : 0] y6_z2q_eta;
    array[bK2_len[7]] vector[assoc_uses[1,7] == 1 && bK2_len[7] > 0 ? nrow_y_Xq[7] : 0] y7_z2q_eta;
    array[bK2_len[8]] vector[assoc_uses[1,8] == 1 && bK2_len[8] > 0 ? nrow_y_Xq[8] : 0] y8_z2q_eta;
    array[bK2_len[9]] vector[assoc_uses[1,9] == 1 && bK2_len[9] > 0 ? nrow_y_Xq[9] : 0] y9_z2q_eta;
    array[bK2_len[10]] vector[assoc_uses[1,10] == 1 && bK2_len[10] > 0 ? nrow_y_Xq[10] : 0] y10_z2q_eta;
    array[bK2_len[11]] vector[assoc_uses[1,11] == 1 && bK2_len[11] > 0 ? nrow_y_Xq[11] : 0] y11_z2q_eta;
    array[bK2_len[12]] vector[assoc_uses[1,12] == 1 && bK2_len[12] > 0 ? nrow_y_Xq[12] : 0] y12_z2q_eta;
    array[bK2_len[13]] vector[assoc_uses[1,13] == 1 && bK2_len[13] > 0 ? nrow_y_Xq[13] : 0] y13_z2q_eta;
    array[bK2_len[14]] vector[assoc_uses[1,14] == 1 && bK2_len[14] > 0 ? nrow_y_Xq[14] : 0] y14_z2q_eta;
    array[bK2_len[15]] vector[assoc_uses[1,15] == 1 && bK2_len[15] > 0 ? nrow_y_Xq[15] : 0] y15_z2q_eta;
    array[bK2_len[16]] vector[assoc_uses[1,16] == 1 && bK2_len[16] > 0 ? nrow_y_Xq[16] : 0] y16_z2q_eta;
    array[bK2_len[17]] vector[assoc_uses[1,17] == 1 && bK2_len[17] > 0 ? nrow_y_Xq[17] : 0] y17_z2q_eta;
    array[bK2_len[18]] vector[assoc_uses[1,18] == 1 && bK2_len[18] > 0 ? nrow_y_Xq[18] : 0] y18_z2q_eta;
    array[bK2_len[19]] vector[assoc_uses[1,19] == 1 && bK2_len[19] > 0 ? nrow_y_Xq[19] : 0] y19_z2q_eta;
    array[bK2_len[20]] vector[assoc_uses[1,20] == 1 && bK2_len[20] > 0 ? nrow_y_Xq[20] : 0] y20_z2q_eta;
    array[assoc_uses[1,1] == 1 && bK2_len[1] > 0 ? nrow_y_Xq[1] : 0] int<lower=0> y1_z2q_id_eta;
    array[assoc_uses[1,2] == 1 && bK2_len[2] > 0 ? nrow_y_Xq[2] : 0] int<lower=0> y2_z2q_id_eta;
    array[assoc_uses[1,3] == 1 && bK2_len[3] > 0 ? nrow_y_Xq[3] : 0] int<lower=0> y3_z2q_id_eta;
    array[assoc_uses[1,4] == 1 && bK2_len[4] > 0 ? nrow_y_Xq[4] : 0] int<lower=0> y4_z2q_id_eta;
    array[assoc_uses[1,5] == 1 && bK2_len[5] > 0 ? nrow_y_Xq[5] : 0] int<lower=0> y5_z2q_id_eta;
    array[assoc_uses[1,6] == 1 && bK2_len[6] > 0 ? nrow_y_Xq[6] : 0] int<lower=0> y6_z2q_id_eta;
    array[assoc_uses[1,7] == 1 && bK2_len[7] > 0 ? nrow_y_Xq[7] : 0] int<lower=0> y7_z2q_id_eta;
    array[assoc_uses[1,8] == 1 && bK2_len[8] > 0 ? nrow_y_Xq[8] : 0] int<lower=0> y8_z2q_id_eta;
    array[assoc_uses[1,9] == 1 && bK2_len[9] > 0 ? nrow_y_Xq[9] : 0] int<lower=0> y9_z2q_id_eta;
    array[assoc_uses[1,10] == 1 && bK2_len[10] > 0 ? nrow_y_Xq[10] : 0] int<lower=0> y10_z2q_id_eta;
    array[assoc_uses[1,11] == 1 && bK2_len[11] > 0 ? nrow_y_Xq[11] : 0] int<lower=0> y11_z2q_id_eta;
    array[assoc_uses[1,12] == 1 && bK2_len[12] > 0 ? nrow_y_Xq[12] : 0] int<lower=0> y12_z2q_id_eta;
    array[assoc_uses[1,13] == 1 && bK2_len[13] > 0 ? nrow_y_Xq[13] : 0] int<lower=0> y13_z2q_id_eta;
    array[assoc_uses[1,14] == 1 && bK2_len[14] > 0 ? nrow_y_Xq[14] : 0] int<lower=0> y14_z2q_id_eta;
    array[assoc_uses[1,15] == 1 && bK2_len[15] > 0 ? nrow_y_Xq[15] : 0] int<lower=0> y15_z2q_id_eta;
    array[assoc_uses[1,16] == 1 && bK2_len[16] > 0 ? nrow_y_Xq[16] : 0] int<lower=0> y16_z2q_id_eta;
    array[assoc_uses[1,17] == 1 && bK2_len[17] > 0 ? nrow_y_Xq[17] : 0] int<lower=0> y17_z2q_id_eta;
    array[assoc_uses[1,18] == 1 && bK2_len[18] > 0 ? nrow_y_Xq[18] : 0] int<lower=0> y18_z2q_id_eta;
    array[assoc_uses[1,19] == 1 && bK2_len[19] > 0 ? nrow_y_Xq[19] : 0] int<lower=0> y19_z2q_id_eta;
    array[assoc_uses[1,20] == 1 && bK2_len[20] > 0 ? nrow_y_Xq[20] : 0] int<lower=0> y20_z2q_id_eta;

  //---- data for calculating derivative of eta in GK quadrature

    // fe design matrix at quadpoints
    matrix[assoc_uses[2,1] == 1 ? nrow_y_Xq[1] : 0, yK[1]] y1_xq_eps;
    matrix[assoc_uses[2,2] == 1 ? nrow_y_Xq[2] : 0, yK[2]] y2_xq_eps;
    matrix[assoc_uses[2,3] == 1 ? nrow_y_Xq[3] : 0, yK[3]] y3_xq_eps;
    matrix[assoc_uses[2,4] == 1 ? nrow_y_Xq[4] : 0, yK[4]] y4_xq_eps;
    matrix[assoc_uses[2,5] == 1 ? nrow_y_Xq[5] : 0, yK[5]] y5_xq_eps;
    matrix[assoc_uses[2,6] == 1 ? nrow_y_Xq[6] : 0, yK[6]] y6_xq_eps;
    matrix[assoc_uses[2,7] == 1 ? nrow_y_Xq[7] : 0, yK[7]] y7_xq_eps;
    matrix[assoc_uses[2,8] == 1 ? nrow_y_Xq[8] : 0, yK[8]] y8_xq_eps;
    matrix[assoc_uses[2,9] == 1 ? nrow_y_Xq[9] : 0, yK[9]] y9_xq_eps;
    matrix[assoc_uses[2,10] == 1 ? nrow_y_Xq[10] : 0, yK[10]] y10_xq_eps;
    matrix[assoc_uses[2,11] == 1 ? nrow_y_Xq[11] : 0, yK[11]] y11_xq_eps;
    matrix[assoc_uses[2,12] == 1 ? nrow_y_Xq[12] : 0, yK[12]] y12_xq_eps;
    matrix[assoc_uses[2,13] == 1 ? nrow_y_Xq[13] : 0, yK[13]] y13_xq_eps;
    matrix[assoc_uses[2,14] == 1 ? nrow_y_Xq[14] : 0, yK[14]] y14_xq_eps;
    matrix[assoc_uses[2,15] == 1 ? nrow_y_Xq[15] : 0, yK[15]] y15_xq_eps;
    matrix[assoc_uses[2,16] == 1 ? nrow_y_Xq[16] : 0, yK[16]] y16_xq_eps;
    matrix[assoc_uses[2,17] == 1 ? nrow_y_Xq[17] : 0, yK[17]] y17_xq_eps;
    matrix[assoc_uses[2,18] == 1 ? nrow_y_Xq[18] : 0, yK[18]] y18_xq_eps;
    matrix[assoc_uses[2,19] == 1 ? nrow_y_Xq[19] : 0, yK[19]] y19_xq_eps;
    matrix[assoc_uses[2,20] == 1 ? nrow_y_Xq[20] : 0, yK[20]] y20_xq_eps;
    
    // offset values at quadpoints
    vector[has_offset[1] && assoc_uses[2,1] == 1 ? nrow_y_Xq[1] : 0] y1_offset_eps;
    vector[has_offset[2] && assoc_uses[2,2] == 1 ? nrow_y_Xq[2] : 0] y2_offset_eps;
    vector[has_offset[3] && assoc_uses[2,3] == 1 ? nrow_y_Xq[3] : 0] y3_offset_eps;
    vector[has_offset[4] && assoc_uses[2,4] == 1 ? nrow_y_Xq[4] : 0] y4_offset_eps;
    vector[has_offset[5] && assoc_uses[2,5] == 1 ? nrow_y_Xq[5] : 0] y5_offset_eps;
    vector[has_offset[6] && assoc_uses[2,6] == 1 ? nrow_y_Xq[6] : 0] y6_offset_eps;
    vector[has_offset[7] && assoc_uses[2,7] == 1 ? nrow_y_Xq[7] : 0] y7_offset_eps;
    vector[has_offset[8] && assoc_uses[2,8] == 1 ? nrow_y_Xq[8] : 0] y8_offset_eps;
    vector[has_offset[9] && assoc_uses[2,9] == 1 ? nrow_y_Xq[9] : 0] y9_offset_eps;
    vector[has_offset[10] && assoc_uses[2,10] == 1 ? nrow_y_Xq[10] : 0] y10_offset_eps;
    vector[has_offset[11] && assoc_uses[2,11] == 1 ? nrow_y_Xq[11] : 0] y11_offset_eps;
    vector[has_offset[12] && assoc_uses[2,12] == 1 ? nrow_y_Xq[12] : 0] y12_offset_eps;
    vector[has_offset[13] && assoc_uses[2,13] == 1 ? nrow_y_Xq[13] : 0] y13_offset_eps;
    vector[has_offset[14] && assoc_uses[2,14] == 1 ? nrow_y_Xq[14] : 0] y14_offset_eps;
    vector[has_offset[15] && assoc_uses[2,15] == 1 ? nrow_y_Xq[15] : 0] y15_offset_eps;
    vector[has_offset[16] && assoc_uses[2,16] == 1 ? nrow_y_Xq[16] : 0] y16_offset_eps;
    vector[has_offset[17] && assoc_uses[2,17] == 1 ? nrow_y_Xq[17] : 0] y17_offset_eps;
    vector[has_offset[18] && assoc_uses[2,18] == 1 ? nrow_y_Xq[18] : 0] y18_offset_eps;
    vector[has_offset[19] && assoc_uses[2,19] == 1 ? nrow_y_Xq[19] : 0] y19_offset_eps;
    vector[has_offset[20] && assoc_uses[2,20] == 1 ? nrow_y_Xq[20] : 0] y20_offset_eps;

    // re design matrix at quadpoints, group factor 1
    array[bK1_len[1]] vector[assoc_uses[2,1] == 1 && bK1_len[1] > 0 ? nrow_y_Xq[1] : 0] y1_z1q_eps;
    array[bK1_len[2]] vector[assoc_uses[2,2] == 1 && bK1_len[2] > 0 ? nrow_y_Xq[2] : 0] y2_z1q_eps;
    array[bK1_len[3]] vector[assoc_uses[2,3] == 1 && bK1_len[3] > 0 ? nrow_y_Xq[3] : 0] y3_z1q_eps;
    array[bK1_len[4]] vector[assoc_uses[2,4] == 1 && bK1_len[4] > 0 ? nrow_y_Xq[4] : 0] y4_z1q_eps;
    array[bK1_len[5]] vector[assoc_uses[2,5] == 1 && bK1_len[5] > 0 ? nrow_y_Xq[5] : 0] y5_z1q_eps;
    array[bK1_len[6]] vector[assoc_uses[2,6] == 1 && bK1_len[6] > 0 ? nrow_y_Xq[6] : 0] y6_z1q_eps;
    array[bK1_len[7]] vector[assoc_uses[2,7] == 1 && bK1_len[7] > 0 ? nrow_y_Xq[7] : 0] y7_z1q_eps;
    array[bK1_len[8]] vector[assoc_uses[2,8] == 1 && bK1_len[8] > 0 ? nrow_y_Xq[8] : 0] y8_z1q_eps;
    array[bK1_len[9]] vector[assoc_uses[2,9] == 1 && bK1_len[9] > 0 ? nrow_y_Xq[9] : 0] y9_z1q_eps;
    array[bK1_len[10]] vector[assoc_uses[2,10] == 1 && bK1_len[10] > 0 ? nrow_y_Xq[10] : 0] y10_z1q_eps;
    array[bK1_len[11]] vector[assoc_uses[2,11] == 1 && bK1_len[11] > 0 ? nrow_y_Xq[11] : 0] y11_z1q_eps;
    array[bK1_len[12]] vector[assoc_uses[2,12] == 1 && bK1_len[12] > 0 ? nrow_y_Xq[12] : 0] y12_z1q_eps;
    array[bK1_len[13]] vector[assoc_uses[2,13] == 1 && bK1_len[13] > 0 ? nrow_y_Xq[13] : 0] y13_z1q_eps;
    array[bK1_len[14]] vector[assoc_uses[2,14] == 1 && bK1_len[14] > 0 ? nrow_y_Xq[14] : 0] y14_z1q_eps;
    array[bK1_len[15]] vector[assoc_uses[2,15] == 1 && bK1_len[15] > 0 ? nrow_y_Xq[15] : 0] y15_z1q_eps;
    array[bK1_len[16]] vector[assoc_uses[2,16] == 1 && bK1_len[16] > 0 ? nrow_y_Xq[16] : 0] y16_z1q_eps;
    array[bK1_len[17]] vector[assoc_uses[2,17] == 1 && bK1_len[17] > 0 ? nrow_y_Xq[17] : 0] y17_z1q_eps;
    array[bK1_len[18]] vector[assoc_uses[2,18] == 1 && bK1_len[18] > 0 ? nrow_y_Xq[18] : 0] y18_z1q_eps;
    array[bK1_len[19]] vector[assoc_uses[2,19] == 1 && bK1_len[19] > 0 ? nrow_y_Xq[19] : 0] y19_z1q_eps;
    array[bK1_len[20]] vector[assoc_uses[2,20] == 1 && bK1_len[20] > 0 ? nrow_y_Xq[20] : 0] y20_z1q_eps;
    array[assoc_uses[2,1] == 1 && bK1_len[1] > 0 ? nrow_y_Xq[1] : 0] int<lower=0> y1_z1q_id_eps;
    array[assoc_uses[2,2] == 1 && bK1_len[2] > 0 ? nrow_y_Xq[2] : 0] int<lower=0> y2_z1q_id_eps;
    array[assoc_uses[2,3] == 1 && bK1_len[3] > 0 ? nrow_y_Xq[3] : 0] int<lower=0> y3_z1q_id_eps;
    array[assoc_uses[2,4] == 1 && bK1_len[4] > 0 ? nrow_y_Xq[4] : 0] int<lower=0> y4_z1q_id_eps;
    array[assoc_uses[2,5] == 1 && bK1_len[5] > 0 ? nrow_y_Xq[5] : 0] int<lower=0> y5_z1q_id_eps;
    array[assoc_uses[2,6] == 1 && bK1_len[6] > 0 ? nrow_y_Xq[6] : 0] int<lower=0> y6_z1q_id_eps;
    array[assoc_uses[2,7] == 1 && bK1_len[7] > 0 ? nrow_y_Xq[7] : 0] int<lower=0> y7_z1q_id_eps;
    array[assoc_uses[2,8] == 1 && bK1_len[8] > 0 ? nrow_y_Xq[8] : 0] int<lower=0> y8_z1q_id_eps;
    array[assoc_uses[2,9] == 1 && bK1_len[9] > 0 ? nrow_y_Xq[9] : 0] int<lower=0> y9_z1q_id_eps;
    array[assoc_uses[2,10] == 1 && bK1_len[10] > 0 ? nrow_y_Xq[10] : 0] int<lower=0> y10_z1q_id_eps;
    array[assoc_uses[2,11] == 1 && bK1_len[11] > 0 ? nrow_y_Xq[11] : 0] int<lower=0> y11_z1q_id_eps;
    array[assoc_uses[2,12] == 1 && bK1_len[12] > 0 ? nrow_y_Xq[12] : 0] int<lower=0> y12_z1q_id_eps;
    array[assoc_uses[2,13] == 1 && bK1_len[13] > 0 ? nrow_y_Xq[13] : 0] int<lower=0> y13_z1q_id_eps;
    array[assoc_uses[2,14] == 1 && bK1_len[14] > 0 ? nrow_y_Xq[14] : 0] int<lower=0> y14_z1q_id_eps;
    array[assoc_uses[2,15] == 1 && bK1_len[15] > 0 ? nrow_y_Xq[15] : 0] int<lower=0> y15_z1q_id_eps;
    array[assoc_uses[2,16] == 1 && bK1_len[16] > 0 ? nrow_y_Xq[16] : 0] int<lower=0> y16_z1q_id_eps;
    array[assoc_uses[2,17] == 1 && bK1_len[17] > 0 ? nrow_y_Xq[17] : 0] int<lower=0> y17_z1q_id_eps;
    array[assoc_uses[2,18] == 1 && bK1_len[18] > 0 ? nrow_y_Xq[18] : 0] int<lower=0> y18_z1q_id_eps;
    array[assoc_uses[2,19] == 1 && bK1_len[19] > 0 ? nrow_y_Xq[19] : 0] int<lower=0> y19_z1q_id_eps;
    array[assoc_uses[2,20] == 1 && bK1_len[20] > 0 ? nrow_y_Xq[20] : 0] int<lower=0> y20_z1q_id_eps;

    // re design matrix at quadpoints, group factor 2
    array[bK2_len[1]] vector[assoc_uses[2,1] == 1 && bK2_len[1] > 0 ? nrow_y_Xq[1] : 0] y1_z2q_eps;
    array[bK2_len[2]] vector[assoc_uses[2,2] == 1 && bK2_len[2] > 0 ? nrow_y_Xq[2] : 0] y2_z2q_eps;
    array[bK2_len[3]] vector[assoc_uses[2,3] == 1 && bK2_len[3] > 0 ? nrow_y_Xq[3] : 0] y3_z2q_eps;
    array[bK2_len[4]] vector[assoc_uses[2,4] == 1 && bK2_len[4] > 0 ? nrow_y_Xq[4] : 0] y4_z2q_eps;
    array[bK2_len[5]] vector[assoc_uses[2,5] == 1 && bK2_len[5] > 0 ? nrow_y_Xq[5] : 0] y5_z2q_eps;
    array[bK2_len[6]] vector[assoc_uses[2,6] == 1 && bK2_len[6] > 0 ? nrow_y_Xq[6] : 0] y6_z2q_eps;
    array[bK2_len[7]] vector[assoc_uses[2,7] == 1 && bK2_len[7] > 0 ? nrow_y_Xq[7] : 0] y7_z2q_eps;
    array[bK2_len[8]] vector[assoc_uses[2,8] == 1 && bK2_len[8] > 0 ? nrow_y_Xq[8] : 0] y8_z2q_eps;
    array[bK2_len[9]] vector[assoc_uses[2,9] == 1 && bK2_len[9] > 0 ? nrow_y_Xq[9] : 0] y9_z2q_eps;
    array[bK2_len[10]] vector[assoc_uses[2,10] == 1 && bK2_len[10] > 0 ? nrow_y_Xq[10] : 0] y10_z2q_eps;
    array[bK2_len[11]] vector[assoc_uses[2,11] == 1 && bK2_len[11] > 0 ? nrow_y_Xq[11] : 0] y11_z2q_eps;
    array[bK2_len[12]] vector[assoc_uses[2,12] == 1 && bK2_len[12] > 0 ? nrow_y_Xq[12] : 0] y12_z2q_eps;
    array[bK2_len[13]] vector[assoc_uses[2,13] == 1 && bK2_len[13] > 0 ? nrow_y_Xq[13] : 0] y13_z2q_eps;
    array[bK2_len[14]] vector[assoc_uses[2,14] == 1 && bK2_len[14] > 0 ? nrow_y_Xq[14] : 0] y14_z2q_eps;
    array[bK2_len[15]] vector[assoc_uses[2,15] == 1 && bK2_len[15] > 0 ? nrow_y_Xq[15] : 0] y15_z2q_eps;
    array[bK2_len[16]] vector[assoc_uses[2,16] == 1 && bK2_len[16] > 0 ? nrow_y_Xq[16] : 0] y16_z2q_eps;
    array[bK2_len[17]] vector[assoc_uses[2,17] == 1 && bK2_len[17] > 0 ? nrow_y_Xq[17] : 0] y17_z2q_eps;
    array[bK2_len[18]] vector[assoc_uses[2,18] == 1 && bK2_len[18] > 0 ? nrow_y_Xq[18] : 0] y18_z2q_eps;
    array[bK2_len[19]] vector[assoc_uses[2,19] == 1 && bK2_len[19] > 0 ? nrow_y_Xq[19] : 0] y19_z2q_eps;
    array[bK2_len[20]] vector[assoc_uses[2,20] == 1 && bK2_len[20] > 0 ? nrow_y_Xq[20] : 0] y20_z2q_eps;
    array[assoc_uses[2,1] == 1 && bK2_len[1] > 0 ? nrow_y_Xq[1] : 0] int<lower=0> y1_z2q_id_eps;
    array[assoc_uses[2,2] == 1 && bK2_len[2] > 0 ? nrow_y_Xq[2] : 0] int<lower=0> y2_z2q_id_eps;
    array[assoc_uses[2,3] == 1 && bK2_len[3] > 0 ? nrow_y_Xq[3] : 0] int<lower=0> y3_z2q_id_eps;
    array[assoc_uses[2,4] == 1 && bK2_len[4] > 0 ? nrow_y_Xq[4] : 0] int<lower=0> y4_z2q_id_eps;
    array[assoc_uses[2,5] == 1 && bK2_len[5] > 0 ? nrow_y_Xq[5] : 0] int<lower=0> y5_z2q_id_eps;
    array[assoc_uses[2,6] == 1 && bK2_len[6] > 0 ? nrow_y_Xq[6] : 0] int<lower=0> y6_z2q_id_eps;
    array[assoc_uses[2,7] == 1 && bK2_len[7] > 0 ? nrow_y_Xq[7] : 0] int<lower=0> y7_z2q_id_eps;
    array[assoc_uses[2,8] == 1 && bK2_len[8] > 0 ? nrow_y_Xq[8] : 0] int<lower=0> y8_z2q_id_eps;
    array[assoc_uses[2,9] == 1 && bK2_len[9] > 0 ? nrow_y_Xq[9] : 0] int<lower=0> y9_z2q_id_eps;
    array[assoc_uses[2,10] == 1 && bK2_len[10] > 0 ? nrow_y_Xq[10] : 0] int<lower=0> y10_z2q_id_eps;
    array[assoc_uses[2,11] == 1 && bK2_len[11] > 0 ? nrow_y_Xq[11] : 0] int<lower=0> y11_z2q_id_eps;
    array[assoc_uses[2,12] == 1 && bK2_len[12] > 0 ? nrow_y_Xq[12] : 0] int<lower=0> y12_z2q_id_eps;
    array[assoc_uses[2,13] == 1 && bK2_len[13] > 0 ? nrow_y_Xq[13] : 0] int<lower=0> y13_z2q_id_eps;
    array[assoc_uses[2,14] == 1 && bK2_len[14] > 0 ? nrow_y_Xq[14] : 0] int<lower=0> y14_z2q_id_eps;
    array[assoc_uses[2,15] == 1 && bK2_len[15] > 0 ? nrow_y_Xq[15] : 0] int<lower=0> y15_z2q_id_eps;
    array[assoc_uses[2,16] == 1 && bK2_len[16] > 0 ? nrow_y_Xq[16] : 0] int<lower=0> y16_z2q_id_eps;
    array[assoc_uses[2,17] == 1 && bK2_len[17] > 0 ? nrow_y_Xq[17] : 0] int<lower=0> y17_z2q_id_eps;
    array[assoc_uses[2,18] == 1 && bK2_len[18] > 0 ? nrow_y_Xq[18] : 0] int<lower=0> y18_z2q_id_eps;
    array[assoc_uses[2,19] == 1 && bK2_len[19] > 0 ? nrow_y_Xq[19] : 0] int<lower=0> y19_z2q_id_eps;
    array[assoc_uses[2,20] == 1 && bK2_len[20] > 0 ? nrow_y_Xq[20] : 0] int<lower=0> y20_z2q_id_eps;

  //---- data for calculating integral of eta in GK quadrature

    // num. of nodes for GK quadrature for area under marker trajectory
    int<lower=0> auc_qnodes;
    int<lower=0> nrow_y_Xq_auc; // num. rows in long. predictor matrix at auc quadpoints
    vector[sum(assoc_uses[3,]) > 0 ? nrow_y_Xq_auc : 0] auc_qwts;

    // fe design matrix at quadpoints
    matrix[assoc_uses[3,1] == 1 ? nrow_y_Xq_auc : 0, yK[1]] y1_xq_auc;
    matrix[assoc_uses[3,2] == 1 ? nrow_y_Xq_auc : 0, yK[2]] y2_xq_auc;
    matrix[assoc_uses[3,3] == 1 ? nrow_y_Xq_auc : 0, yK[3]] y3_xq_auc;
    matrix[assoc_uses[3,4] == 1 ? nrow_y_Xq_auc : 0, yK[4]] y4_xq_auc;
    matrix[assoc_uses[3,5] == 1 ? nrow_y_Xq_auc : 0, yK[5]] y5_xq_auc;
    matrix[assoc_uses[3,6] == 1 ? nrow_y_Xq_auc : 0, yK[6]] y6_xq_auc;
    matrix[assoc_uses[3,7] == 1 ? nrow_y_Xq_auc : 0, yK[7]] y7_xq_auc;
    matrix[assoc_uses[3,8] == 1 ? nrow_y_Xq_auc : 0, yK[8]] y8_xq_auc;
    matrix[assoc_uses[3,9] == 1 ? nrow_y_Xq_auc : 0, yK[9]] y9_xq_auc;
    matrix[assoc_uses[3,10] == 1 ? nrow_y_Xq_auc : 0, yK[10]] y10_xq_auc;
    matrix[assoc_uses[3,11] == 1 ? nrow_y_Xq_auc : 0, yK[11]] y11_xq_auc;
    matrix[assoc_uses[3,12] == 1 ? nrow_y_Xq_auc : 0, yK[12]] y12_xq_auc;
    matrix[assoc_uses[3,13] == 1 ? nrow_y_Xq_auc : 0, yK[13]] y13_xq_auc;
    matrix[assoc_uses[3,14] == 1 ? nrow_y_Xq_auc : 0, yK[14]] y14_xq_auc;
    matrix[assoc_uses[3,15] == 1 ? nrow_y_Xq_auc : 0, yK[15]] y15_xq_auc;
    matrix[assoc_uses[3,16] == 1 ? nrow_y_Xq_auc : 0, yK[16]] y16_xq_auc;
    matrix[assoc_uses[3,17] == 1 ? nrow_y_Xq_auc : 0, yK[17]] y17_xq_auc;
    matrix[assoc_uses[3,18] == 1 ? nrow_y_Xq_auc : 0, yK[18]] y18_xq_auc;
    matrix[assoc_uses[3,19] == 1 ? nrow_y_Xq_auc : 0, yK[19]] y19_xq_auc;
    matrix[assoc_uses[3,20] == 1 ? nrow_y_Xq_auc : 0, yK[20]] y20_xq_auc;
    
    // offset values at quadpoints
    vector[has_offset[1] && assoc_uses[3,1] == 1 ? nrow_y_Xq_auc : 0] y1_offset_auc;
    vector[has_offset[2] && assoc_uses[3,2] == 1 ? nrow_y_Xq_auc : 0] y2_offset_auc;
    vector[has_offset[3] && assoc_uses[3,3] == 1 ? nrow_y_Xq_auc : 0] y3_offset_auc;
    vector[has_offset[4] && assoc_uses[3,4] == 1 ? nrow_y_Xq_auc : 0] y4_offset_auc;
    vector[has_offset[5] && assoc_uses[3,5] == 1 ? nrow_y_Xq_auc : 0] y5_offset_auc;
    vector[has_offset[6] && assoc_uses[3,6] == 1 ? nrow_y_Xq_auc : 0] y6_offset_auc;
    vector[has_offset[7] && assoc_uses[3,7] == 1 ? nrow_y_Xq_auc : 0] y7_offset_auc;
    vector[has_offset[8] && assoc_uses[3,8] == 1 ? nrow_y_Xq_auc : 0] y8_offset_auc;
    vector[has_offset[9] && assoc_uses[3,9] == 1 ? nrow_y_Xq_auc : 0] y9_offset_auc;
    vector[has_offset[10] && assoc_uses[3,10] == 1 ? nrow_y_Xq_auc : 0] y10_offset_auc;
    vector[has_offset[11] && assoc_uses[3,11] == 1 ? nrow_y_Xq_auc : 0] y11_offset_auc;
    vector[has_offset[12] && assoc_uses[3,12] == 1 ? nrow_y_Xq_auc : 0] y12_offset_auc;
    vector[has_offset[13] && assoc_uses[3,13] == 1 ? nrow_y_Xq_auc : 0] y13_offset_auc;
    vector[has_offset[14] && assoc_uses[3,14] == 1 ? nrow_y_Xq_auc : 0] y14_offset_auc;
    vector[has_offset[15] && assoc_uses[3,15] == 1 ? nrow_y_Xq_auc : 0] y15_offset_auc;
    vector[has_offset[16] && assoc_uses[3,16] == 1 ? nrow_y_Xq_auc : 0] y16_offset_auc;
    vector[has_offset[17] && assoc_uses[3,17] == 1 ? nrow_y_Xq_auc : 0] y17_offset_auc;
    vector[has_offset[18] && assoc_uses[3,18] == 1 ? nrow_y_Xq_auc : 0] y18_offset_auc;
    vector[has_offset[19] && assoc_uses[3,19] == 1 ? nrow_y_Xq_auc : 0] y19_offset_auc;
    vector[has_offset[20] && assoc_uses[3,20] == 1 ? nrow_y_Xq_auc : 0] y20_offset_auc;

    // re design matrix at quadpoints, group factor 1
    array[bK1_len[1]] vector[assoc_uses[3,1] == 1 && bK1_len[1] > 0 ? nrow_y_Xq_auc : 0] y1_z1q_auc;
    array[bK1_len[2]] vector[assoc_uses[3,2] == 1 && bK1_len[2] > 0 ? nrow_y_Xq_auc : 0] y2_z1q_auc;
    array[bK1_len[3]] vector[assoc_uses[3,3] == 1 && bK1_len[3] > 0 ? nrow_y_Xq_auc : 0] y3_z1q_auc;
    array[bK1_len[4]] vector[assoc_uses[3,4] == 1 && bK1_len[4] > 0 ? nrow_y_Xq_auc : 0] y4_z1q_auc;
    array[bK1_len[5]] vector[assoc_uses[3,5] == 1 && bK1_len[5] > 0 ? nrow_y_Xq_auc : 0] y5_z1q_auc;
    array[bK1_len[6]] vector[assoc_uses[3,6] == 1 && bK1_len[6] > 0 ? nrow_y_Xq_auc : 0] y6_z1q_auc;
    array[bK1_len[7]] vector[assoc_uses[3,7] == 1 && bK1_len[7] > 0 ? nrow_y_Xq_auc : 0] y7_z1q_auc;
    array[bK1_len[8]] vector[assoc_uses[3,8] == 1 && bK1_len[8] > 0 ? nrow_y_Xq_auc : 0] y8_z1q_auc;
    array[bK1_len[9]] vector[assoc_uses[3,9] == 1 && bK1_len[9] > 0 ? nrow_y_Xq_auc : 0] y9_z1q_auc;
    array[bK1_len[10]] vector[assoc_uses[3,10] == 1 && bK1_len[10] > 0 ? nrow_y_Xq_auc : 0] y10_z1q_auc;
    array[bK1_len[11]] vector[assoc_uses[3,11] == 1 && bK1_len[11] > 0 ? nrow_y_Xq_auc : 0] y11_z1q_auc;
    array[bK1_len[12]] vector[assoc_uses[3,12] == 1 && bK1_len[12] > 0 ? nrow_y_Xq_auc : 0] y12_z1q_auc;
    array[bK1_len[13]] vector[assoc_uses[3,13] == 1 && bK1_len[13] > 0 ? nrow_y_Xq_auc : 0] y13_z1q_auc;
    array[bK1_len[14]] vector[assoc_uses[3,14] == 1 && bK1_len[14] > 0 ? nrow_y_Xq_auc : 0] y14_z1q_auc;
    array[bK1_len[15]] vector[assoc_uses[3,15] == 1 && bK1_len[15] > 0 ? nrow_y_Xq_auc : 0] y15_z1q_auc;
    array[bK1_len[16]] vector[assoc_uses[3,16] == 1 && bK1_len[16] > 0 ? nrow_y_Xq_auc : 0] y16_z1q_auc;
    array[bK1_len[17]] vector[assoc_uses[3,17] == 1 && bK1_len[17] > 0 ? nrow_y_Xq_auc : 0] y17_z1q_auc;
    array[bK1_len[18]] vector[assoc_uses[3,18] == 1 && bK1_len[18] > 0 ? nrow_y_Xq_auc : 0] y18_z1q_auc;
    array[bK1_len[19]] vector[assoc_uses[3,19] == 1 && bK1_len[19] > 0 ? nrow_y_Xq_auc : 0] y19_z1q_auc;
    array[bK1_len[20]] vector[assoc_uses[3,20] == 1 && bK1_len[20] > 0 ? nrow_y_Xq_auc : 0] y20_z1q_auc;
    array[assoc_uses[3,1] == 1 && bK1_len[1] > 0 ? nrow_y_Xq_auc : 0] int<lower=0> y1_z1q_id_auc;
    array[assoc_uses[3,2] == 1 && bK1_len[2] > 0 ? nrow_y_Xq_auc : 0] int<lower=0> y2_z1q_id_auc;
    array[assoc_uses[3,3] == 1 && bK1_len[3] > 0 ? nrow_y_Xq_auc : 0] int<lower=0> y3_z1q_id_auc;
    array[assoc_uses[3,4] == 1 && bK1_len[4] > 0 ? nrow_y_Xq_auc : 0] int<lower=0> y4_z1q_id_auc;
    array[assoc_uses[3,5] == 1 && bK1_len[5] > 0 ? nrow_y_Xq_auc : 0] int<lower=0> y5_z1q_id_auc;
    array[assoc_uses[3,6] == 1 && bK1_len[6] > 0 ? nrow_y_Xq_auc : 0] int<lower=0> y6_z1q_id_auc;
    array[assoc_uses[3,7] == 1 && bK1_len[7] > 0 ? nrow_y_Xq_auc : 0] int<lower=0> y7_z1q_id_auc;
    array[assoc_uses[3,8] == 1 && bK1_len[8] > 0 ? nrow_y_Xq_auc : 0] int<lower=0> y8_z1q_id_auc;
    array[assoc_uses[3,9] == 1 && bK1_len[9] > 0 ? nrow_y_Xq_auc : 0] int<lower=0> y9_z1q_id_auc;
    array[assoc_uses[3,10] == 1 && bK1_len[10] > 0 ? nrow_y_Xq_auc : 0] int<lower=0> y10_z1q_id_auc;
    array[assoc_uses[3,11] == 1 && bK1_len[11] > 0 ? nrow_y_Xq_auc : 0] int<lower=0> y11_z1q_id_auc;
    array[assoc_uses[3,12] == 1 && bK1_len[12] > 0 ? nrow_y_Xq_auc : 0] int<lower=0> y12_z1q_id_auc;
    array[assoc_uses[3,13] == 1 && bK1_len[13] > 0 ? nrow_y_Xq_auc : 0] int<lower=0> y13_z1q_id_auc;
    array[assoc_uses[3,14] == 1 && bK1_len[14] > 0 ? nrow_y_Xq_auc : 0] int<lower=0> y14_z1q_id_auc;
    array[assoc_uses[3,15] == 1 && bK1_len[15] > 0 ? nrow_y_Xq_auc : 0] int<lower=0> y15_z1q_id_auc;
    array[assoc_uses[3,16] == 1 && bK1_len[16] > 0 ? nrow_y_Xq_auc : 0] int<lower=0> y16_z1q_id_auc;
    array[assoc_uses[3,17] == 1 && bK1_len[17] > 0 ? nrow_y_Xq_auc : 0] int<lower=0> y17_z1q_id_auc;
    array[assoc_uses[3,18] == 1 && bK1_len[18] > 0 ? nrow_y_Xq_auc : 0] int<lower=0> y18_z1q_id_auc;
    array[assoc_uses[3,19] == 1 && bK1_len[19] > 0 ? nrow_y_Xq_auc : 0] int<lower=0> y19_z1q_id_auc;
    array[assoc_uses[3,20] == 1 && bK1_len[20] > 0 ? nrow_y_Xq_auc : 0] int<lower=0> y20_z1q_id_auc;
    
    // re design matrix at quadpoints, group factor 2
    array[bK2_len[1]] vector[assoc_uses[3,1] == 1 && bK2_len[1] > 0 ? nrow_y_Xq_auc : 0] y1_z2q_auc;
    array[bK2_len[2]] vector[assoc_uses[3,2] == 1 && bK2_len[2] > 0 ? nrow_y_Xq_auc : 0] y2_z2q_auc;
    array[bK2_len[3]] vector[assoc_uses[3,3] == 1 && bK2_len[3] > 0 ? nrow_y_Xq_auc : 0] y3_z2q_auc;
    array[bK2_len[4]] vector[assoc_uses[3,4] == 1 && bK2_len[4] > 0 ? nrow_y_Xq_auc : 0] y4_z2q_auc;
    array[bK2_len[5]] vector[assoc_uses[3,5] == 1 && bK2_len[5] > 0 ? nrow_y_Xq_auc : 0] y5_z2q_auc;
    array[bK2_len[6]] vector[assoc_uses[3,6] == 1 && bK2_len[6] > 0 ? nrow_y_Xq_auc : 0] y6_z2q_auc;
    array[bK2_len[7]] vector[assoc_uses[3,7] == 1 && bK2_len[7] > 0 ? nrow_y_Xq_auc : 0] y7_z2q_auc;
    array[bK2_len[8]] vector[assoc_uses[3,8] == 1 && bK2_len[8] > 0 ? nrow_y_Xq_auc : 0] y8_z2q_auc;
    array[bK2_len[9]] vector[assoc_uses[3,9] == 1 && bK2_len[9] > 0 ? nrow_y_Xq_auc : 0] y9_z2q_auc;
    array[bK2_len[10]] vector[assoc_uses[3,10] == 1 && bK2_len[10] > 0 ? nrow_y_Xq_auc : 0] y10_z2q_auc;
    array[bK2_len[11]] vector[assoc_uses[3,11] == 1 && bK2_len[11] > 0 ? nrow_y_Xq_auc : 0] y11_z2q_auc;
    array[bK2_len[12]] vector[assoc_uses[3,12] == 1 && bK2_len[12] > 0 ? nrow_y_Xq_auc : 0] y12_z2q_auc;
    array[bK2_len[13]] vector[assoc_uses[3,13] == 1 && bK2_len[13] > 0 ? nrow_y_Xq_auc : 0] y13_z2q_auc;
    array[bK2_len[14]] vector[assoc_uses[3,14] == 1 && bK2_len[14] > 0 ? nrow_y_Xq_auc : 0] y14_z2q_auc;
    array[bK2_len[15]] vector[assoc_uses[3,15] == 1 && bK2_len[15] > 0 ? nrow_y_Xq_auc : 0] y15_z2q_auc;
    array[bK2_len[16]] vector[assoc_uses[3,16] == 1 && bK2_len[16] > 0 ? nrow_y_Xq_auc : 0] y16_z2q_auc;
    array[bK2_len[17]] vector[assoc_uses[3,17] == 1 && bK2_len[17] > 0 ? nrow_y_Xq_auc : 0] y17_z2q_auc;
    array[bK2_len[18]] vector[assoc_uses[3,18] == 1 && bK2_len[18] > 0 ? nrow_y_Xq_auc : 0] y18_z2q_auc;
    array[bK2_len[19]] vector[assoc_uses[3,19] == 1 && bK2_len[19] > 0 ? nrow_y_Xq_auc : 0] y19_z2q_auc;
    array[bK2_len[20]] vector[assoc_uses[3,20] == 1 && bK2_len[20] > 0 ? nrow_y_Xq_auc : 0] y20_z2q_auc;
    array[assoc_uses[3,1] == 1 && bK2_len[1] > 0 ? nrow_y_Xq_auc : 0] int<lower=0> y1_z2q_id_auc;
    array[assoc_uses[3,2] == 1 && bK2_len[2] > 0 ? nrow_y_Xq_auc : 0] int<lower=0> y2_z2q_id_auc;
    array[assoc_uses[3,3] == 1 && bK2_len[3] > 0 ? nrow_y_Xq_auc : 0] int<lower=0> y3_z2q_id_auc;
    array[assoc_uses[3,4] == 1 && bK2_len[4] > 0 ? nrow_y_Xq_auc : 0] int<lower=0> y4_z2q_id_auc;
    array[assoc_uses[3,5] == 1 && bK2_len[5] > 0 ? nrow_y_Xq_auc : 0] int<lower=0> y5_z2q_id_auc;
    array[assoc_uses[3,6] == 1 && bK2_len[6] > 0 ? nrow_y_Xq_auc : 0] int<lower=0> y6_z2q_id_auc;
    array[assoc_uses[3,7] == 1 && bK2_len[7] > 0 ? nrow_y_Xq_auc : 0] int<lower=0> y7_z2q_id_auc;
    array[assoc_uses[3,8] == 1 && bK2_len[8] > 0 ? nrow_y_Xq_auc : 0] int<lower=0> y8_z2q_id_auc;
    array[assoc_uses[3,9] == 1 && bK2_len[9] > 0 ? nrow_y_Xq_auc : 0] int<lower=0> y9_z2q_id_auc;
    array[assoc_uses[3,10] == 1 && bK2_len[10] > 0 ? nrow_y_Xq_auc : 0] int<lower=0> y10_z2q_id_auc;
    array[assoc_uses[3,11] == 1 && bK2_len[11] > 0 ? nrow_y_Xq_auc : 0] int<lower=0> y11_z2q_id_auc;
    array[assoc_uses[3,12] == 1 && bK2_len[12] > 0 ? nrow_y_Xq_auc : 0] int<lower=0> y12_z2q_id_auc;
    array[assoc_uses[3,13] == 1 && bK2_len[13] > 0 ? nrow_y_Xq_auc : 0] int<lower=0> y13_z2q_id_auc;
    array[assoc_uses[3,14] == 1 && bK2_len[14] > 0 ? nrow_y_Xq_auc : 0] int<lower=0> y14_z2q_id_auc;
    array[assoc_uses[3,15] == 1 && bK2_len[15] > 0 ? nrow_y_Xq_auc : 0] int<lower=0> y15_z2q_id_auc;
    array[assoc_uses[3,16] == 1 && bK2_len[16] > 0 ? nrow_y_Xq_auc : 0] int<lower=0> y16_z2q_id_auc;
    array[assoc_uses[3,17] == 1 && bK2_len[17] > 0 ? nrow_y_Xq_auc : 0] int<lower=0> y17_z2q_id_auc;
    array[assoc_uses[3,18] == 1 && bK2_len[18] > 0 ? nrow_y_Xq_auc : 0] int<lower=0> y18_z2q_id_auc;
    array[assoc_uses[3,19] == 1 && bK2_len[19] > 0 ? nrow_y_Xq_auc : 0] int<lower=0> y19_z2q_id_auc;
    array[assoc_uses[3,20] == 1 && bK2_len[20] > 0 ? nrow_y_Xq_auc : 0] int<lower=0> y20_z2q_id_auc;
  //---- data for calculating assoc*data interactions in GK quadrature

    // num assoc pars used in {ev/es/mv/ms}*data interactions
    array[M*4] int<lower=0,upper=a_K> a_K_data;

    // design matrix for interacting with ev/es/mv/ms at quadpoints
    matrix[sum(nrow_y_Xq[1:M]), sum(a_K_data)] y_Xq_data;

    // indexing specifying the rows of y_Xq_data that correspond to
    // each submodel
    array[20,2] int<lower=0> idx_q;

  //---- data for combining lower level units clustered within patients

    array[M] int<lower=0,upper=1> has_grp; // 1 = has clustering below patient level
    int<lower=0,upper=4> grp_assoc; // 1=sum, 2=mean, 3=min, 4=max
    array[nrow_e_Xq,2] int<lower=0> grp_idx;
