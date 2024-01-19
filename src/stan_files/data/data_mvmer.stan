  // population level data
  array[resp_type[1] == 2 ? yNobs[1] : 0] int<lower=0> yInt1; // integer responses
  array[resp_type[2] == 2 ? yNobs[2] : 0] int<lower=0> yInt2;
  array[resp_type[3] == 2 ? yNobs[3] : 0] int<lower=0> yInt3;
  array[resp_type[4] == 2 ? yNobs[4] : 0] int<lower=0> yInt4;
  array[resp_type[5] == 2 ? yNobs[5] : 0] int<lower=0> yInt5;
  array[resp_type[6] == 2 ? yNobs[6] : 0] int<lower=0> yInt6;
  array[resp_type[7] == 2 ? yNobs[7] : 0] int<lower=0> yInt7;
  array[resp_type[8] == 2 ? yNobs[8] : 0] int<lower=0> yInt8;
  array[resp_type[9] == 2 ? yNobs[9] : 0] int<lower=0> yInt9;
  array[resp_type[10] == 2 ? yNobs[10] : 0] int<lower=0> yInt10;
  array[resp_type[11] == 2 ? yNobs[11] : 0] int<lower=0> yInt11;
  array[resp_type[12] == 2 ? yNobs[12] : 0] int<lower=0> yInt12;
  array[resp_type[13] == 2 ? yNobs[13] : 0] int<lower=0> yInt13;
  array[resp_type[14] == 2 ? yNobs[14] : 0] int<lower=0> yInt14;
  array[resp_type[15] == 2 ? yNobs[15] : 0] int<lower=0> yInt15;
  array[resp_type[16] == 2 ? yNobs[16] : 0] int<lower=0> yInt16;
  array[resp_type[17] == 2 ? yNobs[17] : 0] int<lower=0> yInt17;
  array[resp_type[18] == 2 ? yNobs[18] : 0] int<lower=0> yInt18;
  array[resp_type[19] == 2 ? yNobs[19] : 0] int<lower=0> yInt19;
  array[resp_type[20] == 2 ? yNobs[20] : 0] int<lower=0> yInt20;
  vector[resp_type[1] == 1 ? yNobs[1] : 0] yReal1; // real responses
  vector[resp_type[2] == 1 ? yNobs[2] : 0] yReal2;
  vector[resp_type[3] == 1 ? yNobs[3] : 0] yReal3;
  vector[resp_type[4] == 1 ? yNobs[4] : 0] yReal4;
  vector[resp_type[5] == 1 ? yNobs[5] : 0] yReal5;
  vector[resp_type[6] == 1 ? yNobs[6] : 0] yReal6;
  vector[resp_type[7] == 1 ? yNobs[7] : 0] yReal7;
  vector[resp_type[8] == 1 ? yNobs[8] : 0] yReal8;
  vector[resp_type[9] == 1 ? yNobs[9] : 0] yReal9;
  vector[resp_type[10] == 1 ? yNobs[10] : 0] yReal10;
  vector[resp_type[11] == 1 ? yNobs[11] : 0] yReal11;
  vector[resp_type[12] == 1 ? yNobs[12] : 0] yReal12;
  vector[resp_type[13] == 1 ? yNobs[13] : 0] yReal13;
  vector[resp_type[14] == 1 ? yNobs[14] : 0] yReal14;
  vector[resp_type[15] == 1 ? yNobs[15] : 0] yReal15;
  vector[resp_type[16] == 1 ? yNobs[16] : 0] yReal16;
  vector[resp_type[17] == 1 ? yNobs[17] : 0] yReal17;
  vector[resp_type[18] == 1 ? yNobs[18] : 0] yReal18;
  vector[resp_type[19] == 1 ? yNobs[19] : 0] yReal19;
  vector[resp_type[20] == 1 ? yNobs[20] : 0] yReal20;
  matrix[yNeta[1],yK[1]] yX1; // fe design matrix
  matrix[yNeta[2],yK[2]] yX2;
  matrix[yNeta[3],yK[3]] yX3;
  matrix[yNeta[4],yK[4]] yX4;
  matrix[yNeta[5],yK[5]] yX5;
  matrix[yNeta[6],yK[6]] yX6;
  matrix[yNeta[7],yK[7]] yX7;
  matrix[yNeta[8],yK[8]] yX8;
  matrix[yNeta[9],yK[9]] yX9;
  matrix[yNeta[10],yK[10]] yX10;
  matrix[yNeta[11],yK[11]] yX11;
  matrix[yNeta[12],yK[12]] yX12;
  matrix[yNeta[13],yK[13]] yX13;
  matrix[yNeta[14],yK[14]] yX14;
  matrix[yNeta[15],yK[15]] yX15;
  matrix[yNeta[16],yK[16]] yX16;
  matrix[yNeta[17],yK[17]] yX17;
  matrix[yNeta[18],yK[18]] yX18;
  matrix[yNeta[19],yK[19]] yX19;
  matrix[yNeta[20],yK[20]] yX20;
  vector[yK[1]] yXbar1; // predictor means
  vector[yK[2]] yXbar2;
  vector[yK[3]] yXbar3;
  vector[yK[4]] yXbar4;
  vector[yK[5]] yXbar5;
  vector[yK[6]] yXbar6;
  vector[yK[7]] yXbar7;
  vector[yK[8]] yXbar8;
  vector[yK[9]] yXbar9;
  vector[yK[10]] yXbar10;
  vector[yK[11]] yXbar11;
  vector[yK[12]] yXbar12;
  vector[yK[13]] yXbar13;
  vector[yK[14]] yXbar14;
  vector[yK[15]] yXbar15;
  vector[yK[16]] yXbar16;
  vector[yK[17]] yXbar17;
  vector[yK[18]] yXbar18;
  vector[yK[19]] yXbar19;
  vector[yK[20]] yXbar20;

  // family and link (determined by 'append_mvmer_famlink' R function)
  // 1 = gaussian
  // 2 = gamma
  // 3 = inverse gaussian
  // 4 = bernoulli
  // 5 = binomial (n>1)
  // 6 = poisson
  // 7 = negative binomial
  array[M] int<lower=0> family;
  array[M] int<lower=0> link; // varies by family

  // group level data, group factor 1
  array[bK1_len[1]] vector[bK1_len[1] > 0 ? yNeta[1] : 0] y1_Z1; // re design matrix
  array[bK1_len[2]] vector[bK1_len[2] > 0 ? yNeta[2] : 0] y2_Z1;
  array[bK1_len[3]] vector[bK1_len[3] > 0 ? yNeta[3] : 0] y3_Z1;
  array[bK1_len[4]] vector[bK1_len[4] > 0 ? yNeta[4] : 0] y4_Z1;
  array[bK1_len[5]] vector[bK1_len[5] > 0 ? yNeta[5] : 0] y5_Z1;
  array[bK1_len[6]] vector[bK1_len[6] > 0 ? yNeta[6] : 0] y6_Z1;
  array[bK1_len[7]] vector[bK1_len[7] > 0 ? yNeta[7] : 0] y7_Z1;
  array[bK1_len[8]] vector[bK1_len[8] > 0 ? yNeta[8] : 0] y8_Z1;
  array[bK1_len[9]] vector[bK1_len[9] > 0 ? yNeta[9] : 0] y9_Z1;
  array[bK1_len[10]] vector[bK1_len[10] > 0 ? yNeta[10] : 0] y10_Z1;
  array[bK1_len[11]] vector[bK1_len[11] > 0 ? yNeta[11] : 0] y11_Z1;
  array[bK1_len[12]] vector[bK1_len[12] > 0 ? yNeta[12] : 0] y12_Z1;
  array[bK1_len[13]] vector[bK1_len[13] > 0 ? yNeta[13] : 0] y13_Z1;
  array[bK1_len[14]] vector[bK1_len[14] > 0 ? yNeta[14] : 0] y14_Z1;
  array[bK1_len[15]] vector[bK1_len[15] > 0 ? yNeta[15] : 0] y15_Z1;
  array[bK1_len[16]] vector[bK1_len[16] > 0 ? yNeta[16] : 0] y16_Z1;
  array[bK1_len[17]] vector[bK1_len[17] > 0 ? yNeta[17] : 0] y17_Z1;
  array[bK1_len[18]] vector[bK1_len[18] > 0 ? yNeta[18] : 0] y18_Z1;
  array[bK1_len[19]] vector[bK1_len[19] > 0 ? yNeta[19] : 0] y19_Z1;
  array[bK1_len[20]] vector[bK1_len[20] > 0 ? yNeta[20] : 0] y20_Z1;
  array[bK1_len[1] > 0 ? yNeta[1] : 0] int<lower=0> y1_Z1_id; // group indexing for y1_Z1
  array[bK1_len[2] > 0 ? yNeta[2] : 0] int<lower=0> y2_Z1_id; // group indexing for y2_Z1
  array[bK1_len[3] > 0 ? yNeta[3] : 0] int<lower=0> y3_Z1_id; // group indexing for y3_Z1
  array[bK1_len[4] > 0 ? yNeta[4] : 0] int<lower=0> y4_Z1_id; // group indexing for y4_Z1
  array[bK1_len[5] > 0 ? yNeta[5] : 0] int<lower=0> y5_Z1_id; // group indexing for y5_Z1
  array[bK1_len[6] > 0 ? yNeta[6] : 0] int<lower=0> y6_Z1_id; // group indexing for y6_Z1
  array[bK1_len[7] > 0 ? yNeta[7] : 0] int<lower=0> y7_Z1_id; // group indexing for y7_Z1
  array[bK1_len[8] > 0 ? yNeta[8] : 0] int<lower=0> y8_Z1_id; // group indexing for y8_Z1
  array[bK1_len[9] > 0 ? yNeta[9] : 0] int<lower=0> y9_Z1_id; // group indexing for y9_Z1
  array[bK1_len[10] > 0 ? yNeta[10] : 0] int<lower=0> y10_Z1_id; // group indexing for y10_Z1
  array[bK1_len[11] > 0 ? yNeta[11] : 0] int<lower=0> y11_Z1_id; // group indexing for y11_Z1
  array[bK1_len[12] > 0 ? yNeta[12] : 0] int<lower=0> y12_Z1_id; // group indexing for y12_Z1
  array[bK1_len[13] > 0 ? yNeta[13] : 0] int<lower=0> y13_Z1_id; // group indexing for y13_Z1
  array[bK1_len[14] > 0 ? yNeta[14] : 0] int<lower=0> y14_Z1_id; // group indexing for y14_Z1
  array[bK1_len[15] > 0 ? yNeta[15] : 0] int<lower=0> y15_Z1_id; // group indexing for y15_Z1
  array[bK1_len[16] > 0 ? yNeta[16] : 0] int<lower=0> y16_Z1_id; // group indexing for y16_Z1
  array[bK1_len[17] > 0 ? yNeta[17] : 0] int<lower=0> y17_Z1_id; // group indexing for y17_Z1
  array[bK1_len[18] > 0 ? yNeta[18] : 0] int<lower=0> y18_Z1_id; // group indexing for y18_Z1
  array[bK1_len[19] > 0 ? yNeta[19] : 0] int<lower=0> y19_Z1_id; // group indexing for y19_Z1
  array[bK1_len[20] > 0 ? yNeta[20] : 0] int<lower=0> y20_Z1_id; // group indexing for y20_Z1

  // group level data, group factor 2
  array[bK2_len[1]] vector[bK2_len[1] > 0 ? yNeta[1] : 0] y1_Z2; // re design matrix
  array[bK2_len[2]] vector[bK2_len[2] > 0 ? yNeta[2] : 0] y2_Z2;
  array[bK2_len[3]] vector[bK2_len[3] > 0 ? yNeta[3] : 0] y3_Z2;
  array[bK2_len[4]] vector[bK2_len[4] > 0 ? yNeta[4] : 0] y4_Z2;
  array[bK2_len[5]] vector[bK2_len[5] > 0 ? yNeta[5] : 0] y5_Z2;
  array[bK2_len[6]] vector[bK2_len[6] > 0 ? yNeta[6] : 0] y6_Z2;
  array[bK2_len[7]] vector[bK2_len[7] > 0 ? yNeta[7] : 0] y7_Z2;
  array[bK2_len[8]] vector[bK2_len[8] > 0 ? yNeta[8] : 0] y8_Z2;
  array[bK2_len[9]] vector[bK2_len[9] > 0 ? yNeta[9] : 0] y9_Z2;
  array[bK2_len[10]] vector[bK2_len[10] > 0 ? yNeta[10] : 0] y10_Z2;
  array[bK2_len[11]] vector[bK2_len[11] > 0 ? yNeta[11] : 0] y11_Z2;
  array[bK2_len[12]] vector[bK2_len[12] > 0 ? yNeta[12] : 0] y12_Z2;
  array[bK2_len[13]] vector[bK2_len[13] > 0 ? yNeta[13] : 0] y13_Z2;
  array[bK2_len[14]] vector[bK2_len[14] > 0 ? yNeta[14] : 0] y14_Z2;
  array[bK2_len[15]] vector[bK2_len[15] > 0 ? yNeta[15] : 0] y15_Z2;
  array[bK2_len[16]] vector[bK2_len[16] > 0 ? yNeta[16] : 0] y16_Z2;
  array[bK2_len[17]] vector[bK2_len[17] > 0 ? yNeta[17] : 0] y17_Z2;
  array[bK2_len[18]] vector[bK2_len[18] > 0 ? yNeta[18] : 0] y18_Z2;
  array[bK2_len[19]] vector[bK2_len[19] > 0 ? yNeta[19] : 0] y19_Z2;
  array[bK2_len[20]] vector[bK2_len[20] > 0 ? yNeta[20] : 0] y20_Z2;
  array[bK2_len[1] > 0 ? yNeta[1] : 0] int<lower=0> y1_Z2_id; // group indexing for y1_Z2
  array[bK2_len[2] > 0 ? yNeta[2] : 0] int<lower=0> y2_Z2_id; // group indexing for y2_Z2
  array[bK2_len[3] > 0 ? yNeta[3] : 0] int<lower=0> y3_Z2_id; // group indexing for y3_Z2
  array[bK2_len[4] > 0 ? yNeta[4] : 0] int<lower=0> y4_Z2_id; // group indexing for y4_Z2
  array[bK2_len[5] > 0 ? yNeta[5] : 0] int<lower=0> y5_Z2_id; // group indexing for y5_Z2
  array[bK2_len[6] > 0 ? yNeta[6] : 0] int<lower=0> y6_Z2_id; // group indexing for y6_Z2
  array[bK2_len[7] > 0 ? yNeta[7] : 0] int<lower=0> y7_Z2_id; // group indexing for y7_Z2
  array[bK2_len[8] > 0 ? yNeta[8] : 0] int<lower=0> y8_Z2_id; // group indexing for y8_Z2
  array[bK2_len[9] > 0 ? yNeta[9] : 0] int<lower=0> y9_Z2_id; // group indexing for y9_Z2
  array[bK2_len[10] > 0 ? yNeta[10] : 0] int<lower=0> y10_Z2_id; // group indexing for y10_Z2
  array[bK2_len[11] > 0 ? yNeta[11] : 0] int<lower=0> y11_Z2_id; // group indexing for y11_Z2
  array[bK2_len[12] > 0 ? yNeta[12] : 0] int<lower=0> y12_Z2_id; // group indexing for y12_Z2
  array[bK2_len[13] > 0 ? yNeta[13] : 0] int<lower=0> y13_Z2_id; // group indexing for y13_Z2
  array[bK2_len[14] > 0 ? yNeta[14] : 0] int<lower=0> y14_Z2_id; // group indexing for y14_Z2
  array[bK2_len[15] > 0 ? yNeta[15] : 0] int<lower=0> y15_Z2_id; // group indexing for y15_Z2
  array[bK2_len[16] > 0 ? yNeta[16] : 0] int<lower=0> y16_Z2_id; // group indexing for y16_Z2
  array[bK2_len[17] > 0 ? yNeta[17] : 0] int<lower=0> y17_Z2_id; // group indexing for y17_Z2
  array[bK2_len[18] > 0 ? yNeta[18] : 0] int<lower=0> y18_Z2_id; // group indexing for y18_Z2
  array[bK2_len[19] > 0 ? yNeta[19] : 0] int<lower=0> y19_Z2_id; // group indexing for y19_Z2
  array[bK2_len[20] > 0 ? yNeta[20] : 0] int<lower=0> y20_Z2_id; // group indexing for y20_Z2

  // prior family: 0 = none, 1 = normal, 2 = student_t, 3 = hs, 4 = hs_plus,
  //   5 = laplace, 6 = lasso, 7 = product_normal
  array[3] int<lower=0,upper=7> y_prior_dist;
  array[M] int<lower=0,upper=2> y_prior_dist_for_intercept;

  // prior family: 0 = none, 1 = normal, 2 = student_t, 3 = exponential
  array[M] int<lower=0,upper=3> y_prior_dist_for_aux;

  // prior family: 1 = decov, 2 = lkj
  int<lower=1,upper=2> prior_dist_for_cov;

  // flag indicating whether to draw from the prior
  int<lower=0,upper=1> prior_PD;  // 1 = yes
  
  // offset
  array[3] int<lower=0,upper=1> has_offset;  // 0 = No, 1 = Yes
  vector[has_offset[1] ? yNeta[1] : 0] y1_offset;
  vector[has_offset[2] ? yNeta[2] : 0] y2_offset;
  vector[has_offset[3] ? yNeta[3] : 0] y3_offset;
  vector[has_offset[4] ? yNeta[4] : 0] y4_offset;
  vector[has_offset[5] ? yNeta[5] : 0] y5_offset;
  vector[has_offset[6] ? yNeta[6] : 0] y6_offset;
  vector[has_offset[7] ? yNeta[7] : 0] y7_offset;
  vector[has_offset[8] ? yNeta[8] : 0] y8_offset;
  vector[has_offset[9] ? yNeta[9] : 0] y9_offset;
  vector[has_offset[10] ? yNeta[10] : 0] y10_offset;
  vector[has_offset[11] ? yNeta[11] : 0] y11_offset;
  vector[has_offset[12] ? yNeta[12] : 0] y12_offset;
  vector[has_offset[13] ? yNeta[13] : 0] y13_offset;
  vector[has_offset[14] ? yNeta[14] : 0] y14_offset;
  vector[has_offset[15] ? yNeta[15] : 0] y15_offset;
  vector[has_offset[16] ? yNeta[16] : 0] y16_offset;
  vector[has_offset[17] ? yNeta[17] : 0] y17_offset;
  vector[has_offset[18] ? yNeta[18] : 0] y18_offset;
  vector[has_offset[19] ? yNeta[19] : 0] y19_offset;
  vector[has_offset[20] ? yNeta[20] : 0] y20_offset;
