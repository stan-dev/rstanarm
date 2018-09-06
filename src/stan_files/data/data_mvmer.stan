  // population level data
  int<lower=0> yInt1[resp_type[1] == 2 ? yNobs[1] : 0]; // integer responses
  int<lower=0> yInt2[resp_type[2] == 2 ? yNobs[2] : 0];
  int<lower=0> yInt3[resp_type[3] == 2 ? yNobs[3] : 0];
  vector[resp_type[1] == 1 ? yNobs[1] : 0] yReal1; // real responses
  vector[resp_type[2] == 1 ? yNobs[2] : 0] yReal2;
  vector[resp_type[3] == 1 ? yNobs[3] : 0] yReal3;
  matrix[yNeta[1],yK[1]] yX1; // fe design matrix
  matrix[yNeta[2],yK[2]] yX2;
  matrix[yNeta[3],yK[3]] yX3;
  vector[yK[1]] yXbar1; // predictor means
  vector[yK[2]] yXbar2;
  vector[yK[3]] yXbar3;

  // family and link (determined by 'append_mvmer_famlink' R function)
  // 1 = gaussian
  // 2 = gamma
  // 3 = inverse gaussian
  // 4 = bernoulli
  // 5 = binomial (n>1)
  // 6 = poisson
  // 7 = negative binomial
  int<lower=0> family[M];
  int<lower=0> link[M]; // varies by family

  // group level data, group factor 1
  vector[bK1_len[1] > 0 ? yNeta[1] : 0] y1_Z1[bK1_len[1]]; // re design matrix
  vector[bK1_len[2] > 0 ? yNeta[2] : 0] y2_Z1[bK1_len[2]];
  vector[bK1_len[3] > 0 ? yNeta[3] : 0] y3_Z1[bK1_len[3]];
  int<lower=0> y1_Z1_id[bK1_len[1] > 0 ? yNeta[1] : 0]; // group indexing for y1_Z1
  int<lower=0> y2_Z1_id[bK1_len[2] > 0 ? yNeta[2] : 0]; // group indexing for y2_Z1
  int<lower=0> y3_Z1_id[bK1_len[3] > 0 ? yNeta[3] : 0]; // group indexing for y3_Z1

  // group level data, group factor 2
  vector[bK2_len[1] > 0 ? yNeta[1] : 0] y1_Z2[bK2_len[1]]; // re design matrix
  vector[bK2_len[2] > 0 ? yNeta[2] : 0] y2_Z2[bK2_len[2]];
  vector[bK2_len[3] > 0 ? yNeta[3] : 0] y3_Z2[bK2_len[3]];
  int<lower=0> y1_Z2_id[bK2_len[1] > 0 ? yNeta[1] : 0]; // group indexing for y1_Z2
  int<lower=0> y2_Z2_id[bK2_len[2] > 0 ? yNeta[2] : 0]; // group indexing for y2_Z2
  int<lower=0> y3_Z2_id[bK2_len[3] > 0 ? yNeta[3] : 0]; // group indexing for y3_Z2

  // prior family: 0 = none, 1 = normal, 2 = student_t, 3 = hs, 4 = hs_plus,
  //   5 = laplace, 6 = lasso, 7 = product_normal
  int<lower=0,upper=7> y_prior_dist[3];
  int<lower=0,upper=2> y_prior_dist_for_intercept[M];

  // prior family: 0 = none, 1 = normal, 2 = student_t, 3 = exponential
  int<lower=0,upper=3> y_prior_dist_for_aux[M];

  // prior family: 1 = decov, 2 = lkj
  int<lower=1,upper=2> prior_dist_for_cov;

  // flag indicating whether to draw from the prior
  int<lower=0,upper=1> prior_PD;  // 1 = yes
