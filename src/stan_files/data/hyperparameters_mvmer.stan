  // hyperparameter values are set to 0 if there is no prior

  // coefficients
  vector[yK[1]] y_prior_mean1;
  vector[yK[2]] y_prior_mean2;
  vector[yK[3]] y_prior_mean3;
  vector[yK[4]] y_prior_mean4;
  vector[yK[5]] y_prior_mean5;
  vector[yK[6]] y_prior_mean6;
  vector[yK[7]] y_prior_mean7;
  vector[yK[8]] y_prior_mean8;
  vector[yK[9]] y_prior_mean9;
  vector[yK[10]] y_prior_mean10;
  vector[yK[11]] y_prior_mean11;
  vector[yK[12]] y_prior_mean12;
  vector[yK[13]] y_prior_mean13;
  vector[yK[14]] y_prior_mean14;
  vector[yK[15]] y_prior_mean15;
  vector[yK[16]] y_prior_mean16;
  vector[yK[17]] y_prior_mean17;
  vector[yK[18]] y_prior_mean18;
  vector[yK[19]] y_prior_mean19;
  vector[yK[20]] y_prior_mean20;
  vector<lower=0>[yK[1]] y_prior_scale1;
  vector<lower=0>[yK[2]] y_prior_scale2;
  vector<lower=0>[yK[3]] y_prior_scale3;
  vector<lower=0>[yK[4]] y_prior_scale4;
  vector<lower=0>[yK[5]] y_prior_scale5;
  vector<lower=0>[yK[6]] y_prior_scale6;
  vector<lower=0>[yK[7]] y_prior_scale7;
  vector<lower=0>[yK[8]] y_prior_scale8;
  vector<lower=0>[yK[9]] y_prior_scale9;
  vector<lower=0>[yK[10]] y_prior_scale10;
  vector<lower=0>[yK[11]] y_prior_scale11;
  vector<lower=0>[yK[12]] y_prior_scale12;
  vector<lower=0>[yK[13]] y_prior_scale13;
  vector<lower=0>[yK[14]] y_prior_scale14;
  vector<lower=0>[yK[15]] y_prior_scale15;
  vector<lower=0>[yK[16]] y_prior_scale16;
  vector<lower=0>[yK[17]] y_prior_scale17;
  vector<lower=0>[yK[18]] y_prior_scale18;
  vector<lower=0>[yK[19]] y_prior_scale19;
  vector<lower=0>[yK[20]] y_prior_scale20;
  vector<lower=0>[yK[1]] y_prior_df1;
  vector<lower=0>[yK[2]] y_prior_df2;
  vector<lower=0>[yK[3]] y_prior_df3;
  vector<lower=0>[yK[4]] y_prior_df4;
  vector<lower=0>[yK[5]] y_prior_df5;
  vector<lower=0>[yK[6]] y_prior_df6;
  vector<lower=0>[yK[7]] y_prior_df7;
  vector<lower=0>[yK[8]] y_prior_df8;
  vector<lower=0>[yK[9]] y_prior_df9;
  vector<lower=0>[yK[10]] y_prior_df10;
  vector<lower=0>[yK[11]] y_prior_df11;
  vector<lower=0>[yK[12]] y_prior_df12;
  vector<lower=0>[yK[13]] y_prior_df13;
  vector<lower=0>[yK[14]] y_prior_df14;
  vector<lower=0>[yK[15]] y_prior_df15;
  vector<lower=0>[yK[16]] y_prior_df16;
  vector<lower=0>[yK[17]] y_prior_df17;
  vector<lower=0>[yK[18]] y_prior_df18;
  vector<lower=0>[yK[19]] y_prior_df19;
  vector<lower=0>[yK[20]] y_prior_df20;
  vector<lower=0>[M] y_global_prior_df;    // for hs priors only
  vector<lower=0>[M] y_global_prior_scale; // for hs priors only
  vector<lower=0>[M] y_slab_df;            // for hs priors only
  vector<lower=0>[M] y_slab_scale;         // for hs priors only

  // intercepts
  vector[M] y_prior_mean_for_intercept;
  vector<lower=0>[M] y_prior_scale_for_intercept;
  vector<lower=0>[M] y_prior_df_for_intercept;

  // auxiliary params
  vector<lower=0>[M] y_prior_mean_for_aux;
  vector<lower=0>[M] y_prior_scale_for_aux;
  vector<lower=0>[M] y_prior_df_for_aux;

  // decov prior stuff
  int<lower=0> len_concentration;
  int<lower=0> len_regularization;
  vector<lower=0>[t] b_prior_shape;
  vector<lower=0>[t] b_prior_scale;
  array[len_concentration]  real<lower=0> b_prior_concentration;
  array[len_regularization] real<lower=0> b_prior_regularization;

  // lkj prior stuff
  vector<lower=0>[bK1] b1_prior_scale;
  vector<lower=0>[bK2] b2_prior_scale;
  vector<lower=0>[bK1] b1_prior_df;
  vector<lower=0>[bK2] b2_prior_df;
  real<lower=0> b1_prior_regularization;
  real<lower=0> b2_prior_regularization;
