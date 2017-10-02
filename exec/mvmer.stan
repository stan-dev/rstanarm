#include "Columbia_copyright.stan"
#include "Brilleman_copyright.stan"
#include "license.stan" // GPL3+

// Multivariate GLM with group-specific terms
functions {
  #include "common_functions.stan"
  #include "bernoulli_likelihoods.stan"
  #include "binomial_likelihoods.stan"
  #include "continuous_likelihoods.stan"
  #include "count_likelihoods.stan"
  #include "jm_functions.stan"
  
  /** 
  * Return the lower/upper bound for the specified intercept type
  *
  * @param intercept_type An integer specifying the type of intercept; 
  *   0=no intercept, 1=unbounded, 2=lower bounded, 3=upper bounded
  * @return A real, corresponding to the lower bound
  */
  real lb(int intercept_type) {
    real lb;
    if (intercept_type == 2) lb = 0;
    else lb = -Inf;
    return lb;
  }
  real ub(int intercept_type) {
    real ub;
    if (intercept_type == 3) ub = 0;
    else ub = Inf;
    return ub;
  } 

  /** 
  * Scale the auxiliary parameter based on prior information
  *
  * @param aux_unscaled A real, the unscaled auxiliary parameter
  * @param prior_dist Integer, the type of prior distribution
  * @param prior_mean,prior_scale Real scalars, the mean and scale
  *   of the prior distribution
  * @return A real, corresponding to the scaled auxiliary parameter
  */  
  real make_aux(real aux_unscaled, int prior_dist, 
                real prior_mean, real prior_scale) {
    real aux;
    if (prior_dist == 0) // none
      aux = aux_unscaled;
    else {
      aux = prior_scale * aux_unscaled;
      if (prior_dist <= 2) // normal or student_t
        aux = aux + prior_mean;
    }
    return aux;
  }

  /** 
  * Scale the primitive population level parameters based on prior information
  *
  * @param z_beta A vector of primitive parameters
  * @param prior_dist Integer, the type of prior distribution
  * @param prior_mean,prior_scale Vectors of mean and scale parameters
  *   for the prior distributions
  * @return A vector containing the population level parameters (coefficients)
  */ 
  vector make_beta(vector z_beta, int prior_dist, vector prior_mean, 
                   vector prior_scale, vector prior_df, real global_prior_scale,
                   real[] global, vector[] local, real[] ool, vector[] mix, 
                   vector aux, int family) {
    if (prior_dist == 0) beta = z_beta;
    else if (prior_dist == 1) beta = z_beta .* prior_scale + prior_mean;
    else if (prior_dist == 2) for (k in 1:rows(prior_mean)) {
      beta[k] = CFt(z_beta[k], prior_df[k]) * prior_scale[k] + prior_mean[k];
    }
    else if (prior_dist == 3) {
      if (family == 1) // don't need is_continuous since family == 1 is gaussian in mvmer
        beta = hs_prior(z_beta, global, local, global_prior_scale, aux[1]);
      else 
        beta = hs_prior(z_beta, global, local, global_prior_scale, 1.0);
    }
    else if (prior_dist == 4) {
  	  if (family == 1) // don't need is_continuous since family == 1 is gaussian in mvmer
        beta = hsplus_prior(z_beta, global, local, global_prior_scale, aux[1]);
      else 
        beta = hsplus_prior(z_beta, global, local, global_prior_scale, 1.0);
    }
    else if (prior_dist == 5) // laplace
      beta = prior_mean + prior_scale .* sqrt(2 * mix[1]) .* z_beta;
    else if (prior_dist == 6) // lasso
      beta = prior_mean + ool * prior_scale .* sqrt(2 * mix[1]) .* z_beta;  
                     
  }

  /** 
  * Increment the target with the log-likelihood for the glmer submodel
  *
  * @param z_beta A vector of primitive parameters
  * @param prior_dist Integer, the type of prior distribution
  * @param prior_mean,prior_scale Vectors of mean and scale parameters
  *   for the prior distributions
  * @return A vector containing the population level parameters (coefficients)
  */   
  void glm_lp(vector y_real, int[] y_integer, vector eta, real aux, 
              int family, int link, real sum_log_y, vector sqrt_y) {
    if (family == 1) {  // gaussian
      if (link == 1) target += normal_lpdf(y_real | eta, aux);
      else if (link == 2) target += lognormal_lpdf(y_real | eta, aux);
      else target += normal_lpdf(y_real | divide_real_by_vector(1, eta), aux);
    }
    else if (family == 2) {  // gamma
      target += GammaReg(y_real, eta, aux, link, sum_log_y);
    }
    else if (family == 3) {  // inverse gaussian 
      target += inv_gaussian(y_real, linkinv_inv_gaussian(eta, link), 
                             aux, sum_log_y, sqrt_y);
    }
    else if (family == 4) {  // bernoulli
      if (link == 1) target += bernoulli_logit_lpmf(y_integer | eta);
      else target += bernoulli_lpmf(y_integer | linkinv_bern(eta, link));
    }
    else if (family == 5) {  // binomial
      // binomial with num trials > 1 has been removed	  
    }
    else if (family == 6 || family == 8) {  // poisson or poisson-gamma
      if (link == 1) target += poisson_log_lpmf(y_integer | eta);
      else target += poisson_lpmf(y_integer | linkinv_count(eta, link));
    }
    else if (family == 7) {  // negative binomial
	    if (link == 1) target += neg_binomial_2_log_lpmf(y_integer | eta, aux);
      else target += neg_binomial_2_lpmf(y_integer | linkinv_count(eta, link), aux);
    }    
  }  


 
}
data {
  
  // dimensions, data
  int<lower=1,upper=3> M; // num submodels with data (limit of 3)
  int<lower=0,upper=2> resp_type[3]; // 1=real,2=integer,0=none
  int<lower=0,upper=3> intercept_type[3]; // 1=unbounded,2=lob,3=upb,0=none
  int<lower=0,upper=1> has_aux[3]; // has auxiliary param
  int<lower=0> yNobs[3]; // num observations
  int<lower=0> yNeta[3]; // required length of eta
  int<lower=0> yK[3]; // num predictors
  int<lower=0> yInt1[resp_type[1] == 2 ? yNobs[1] : 0]; // integer responses
  int<lower=0> yInt2[resp_type[2] == 2 ? yNobs[2] : 0];
  int<lower=0> yInt2[resp_type[3] == 2 ? yNobs[3] : 0];
  vector[resp_type[1] == 1 ? yNobs[1] : 0] yReal1; // real responses
  vector[resp_type[2] == 1 ? yNobs[2] : 0] yReal2;  
  vector[resp_type[3] == 1 ? yNobs[3] : 0] yReal3;  
  matrix[yN[1],yK[1]] yX1; // fe design matrix
  matrix[yN[2],yK[2]] yX2; 
  matrix[yN[3],yK[3]] yX3; 
  vector[yK[1]] yXbar1; // predictor means
  vector[yK[2]] yXbar2;
  vector[yK[3]] yXbar3;
  int<lower=0> bK1; // total num params for group factor 1
  int<lower=0> bN1; // num groups for group factor 1
  int<lower=0> bK2; // total num params for group factor 2
  int<lower=0> bN2; // num groups for group factor 2
  int<lower=0> bK1_len[3]; // num params in each submodel for group factor 1
  int<lower=0> bK1_beg[3]; // num params in each submodel for group factor 1
  int<lower=0> bK1_end[3]; // num params in each submodel for group factor 1
  int<lower=0> yP2[3]; // num params in each submodel for group factor 2
 
  // family and link
  int<lower=0> family[M];
  int<lower=0> link[M]; // varies by family 
  
  // prior family: 0 = none, 1 = normal, 2 = student_t, 3 = hs, 4 = hs_plus, 
  //   5 = laplace, 6 = lasso, 7 = product_normal
  int<lower=0,upper=7> y_prior_dist[M];
  int<lower=0,upper=2> y_prior_dist_for_intercept[M]; 
  
  // prior family: 0 = none, 1 = normal, 2 = student_t, 3 = exponential
  int<lower=0,upper=3> y_prior_dist_for_aux[M];

  // hyperparameters, values are set to 0 if there is no prior
  vector<lower=0>[yK[1]] y_prior_scale1;
  vector<lower=0>[yK[2]] y_prior_scale2;
  vector<lower=0>[yK[3]] y_prior_scale3;
  vector<lower=0>[M] y_prior_scale_for_intercept;
  vector<lower=0>[M] y_prior_scale_for_aux;
  vector[yK[1]] y_prior_mean1;
  vector[yK[2]] y_prior_mean2;
  vector[yK[3]] y_prior_mean3;
  vector[M] y_prior_mean_for_intercept;
  vector<lower=0>[M] y_prior_mean_for_aux;
  vector<lower=0>[yK[1]] y_prior_df1;
  vector<lower=0>[yK[2]] y_prior_df2;
  vector<lower=0>[yK[3]] y_prior_df3;
  vector<lower=0>[M] y_prior_df_for_intercept;
  vector<lower=0>[M] y_prior_df_for_aux;
  vector<lower=0>[M] y_global_prior_df;    // for hs priors only 
  vector<lower=0>[M] y_global_prior_scale; // for hs priors only
}
transformed data {   
  int sum_has_aux = sum(has_aux);
  int<lower=0> yHs1 = get_nvars_for_hs(y_prior_dist[1]);
  int<lower=0> yHs2 = get_nvars_for_hs(y_prior_dist[2]);
  int<lower=0> yHs3 = get_nvars_for_hs(y_prior_dist[3]);
}
parameters {
  // intercepts
  real<lower=lb(intercept_type[1],upper=ub(intercept_type[1]))>
    yGamma1[intercept_type[1] > 0];
  real<lower=lb(intercept_type[2],upper=ub(intercept_type[2]))> 
    yGamma2[intercept_type[2] > 0]; 
  real<lower=lb(intercept_type[3],upper=ub(intercept_type[3]))> 
    yGamma3[intercept_type[3] > 0];
  
  // population level primitive params  
  vector[yK[1]] z_yBeta1; 
  vector[yK[2]] z_yBeta2;
  vector[yK[3]] z_yBeta3;

  // group level params for first grouping factor
    // group-level sds
    vector<lower=0>[bK1] bSd1; 
    // unscaled group-level params 
    vector[bK1 == 1 ? bN1 : 0] z_bVec1; 
    matrix[bK1 >  1 ? bK1 : 0, bK1 >  1 ? bN1 : 0] z_bMat1; 
    // cholesky factor of corr matrix (if > 1 random effect)
    cholesky_factor_corr[bK1 > 1 ? bK1 : 0] bCholesky1;  
  
  // group level params for second grouping factor
    // group-level sds
    vector<lower=0>[bK2] bSd2; 
    // unscaled group-level params 
    vector[bK2 == 1 ? bN2 : 0] z_bVec2; 
    matrix[bK2 >  1 ? bK2 : 0, bK2 >  1 ? bN2 : 0] z_bMat2;
    // cholesky factor of corr matrix (if > 1 random effect)
    cholesky_factor_corr[bK2 > 1 ? bK2 : 0] bCholesky2;  
  
  // auxiliary params, interpretation depends on family
  vector<lower=0>[has_aux[1]] yAux1_unscaled; 
  vector<lower=0>[has_aux[2]] yAux2_unscaled; 
  vector<lower=0>[has_aux[3]] yAux3_unscaled; 
  
  // params for priors
  real<lower=0> yGlobal1[yHs1];
  real<lower=0> yGlobal2[yHs2];
  real<lower=0> yGlobal3[yHs3];
  vector<lower=0>[yK[1]] yLocal1[yHs1];
  vector<lower=0>[yK[2]] yLocal2[yHs2];
  vector<lower=0>[yK[3]] yLocal3[yHs3];
  real<lower=0> yOol1[y_prior_dist[1] == 6]; // one_over_lambda
  real<lower=0> yOol2[y_prior_dist[2] == 6];
  real<lower=0> yOol3[y_prior_dist[3] == 6];
  vector<lower=0>[yK[1]] yMix1[y_prior_dist[1] == 5 || y_prior_dist[1] == 6];
  vector<lower=0>[yK[2]] yMix2[y_prior_dist[2] == 5 || y_prior_dist[2] == 6];
  vector<lower=0>[yK[3]] yMix3[y_prior_dist[3] == 5 || y_prior_dist[3] == 6];
}
transformed parameters { 
  // group-level params for first grouping factor
  vector[bK1 == 1 ? bN1 : 0] bVec1 = bSd1[1] * (z_bVec1); 
  matrix[bK1 >  1 ? bN1 : 0, bK1 >  1 ? bK1 : 0] 
    bMat1 = (diag_pre_multiply(bSd1, bCholesky1) * z_bMat1); 
 
  // group level params for second grouping factor
  vector[bK2 == 1 ? bN2 : 0] bVec2 = bSd2[1] * (z_bVec2); 
  matrix[bK2 >  1 ? bN2 : 0, bK2 >  1 ? bK2 : 0] 
    bMat2 = (diag_pre_multiply(bSd2, bCholesky2) * z_bMat2); 

  // population level params, auxiliary params
  vector[yK[1]] yBeta1;
  vector[yK[2]] yBeta2;
  vector[yK[3]] yBeta3;
  vector[has_aux[1]] yAux1;  
  vector[has_aux[2]] yAux2;  
  vector[has_aux[3]] yAux3;  
  
  if (has_aux[1] == 1)
    yAux1 = make_aux(yAux1_unscaled, y_prior_dist_for_aux[1], 
                     y_prior_mean_for_aux[1], y_prior_scale_for_aux[1]);
  if (yK[1] > 0)
    yBeta1 = make_beta(z_yBeta1, y_prior_dist[1], y_prior_mean1, 
                       y_prior_scale1, y_prior_df1, y_global_prior_scale[1],  
                       yGlobal1, yLocal1, yOol1, yMix1, yAux1, family[1]);
  if (M > 1) {
    if (has_aux[2] == 1)
      yAux2 = make_aux(yAux2_unscaled, y_prior_dist_for_aux[2], 
                       y_prior_mean_for_aux[2], y_prior_scale_for_aux[2]);
    if (yK[2] > 0)
      yBeta2 = make_beta(z_yBeta2, y_prior_dist[2], y_prior_mean2, 
                         y_prior_scale2, y_prior_df2, y_global_prior_scale[2],  
                         yGlobal2, yLocal2, yOol2, yMix2, yAux2, family[2]);
  } 
  if (M > 2) {
    if (has_aux[3] == 1)
      yAux3 = make_aux(yAux3_unscaled, y_prior_dist_for_aux[3], 
                       y_prior_mean_for_aux[3], y_prior_scale_for_aux[3]);
    if (yK[3] > 0)
      yBeta3 = make_beta(z_yBeta3, y_prior_dist[3], y_prior_mean3, 
                         y_prior_scale3, y_prior_df3, y_global_prior_scale[3],  
                         yGlobal3, yLocal3, yOol3, yMix3, yAux3, family[3]);
  }  
}
model {
  vector[yN[1]] yEta1; // linear predictor
  vector[yN[2]] yEta2;
  vector[yN[3]] yEta3;
  
  // Linear predictor for submodel 1
  if (yK[1] > 0) yEta1 = yX1 * yBeta1;
  else yEta1 = rep(0.0, yN[1]);
  if (intercept_type[1] == 1) yEta1 = yEta1 + yGamma1;
  else if (intercept_type[1] == 2) yEta1 = yEta1 + yGamma1 - max(yEta1);
  else if (intercept_type[1] == 3) yEta1 = yEta1 + yGamma1 - min(yEta1);
  if (bK1_len[1] > 0) { // includes group factor 1
    if (bK1 == 1) { // one group level param
      for (n in 1:yN[1])
        yEta1[n] = yEta1[n] + (bVec1[y1_b1_idx[n]]) * y1_Z1[1,n];
    }
    else if (bK1 > 1) {// more than one group level param
      for (k in bK1_beg[1]:bK1_end[1])
        for (n in 1:yN[1])
          yEta1[n] = yEta1[n] + (bMat1[y1_b1_idx[n],k]) * y1_Z1[k,n];
    }
  }
  if (bK2_len[1] > 0) { // includes group factor 2
    if (bK2 == 1) { // one group level param
      for (n in 1:yN[1])
        yEta1[n] = yEta1[n] + (bVec2[y1_b2_idx[n]]) * y1_Z2[1,n];
    }
    else if (bK2 > 1) { // more than one group level param
      for (k in bK1_beg[1]:bK1_end[1])
        for (n in 1:yN[1])
          yEta1[n] = yEta1[n] + (bMat2[y1_b2_idx[n],k]) * y1_Z2[k,n];
    }
  }
  
  // Linear predictor for submodel 2
  if (M > 1) {
    if (yK[2] > 0) yEta2 = yX2 * yBeta2;
    else yEta2 = rep(0.0, yN[2]);
    if (intercept_type[2] == 1) yEta2 = yEta2 + yGamma2;
    else if (intercept_type[2] == 2) yEta2 = yEta2 + yGamma2 - max(yEta2);
    else if (intercept_type[2] == 3) yEta2 = yEta2 + yGamma2 - min(yEta2);
    if (b1_len[2] > 0) { // includes group factor 1
      if (bK1 == 1) { // one group level param
        for (n in 1:yN[2])
          yEta2[n] = yEta2[n] + (bVec1[y2_b1_idx[n]]) * y2_Z1[1,n];
      }
      else if (bK1 > 1) { // more than one group level param
        int k_shift = b1_len[1];
        for (k in b1_beg[2]:b1_end[2])
          for (n in 1:yN[2])
            yEta2[n] = yEta2[n] + (bMat1[y2_b1_idx[n],k]) * y2_Z1[k-k_shift,n];
      }
    }
    if (b2_len[2] > 0) { // includes group factor 2
      if (bK2 == 1) // one group level param
        for (n in 1:yN[2])
          yEta2[n] = yEta2[n] + (bVec2[y2_b2_idx[n]]) * y2_Z2[1,n]
      else if (bK2 > 1) // more than one group level param
        int k_shift = b2_len[1];
        for (k in b2_beg[2]:b2_end[2])
          for (n in 1:yN[1])
            yEta1[n] = yEta1[n] + (bMat2[y2_b2_idx[n],k]) * y2_Z2[k-k_shift,n]
    }
  }
  
  // Linear predictor for submodel 3
  if (M > 2) {
    if (yK[3] > 0) yEta3 = yX3 * yBeta3;
    else yEta3 = rep(0.0, yN[3]);
    if (intercept_type[3] == 1) yEta3 = yEta3 + yGamma3;
    else if (intercept_type[3] == 2) yEta3 = yEta3 + yGamma3 - max(yEta3);
    else if (intercept_type[3] == 3) yEta3 = yEta3 + yGamma3 - min(yEta3);
    if (b1_len[3] > 0) { // includes group factor 1
      if (bK1 == 1) { // one group level param
        for (n in 1:yN[3])
          yEta3[n] = yEta3[n] + (bVec1[y3_b1_idx[n]]) * y3_Z1[1,n];
      }
      else if (bK1 > 1) { // more than one group level param
        int k_shift = sum(b1_len[1:2]);
        for (k in b1_beg[3]:b1_end[3])
          for (n in 1:yN[3])
            yEta3[n] = yEta3[n] + (bMat1[y3_b1_idx[n],k]) * y3_Z1[k-k_shift,n];
      }
    }
    if (b2_len[3] > 0) { // includes group factor 2
      if (bK2 == 1) // one group level param
        for (n in 1:yN[3])
          yEta3[n] = yEta3[n] + (bVec2[y3_b2_idx[n]]) * y3_Z2[1,n]
      else if (bK2 > 1) // more than one group level param
        int k_shift = sum(b2_len[1:2]);
        for (k in b2_beg[3]:b2_end[3])
          for (n in 1:yN[3])
            yEta3[n] = yEta3[n] + (bMat2[y3_b2_idx[n],k]) * y3_Z2[k-k_shift,n]
    }    
  }
  
  // Log-likelihoods
  if (prior_pd == 0) {
    glm_lp(yReal1, yInt1, yEta1, yAux1, family[1], link[1], sum_log_y1, sqrt_y1);
    if (M > 1)
      glm_lp(yReal2, yInt2, yEta2, yAux2, family[2], link[2], sum_log_y2, sqrt_y2);
    if (M > 2)
      glm_lp(yReal3, yInt3, yEta3, yAux3, family[3], link[3], sum_log_y3, sqrt_y3);
  }
  
  // Log-priors, auxiliary params	
  if (has_aux[1] == 1) 
    aux_lp(yAux1_unscaled, y_prior_dist_for_aux[1], 
           y_prior_scale_for_aux[1], y_prior_df_for_aux[1]);
  if (M > 1 && has_aux[2] == 1)
    aux_lp(yAux2_unscaled, y_prior_dist_for_aux[2], 
           y_prior_scale_for_aux[2], y_prior_df_for_aux[2]);
  if (M > 2 && has_aux[3] == 1)
    aux_lp(yAux3_unscaled, y_prior_dist_for_aux[3], 
           y_prior_scale_for_aux[3], y_prior_df_for_aux[3]);

  // Log priors, intercepts
  if (intercept_type[1] > 0) 
    gamma_lp(yGamma1[1], y_prior_dist_for_intercept[1], y_prior_mean_for_intercept[1], 
             y_prior_scale_for_intercept[1], y_prior_df_for_intercept)[1];  
  if (M > 1 && intercept_type[2] > 0) 
    gamma_lp(yGamma1[2], y_prior_dist_for_intercept[2], y_prior_mean_for_intercept[2], 
             y_prior_scale_for_intercept[2], y_prior_df_for_intercept)[2];  
  if (M > 2 && intercept_type[3] > 0) 
    gamma_lp(yGamma1[3], y_prior_dist_for_intercept[3], y_prior_mean_for_intercept[3], 
             y_prior_scale_for_intercept[3], y_prior_df_for_intercept)[3];  

  // Log priors, population level params
  if (yK[1] > 0)
    beta_lp(z_yBeta1, y_prior_dist[1], y_prior_scale1, y_prior_df1, 
            y_global_prior_df[1], yLocal1, yGlobal1, yMix1, yOol1)
  if (M > 1 && yK[2] > 0)
    beta_lp(z_yBeta2, y_prior_dist[2], y_prior_scale2, y_prior_df2, 
            y_global_prior_df[2], yLocal2, yGlobal2, yMix2, yOol2)
  if (M > 2 && yK[3] > 0)
    beta_lp(z_yBeta3, y_prior_dist[3], y_prior_scale3, y_prior_df3, 
            y_global_prior_df[3], yLocal3, yGlobal3, yMix3, yOol3)
  
  // Log priors, group level terms
  if (bK1 > 0) // sds for group factor 1
    target += student_t_lpdf(bSd1 | 3, 0, b1_prior_scale);
  if (bK2 > 0) // sds for group factor 2
    target += student_t_lpdf(bSd2 | 3, 0, b2_prior_scale);
  if (bK1 == 1) // primitive group level params for group factor 1
    target += normal_lpdf(z_bVec1 | 0, 1); 
  if (bK2 == 1) // primitive group level params for group factor 2
    target += normal_lpdf(z_bVec2 | 0, 1); 
  if (bK1 > 1) {
    // primitive group level params for group factor 1
    target += normal_lpdf(to_vector(z_bMat1) | 0, 1); 
    // corr matrix for group factor 1 
    target += lkj_corr_cholesky_lpdf(bCholesky1 | 1);
  }  
  if (bK2 > 1) {
    // primitive group level params for group factor 2
    target += normal_lpdf(to_vector(z_bMat2) | 0, 1); 
    // corr matrix for group factor 2
    target += lkj_corr_cholesky_lpdf(bCholesky2 | 1);
  }
}
generated quantities {
  corr_matrix[bK1 > 1 ? bK1 : 0] 
    bCor1 = multiply_lower_tri_self_transpose(bCholesky1); 
  corr_matrix[bK2 > 1 ? bK2 : 0] 
    bCor2 = multiply_lower_tri_self_transpose(bCholesky2); 

  #include "gen_quantities_mvmer.stan"  // defines alpha, mean_PPD
}
