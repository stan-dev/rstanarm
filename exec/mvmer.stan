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
 
  // family and link
  int<lower=0> family[M];
  int<lower=0> link[M]; // varies by family 
  
  // prior family: 0 = none, 1 = normal, 2 = student_t, 3 = hs, 4 = hs_plus, 
  //   5 = laplace, 6 = lasso, 7 = product_normal
  int<lower=0,upper=7> prior_dist[M];
  int<lower=0,upper=2> prior_dist_for_intercept[M]; 
  
  // prior family: 0 = none, 1 = normal, 2 = student_t, 3 = exponential
  int<lower=0,upper=3> prior_dist_for_aux[M];

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
  int sum_has_aux;
  sum_has_aux = sum(has_aux);
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
    vector[bK1 == 1 ? bN1 : 0] 
      z_bVec1; 
    matrix[bK1 >  1 ? bK1 : 0, bK1 >  1 ? bN1 : 0] 
      z_bMat1; 
    // cholesky factor of corr matrix (if > 1 random effect)
    cholesky_factor_corr[bK1 > 1 ? bK1 : 0] bCholesky1;  
  
  // group level params for second grouping factor
    // group-level sds
    vector<lower=0>[bK2] bSd2; 
    // unscaled group-level params 
    vector[bK2 == 1 ? bN2 : 0] 
      z_bVec2; 
    matrix[bK2 >  1 ? bK2 : 0, bK2 >  1 ? bN2 : 0] 
      z_bMat2;
    // cholesky factor of corr matrix (if > 1 random effect)
    cholesky_factor_corr[bK2 > 1 ? bK2 : 0] bCholesky2;  
  
  // auxiliary params, interpretation depends on family
  vector<lower=0>[sum_has_aux] aux_unscaled; 
  
  // params for priors
  real<lower=0> global[len_global];
  vector<lower=0>[len_local2] local2[(len_local2 > 0) ? 2 : 0];
  vector<lower=0>[len_local4] local4[(len_local4 > 0) ? 4 : 0];
  vector<lower=0>[len_mix] mix[(len_mix > 0)];
  real<lower=0> ool[len_ool]; // one_over_lambda
  vector<lower=0>[len_noise] noise[(len_noise > 0)];
  
}
transformed parameters { 
  // group-level params for first grouping factor
  vector[bK1 == 1 ? bN1 : 0] 
    bVec1 = bSd1[1] * (z_bVec1); 
  matrix[bK1 >  1 ? bN1 : 0, bK1 >  1 ? bK1 : 0] 
    bMat1 = (diag_pre_multiply(bSd1, bCholesky1) * z_bMat1); 
 
  // group level params for second grouping factor
  vector[bK2 == 1 ? bN2 : 0] 
    bVec2 = bSd2[1] * (z_bVec2); 
  matrix[bK2 >  1 ? bN2 : 0, bK2 >  1 ? bK2 : 0] 
    bMat2 = (diag_pre_multiply(bSd2, bCholesky2) * z_bMat2); 

  // population level params, auxiliary params
  vector[yK[1]] yBeta1;
  vector[yK[2]] yBeta2;
  vector[yK[3]] yBeta3;
  vector[sum_has_aux] aux;  

  // auxiliary parameters
  if (sum_has_aux > 0) {
    int aux_mark = 1;
    for (m in 1:M) {
      if (has_aux[m] == 1) {
        if (prior_dist_for_aux[m] == 0) // none
          aux[aux_mark] = aux_unscaled[aux_mark];
        else {
          aux[aux_mark] = prior_scale_for_aux[m] * aux_unscaled[aux_mark];
          if (prior_dist_for_aux[m] <= 2) // normal or student_t
            aux[aux_mark] = aux[aux_mark] + prior_mean_for_aux[m];
        }
        aux_mark = aux_mark + 1;
      }
    }
  }
  
  // population level params
  yBeta1 = make_beta(y_prior_dist[1], y_prior_mean1, 
                     y_prior_scale1, y_prior_df1, y_global1, y_local1,
                     y_global_prior_scale[1], y_global_prior_df[1], y_ool1, y_mix1)
  if (M >= 2) {
    yBeta1 = make_beta(y_prior_dist[2], y_prior_mean2, 
                       y_prior_scale2, y_prior_df2, y_global2, y_local2,
                       y_global_prior_scale[2], y_global_prior_df[2], y_ool2, y_mix2)
  }
  if (M >= 3) {
    yBeta1 = make_beta(y_prior_dist[3], y_prior_mean3, 
                       y_prior_scale3, y_prior_df3, y_global3, y_local3,
                       y_global_prior_scale[3], y_global_prior_df[3], y_ool3, y_mix3)
  }  
  
  if (prior_special_case == 1) { 
    // same prior for all submodels (none, normal or student t)
    if      (prior_dist[1] == 0) beta = z_beta;
    else if (prior_dist[1] == 1) beta = z_beta .* prior_scale + prior_mean;
    else if (prior_dist[1] == 2) for (k in 1:K) {
      beta[k] = CFt(z_beta[k], prior_df[k]) * prior_scale[k] + prior_mean[k];
    }    
  }
  else {  
    // different prior for each submodel
  for (m in 1:M) {
    int K1 = idx_K[m,1];  // indexing for beta vector
    int K2 = idx_K[m,2];  // indexing for beta vector
    if      (prior_dist[m] == 0) beta[K1:K2] = z_beta[K1:K2];
    else if (prior_dist[m] == 1) beta[K1:K2] = z_beta[K1:K2] .* prior_scale[K1:K2] + prior_mean[K1:K2];
    else if (prior_dist[m] == 2) for (k in K1:K2) {
      beta[k] = CFt(z_beta[k], prior_df[k]) * prior_scale[k] + prior_mean[k];
    }
    else if (prior_dist[m] == 3) {
	  int G1 = idx_global[m,1];  // indexing for global params
	  int G2 = idx_global[m,2];  // indexing for global params
	  int L1 = idx_local2[m,1];  // indexing for local params
	  int L2 = idx_local2[m,2];  // indexing for local params	  
      if (family[m] == 1) { // don't need is_continuous since family == 1 is gaussian in mvmer
		    int aux_mark = sum(has_aux[1:m]);
        beta[K1:K2] = hs_prior(z_beta[K1:K2], global[G1:G2], local2[,L1:L2], global_prior_scale[m], aux[aux_mark]);
	  }
      else beta[K1:K2] = hs_prior(z_beta[K1:K2], global[G1:G2], local2[,L1:L2], global_prior_scale[m], 1);
    }
    else if (prior_dist[m] == 4) {
      int G1 = idx_global[m,1];  // indexing for global params
	  int G2 = idx_global[m,2];  // indexing for global params
	  int L1 = idx_local4[m,1];  // indexing for local params
	  int L2 = idx_local4[m,2];  // indexing for local params
  	  if (family[m] == 1) { // don't need is_continuous since family == 1 is gaussian in mvmer
	  	int aux_mark = sum(has_aux[1:m]);
        beta[K1:K2] = hsplus_prior(z_beta[K1:K2], global[G1:G2], local4[,L1:L2], global_prior_scale[m], aux[aux_mark]);
	  }
      else beta[K1:K2] = hsplus_prior(z_beta[K1:K2], global[G1:G2], local4[,L1:L2], global_prior_scale[m], 1);
    }
    else if (prior_dist[m] == 5) // laplace
      beta[K1:K2] = prior_mean[K1:K2] + prior_scale[K1:K2] .* sqrt(2 * mix[1][idx_mix[m,1]:idx_mix[m,2]]) .* z_beta[K1:K2];
    else if (prior_dist[m] == 6) // lasso
      beta[K1:K2] = prior_mean[K1:K2] + ool[idx_ool[m]] * prior_scale[K1:K2] .* sqrt(2 * mix[1][idx_mix[m,1]:idx_mix[m,2]]) .* z_beta[K1:K2];  
  }  

}
model {
  int aux_mark = 1; 
  vector[yN[1]] yEta1; // linear predictor
  vector[yN[2]] yEta2;
  vector[yN[3]] yEta3;
  if (K > 0) {
    if (dense_X) eta = X[1] * beta;
    else eta = csr_matrix_times_vector(N, K, w_X, v_X, u_X, beta);
  }

  if (t > 0) {
    #include "eta_add_Zb.stan"
  } 
  // Log-likelihoods
  for (m in 1:M) {
    vector[NM[m]] eta_tmp;	            // eta for just one submodel 
    eta_tmp = eta[idx[m,1]:idx[m,2]];   // eta for just one submodel 
    #include "eta_intercept_mvmer.stan"	// adds intercept or shifts eta
    if (family[m] == 8) {  // poisson-gamma mixture
	  #include "eta_add_noise_mvmer.stan"
    }    
    #include "mvmer_lp.stan" // increments target with long log-liks
    if (has_aux[m] == 1) aux_mark = aux_mark + 1;
  }
  // Log-priors	
  for (m in 1:M) { 
    if (has_aux[m] == 1) {
      aux_mark = sum(has_aux[1:m]);
      aux_lp(aux_unscaled[aux_mark], prior_dist_for_aux[m], prior_scale_for_aux[m], prior_df_for_aux[m])
    } 
  }
  if (t > 0) decov_lp(z_b, z_T, rho, zeta, tau, 
                      regularization, delta, shape, t, p);
}
generated quantities {
  #include "gen_quantities_mvmer.stan"  // defines alpha, mean_PPD
}
