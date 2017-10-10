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
  * Get the indices corresponding to the lower tri of a square matrix
  *
  * @param dim The number of rows in the square matrix
  * @return A vector of indices
  */ 
  int[] lower_tri_indices(int dim) {
    int indices[dim + choose(dim, 2)];
    int mark = 1;
    for (r in 1:dim) {
      for (c in r:dim) {
        indices[mark] = (r - 1) * dim + c;
        mark = mark + 1;
      }
    }
    return indices;
  }
  
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
    else lb = negative_infinity();
    return lb;
  }
  real ub(int intercept_type) {
    real ub;
    if (intercept_type == 3) ub = 0;
    else ub = positive_infinity();
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
                   real[] aux, int family) {
    vector[rows(z_beta)] beta;                 
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
      beta = prior_mean + ool[1] * prior_scale .* sqrt(2 * mix[1]) .* z_beta;  
    return beta;                 
  }

  /** 
  * Evaluate the linear predictor for the glmer submodel
  *
  * @param X Design matrix for fe
  * @param Z1 Design matrix for re, for first grouping factor
  * @param Z2 Design matrix for re, for second grouping factor
  * @param Z1_id Group indexing for Z1
  * @param Z2_id Group indexing for Z2
  * @param gamma The intercept parameter
  * @param beta Vector of population level parameters
  * @param b1Mat Matrix of group level params for first grouping factor
  * @param b2Mat Matrix of group level params for second grouping factor
  * @param b1Mat_colshift,b2Mat_colshift Number of columns in b1Mat/b2Mat
  *   that correpond to group level params from prior glmer submodels
  * @param intercept_type The type of intercept parameter (0 = none,
  *   1 = unbounded, 2 = lower bound, 3 = upper bound)
  * @return A vector containing the linear predictor for the glmer submodel
  */  
  vector evaluate_eta(matrix X, vector[] Z1, vector[] Z2, int[] Z1_id, int[] Z2_id, 
                      real[] gamma, vector beta, matrix b1Mat, matrix b2Mat, 
                      int b1Mat_colshift, int b2Mat_colshift, 
                      int intercept_type) {
    int N = rows(X);    // num rows in design matrix
    int K = rows(beta); // num predictors
    int p1 = size(Z1);  // num group level params for group factor 1
    int p2 = size(Z2);  // num group level params for group factor 2
    vector[N] eta;
    
    if (K > 0) eta = X * beta;
    else eta = rep_vector(0.0, N);
    
    if (intercept_type > 0) { // submodel has an intercept
      if (intercept_type == 1) eta = eta + gamma[1];
      else if (intercept_type == 2) eta = eta + gamma[1] - max(eta);
      else if (intercept_type == 3) eta = eta + gamma[1] - min(eta);
    }
    
    if (p1 > 0) { // submodel includes group factor 1
      for (k in 1:p1)
        for (n in 1:N)
          eta[n] = eta[n] + (b1Mat[Z1_id[n], k+b1Mat_colshift]) * Z1[k,n];
    }
    if (p2 > 0) { // submodel includes group factor 2
      for (k in 1:p2)
        for (n in 1:N)
          eta[n] = eta[n] + (b2Mat[Z2_id[n], k+b2Mat_colshift]) * Z2[k,n];
    }
    
    return eta;
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
  void glm_lp(vector y_real, int[] y_integer, vector eta, real[] aux, 
              int family, int link, real sum_log_y, vector sqrt_y, vector log_y) {
    if (family == 1) {  // gaussian
      if (link == 1) target += normal_lpdf(y_real | eta, aux[1]);
      else if (link == 2) target += lognormal_lpdf(y_real | eta, aux[1]);
      else target += normal_lpdf(y_real | divide_real_by_vector(1, eta), aux[1]);
    }
    else if (family == 2) {  // gamma
      target += GammaReg(y_real, eta, aux[1], link, sum_log_y);
    }
    else if (family == 3) {  // inverse gaussian 
      target += inv_gaussian(y_real, linkinv_inv_gaussian(eta, link), 
                             aux[1], sum_log_y, sqrt_y);
    }
    else if (family == 4) {  // bernoulli
      if (link == 1) target += bernoulli_logit_lpmf(y_integer | eta);
      else target += bernoulli_lpmf(y_integer | linkinv_bern(eta, link));
    }
    else if (family == 5) {  // binomial
      reject("Binomial with >1 trials not allowed.");
    }
    else if (family == 6 || family == 8) {  // poisson or poisson-gamma
      if (link == 1) target += poisson_log_lpmf(y_integer | eta);
      else target += poisson_lpmf(y_integer | linkinv_count(eta, link));
    }
    else if (family == 7) {  // negative binomial
	    if (link == 1) target += neg_binomial_2_log_lpmf(y_integer | eta, aux[1]);
      else target += neg_binomial_2_lpmf(y_integer | linkinv_count(eta, link), aux[1]);
    }
    else reject("Invalid family.");
  }  

  /** 
  * Create group-specific coefficients, see section 2 of
  * https://cran.r-project.org/web/packages/lme4/vignettes/lmer.pdf
  *
  * @param z_b Vector whose elements are iid normal(0,sigma) a priori
  * @param theta Vector with covariance parameters as defined in lme4
  * @param p An integer array with the number variables on the LHS of each |
  * @param l An integer array with the number of levels for the factor(s) on 
  *   the RHS of each |
  * @param i The index of the grouping factor for which you want to return
  *   the group-specific coefficients for
  * @return An array of group-specific coefficients for grouping factor i
  */
  matrix make_b_matrix(vector z_b, vector theta_L, int[] p, int[] l, int i) {
    matrix[p[i],l[i]] b_matrix;
    int nc = p[i];
    int b_mark = 1;
    int theta_L_mark = 1;
    if (i > 1) {
      for (j in 1:(i-1)) {
        theta_L_mark = theta_L_mark + p[j] + choose(p[j], 2);
        b_mark = b_mark + p[j] * l[j];
      }
    }
    if (nc == 1) {
      real theta_L_start = theta_L[theta_L_mark];
      for (s in b_mark:(b_mark + l[i] - 1)) 
        b_matrix[nc,s] = theta_L_start * z_b[s];
    }
    else {
      matrix[nc,nc] T_i = rep_matrix(0, nc, nc);
      for (c in 1:nc) {
        T_i[c,c] = theta_L[theta_L_mark];
        theta_L_mark = theta_L_mark + 1;
        for(r in (c+1):nc) {
          T_i[r,c] = theta_L[theta_L_mark];
          theta_L_mark = theta_L_mark + 1;
        }
      }
      for (j in 1:l[i]) {
        vector[nc] temp = T_i * segment(z_b, b_mark, nc);
        b_matrix[,j] = temp;
        b_mark = b_mark + nc;
      }
    }
    return b_matrix';
  }

  /** 
  * Evaluate the mean of the posterior predictive distribution
  *
  * @param mu Vector containing the mean of the posterior predictive 
  *   distribution for each observation (ie. the linear predictor after
  *   applying the inverse link function).
  * @param real The auxiliary parameter for the glmer submodel. This will be
  *   an empty array if the submodel does not have an auxiliary parameter
  * @param family An integer specifying the family
  * @return A real, the mean of the posterior predictive distribution
  */ 
  real mean_PPD_rng(vector mu, real[] aux, int family) {
    int N = rows(mu);
    real mean_PPD = 0;	
    if (family == 1) { // gaussian
      for (n in 1:N) 
        mean_PPD = mean_PPD + normal_rng(mu[n], aux[1]);
    }
    else if (family == 2) {  // gamma
      for (n in 1:N) 
        mean_PPD = mean_PPD + gamma_rng(aux[1], aux[1] / mu[n]);
    }
    else if (family == 3) {  // inverse gaussian
      for (n in 1:N) 
        mean_PPD = mean_PPD + inv_gaussian_rng(mu[n], aux[1]);
    }
    else if (family == 4) {  // bernoulli
      for (n in 1:N) 
        mean_PPD = mean_PPD + bernoulli_rng(mu[n]);
    } 
    else if (family == 5) {  // binomial
      reject("Binomial with >1 trials not allowed.");
    }
    else if (family == 6 || family == 8) { 
      real poisson_max = pow(2.0, 30.0);
      for (n in 1:N) {  // poisson or poisson-gamma
        if (mu[n] < poisson_max) 
          mean_PPD = mean_PPD + poisson_rng(mu[n]);
        else 
          mean_PPD = mean_PPD + normal_rng(mu[n], sqrt(mu[n]));
      }
    }
    else if (family == 7) {
      real poisson_max = pow(2.0, 30.0);
      for (n in 1:N) {  // negative binomial
        real gamma_temp;
        if (is_inf(aux[1])) 
          gamma_temp = mu[n];
        else 
          gamma_temp = gamma_rng(aux[1], aux[1] / mu[n]);
        if (gamma_temp < poisson_max) 
          mean_PPD = mean_PPD + poisson_rng(gamma_temp);
        else 
          mean_PPD = mean_PPD + normal_rng(gamma_temp, sqrt(gamma_temp));
      }
    }		  
    return mean_PPD / N;
  }
  
}
data {
  
  // population level dimensions and data
  int<lower=1,upper=3> M; // num submodels with data (limit of 3)
  int<lower=0,upper=1> has_aux[3]; // has auxiliary param
  int<lower=0,upper=2> resp_type[3]; // 1=real,2=integer,0=none
  int<lower=0,upper=3> intercept_type[3]; // 1=unbounded,2=lob,3=upb,0=none
  int<lower=0> yNobs[3]; // num observations
  int<lower=0> yNeta[3]; // required length of eta
  int<lower=0> yK[3]; // num predictors
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
  
  // prior family: 1 = decov, 2 = lkj
  int<lower=1,upper=2> prior_dist_for_cov;
  
  // group level dimensions and data, group factor 1
  int<lower=0> bN1; // num groups
  int<lower=0> bK1; // total num params
  int<lower=0> bK1_len[3]; // num params in each submodel
  int<lower=0> bK1_idx[3,2]; // beg/end index for group params
  vector[bK1_len[1] > 0 ? yNeta[1] : 0] y1_Z1[bK1_len[1]]; // re design matrix
  vector[bK1_len[2] > 0 ? yNeta[2] : 0] y2_Z1[bK1_len[2]];
  vector[bK1_len[3] > 0 ? yNeta[3] : 0] y3_Z1[bK1_len[3]];
  int<lower=0> y1_Z1_id[bK1_len[1] > 0 ? yNeta[1] : 0]; // group indexing for y1_Z1
  int<lower=0> y2_Z1_id[bK1_len[2] > 0 ? yNeta[2] : 0]; // group indexing for y2_Z1
  int<lower=0> y3_Z1_id[bK1_len[3] > 0 ? yNeta[3] : 0]; // group indexing for y3_Z1
  
	// group level dimensions and data, group factor 2
	int<lower=0> bN2; // num groups
  int<lower=0> bK2; // total num params
  int<lower=0> bK2_len[3]; // num params in each submodel
  int<lower=0> bK2_idx[3,2]; // beg/end index for group params
  vector[bK2_len[1] > 0 ? yNeta[1] : 0] y1_Z2[bK2_len[1]]; // re design matrix
  vector[bK2_len[2] > 0 ? yNeta[2] : 0] y2_Z2[bK2_len[2]];
  vector[bK2_len[3] > 0 ? yNeta[3] : 0] y3_Z2[bK2_len[3]];
  int<lower=0> y1_Z2_id[bK2_len[1] > 0 ? yNeta[1] : 0]; // group indexing for y1_Z2
  int<lower=0> y2_Z2_id[bK2_len[2] > 0 ? yNeta[2] : 0]; // group indexing for y2_Z2
  int<lower=0> y3_Z2_id[bK2_len[3] > 0 ? yNeta[3] : 0]; // group indexing for y3_Z2

  // family and link
  int<lower=0> family[M];
  int<lower=0> link[M]; // varies by family 
  
  // prior family: 0 = none, 1 = normal, 2 = student_t, 3 = hs, 4 = hs_plus, 
  //   5 = laplace, 6 = lasso, 7 = product_normal
  int<lower=0,upper=7> y_prior_dist[3];
  int<lower=0,upper=2> y_prior_dist_for_intercept[M]; 
  
  // prior family: 0 = none, 1 = normal, 2 = student_t, 3 = exponential
  int<lower=0,upper=3> y_prior_dist_for_aux[M];
  

  // hyperparameters, values are set to 0 if there is no prior
  
    // coefficients
    vector[yK[1]] y_prior_mean1;
    vector[yK[2]] y_prior_mean2;
    vector[yK[3]] y_prior_mean3;
    vector<lower=0>[yK[1]] y_prior_scale1;
    vector<lower=0>[yK[2]] y_prior_scale2;
    vector<lower=0>[yK[3]] y_prior_scale3;
    vector<lower=0>[yK[1]] y_prior_df1;
    vector<lower=0>[yK[2]] y_prior_df2;
    vector<lower=0>[yK[3]] y_prior_df3;
    vector<lower=0>[M] y_global_prior_df;    // for hs priors only 
    vector<lower=0>[M] y_global_prior_scale; // for hs priors only
  
    // intercepts  
    vector[M] y_prior_mean_for_intercept;
    vector<lower=0>[M] y_prior_scale_for_intercept;
    vector<lower=0>[M] y_prior_df_for_intercept;
    
    // auxiliary params
    vector<lower=0>[M] y_prior_mean_for_aux;
    vector<lower=0>[M] y_prior_scale_for_aux;
    vector<lower=0>[M] y_prior_df_for_aux;
    
    // group level effects, decov prior
    int<lower=0> t;    // num. terms (maybe 0) with a | in the glmer formula
    int<lower=1> p[t]; // num. variables on the LHS of each |
    int<lower=1> l[t]; // num. levels for the factor(s) on the RHS of each |
    int<lower=0> q;    // conceptually equals \sum_{i=1}^t p_i \times l_i
    int<lower=0> len_theta_L; // length of the theta_L vector
    int<lower=0> len_concentration;
    int<lower=0> len_regularization;
    vector<lower=0>[t] b_prior_shape; 
    vector<lower=0>[t] b_prior_scale;
    real<lower=0> b_prior_concentration[len_concentration];
    real<lower=0> b_prior_regularization[len_regularization];
    
    // group level effects, lkj prior
    vector<lower=0>[bK1] b1_prior_scale;
    vector<lower=0>[bK2] b2_prior_scale;
    vector<lower=0>[bK1] b1_prior_df;
    vector<lower=0>[bK2] b2_prior_df;
    real<lower=0> b1_prior_regularization;
    real<lower=0> b2_prior_regularization;
    
    // flag indicating whether to draw from the prior
    int<lower=0,upper=1> prior_PD;  // 1 = yes
}
transformed data {
  // dimensions for hs priors
  int<lower=0> yHs1 = get_nvars_for_hs(M > 0 ? y_prior_dist[1] : 0);
  int<lower=0> yHs2 = get_nvars_for_hs(M > 1 ? y_prior_dist[2] : 0);
  int<lower=0> yHs3 = get_nvars_for_hs(M > 2 ? y_prior_dist[3] : 0);
 
	// data for decov prior
  int<lower=0> len_z_T = 0;
  int<lower=0> len_var_group = sum(p) * (t > 0);
  int<lower=0> len_rho = sum(p) - t;
  int<lower=1> pos = 1;
  real<lower=0> delta[len_concentration]; 
  
  // data for lkj prior
  int bCov1_idx[prior_dist_for_cov == 2 ? (bK1 + choose(bK1, 2)) : 0];
  int bCov2_idx[prior_dist_for_cov == 2 ? (bK2 + choose(bK2, 2)) : 0];
  
  // transformations of data
  real sum_log_y1 = M > 0 && (family[1] == 2 || family[1] == 3) ? 
    sum(log(yReal1)) : not_a_number();
  real sum_log_y2 = M > 1 && (family[2] == 2 || family[2] == 3) ? 
    sum(log(yReal2)) : not_a_number();
  real sum_log_y3 = M > 2 && (family[3] == 2 || family[3] == 3) ? 
    sum(log(yReal3)) : not_a_number();
  vector[M > 0 && family[1] == 3 ? yNobs[1] : 0] sqrt_y1;
  vector[M > 1 && family[2] == 3 ? yNobs[2] : 0] sqrt_y2;
  vector[M > 2 && family[3] == 3 ? yNobs[3] : 0] sqrt_y3;
  vector[M > 0 && family[1] == 3 ? yNobs[1] : 0] log_y1;
  vector[M > 1 && family[2] == 3 ? yNobs[2] : 0] log_y2;
  vector[M > 2 && family[3] == 3 ? yNobs[3] : 0] log_y3;
  if (M > 0 && family[1] == 3) { 
		sqrt_y1 = sqrt(yReal1); 
		log_y1 = log(yReal1);
	} 
  if (M > 1 && family[2] == 3) { 
		sqrt_y2 = sqrt(yReal2); 
		log_y2 = log(yReal2); 
	}  
  if (M > 2 && family[3] == 3) { 
		sqrt_y3 = sqrt(yReal3); 
		log_y3 = log(yReal3); 
	} 

  // data for decov prior
  if (prior_dist_for_cov == 1) {
    for (i in 1:t) {
      if (p[i] > 1) {
        for (j in 1:p[i]) {
          delta[pos] = b_prior_concentration[j];
          pos = pos + 1;
        }
      }
      for (j in 3:p[i]) len_z_T = len_z_T + p[i] - 1;
    }	  
  }
  
  // data for lkj prior
  if (prior_dist_for_cov == 2) {
    if (bK1 > 0) 
      bCov1_idx = lower_tri_indices(bK1);
    if (bK2 > 0) 
      bCov2_idx = lower_tri_indices(bK2); 
  }
}
parameters {
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
  real<lower=0> yOol1[y_prior_dist[1] == 6]; // one_over_lambda
  real<lower=0> yOol2[y_prior_dist[2] == 6];
  real<lower=0> yOol3[y_prior_dist[3] == 6];
  vector<lower=0>[yK[1]] yMix1[y_prior_dist[1] == 5 || y_prior_dist[1] == 6];
  vector<lower=0>[yK[2]] yMix2[y_prior_dist[2] == 5 || y_prior_dist[2] == 6];
  vector<lower=0>[yK[3]] yMix3[y_prior_dist[3] == 5 || y_prior_dist[3] == 6];
}
transformed parameters { 
  vector[yK[1]] yBeta1; // population level params
  vector[yK[2]] yBeta2;
  vector[yK[3]] yBeta3;
  real yAux1[has_aux[1]]; // auxiliary params 
  real yAux2[has_aux[2]];  
  real yAux3[has_aux[3]];
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
                       yGlobal1, yLocal1, yOol1, yMix1, yAux1, family[1]);
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
                         yGlobal2, yLocal2, yOol2, yMix2, yAux2, family[2]);
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
                         yGlobal3, yLocal3, yOol3, yMix3, yAux3, family[3]);
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
  
}
model {
  vector[yNeta[1]] yEta1; // linear predictor
  vector[yNeta[2]] yEta2;
  vector[yNeta[3]] yEta3;
  int bMat1_colshift = 0; // column shift in bMat1
  int bMat2_colshift = 0; // column shift in bMat2

  // Linear predictor for submodel 1
  yEta1 = evaluate_eta(yX1, y1_Z1, y1_Z2, y1_Z1_id, y1_Z2_id, yGamma1, yBeta1, 
                       bMat1, bMat2, bMat1_colshift, bMat2_colshift, intercept_type[1]);

  // Linear predictor for submodel 2
  if (M > 1) {
    bMat1_colshift = bMat1_colshift + bK1_len[1];
    bMat2_colshift = bMat2_colshift + bK2_len[1];
    yEta2 = evaluate_eta(yX2, y2_Z1, y2_Z2, y2_Z1_id, y2_Z2_id, yGamma2, yBeta2, 
                         bMat1, bMat2, bMat1_colshift, bMat2_colshift, intercept_type[2]);
  }
  
  // Linear predictor for submodel 3
  if (M > 2) {
    bMat1_colshift = bMat1_colshift + bK1_len[2];
    bMat2_colshift = bMat2_colshift + bK2_len[2];
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
            y_global_prior_df[1], yLocal1, yGlobal1, yMix1, yOol1)
  if (M > 1 && yK[2] > 0)
    beta_lp(z_yBeta2, y_prior_dist[2], y_prior_scale2, y_prior_df2, 
            y_global_prior_df[2], yLocal2, yGlobal2, yMix2, yOol2)
  if (M > 2 && yK[3] > 0)
    beta_lp(z_yBeta3, y_prior_dist[3], y_prior_scale3, y_prior_df3, 
            y_global_prior_df[3], yLocal3, yGlobal3, yMix3, yOol3)
  
  // Log priors, group level terms
  if (prior_dist_for_cov == 1) { // decov
    decov_lp(z_b, z_T, rho, zeta, tau, b_prior_regularization, 
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

}
generated quantities {
  real mean_PPD[M];
  real yAlpha1[intercept_type[1] > 0];
  real yAlpha2[intercept_type[2] > 0];
  real yAlpha3[intercept_type[3] > 0];
  vector[prior_dist_for_cov == 2 && bK1 > 0 ? size(bCov1_idx) : 0] bCov1;
  vector[prior_dist_for_cov == 2 && bK2 > 0 ? size(bCov2_idx) : 0] bCov2; 
  
  // Evaluate mean_PPD
  {
    int bMat1_colshift = 0; // column shift in bMat1
    int bMat2_colshift = 0; // column shift in bMat2

    // Linear predictor for submodel 1
    if (M > 0) {
      vector[yNeta[1]] yEta1; // linear predictor
      yEta1 = evaluate_eta(yX1, y1_Z1, y1_Z2, y1_Z1_id, y1_Z2_id, yGamma1, yBeta1, 
                           bMat1, bMat2, bMat1_colshift, bMat2_colshift, intercept_type[1]);
      yEta1 = evaluate_mu(yEta1, family[1], link[1]);
      mean_PPD[1] = mean_PPD_rng(yEta1, yAux1, family[1]);
    }
    
    // Linear predictor for submodel 2
    if (M > 1) {
      vector[yNeta[2]] yEta2;
      bMat1_colshift = bMat1_colshift + bK1_len[1];
      bMat2_colshift = bMat2_colshift + bK2_len[1];
      yEta2 = evaluate_eta(yX2, y2_Z1, y2_Z2, y2_Z1_id, y2_Z2_id, yGamma2, yBeta2, 
                           bMat1, bMat2, bMat1_colshift, bMat2_colshift, intercept_type[2]);
      yEta2 = evaluate_mu(yEta2, family[2], link[2]);
      mean_PPD[2] = mean_PPD_rng(yEta2, yAux2, family[2]);
    }
    
    // Linear predictor for submodel 3
    if (M > 2) {
      vector[yNeta[3]] yEta3;
      bMat1_colshift = bMat1_colshift + bK1_len[2];
      bMat2_colshift = bMat2_colshift + bK2_len[2];
      yEta3 = evaluate_eta(yX3, y3_Z1, y3_Z2, y3_Z1_id, y3_Z2_id, yGamma3, yBeta3, 
                           bMat1, bMat2, bMat1_colshift, bMat2_colshift, intercept_type[3]);
      yEta3 = evaluate_mu(yEta3, family[3], link[3]);
      mean_PPD[3] = mean_PPD_rng(yEta3, yAux3, family[3]);
    } 
  }  
  
  // Transform intercept parameters
	if (intercept_type[1] > 0) 
    yAlpha1[1] = yGamma1[1] - dot_product(yXbar1, yBeta1);
  if (M > 1 && intercept_type[2] > 0) 
    yAlpha2[1] = yGamma2[1] - dot_product(yXbar2, yBeta2);
  if (M > 2 && intercept_type[3] > 0) 
    yAlpha3[1] = yGamma3[1] - dot_product(yXbar3, yBeta3);
	
	// Transform variance-covariance matrices
	
	  // Grouping factor 1
  	if (prior_dist_for_cov == 2 && bK1 == 1) {
  	  bCov1[1] = bSd1[1] * bSd1[1];
  	} 
  	else if (prior_dist_for_cov == 2 && bK1 > 1) {
  	  bCov1 = to_vector(quad_form_diag(
  	    multiply_lower_tri_self_transpose(bCholesky1), bSd1))[bCov1_idx];
    }
    
    // Grouping factor 2
    if (prior_dist_for_cov == 2 && bK2 == 1) {
  	  bCov2[1] = bSd2[1] * bSd2[1];
  	} 
  	else if (prior_dist_for_cov == 2 && bK2 > 1) {
    	bCov2 = to_vector(quad_form_diag(
    	   multiply_lower_tri_self_transpose(bCholesky2), bSd2))[bCov2_idx]; 
    }
}
