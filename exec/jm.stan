#include "Columbia_copyright.stan"
#include "Monash_copyright.stan"
#include "license.stan" // GPL3+

// Shared parameter joint model
functions {
  #include "common_functions.stan"
  #include "bernoulli_likelihoods.stan"
  #include "binomial_likelihoods.stan"
  #include "continuous_likelihoods.stan"
  #include "count_likelihoods.stan"
  #include "jm_functions.stan"
}
data {
  
  // dimensions
  int<lower=1> M; // num. of long. submodels
  int<lower=0> Npat; // num. individuals (equal to l[1] - 1)
  int<lower=0> y_N[M]; // num. of obs. in each long. submodel
  int<lower=0> sum_y_N; // total num. of obs. across all long submodels
  int<lower=0> sum_y_real_N; // total num. of obs. across all long submodels with real outcomes
  int<lower=0> sum_y_int_N; // total num. of obs. across all long submodels with integer outcomes
  int<lower=0,upper=sum_y_N> y_beg[M]; // index of first obs. for each submodel
  int<lower=0,upper=sum_y_N> y_end[M]; // index of last obs. for each submodel
  int<lower=0,upper=sum_y_real_N> y_real_beg[M]; // index of first obs. for each submodel
  int<lower=0,upper=sum_y_real_N> y_real_end[M]; // index of last obs. for each submodel
  int<lower=0,upper=sum_y_int_N> y_int_beg[M]; // index of first obs. for each submodel
  int<lower=0,upper=sum_y_int_N> y_int_end[M]; // index of last obs. for each submodel
  int<lower=0> y_N01[M,2];  // num. of bernoulli 0/1 observations in each long. submodel
  int<lower=0> y_K[M];  // num. of predictors in each long. submodel
  int<lower=0> sum_y_K; // total num. of predictors across all long submodels
  int<lower=0,upper=sum_y_K> y_K_beg[M];  // index of first parameter for each submodel
  int<lower=0,upper=sum_y_K> y_K_end[M];  // index of last parameter for each submodel
  int<lower=0> e_K; // num. of predictors in event submodel
  int<lower=0> a_K; // num. of association parameters
  int<lower=0> quadnodes;  // num. of nodes for Gauss-Kronrod quadrature 
  int<lower=0> Npat_times_quadnodes;
  int<lower=0,upper=M> sum_y_has_intercept; // num. submodels w/ intercept
  int<lower=0,upper=M> sum_y_has_intercept_unbound; // num. submodels w/ unbounded intercept
  int<lower=0,upper=M> sum_y_has_intercept_lobound; // num. submodels w/ lower bounded intercept
  int<lower=0,upper=M> sum_y_has_intercept_upbound; // num. submodels w/ upper bounded intercept
  int<lower=0,upper=M> sum_y_has_aux; // num. submodels w/ auxiliary term
  int<lower=0,upper=1> has_weights;  // 1 = Yes
  
  // data for longitudinal submodel(s)
  int<lower=1> family[M];  // family
  int<lower=1> link[M];  // link function, varies by family
  int<lower=0,upper=1> y_has_intercept[M];  // 1 = yes
  int<lower=0,upper=1> y_has_intercept_unbound[M];  // intercept unbounded
  int<lower=0,upper=1> y_has_intercept_lobound[M];  // intercept lower bound at 0
  int<lower=0,upper=1> y_has_intercept_upbound[M];  // intercept upper bound at 0
  //int<lower=0,upper=1> y_has_offset;  // 1 = Yes
  int<lower=0,upper=1> y_has_aux[M];  // 1 = Yes
  vector[sum_y_real_N] y_real;  // outcome vector, reals                
  int<lower=0> y_int[sum_y_int_N];  // outcome vector, integers                
  vector[sum_y_K] y_xbar; // predictor means
  matrix[sum_y_N,sum_y_K] y_X;  // predictor matrix, centred          
  vector[sum_y_N] y_weights;  // weights, set to zero if not used
  int<lower=0> trials[sum_y_N];  // num. binomial trials, set to zero if not used
  //vector[sum_y_N*y_has_offset] y_offset; 
  int<lower=0> num_non_zero;   // number of non-zero elements in the Z matrix
  vector[num_non_zero] w;  // non-zero elements in the implicit Z matrix
  int<lower=0> v[num_non_zero]; // column indices for w  
  int<lower=0> u[(sum_y_N+1)];  // where the non-zeros start in each row 
  
  // data for event submodel
  int<lower=0,upper=3> basehaz_type;  // 1 = Weibull, 2 = B-splines, 3 = piecewise
  int<lower=0> basehaz_df;  // df for baseline hazard
  int<lower=0,upper=1> e_has_intercept;  // 1 = yes
  int<lower=0> nrow_y_Xq;     // num. rows in long. predictor matrix at quad points
  int<lower=0> nrow_e_Xq;   // num. rows in event predictor matrix at quad points
  matrix[nrow_e_Xq,e_K] e_Xq;         // predictor matrix (event submodel) at quadpoints, centred
  vector[nrow_e_Xq] e_times;          // event times and unstandardised quadrature points
  matrix[nrow_e_Xq,basehaz_df] e_basehaz_X; // design matrix (basis terms) for baseline hazard
  vector[nrow_e_Xq] e_d; // event indicator, followed by dummy indicator for quadpoints
  vector[e_K] e_xbar;   // predictor means (event submodel)
  vector[Npat] e_weights;           // weights, set to zero if not used
  vector[Npat_times_quadnodes] e_weights_rep;   // repeated weights, set to zero if not used
  vector[Npat_times_quadnodes] quadweight;
    
  // data for association structure
  int<lower=0,upper=1> assoc;           // 0 = no assoc structure, 1 = any assoc structure
  int<lower=0,upper=1> assoc_uses[8];   // which components required to build association terms
  int<lower=0,upper=1> has_assoc[18,M]; // which association terms does each submodel use
  int<lower=0> sum_size_which_b;           // num. of shared random effects
  int<lower=0> size_which_b[M];            // num. of shared random effects for each long submodel
  int<lower=1> which_b_zindex[sum_size_which_b];  // which random effects are shared for each long submodel
  int<lower=0> sum_size_which_coef;        // num. of shared random effects incl fixed component
  int<lower=0> size_which_coef[M];         // num. of shared random effects incl fixed component for each long submodel
  int<lower=1> which_coef_zindex[sum_size_which_coef];  // which random effects are shared incl fixed component for each long submodel
  int<lower=1> which_coef_xindex[sum_size_which_coef];  // which fixed effects are shared
  int<lower=0,upper=a_K> sum_a_K_data;  // total num pars used in assoc*data interactions
  int<lower=0,upper=sum_a_K_data> a_K_data[M*4]; // num pars used in assoc*data interactions, by submodel and by ev/es/mv/ms interactions
  int<lower=0> sum_size_which_interactions; // total num pars used in assoc*assoc interactions
  int<lower=0,upper=sum_size_which_interactions> size_which_interactions[M*4]; // num pars used in assoc*data interactions, by submodel and by evev/evmv/mvev/mvmv interactions
  int<lower=1> which_interactions[sum_size_which_interactions];  // which terms to interact with

  // data for calculating eta in GK quadrature
  matrix[M*nrow_y_Xq*(assoc_uses[1]>0),sum_y_K] y_Xq_eta; // predictor matrix (long submodel) at quadpoints, centred     
  int<lower=0> num_non_zero_Zq_eta;    // number of non-zero elements in the Z matrix (at quadpoints)
  vector[num_non_zero_Zq_eta] w_Zq_eta;  // non-zero elements in the implicit Z matrix (at quadpoints)
  int<lower=0> v_Zq_eta[num_non_zero_Zq_eta]; // column indices for w (at quadpoints)
  int<lower=0> u_Zq_eta[(M*nrow_y_Xq*(assoc_uses[1]>0) + 1)]; // where the non-zeros start in each row (at quadpoints)

  // data for calculating slope in GK quadrature
  real<lower=0> eps;  // time shift used for numerically calculating derivative
  matrix[M*nrow_y_Xq*(assoc_uses[2]>0),sum_y_K] 
    y_Xq_eps; // predictor matrix (long submodel) at quadpoints plus time shift of epsilon              
  int<lower=0> num_non_zero_Zq_eps;        // number of non-zero elements in the Zq_eps matrix (at quadpoints plus time shift of epsilon)
  vector[num_non_zero_Zq_eps] w_Zq_eps;    // non-zero elements in the implicit Zq_eps matrix (at quadpoints plus time shift of epsilon)
  int<lower=0> v_Zq_eps[num_non_zero_Zq_eps]; // column indices for w (at quadpoints plus time shift of epsilon)
  int<lower=0> u_Zq_eps[(M*nrow_y_Xq*(assoc_uses[2]>0) + 1)]; 
    // where the non-zeros start in each row (at quadpoints plus time shift of epsilon)

  // data for calculating lag in GK quadrature
  matrix[M*nrow_y_Xq*(assoc_uses[3]>0),sum_y_K] 
    y_Xq_lag; // predictor matrix (long submodel) at lagged quadpoints            
  int<lower=0> num_non_zero_Zq_lag;        // number of non-zero elements in the Zq_lag matrix (at lagged quadpoints)
  vector[num_non_zero_Zq_lag] w_Zq_lag;    // non-zero elements in the implicit Zq_lag matrix (at lagged quadpointsn)
  int<lower=0> v_Zq_lag[num_non_zero_Zq_lag]; // column indices for w (at lagged quadpoints)
  int<lower=0> u_Zq_lag[(M*nrow_y_Xq*(assoc_uses[3]>0) + 1)]; 
    // where the non-zeros start in each row (at lagged quadpoints)

  // data for calculating auc in GK quadrature
  int<lower=0> nrow_y_Xq_auc;     // num. rows in long. predictor matrix at auc quad points
  int<lower=0> auc_quadnodes;              // num. of nodes for Gauss-Kronrod quadrature for area under marker trajectory 
  vector[nrow_y_Xq_auc*(assoc_uses[4]>0)] auc_quadweights;
  matrix[M*nrow_y_Xq_auc*(assoc_uses[4]>0),sum_y_K] 
    y_Xq_auc; // predictor matrix (long submodel) at auc quadpoints            
  int<lower=0> num_non_zero_Zq_auc;        // number of non-zero elements in the Zq_lag matrix (at auc quadpoints)
  vector[num_non_zero_Zq_auc] w_Zq_auc;    // non-zero elements in the implicit Zq_lag matrix (at auc quadpointsn)
  int<lower=0> v_Zq_auc[num_non_zero_Zq_auc]; // column indices for w (at auc quadpoints)
  int<lower=0> u_Zq_auc[(M*nrow_y_Xq_auc*(assoc_uses[4]>0) + 1)]; 
    // where the non-zeros start in each row (at auc quadpoints)

  // data for calculating interactions in GK quadrature
  vector[nrow_e_Xq] y_Xq_data[sum_a_K_data]; // design matrix for interacting with ev/es/mv/ms at quadpoints            
  // data for random effects model
  // declares t, p[t], l[t], q, len_theta_L, shape, scale, {len_}concentration, {len_}regularization
  #include "glmer_stuff.stan"  
  int<lower=1,upper=t> t_i;     // index of grouping factor corresponding to patient-level
  int<lower=0> pmat[t,M];       // num. random effects for each grouping factor (t) in each submodel (M)
  int<lower=0> qmat[t,M];       // = l * pmat --> num. random coefs for each grouping factor in each submodel
  int<lower=0> q1[t];           // = l * p --> num. random coefs for each grouping factor
  int<lower=0> q2[M];           // num. random coefs for each submodel
  int<lower=0> len_b;           // length of the b vector

  // prior family: 0 = none, 1 = normal, 2 = student_t, 3 = hs, 4 = hs_plus, 
  //   5 = laplace, 6 = lasso, 7 = product_normal
  int<lower=0,upper=7> y_prior_dist;
  int<lower=0,upper=7> e_prior_dist;
  int<lower=0,upper=7> a_prior_dist;
  int<lower=0,upper=2> y_prior_dist_for_intercept[M];  
  int<lower=0,upper=2> e_prior_dist_for_intercept;
  
  // prior family: 0 = none, 1 = normal, 2 = student_t, 3 = exponential
  int<lower=0,upper=3> y_prior_dist_for_aux[M];  
  int<lower=0,upper=3> e_prior_dist_for_aux;
    
  // hyperparameter values are set to 0 if there is no prior
  vector[sum_y_K]             y_prior_mean;
  vector[M]                   y_prior_mean_for_intercept;
  vector[M]                   y_prior_mean_for_aux;
  vector[e_K]                 e_prior_mean;
  real                        e_prior_mean_for_intercept;
  vector[basehaz_df]          e_prior_mean_for_aux;
  vector[a_K]                 a_prior_mean;
  vector<lower=0>[sum_y_K]    y_prior_scale;
  vector<lower=0>[M]          y_prior_scale_for_intercept;
  vector<lower=0>[M]          y_prior_scale_for_aux;
  vector<lower=0>[e_K]        e_prior_scale;
  real<lower=0>               e_prior_scale_for_intercept;
  vector<lower=0>[basehaz_df] e_prior_scale_for_aux;
  vector<lower=0>[a_K]        a_prior_scale;
  vector<lower=0>[sum_y_K]    y_prior_df;
  vector<lower=0>[M]          y_prior_df_for_intercept;
  vector<lower=0>[M]          y_prior_df_for_aux;
  vector<lower=0>[e_K]        e_prior_df;
  real<lower=0>               e_prior_df_for_intercept; 
  vector<lower=0>[basehaz_df] e_prior_df_for_aux; 
  vector<lower=0>[a_K]        a_prior_df;
  real<lower=0> y_global_prior_scale; // for hs priors only
  real<lower=0> e_global_prior_scale;
  real<lower=0> a_global_prior_scale; 
  real<lower=0> y_global_prior_df;    
  real<lower=0> e_global_prior_df;  
  real<lower=0> a_global_prior_df;    
  int<lower=2>  y_num_normals[y_prior_dist == 7 ? sum_y_K : 0];
  int<lower=2>  e_num_normals[e_prior_dist == 7 ? e_K : 0];
  int<lower=2>  a_num_normals[a_prior_dist == 7 ? a_K : 0];
  
  // flag indicating whether to draw from the prior
  int<lower=0,upper=1> prior_PD;  // 1 = yes
  
  // flag indiciating whether to accumulate lp for each submodel
  int<lower=0,upper=1> long_lp;  // 1 = yes
  int<lower=0,upper=1> event_lp;  // 1 = yes

}
transformed data {
  real poisson_max = pow(2.0, 30.0);
  real sum_log_y[M];
  vector[sum_y_real_N] sqrt_y = rep_vector(not_a_number(), sum_y_real_N);
  vector[sum_y_real_N] log_y  = rep_vector(not_a_number(), sum_y_real_N);
  int<lower=0> y_hs = get_nvars_for_hs(y_prior_dist);
  int<lower=0> e_hs = get_nvars_for_hs(e_prior_dist);                 
  int<lower=0> a_hs = get_nvars_for_hs(a_prior_dist);                 
  int<lower=0> len_z_T = 0;
  int<lower=0> len_var_group = sum(p) * (t > 0);
  int<lower=0> len_rho = sum(p) - t;
  real<lower=0> delta[len_concentration];
  int<lower=1> pos = 1;
 
  // calculate transformations of outcome
  for (m in 1:M) {
    sum_log_y[m] = not_a_number();
    if (family[m] == 2 || family[m] == 3) {
      sum_log_y[m] = sum(log(y_real[y_real_beg[m]:y_real_end[m]]));
    }
    if (family[m] == 3) {
		  sqrt_y[y_real_beg[m]:y_real_end[m]] = sqrt(y_real[y_real_beg[m]:y_real_end[m]]);
		  log_y[y_real_beg[m]:y_real_end[m]] = log(y_real[y_real_beg[m]:y_real_end[m]]);
    }
  }
  
  // prior for covariance
  for (i in 1:t) {
    if (p[i] > 1) {
      for (j in 1:p[i]) {
        delta[pos] = concentration[j];
        pos = pos + 1;
      }
    }
    for (j in 3:p[i]) len_z_T = len_z_T + p[i] - 1;
  }
  
}
parameters {
  // intercepts
  real          y_gamma_unbound[sum_y_has_intercept_unbound];
  real<lower=0> y_gamma_lobound[sum_y_has_intercept_lobound];  
  real<upper=0> y_gamma_upbound[sum_y_has_intercept_upbound];  
  real          e_gamma[e_has_intercept];

  // primative coefficients
  vector[y_prior_dist == 7 ? sum(y_num_normals) : sum_y_K] y_z_beta;
  vector[e_prior_dist == 7 ? sum(e_num_normals) : e_K]     e_z_beta;
  vector[a_prior_dist == 7 ? sum(a_num_normals) : a_K]     a_z_beta;
  
  // unscaled auxiliary parameters (for longitudinal submodels)
  vector<lower=0>[sum_y_has_aux] y_aux_unscaled; # interpretation depends on family!
  #vector<lower=0>[sum_y_noise_N] y_noise; // do not store this

  // unscaled weibull shape parameter, unscaled bs coefs on log basehaz, or unscaled coefs for piecewise constant basehaz
  vector<lower=(basehaz_type == 1 ? 0 : negative_infinity())>[basehaz_df] e_aux_unscaled;       

  // parameters for random effects model
  vector[len_b] z_b;
  vector[len_z_T] z_T;
  vector<lower=0,upper=1>[len_rho] rho;
  vector<lower=0>[len_concentration] zeta;
  vector<lower=0>[t] tau;  

  // parameters for priors
  real<lower=0> y_global[y_hs];
  real<lower=0> e_global[e_hs];
  real<lower=0> a_global[a_hs];
  vector<lower=0>[(y_hs>0)*sum_y_K] y_local[y_hs];
  vector<lower=0>[(e_hs>0)*e_K]     e_local[e_hs];
  vector<lower=0>[(a_hs>0)*a_K]     a_local[a_hs];
  vector<lower=0>[sum_y_K] y_S[y_prior_dist == 5 || y_prior_dist == 6];
  vector<lower=0>[e_K]     e_S[e_prior_dist == 5 || e_prior_dist == 6];
  vector<lower=0>[a_K]     a_S[a_prior_dist == 5 || a_prior_dist == 6];
  real<lower=0> y_one_over_lambda[y_prior_dist == 6];  
  real<lower=0> e_one_over_lambda[e_prior_dist == 6];  
  real<lower=0> a_one_over_lambda[a_prior_dist == 6];  
  
}
transformed parameters {
  // parameters for longitudinal submodel(s)
  vector[sum_y_K] y_beta;
  vector[sum_y_has_aux] y_aux;
  
  // parameters for event submodel
  vector[e_K] e_beta; 
  vector[basehaz_df] e_aux;       
  
  // parameters for GK quadrature  
  vector[nrow_e_Xq] e_eta_q;      // linear predictor (event submodel) evaluated at quadpoints
  vector[nrow_e_Xq] log_basehaz;      // baseline hazard evaluated at quadpoints
  vector[nrow_e_Xq] ll_haz_q;     // log hazard contribution to the log likelihood for the event model at event time and quad points
  vector[Npat] ll_haz_eventtime;  // log hazard contribution to the log likelihood for the event model AT the event time only
  vector[Npat_times_quadnodes] ll_haz_quadtime;    // log hazard for the event model AT the quadrature points only
  vector[Npat_times_quadnodes] ll_surv_eventtime;  // log survival contribution to the log likelihood for the event model AT the event time
  real sum_ll_haz_eventtime;
  real sum_ll_surv_eventtime;
  real ll_event;                                   // log likelihood for the event model    
  
  // parameters for association structure  
  vector[a_K] a_beta;           

  // parameters for random effects model
  vector[len_theta_L] theta_L; 
  vector[len_b] b; 
  vector[len_b] b_by_model; 

  // coefficients
  y_beta = generate_beta(y_z_beta, y_prior_dist, y_prior_mean, y_prior_scale, y_prior_df, y_global,
                         y_local, y_global_prior_scale, y_one_over_lambda, y_S, y_num_normals);  
  e_beta = generate_beta(e_z_beta, e_prior_dist, e_prior_mean, e_prior_scale, e_prior_df, e_global,
                         e_local, e_global_prior_scale, e_one_over_lambda, e_S, e_num_normals);  
  a_beta = generate_beta(a_z_beta, a_prior_dist, a_prior_mean, a_prior_scale, a_prior_df, a_global,
                         a_local, a_global_prior_scale, a_one_over_lambda, a_S, a_num_normals);  

  // auxiliary parameters (for longitudinal submodels)
  // NB aux is not currently used for scaling in hs or hs_plus priors for JM models
  if (sum_y_has_aux > 0) {
    int mark;
    mark = 1;
    for (m in 1:M) {
      if (y_has_aux[m] == 1) {
        if (y_prior_dist_for_aux[m] == 0) // none
          y_aux[mark] = y_aux_unscaled[mark];
        else {
          y_aux[mark] = y_prior_scale_for_aux[m] * y_aux_unscaled[mark];
          if (y_prior_dist_for_aux[m] <= 2) // normal or student_t
            y_aux[mark] = y_aux[mark] + y_prior_mean_for_aux[m];
        }
        mark = mark + 1;
      }
    }
  }        
      
  // parameters for baseline hazard
  if (e_prior_dist_for_aux == 0) // none
    e_aux = e_aux_unscaled;
  else {
    e_aux = e_prior_scale_for_aux .* e_aux_unscaled;
    if (e_prior_dist_for_aux <= 2) // normal or student_t
      e_aux = e_aux + e_prior_mean_for_aux;
  }

  // parameters for random effects model
  if (t > 0) {
    theta_L = make_theta_L(len_theta_L, p, 1.0, tau, scale, zeta, rho, z_T);
    b = make_b(z_b, theta_L, p, l);
    if (M > 1) b_by_model = reorder_b(b, p, pmat, q1, q2, qmat, l, M);
	  else b_by_model = b;
  }
  
  //---------------
  // GK quadrature
  //---------------
  
  // Event submodel: linear predictor at event and quad times
  if (e_K > 0) e_eta_q = e_Xq * e_beta;
  else e_eta_q = rep_vector(0.0, nrow_e_Xq);
  if (e_has_intercept == 1) {
    e_eta_q = e_eta_q + e_gamma[1];
  }
  else {
    // correction to eta if model has no intercept (because X is centered)
    e_eta_q = e_eta_q + dot_product(e_xbar, e_beta); 
  }

  if (assoc == 1) {

    vector[M*nrow_y_Xq*(assoc_uses[1]>0)]     y_eta_q;     // linear predictor (all long submodels) evaluated at quadpoints
    vector[M*nrow_y_Xq*(assoc_uses[2]>0)]     y_eta_q_eps; // linear predictor (all long submodels) evaluated at quadpoints plus time shift of epsilon
    vector[M*nrow_y_Xq*(assoc_uses[3]>0)]     y_eta_q_lag; // linear predictor (all long submodels) evaluated at lagged quadpoints
    vector[M*nrow_y_Xq_auc*(assoc_uses[4]>0)] y_eta_q_auc; // linear predictor (all long submodels) evaluated at auc quadpoints
    vector[nrow_y_Xq]     y_eta_qwide    [M];    
    vector[nrow_y_Xq]     y_eta_qwide_eps[M];  
    vector[nrow_y_Xq]     y_eta_qwide_lag[M];  
    vector[nrow_y_Xq_auc] y_eta_qwide_auc[M];  
    vector[nrow_y_Xq]     y_qwide    [M];         
    vector[nrow_y_Xq]     y_qwide_eps[M];     
    vector[nrow_y_Xq]     y_qwide_lag[M];     
    vector[nrow_y_Xq_auc] y_qwide_auc[M];

    // mark tracks indexing within a_beta vector, which is the 
    // vector of association parameters
    int mark = 0;
    // mark2 tracks indexing within a_K_data vector, which is the 
    // vector specifying the number of columns used for each possible 
    // type of association term by data interaction
    int mark2 = 0;
    // mark3 tracks indexing within size_which_interactions vector
    int mark3 = 0;
    
    //--------------------------------------
    // Prep work for association structures
    //--------------------------------------

    # The m loop take the linear predictors and/or expected
    # values evaluated in the previous step and, for submodel m, 
    # adds the intercept and stores the result in an element of an 
    # array rather than having all submodels in a single vector

    // Linear predictor
    if (assoc_uses[1] == 1) {
      if (sum_y_K > 0) y_eta_q = y_Xq_eta * y_beta;
      else y_eta_q = rep_vector(0.0, (M*nrow_y_Xq));
      //if (y_has_offset == 1) y_eta_q = y_eta_q + y_offset; # how to handle offset?
      y_eta_q = y_eta_q + csr_matrix_times_vector((M*nrow_y_Xq), len_b, w_Zq_eta, v_Zq_eta, u_Zq_eta, b_by_model);
      for (m in 1:M) {
        y_eta_qwide[m] = add_intercept(
          y_eta_q, m, nrow_y_Xq, y_has_intercept, y_has_intercept_unbound, 
          y_has_intercept_lobound, y_has_intercept_upbound, y_gamma_unbound, 
          y_gamma_lobound, y_gamma_upbound, y_xbar, y_beta, y_K);
      }
    }

    // Linear predictor at time plus epsilon
    if (assoc_uses[2] == 1) {
      if (sum_y_K > 0) y_eta_q_eps = y_Xq_eps * y_beta;
      else y_eta_q_eps = rep_vector(0.0, (M*nrow_y_Xq));
      //if (y_has_offset == 1) y_eta_q_eps = y_eta_q_eps + y_offset; # how to handle offset?
      y_eta_q_eps = y_eta_q_eps + csr_matrix_times_vector((M*nrow_y_Xq), len_b, w_Zq_eps, v_Zq_eps, u_Zq_eps, b_by_model);
      for (m in 1:M) {
        y_eta_qwide_eps[m] = add_intercept(
          y_eta_q_eps, m, nrow_y_Xq, y_has_intercept, y_has_intercept_unbound, 
          y_has_intercept_lobound, y_has_intercept_upbound, y_gamma_unbound, 
          y_gamma_lobound, y_gamma_upbound, y_xbar, y_beta, y_K);
      }
    }

    // Linear predictor at lagged time
    if (assoc_uses[3] == 1) {
      if (sum_y_K > 0) y_eta_q_lag = y_Xq_lag * y_beta;
      else y_eta_q_lag = rep_vector(0.0, (M*nrow_y_Xq));
      //if (y_has_offset == 1) y_eta_q_lag = y_eta_q_lag + y_offset; # how to handle offset?
      y_eta_q_lag = y_eta_q_lag + csr_matrix_times_vector((M*nrow_y_Xq), len_b, w_Zq_lag, v_Zq_lag, u_Zq_lag, b_by_model);
      for (m in 1:M) {
        y_eta_qwide_lag[m] = add_intercept(
          y_eta_q_lag, m, nrow_y_Xq, y_has_intercept, y_has_intercept_unbound, 
          y_has_intercept_lobound, y_has_intercept_upbound, y_gamma_unbound, 
          y_gamma_lobound, y_gamma_upbound, y_xbar, y_beta, y_K);
      }
    }  

    // Linear predictor at auc quadpoints
    if (assoc_uses[4] == 1) {
      if (sum_y_K > 0) y_eta_q_auc = y_Xq_auc * y_beta;
      else y_eta_q_auc = rep_vector(0.0, (M*nrow_y_Xq_auc));
      //if (y_has_offset == 1) y_eta_q_auc = y_eta_q_auc + y_offset; # how to handle offset?
      y_eta_q_auc = y_eta_q_auc + csr_matrix_times_vector((M*nrow_y_Xq_auc), len_b, w_Zq_auc, v_Zq_auc, u_Zq_auc, b_by_model);
      for (m in 1:M) {
        y_eta_qwide_auc[m] = add_intercept(
          y_eta_q_auc, m, nrow_y_Xq_auc, y_has_intercept, y_has_intercept_unbound, 
          y_has_intercept_lobound, y_has_intercept_upbound, y_gamma_unbound, 
          y_gamma_lobound, y_gamma_upbound, y_xbar, y_beta, y_K);
      }
    }   

    // Expected value 
    if (assoc_uses[5] == 1) 
      for (m in 1:M) 
        y_qwide[m] = evaluate_mu(y_eta_qwide[m], family[m], link[m]);
      
    // Expected value at time plus epsilon
    if (assoc_uses[6] == 1)
      for (m in 1:M) 
        y_qwide_eps[m] = evaluate_mu(y_eta_qwide_eps[m], family[m], link[m]);

    // Expected value at lagged time
    if (assoc_uses[7] == 1) 
      for (m in 1:M) 
        y_qwide_lag[m] = evaluate_mu(y_eta_qwide_lag[m], family[m], link[m]);

    // Expected value at auc quadpoints
    if (assoc_uses[8] == 1) 
      for (m in 1:M) 
        y_qwide_auc[m] = evaluate_mu(y_eta_qwide_auc[m], family[m], link[m]);

    //---------------------------------
    // Evaluate association structures
    //---------------------------------
    
    # !!! Be careful that indexing of has_assoc matches stan_jm file
 
    for (m in 1:M) {

      // etavalue and any interactions
      mark2 = mark2 + 1; // count even if assoc type isn't used
      if (has_assoc[1,m] == 1) { # etavalue
        mark = mark + 1;
	      e_eta_q = e_eta_q + a_beta[mark] * y_eta_qwide[m];
      }	
      if (has_assoc[11,m] == 1) { # etavalue*data
  	    int tmp;
  	    int j_shift;
  	    if (mark2 == 1) j_shift = 0;
  	    else j_shift = sum(a_K_data[1:(mark2-1)]);
  	    tmp = a_K_data[mark2];  
        for (j in 1:tmp) {
          int sel;
          sel = j_shift + j;
          mark = mark + 1;
          e_eta_q = e_eta_q + (y_eta_qwide[m] .* y_Xq_data[sel]) * a_beta[mark];
        }
      }
      mark3 = mark3 + 1; // count even if assoc type isn't used
      if (has_assoc[15,m] == 1) { # etavalue*etavalue
        for (j in 1:size_which_interactions[mark3]) { 
          int sel;
    	    int j_shift;
     	    if (mark3 == 1) j_shift = 0;
    	    else j_shift = sum(size_which_interactions[1:(mark3-1)]);
    	    sel = which_interactions[j+j_shift];
  	      mark = mark + 1;
          e_eta_q = e_eta_q + (y_eta_qwide[m] .* y_eta_qwide[sel]) * a_beta[mark];  
       }
      }
      mark3 = mark3 + 1; // count even if assoc type isn't used
      if (has_assoc[16,m] == 1) { # etavalue*muvalue
        for (j in 1:size_which_interactions[mark3]) { 
  	      int sel;
    	    int j_shift;
    	    if (mark3 == 1) j_shift = 0;
    	    else j_shift = sum(size_which_interactions[1:(mark3-1)]);
    	    sel = which_interactions[j+j_shift];
  	      mark = mark + 1;
          e_eta_q = e_eta_q + (y_eta_qwide[m] .* y_qwide[sel]) * a_beta[mark];  
        }
      }
      
      // etaslope and any interactions
      mark2 = mark2 + 1;
      if ((has_assoc[2,m] == 1) || (has_assoc[12,m] == 1)) {
        vector[nrow_y_Xq] dydt_eta_q;
        dydt_eta_q = (y_eta_qwide_eps[m] - y_eta_qwide[m]) / eps;
        if (has_assoc[2,m] == 1) { # etaslope
          mark = mark + 1;
          e_eta_q = e_eta_q + a_beta[mark] * dydt_eta_q;
        }
        if (has_assoc[13,m] == 1) { # etaslope*data
    	    int tmp;
    	    int j_shift;
    	    if (mark2 == 1) j_shift = 0;
    	    else j_shift = sum(a_K_data[1:(mark2-1)]);
    	    tmp = a_K_data[mark2];  
          for (j in 1:tmp) {
            int sel;
            sel = j_shift + j;
            mark = mark + 1;
            e_eta_q = e_eta_q + (dydt_eta_q .* y_Xq_data[sel]) * a_beta[mark];
          }    
        }         
      }
      
      // etalag
      if (has_assoc[3,m] == 1) { # etalag
        mark = mark + 1;
        e_eta_q = e_eta_q + a_beta[mark] * y_eta_qwide_lag[m];          
      }  
      
      // etaauc
      if (has_assoc[4,m] == 1) { # etaauc
        vector[nrow_y_Xq] y_eta_q_auc_tmp;  
        mark = mark + 1;
        for (r in 1:nrow_y_Xq) {
          vector[auc_quadnodes] val_tmp;
          vector[auc_quadnodes] wgt_tmp;
          val_tmp = y_eta_qwide_auc[m,((r-1) * auc_quadnodes + 1):(r * auc_quadnodes)];
          wgt_tmp = auc_quadweights[((r-1) * auc_quadnodes + 1):(r * auc_quadnodes)];
          y_eta_q_auc_tmp[r] = sum(wgt_tmp .* val_tmp);
        }
        e_eta_q = e_eta_q + a_beta[mark] * y_eta_q_auc_tmp;          
      }       
      
      // muvalue and any interactions
      mark2 = mark2 + 1;
      if (has_assoc[5,m] == 1) { # muvalue
        mark = mark + 1;
        e_eta_q = e_eta_q + a_beta[mark] * y_qwide[m]; 
      }
      if (has_assoc[12,m] == 1) { # muvalue*data
  	    int tmp;
  	    int j_shift;
  	    if (mark2 == 1) j_shift = 0;
  	    else j_shift = sum(a_K_data[1:(mark2-1)]);
  	    tmp = a_K_data[mark2];  
        for (j in 1:tmp) {
          int sel;
          sel = j_shift + j;
          mark = mark + 1;
          e_eta_q = e_eta_q + (y_qwide[m] .* y_Xq_data[sel]) * a_beta[mark];
        }      
      } 
      mark3 = mark3 + 1; // count even if assoc type isn't used
      if (has_assoc[17,m] == 1) { # muvalue*etavalue
        for (j in 1:size_which_interactions[mark3]) { 
          int sel;
    	    int j_shift;
     	    if (mark3 == 1) j_shift = 0;
    	    else j_shift = sum(size_which_interactions[1:(mark3-1)]);
    	    sel = which_interactions[j+j_shift];
  	      mark = mark + 1;
          e_eta_q = e_eta_q + (y_qwide[m] .* y_eta_qwide[sel]) * a_beta[mark];  
       }
      }      
      mark3 = mark3 + 1; // count even if assoc type isn't used
      if (has_assoc[18,m] == 1) { # muvalue*muvalue
        for (j in 1:size_which_interactions[mark3]) { 
          int sel;
    	    int j_shift;
     	    if (mark3 == 1) j_shift = 0;
    	    else j_shift = sum(size_which_interactions[1:(mark3-1)]);
    	    sel = which_interactions[j+j_shift];
  	      mark = mark + 1;
          e_eta_q = e_eta_q + (y_qwide[m] .* y_qwide[sel]) * a_beta[mark];  
       }
      }      
      
      // muslope and any interactions
      mark2 = mark2 + 1;
      if ((has_assoc[6,m] == 1) || (has_assoc[14,m] == 1)) {
        vector[nrow_y_Xq] dydt_q;
        dydt_q = (y_qwide_eps[m] - y_qwide[m]) / eps;
        if (has_assoc[6,m] == 1) { # muslope
          mark = mark + 1;
          e_eta_q = e_eta_q + a_beta[mark] * dydt_q;          
        }
        if (has_assoc[14,m] == 1) { # muslope*data
    	    int tmp;
    	    int j_shift;
    	    if (mark2 == 1) j_shift = 0;
    	    else j_shift = sum(a_K_data[1:(mark2-1)]);
    	    tmp = a_K_data[mark2];  
          for (j in 1:tmp) {
            int sel;
            sel = j_shift + j;
            mark = mark + 1;
            e_eta_q = e_eta_q + (dydt_q .* y_Xq_data[sel]) * a_beta[mark];
          }          
        } 
      }
      
      // mulag
      if (has_assoc[7,m] == 1) { # mulag
        mark = mark + 1;
        e_eta_q = e_eta_q + a_beta[mark] * y_qwide_lag[m];          
      }

      // muauc
      if (has_assoc[8,m] == 1) { # muauc
        vector[nrow_y_Xq] y_qwide_auc_tmp;  
        mark = mark + 1;
        for (r in 1:nrow_y_Xq) {
          vector[auc_quadnodes] val_tmp;
          vector[auc_quadnodes] wgt_tmp;
          val_tmp = y_qwide_auc[m,((r-1) * auc_quadnodes + 1):(r * auc_quadnodes)];
          wgt_tmp = auc_quadweights[((r-1) * auc_quadnodes + 1):(r * auc_quadnodes)];
          y_qwide_auc_tmp[r] = sum(wgt_tmp .* val_tmp);
        }
        e_eta_q = e_eta_q + a_beta[mark] * y_qwide_auc_tmp;          
      }  

    }
    
    // shared random effects
  	if (sum_size_which_b > 0) {
  	  int mark_beg;  // used to define segment of a_beta
  	  int mark_end;
  	  matrix[nrow_e_Xq,sum_size_which_b] x_assoc_shared_b;	  
      mark_beg = mark + 1;	  
  	  mark_end = mark + sum_size_which_b;
  	  x_assoc_shared_b = make_x_assoc_shared_b(
  	    b, l, p, pmat, Npat, quadnodes, which_b_zindex,
  	    sum_size_which_b, size_which_b, t_i, M);
  	  e_eta_q = e_eta_q + x_assoc_shared_b * a_beta[mark_beg:mark_end];
  	  mark = mark + sum_size_which_b;
    }	
  	if (sum_size_which_coef > 0) {
  	  int mark_beg;  // used to define segment of a_beta
  	  int mark_end;
  	  matrix[nrow_e_Xq,sum_size_which_coef] x_assoc_shared_coef;	  
      mark_beg = mark + 1;	  
  	  mark_end = mark + sum_size_which_coef;
  	  x_assoc_shared_coef = make_x_assoc_shared_coef(
  	    b, y_beta, y_K, M, t_i, l, p, pmat, Npat, quadnodes,
  	    sum_size_which_coef, size_which_coef,
  	    which_coef_zindex, which_coef_xindex,
  	    y_has_intercept, y_has_intercept_unbound,
  	    y_has_intercept_lobound, y_has_intercept_upbound,
  	    y_gamma_unbound, y_gamma_lobound, y_gamma_upbound);
  	  e_eta_q = e_eta_q + x_assoc_shared_coef * a_beta[mark_beg:mark_end];
  	  mark = mark + sum_size_which_coef;
    }    
  }

  //-----------------------------------
  // Log-likelihood for event submodel
  //-----------------------------------
  
  // Log baseline hazard at event times and unstandardised quadrature points
  if (basehaz_type == 1) 
    log_basehaz = log(e_aux[1]) + e_basehaz_X * (e_aux - 1);
  else log_basehaz = e_basehaz_X * e_aux;	

  // Log hazard at event times and unstandardised quadrature points
  ll_haz_q = e_d .* (log_basehaz + e_eta_q);
					  
  // Partition event times and quadrature points
  ll_haz_eventtime = segment(ll_haz_q, 1, Npat);
  ll_haz_quadtime  = segment(ll_haz_q, (Npat + 1), Npat_times_quadnodes);

  // Log survival contribution to the likelihood (by summing over the 
  // quadrature points to get the approximate integral)
  // NB quadweight already incorporates the (b-a)/2 scaling such that the
  // integral is evaluated over limits (a,b) rather than (-1,+1)
  ll_surv_eventtime = quadweight .* exp(ll_haz_quadtime);        

  // Log likelihood for event model
  if (has_weights == 0) {  # unweighted log likelihood
    sum_ll_haz_eventtime = sum(ll_haz_eventtime);
    sum_ll_surv_eventtime = sum(ll_surv_eventtime);
  } 
  else {  # weighted log likelihood
    sum_ll_haz_eventtime = sum(e_weights .* ll_haz_eventtime);
    sum_ll_surv_eventtime = sum(e_weights_rep .* ll_surv_eventtime);
  }
  ll_event = sum_ll_haz_eventtime - sum_ll_surv_eventtime;				  

}

model {
  int aux_mark;
  int nois_mark;
  vector[sum_y_N] y_eta;                                     

  aux_mark = 1;
  nois_mark = 1;
  
  // Longitudinal submodel(s): regression equations

  if (sum_y_K > 0) y_eta = y_X * y_beta;
  else y_eta = rep_vector(0.0, sum_y_N);
  //if (y_has_offset == 1) y_eta = y_eta + y_offset;

  y_eta = y_eta + csr_matrix_times_vector(sum_y_N, len_b, w, v, u, b_by_model);
  for (m in 1:M) {
    int mark_beg;
    int mark_end;
    vector[y_N[m]] y_eta_tmp;	
    
    y_eta_tmp = y_eta[y_beg[m]:y_end[m]]; 

    if (y_has_intercept[m] == 1) {
      if (y_has_intercept_unbound[m] == 1)
        y_eta_tmp = y_eta_tmp + 
                    y_gamma_unbound[sum(y_has_intercept_unbound[1:m])];
      else if (y_has_intercept_lobound[m] == 1)
        y_eta_tmp = y_eta_tmp - min(y_eta_tmp) + 
                    y_gamma_lobound[sum(y_has_intercept_lobound[1:m])];
      else if (y_has_intercept_upbound[m] == 1)
        y_eta_tmp = y_eta_tmp - max(y_eta_tmp) + 
                    y_gamma_lobound[sum(y_has_intercept_upbound[1:m])];					
    }	
    
    if (m == 1) mark_beg = 1;
    else mark_beg = sum(y_K[1:(m-1)]) + 1;
    mark_end = sum(y_K[1:m]); 
    // correction to eta if model has no intercept (if X is centered)
    y_eta_tmp = y_eta_tmp + 
      dot_product(y_xbar[mark_beg:mark_end], y_beta[mark_beg:mark_end]); 
    
#    if (family[m] == 8) {  # poisson-gamma mixture
#      if      (link[m] == 1) y_eta_tmp = y_eta_tmp + log(y_aux[aux_mark]) + log(y_noise[nois_mark]);
#      else if (link[m] == 2) y_eta_tmp = y_eta_tmp * y_aux[aux_mark] .* y_noise[nois_mark];
#      else                   y_eta_tmp = y_eta_tmp + sqrt(y_aux[aux_mark]) + sqrt_vec(y_noise[nois_mark]);
#    }    

    // Log-likelihood for longitudinal submodel(s)
    if (has_weights == 0 && prior_PD == 0 && long_lp == 1) { # unweighted log-likelihoods
      if (family[m] == 1) {
        if (link[m] == 1)      target += normal_lpdf(y_real[y_real_beg[m]:y_real_end[m]] | y_eta_tmp, y_aux[aux_mark]);
        else if (link[m] == 2) target += lognormal_lpdf(y_real[y_real_beg[m]:y_real_end[m]] | y_eta_tmp, y_aux[aux_mark]);
        else target += normal_lpdf(y_real[y_real_beg[m]:y_real_end[m]] | divide_real_by_vector(1, y_eta_tmp), y_aux[aux_mark]);
      }
      else if (family[m] == 2) {
        target += GammaReg(y_real[y_real_beg[m]:y_real_end[m]], y_eta_tmp, y_aux[aux_mark], link[m], sum_log_y[m]);
      }
      else if (family[m] == 3) {
	      vector[y_N[m]] sqrt_y_tmp;
		    sqrt_y_tmp = sqrt_y[y_real_beg[m]:y_real_end[m]]; 
        target += inv_gaussian(y_real[y_real_beg[m]:y_real_end[m]], 
                               linkinv_inv_gaussian(y_eta_tmp, link[m]), 
                               y_aux[aux_mark], sum_log_y[m], sqrt_y_tmp);
      }
	    else if (family[m] == 4) {
		    vector[y_N01[m,1]] y_eta0_tmp;
		    vector[y_N01[m,2]] y_eta1_tmp;
	      real dummy;  // irrelevant but useful for testing
		    y_eta0_tmp = segment(y_eta_tmp, 1, y_N01[m,1]);
		    y_eta1_tmp = segment(y_eta_tmp, (y_N01[m,1] + 1), y_N01[m,2]);
	      dummy = ll_bern_lp(y_eta0_tmp, y_eta1_tmp, link[m], y_N01[m,]);	  
	    }
	    else if (family[m] == 5) {
		    int trials_tmp[y_N[m]];		
	      real dummy;  // irrelevant but useful for testing
		    trials_tmp = trials[y_beg[m]:y_end[m]];
        dummy = ll_binom_lp(y_int[y_int_beg[m]:y_int_end[m]], trials_tmp, y_eta_tmp, link[m]);	  
	    }
	    else if (family[m] == 6 || family[m] == 8) {
        if (link[m] == 1) target += poisson_log_lpmf(y_int[y_int_beg[m]:y_int_end[m]] | y_eta_tmp);
        else target += poisson_lpmf(y_int[y_int_beg[m]:y_int_end[m]] | linkinv_count(y_eta_tmp, link[m]));
	    }
	    else if (family[m] == 7) {
  	    if (link[m] == 1) target += neg_binomial_2_log_lpmf(y_int[y_int_beg[m]:y_int_end[m]] | y_eta_tmp, y_aux[aux_mark]);
        else target += neg_binomial_2_lpmf(y_int[y_int_beg[m]:y_int_end[m]] | 
                                           linkinv_count(y_eta_tmp, link[m]), y_aux[aux_mark]);
	    }	    
    }    
    else if (prior_PD == 0 && long_lp == 1) { # weighted log-likelihoods
  	  vector[y_N[m]] y_weights_tmp;	  
  	  vector[y_N[m]] summands;
      y_weights_tmp = y_weights[y_beg[m]:y_end[m]];	  
  	  if (family[m] == 1) {
  	    summands = pw_gauss(y_real[y_real_beg[m]:y_real_end[m]], y_eta_tmp, y_aux[aux_mark], link[m]);
  	    target += dot_product(y_weights_tmp, summands);	  	    
  	  }
  	  else if (family[m] == 2) {
  	    summands = pw_gamma(y_real[y_real_beg[m]:y_real_end[m]], y_eta_tmp, y_aux[aux_mark], link[m]);
  	    target += dot_product(y_weights_tmp, summands);	  	    
  	  }
  	  else if (family[m] == 3) {
  	    vector[y_N[m]] log_y_tmp;	  
  	    vector[y_N[m]] sqrt_y_tmp;  	    
    		log_y_tmp = log_y[y_beg[m]:y_end[m]];
    		sqrt_y_tmp = sqrt_y[y_beg[m]:y_end[m]];
  	    summands = pw_inv_gaussian(y_real[y_real_beg[m]:y_real_end[m]], y_eta_tmp, y_aux[aux_mark], 
  		                           link[m], log_y_tmp, sqrt_y_tmp);
  	    target += dot_product(y_weights_tmp, summands);
  	  }
	    else if (family[m] == 4) {
    		vector[y_N01[m,1]] y_weights0_tmp;
    		vector[y_N01[m,2]] y_weights1_tmp;
    		vector[y_N01[m,1]] y_eta0_tmp;
    		vector[y_N01[m,2]] y_eta1_tmp;
    		y_eta0_tmp = segment(y_eta_tmp, 1, y_N01[m,1]);
    		y_eta1_tmp = segment(y_eta_tmp, (y_N01[m,1] + 1), y_N01[m,2]);
    		y_weights0_tmp = segment(y_weights_tmp, 1, y_N01[m,1]);
    		y_weights1_tmp = segment(y_weights_tmp, (y_N01[m,1] + 1), y_N01[m,2]);		
        target += dot_product(y_weights0_tmp, pw_bern(0, y_eta0_tmp, link[m]));
        target += dot_product(y_weights1_tmp, pw_bern(1, y_eta1_tmp, link[m]));
  	  }
  	  else if (family[m] == 5) {
    		int trials_tmp[y_N[m]];		
    		trials_tmp = trials[y_beg[m]:y_end[m]];
        target += dot_product(y_weights_tmp, 
  		                      pw_binom(y_int[y_int_beg[m]:y_int_end[m]], trials_tmp, y_eta_tmp, link[m]));
  	  }
  	  else if (family[m] == 6 || family[m] == 8) {
        target += dot_product(y_weights_tmp, pw_pois(y_int[y_int_beg[m]:y_int_end[m]], y_eta_tmp, link[m]));    		
  	  }
  	  else if (family[m] == 7) {
        target += dot_product(y_weights_tmp, pw_nb(y_int[y_int_beg[m]:y_int_end[m]], y_eta_tmp, y_aux[aux_mark], link[m]));
  	  }  	  
    }
    
    if (y_has_aux[m] == 1) aux_mark = aux_mark + 1;
#    if (y_has_noise[m] == 1)      nois_mark = nois_mark + 1;
  }  
                           
  // Log-likelihood for event submodel  
  // NB weights already incorporated in ll calculation in transformed param block
  if (prior_PD == 0 && event_lp == 1) target += ll_event;  
  
  //--------
  // Priors
  //--------
  
  // Log-priors for coefficients
  beta_lp(y_z_beta, y_prior_dist, y_prior_scale, y_prior_df, 
          y_global_prior_df, y_local, y_global, y_S, y_one_over_lambda)
  beta_lp(e_z_beta, e_prior_dist, e_prior_scale, e_prior_df, 
          e_global_prior_df, e_local, e_global, e_S, e_one_over_lambda)
  beta_lp(a_z_beta, a_prior_dist, a_prior_scale, a_prior_df, 
          a_global_prior_df, a_local, a_global, a_S, a_one_over_lambda)

  // Log-priors for intercept parameters
  if (sum_y_has_intercept > 0) {
    int mark1 = 1; // indexing for unbounded intercepts
    int mark2 = 1; // indexing for lower bounded intercepts
    int mark3 = 1; // indexing for upper bounded intercepts
    for (m in 1:M) {
      // unbounded intercept
      if (y_has_intercept_unbound[m] == 1) {
        gamma_lp(y_gamma_unbound[mark1], y_prior_dist_for_intercept[m], 
                 y_prior_mean_for_intercept[m], y_prior_scale_for_intercept[m], 
                 y_prior_df_for_intercept[m]);
        mark1 = mark1 + 1;
      }
      // lower bounded intercept
      else if (y_has_intercept_lobound[m] == 1) {
        gamma_lp(y_gamma_lobound[mark2], y_prior_dist_for_intercept[m], 
                 y_prior_mean_for_intercept[m], y_prior_scale_for_intercept[m], 
                 y_prior_df_for_intercept[m]);
        mark2 = mark2 + 1;
      }
      // upper bounded intercept
      else if (y_has_intercept_lobound[m] == 1) {
        gamma_lp(y_gamma_lobound[mark3], y_prior_dist_for_intercept[m], 
                 y_prior_mean_for_intercept[m], y_prior_scale_for_intercept[m], 
                 y_prior_df_for_intercept[m]);
        mark3 = mark3 + 1;
      }
    }
  }
  if (e_has_intercept == 1) 
    gamma_lp(e_gamma[1], e_prior_dist_for_intercept, e_prior_mean_for_intercept, 
             e_prior_scale_for_intercept, e_prior_df_for_intercept);  

  // Log-prior for auxiliary parameters
  aux_mark = 1;
  for (m in 1:M) { 
    if (y_has_aux[m] == 1) {
      aux_lp(y_aux_unscaled[aux_mark], y_prior_dist_for_aux[m],
             y_prior_scale_for_aux[m], y_prior_df_for_aux[m])
      aux_mark = aux_mark + 1;
    } 
  }
   
  // Log-prior for baseline hazard parameters
  basehaz_lp(e_aux_unscaled, e_prior_dist_for_aux, 
             e_prior_scale_for_aux, e_prior_df_for_aux)
  
  // Log-prior for random effects model
  if (t > 0) decov_lp(z_b, z_T, rho, zeta, tau, 
                      regularization, delta, shape, t, p);
   
}
