  // prior family: 0 = none, 1 = normal, 2 = student_t, 3 = hs, 4 = hs_plus, 
  //   5 = laplace, 6 = lasso
  int<lower=0,upper=6> a_prior_dist;

  // data for association structure
  int<lower=0> a_K;                     // num. of association parameters
  int<lower=0,upper=1> assoc;           // 0 = no assoc structure, 1 = any assoc structure
  int<lower=0,upper=1> assoc_uses[6];   // which components required to build association terms
  int<lower=0,upper=1> has_assoc[16,M]; // which association terms does each submodel use
  int<lower=0,upper=1> has_clust;    // 1 = has clustering below patient level
  matrix<lower=0,upper=nrow_y_Xq>[nrow_e_Xq,nrow_y_Xq] clust_mat; // design matrix used for summing across lower level clusters for each patient
  int<lower=0> sum_size_which_b;        // num. of shared random effects
  int<lower=0> size_which_b[M];         // num. of shared random effects for each long submodel
  int<lower=1> which_b_zindex[sum_size_which_b]; // which random effects are shared for each long submodel
  int<lower=0> sum_size_which_coef;     // num. of shared random effects incl fixed component
  int<lower=0> size_which_coef[M];      // num. of shared random effects incl fixed component for each long submodel
  int<lower=1> which_coef_zindex[sum_size_which_coef]; // which random effects are shared incl fixed component
  int<lower=1> which_coef_xindex[sum_size_which_coef]; // which fixed effects are shared
  int<lower=0,upper=a_K> sum_a_K_data;  // total num pars used in assoc*data interactions
  int<lower=0,upper=sum_a_K_data> a_K_data[M*4]; // num pars used in assoc*data interactions, by submodel and by ev/es/mv/ms interactions
  int<lower=0> sum_size_which_interactions; // total num pars used in assoc*assoc interactions
  int<lower=0,upper=sum_size_which_interactions> size_which_interactions[M*4]; // num pars used in assoc*assoc interactions, by submodel and by evev/evmv/mvev/mvmv interactions
  int<lower=1> which_interactions[sum_size_which_interactions];  // which terms to interact with

  // data for calculating eta in GK quadrature
  matrix[M*nrow_y_Xq*(assoc_uses[1]>0),K] y_Xq_eta; // predictor matrix (long submodel) at quadpoints, centred     
  int<lower=0> nnz_Zq_eta;    // number of non-zero elements in the Z matrix (at quadpoints)
  vector[nnz_Zq_eta] w_Zq_eta;  // non-zero elements in the implicit Z matrix (at quadpoints)
  int<lower=0> v_Zq_eta[nnz_Zq_eta]; // column indices for w (at quadpoints)
  int<lower=0> u_Zq_eta[(M*nrow_y_Xq*(assoc_uses[1]>0) + 1)]; // where the non-zeros start in each row (at quadpoints)

  // data for calculating slope in GK quadrature
  real<lower=0> eps;  // time shift used for numerically calculating derivative
  matrix[M*nrow_y_Xq*(assoc_uses[2]>0),K] 
    y_Xq_eps; // predictor matrix (long submodel) at quadpoints plus time shift of epsilon              
  int<lower=0> nnz_Zq_eps;        // number of non-zero elements in the Zq_eps matrix (at quadpoints plus time shift of epsilon)
  vector[nnz_Zq_eps] w_Zq_eps;    // non-zero elements in the implicit Zq_eps matrix (at quadpoints plus time shift of epsilon)
  int<lower=0> v_Zq_eps[nnz_Zq_eps]; // column indices for w (at quadpoints plus time shift of epsilon)
  int<lower=0> u_Zq_eps[(M*nrow_y_Xq*(assoc_uses[2]>0) + 1)]; 
    // where the non-zeros start in each row (at quadpoints plus time shift of epsilon)

  // data for calculating auc in GK quadrature
  int<lower=0> nrow_y_Xq_auc;     // num. rows in long. predictor matrix at auc quad points
  int<lower=0> auc_quadnodes;              // num. of nodes for Gauss-Kronrod quadrature for area under marker trajectory 
  vector[nrow_y_Xq_auc*(assoc_uses[3]>0)] auc_quadweights;
  matrix[M*nrow_y_Xq_auc*(assoc_uses[3]>0),K] 
    y_Xq_auc; // predictor matrix (long submodel) at auc quadpoints            
  int<lower=0> nnz_Zq_auc;        // number of non-zero elements in the Zq_lag matrix (at auc quadpoints)
  vector[nnz_Zq_auc] w_Zq_auc;    // non-zero elements in the implicit Zq_lag matrix (at auc quadpointsn)
  int<lower=0> v_Zq_auc[nnz_Zq_auc]; // column indices for w (at auc quadpoints)
  int<lower=0> u_Zq_auc[(M*nrow_y_Xq_auc*(assoc_uses[3]>0) + 1)]; 
    // where the non-zeros start in each row (at auc quadpoints)

  // data for calculating assoc*data interactions in GK quadrature
  vector[nrow_y_Xq] y_Xq_data[sum_a_K_data]; // design matrix for interacting with ev/es/mv/ms at quadpoints 
