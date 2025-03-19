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
    array[6,3] int<lower=0,upper=1> assoc_uses;

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

    array[3] int<lower=0> nrow_y_Xq;     // num. rows in long. predictor matrix at quadpoints

    // fe design matrix at quadpoints
    matrix[assoc_uses[1,1] == 1 ? nrow_y_Xq[1] : 0, yK[1]] y1_xq_eta;
    matrix[assoc_uses[1,2] == 1 ? nrow_y_Xq[2] : 0, yK[2]] y2_xq_eta;
    matrix[assoc_uses[1,3] == 1 ? nrow_y_Xq[3] : 0, yK[3]] y3_xq_eta;
    
    // offset values at quadpoints
    vector[has_offset[1] && assoc_uses[1,1] == 1 ? nrow_y_Xq[1] : 0] y1_offset_eta;
    vector[has_offset[2] && assoc_uses[1,2] == 1 ? nrow_y_Xq[2] : 0] y2_offset_eta;
    vector[has_offset[3] && assoc_uses[1,3] == 1 ? nrow_y_Xq[3] : 0] y3_offset_eta;

    // re design matrix at quadpoints, group factor 1
    array[bK1_len[1]] vector[assoc_uses[1,1] == 1 && bK1_len[1] > 0 ? nrow_y_Xq[1] : 0] y1_z1q_eta;
    array[bK1_len[2]] vector[assoc_uses[1,2] == 1 && bK1_len[2] > 0 ? nrow_y_Xq[2] : 0] y2_z1q_eta;
    array[bK1_len[3]] vector[assoc_uses[1,3] == 1 && bK1_len[3] > 0 ? nrow_y_Xq[3] : 0] y3_z1q_eta;
    array[assoc_uses[1,1] == 1 && bK1_len[1] > 0 ? nrow_y_Xq[1] : 0] int<lower=0> y1_z1q_id_eta;
    array[assoc_uses[1,2] == 1 && bK1_len[2] > 0 ? nrow_y_Xq[2] : 0] int<lower=0> y2_z1q_id_eta;
    array[assoc_uses[1,3] == 1 && bK1_len[3] > 0 ? nrow_y_Xq[3] : 0] int<lower=0> y3_z1q_id_eta;

    // re design matrix at quadpoints, group factor 2
    array[bK2_len[1]] vector[assoc_uses[1,1] == 1 && bK2_len[1] > 0 ? nrow_y_Xq[1] : 0] y1_z2q_eta;
    array[bK2_len[2]] vector[assoc_uses[1,2] == 1 && bK2_len[2] > 0 ? nrow_y_Xq[2] : 0] y2_z2q_eta;
    array[bK2_len[3]] vector[assoc_uses[1,3] == 1 && bK2_len[3] > 0 ? nrow_y_Xq[3] : 0] y3_z2q_eta;
    array[assoc_uses[1,1] == 1 && bK2_len[1] > 0 ? nrow_y_Xq[1] : 0] int<lower=0> y1_z2q_id_eta;
    array[assoc_uses[1,2] == 1 && bK2_len[2] > 0 ? nrow_y_Xq[2] : 0] int<lower=0> y2_z2q_id_eta;
    array[assoc_uses[1,3] == 1 && bK2_len[3] > 0 ? nrow_y_Xq[3] : 0] int<lower=0> y3_z2q_id_eta;

  //---- data for calculating derivative of eta in GK quadrature

    // fe design matrix at quadpoints
    matrix[assoc_uses[2,1] == 1 ? nrow_y_Xq[1] : 0, yK[1]] y1_xq_eps;
    matrix[assoc_uses[2,2] == 1 ? nrow_y_Xq[2] : 0, yK[2]] y2_xq_eps;
    matrix[assoc_uses[2,3] == 1 ? nrow_y_Xq[3] : 0, yK[3]] y3_xq_eps;
    
    // offset values at quadpoints
    vector[has_offset[1] && assoc_uses[2,1] == 1 ? nrow_y_Xq[1] : 0] y1_offset_eps;
    vector[has_offset[2] && assoc_uses[2,2] == 1 ? nrow_y_Xq[2] : 0] y2_offset_eps;
    vector[has_offset[3] && assoc_uses[2,3] == 1 ? nrow_y_Xq[3] : 0] y3_offset_eps;

    // re design matrix at quadpoints, group factor 1
    array[bK1_len[1]] vector[assoc_uses[2,1] == 1 && bK1_len[1] > 0 ? nrow_y_Xq[1] : 0] y1_z1q_eps;
    array[bK1_len[2]] vector[assoc_uses[2,2] == 1 && bK1_len[2] > 0 ? nrow_y_Xq[2] : 0] y2_z1q_eps;
    array[bK1_len[3]] vector[assoc_uses[2,3] == 1 && bK1_len[3] > 0 ? nrow_y_Xq[3] : 0] y3_z1q_eps;
    array[assoc_uses[2,1] == 1 && bK1_len[1] > 0 ? nrow_y_Xq[1] : 0] int<lower=0> y1_z1q_id_eps;
    array[assoc_uses[2,2] == 1 && bK1_len[2] > 0 ? nrow_y_Xq[2] : 0] int<lower=0> y2_z1q_id_eps;
    array[assoc_uses[2,3] == 1 && bK1_len[3] > 0 ? nrow_y_Xq[3] : 0] int<lower=0> y3_z1q_id_eps;

    // re design matrix at quadpoints, group factor 2
    array[bK2_len[1]] vector[assoc_uses[2,1] == 1 && bK2_len[1] > 0 ? nrow_y_Xq[1] : 0] y1_z2q_eps;
    array[bK2_len[2]] vector[assoc_uses[2,2] == 1 && bK2_len[2] > 0 ? nrow_y_Xq[2] : 0] y2_z2q_eps;
    array[bK2_len[3]] vector[assoc_uses[2,3] == 1 && bK2_len[3] > 0 ? nrow_y_Xq[3] : 0] y3_z2q_eps;
    array[assoc_uses[2,1] == 1 && bK2_len[1] > 0 ? nrow_y_Xq[1] : 0] int<lower=0> y1_z2q_id_eps;
    array[assoc_uses[2,2] == 1 && bK2_len[2] > 0 ? nrow_y_Xq[2] : 0] int<lower=0> y2_z2q_id_eps;
    array[assoc_uses[2,3] == 1 && bK2_len[3] > 0 ? nrow_y_Xq[3] : 0] int<lower=0> y3_z2q_id_eps;

  //---- data for calculating integral of eta in GK quadrature

    // num. of nodes for GK quadrature for area under marker trajectory
    int<lower=0> auc_qnodes;
    int<lower=0> nrow_y_Xq_auc; // num. rows in long. predictor matrix at auc quadpoints
    vector[sum(assoc_uses[3,]) > 0 ? nrow_y_Xq_auc : 0] auc_qwts;

    // fe design matrix at quadpoints
    matrix[assoc_uses[3,1] == 1 ? nrow_y_Xq_auc : 0, yK[1]] y1_xq_auc;
    matrix[assoc_uses[3,2] == 1 ? nrow_y_Xq_auc : 0, yK[2]] y2_xq_auc;
    matrix[assoc_uses[3,3] == 1 ? nrow_y_Xq_auc : 0, yK[3]] y3_xq_auc;
    
    // offset values at quadpoints
    vector[has_offset[1] && assoc_uses[3,1] == 1 ? nrow_y_Xq_auc : 0] y1_offset_auc;
    vector[has_offset[2] && assoc_uses[3,2] == 1 ? nrow_y_Xq_auc : 0] y2_offset_auc;
    vector[has_offset[3] && assoc_uses[3,3] == 1 ? nrow_y_Xq_auc : 0] y3_offset_auc;

    // re design matrix at quadpoints, group factor 1
    array[bK1_len[1]] vector[assoc_uses[3,1] == 1 && bK1_len[1] > 0 ? nrow_y_Xq_auc : 0] y1_z1q_auc;
    array[bK1_len[2]] vector[assoc_uses[3,2] == 1 && bK1_len[2] > 0 ? nrow_y_Xq_auc : 0] y2_z1q_auc;
    array[bK1_len[3]] vector[assoc_uses[3,3] == 1 && bK1_len[3] > 0 ? nrow_y_Xq_auc : 0] y3_z1q_auc;
    array[assoc_uses[3,1] == 1 && bK1_len[1] > 0 ? nrow_y_Xq_auc : 0] int<lower=0> y1_z1q_id_auc;
    array[assoc_uses[3,2] == 1 && bK1_len[2] > 0 ? nrow_y_Xq_auc : 0] int<lower=0> y2_z1q_id_auc;
    array[assoc_uses[3,3] == 1 && bK1_len[3] > 0 ? nrow_y_Xq_auc : 0] int<lower=0> y3_z1q_id_auc;

    // re design matrix at quadpoints, group factor 2
    array[bK2_len[1]] vector[assoc_uses[3,1] == 1 && bK2_len[1] > 0 ? nrow_y_Xq_auc : 0] y1_z2q_auc;
    array[bK2_len[2]] vector[assoc_uses[3,2] == 1 && bK2_len[2] > 0 ? nrow_y_Xq_auc : 0] y2_z2q_auc;
    array[bK2_len[3]] vector[assoc_uses[3,3] == 1 && bK2_len[3] > 0 ? nrow_y_Xq_auc : 0] y3_z2q_auc;
    array[assoc_uses[3,1] == 1 && bK2_len[1] > 0 ? nrow_y_Xq_auc : 0] int<lower=0> y1_z2q_id_auc;
    array[assoc_uses[3,2] == 1 && bK2_len[2] > 0 ? nrow_y_Xq_auc : 0] int<lower=0> y2_z2q_id_auc;
    array[assoc_uses[3,3] == 1 && bK2_len[3] > 0 ? nrow_y_Xq_auc : 0] int<lower=0> y3_z2q_id_auc;

  //---- data for calculating assoc*data interactions in GK quadrature

    // num assoc pars used in {ev/es/mv/ms}*data interactions
    array[M*4] int<lower=0,upper=a_K> a_K_data;

    // design matrix for interacting with ev/es/mv/ms at quadpoints
    matrix[sum(nrow_y_Xq[1:M]), sum(a_K_data)] y_Xq_data;

    // indexing specifying the rows of y_Xq_data that correspond to
    // each submodel
    array[3,2] int<lower=0> idx_q;

  //---- data for combining lower level units clustered within patients

    array[M] int<lower=0,upper=1> has_grp; // 1 = has clustering below patient level
    int<lower=0,upper=4> grp_assoc; // 1=sum, 2=mean, 3=min, 4=max
    array[nrow_e_Xq,2] int<lower=0> grp_idx;
