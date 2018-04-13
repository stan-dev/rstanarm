  // prior family: 0 = none, 1 = normal, 2 = student_t, 3 = hs, 4 = hs_plus,
  //   5 = laplace, 6 = lasso
  int<lower=0,upper=6> a_prior_dist;

  //--- dimensions for association structure

    // num. of association parameters
    int<lower=0> a_K;

    // used for centering assoc terms
    vector[a_K] a_xbar;

    // 0 = no assoc structure, 1 = any assoc structure
    int<lower=0,upper=1> assoc;

    // which components are required to build association terms
    int<lower=0,upper=1> assoc_uses[6,3];

    // which association terms does each submodel use
    int<lower=0,upper=1> has_assoc[16,M];

    // num. of shared random effects
    int<lower=0> sum_size_which_b;

    // num. of shared random effects for each long submodel
    int<lower=0> size_which_b[M];

    // which random effects are shared for each long submodel
    int<lower=1> which_b_zindex[sum_size_which_b];

    // num. of shared random effects incl fixed component
    int<lower=0> sum_size_which_coef;

    // num. of shared random effects incl fixed component for each long submodel
    int<lower=0> size_which_coef[M];

    // which random effects are shared incl fixed component
    int<lower=1> which_coef_zindex[sum_size_which_coef];

    // which fixed effects are shared
    int<lower=1> which_coef_xindex[sum_size_which_coef];

    // total num pars used in assoc*assoc interactions
    int<lower=0> sum_size_which_interactions;

    // num pars used in assoc*assoc interactions, by submodel
    //   and by evev/evmv/mvev/mvmv interactions
    int<lower=0,upper=sum_size_which_interactions> size_which_interactions[M*4];

    // which terms to interact with
    int<lower=1> which_interactions[sum_size_which_interactions];

  //---- data for calculating eta in GK quadrature

    int<lower=0> nrow_y_Xq[3];     // num. rows in long. predictor matrix at quadpoints

    // fe design matrix at quadpoints
    matrix[assoc_uses[1,1] == 1 ? nrow_y_Xq[1] : 0, yK[1]] y1_xq_eta;
    matrix[assoc_uses[1,2] == 1 ? nrow_y_Xq[2] : 0, yK[2]] y2_xq_eta;
    matrix[assoc_uses[1,3] == 1 ? nrow_y_Xq[3] : 0, yK[3]] y3_xq_eta;

    // re design matrix at quadpoints, group factor 1
    vector[assoc_uses[1,1] == 1 && bK1_len[1] > 0 ? nrow_y_Xq[1] : 0] y1_z1q_eta[bK1_len[1]];
    vector[assoc_uses[1,2] == 1 && bK1_len[2] > 0 ? nrow_y_Xq[2] : 0] y2_z1q_eta[bK1_len[2]];
    vector[assoc_uses[1,3] == 1 && bK1_len[3] > 0 ? nrow_y_Xq[3] : 0] y3_z1q_eta[bK1_len[3]];
    int<lower=0> y1_z1q_id_eta[assoc_uses[1,1] == 1 && bK1_len[1] > 0 ? nrow_y_Xq[1] : 0];
    int<lower=0> y2_z1q_id_eta[assoc_uses[1,2] == 1 && bK1_len[2] > 0 ? nrow_y_Xq[2] : 0];
    int<lower=0> y3_z1q_id_eta[assoc_uses[1,3] == 1 && bK1_len[3] > 0 ? nrow_y_Xq[3] : 0];

    // re design matrix at quadpoints, group factor 2
    vector[assoc_uses[1,1] == 1 && bK2_len[1] > 0 ? nrow_y_Xq[1] : 0] y1_z2q_eta[bK2_len[1]];
    vector[assoc_uses[1,2] == 1 && bK2_len[2] > 0 ? nrow_y_Xq[2] : 0] y2_z2q_eta[bK2_len[2]];
    vector[assoc_uses[1,3] == 1 && bK2_len[3] > 0 ? nrow_y_Xq[3] : 0] y3_z2q_eta[bK2_len[3]];
    int<lower=0> y1_z2q_id_eta[assoc_uses[1,1] == 1 && bK2_len[1] > 0 ? nrow_y_Xq[1] : 0];
    int<lower=0> y2_z2q_id_eta[assoc_uses[1,2] == 1 && bK2_len[2] > 0 ? nrow_y_Xq[2] : 0];
    int<lower=0> y3_z2q_id_eta[assoc_uses[1,3] == 1 && bK2_len[3] > 0 ? nrow_y_Xq[3] : 0];

  //---- data for calculating derivative of eta in GK quadrature

    // fe design matrix at quadpoints
    matrix[assoc_uses[2,1] == 1 ? nrow_y_Xq[1] : 0, yK[1]] y1_xq_eps;
    matrix[assoc_uses[2,2] == 1 ? nrow_y_Xq[2] : 0, yK[2]] y2_xq_eps;
    matrix[assoc_uses[2,3] == 1 ? nrow_y_Xq[3] : 0, yK[3]] y3_xq_eps;

    // re design matrix at quadpoints, group factor 1
    vector[assoc_uses[2,1] == 1 && bK1_len[1] > 0 ? nrow_y_Xq[1] : 0] y1_z1q_eps[bK1_len[1]];
    vector[assoc_uses[2,2] == 1 && bK1_len[2] > 0 ? nrow_y_Xq[2] : 0] y2_z1q_eps[bK1_len[2]];
    vector[assoc_uses[2,3] == 1 && bK1_len[3] > 0 ? nrow_y_Xq[3] : 0] y3_z1q_eps[bK1_len[3]];
    int<lower=0> y1_z1q_id_eps[assoc_uses[2,1] == 1 && bK1_len[1] > 0 ? nrow_y_Xq[1] : 0];
    int<lower=0> y2_z1q_id_eps[assoc_uses[2,2] == 1 && bK1_len[2] > 0 ? nrow_y_Xq[2] : 0];
    int<lower=0> y3_z1q_id_eps[assoc_uses[2,3] == 1 && bK1_len[3] > 0 ? nrow_y_Xq[3] : 0];

    // re design matrix at quadpoints, group factor 2
    vector[assoc_uses[2,1] == 1 && bK2_len[1] > 0 ? nrow_y_Xq[1] : 0] y1_z2q_eps[bK2_len[1]];
    vector[assoc_uses[2,2] == 1 && bK2_len[2] > 0 ? nrow_y_Xq[2] : 0] y2_z2q_eps[bK2_len[2]];
    vector[assoc_uses[2,3] == 1 && bK2_len[3] > 0 ? nrow_y_Xq[3] : 0] y3_z2q_eps[bK2_len[3]];
    int<lower=0> y1_z2q_id_eps[assoc_uses[2,1] == 1 && bK2_len[1] > 0 ? nrow_y_Xq[1] : 0];
    int<lower=0> y2_z2q_id_eps[assoc_uses[2,2] == 1 && bK2_len[2] > 0 ? nrow_y_Xq[2] : 0];
    int<lower=0> y3_z2q_id_eps[assoc_uses[2,3] == 1 && bK2_len[3] > 0 ? nrow_y_Xq[3] : 0];

  //---- data for calculating integral of eta in GK quadrature

    // num. of nodes for GK quadrature for area under marker trajectory
    int<lower=0> auc_qnodes;
    int<lower=0> nrow_y_Xq_auc; // num. rows in long. predictor matrix at auc quadpoints
    vector[sum(assoc_uses[3,]) > 0 ? nrow_y_Xq_auc : 0] auc_qwts;

    // fe design matrix at quadpoints
    matrix[assoc_uses[3,1] == 1 ? nrow_y_Xq_auc : 0, yK[1]] y1_xq_auc;
    matrix[assoc_uses[3,2] == 1 ? nrow_y_Xq_auc : 0, yK[2]] y2_xq_auc;
    matrix[assoc_uses[3,3] == 1 ? nrow_y_Xq_auc : 0, yK[3]] y3_xq_auc;

    // re design matrix at quadpoints, group factor 1
    vector[assoc_uses[3,1] == 1 && bK1_len[1] > 0 ? nrow_y_Xq_auc : 0] y1_z1q_auc[bK1_len[1]];
    vector[assoc_uses[3,2] == 1 && bK1_len[2] > 0 ? nrow_y_Xq_auc : 0] y2_z1q_auc[bK1_len[2]];
    vector[assoc_uses[3,3] == 1 && bK1_len[3] > 0 ? nrow_y_Xq_auc : 0] y3_z1q_auc[bK1_len[3]];
    int<lower=0> y1_z1q_id_auc[assoc_uses[3,1] == 1 && bK1_len[1] > 0 ? nrow_y_Xq_auc : 0];
    int<lower=0> y2_z1q_id_auc[assoc_uses[3,2] == 1 && bK1_len[2] > 0 ? nrow_y_Xq_auc : 0];
    int<lower=0> y3_z1q_id_auc[assoc_uses[3,3] == 1 && bK1_len[3] > 0 ? nrow_y_Xq_auc : 0];

    // re design matrix at quadpoints, group factor 2
    vector[assoc_uses[3,1] == 1 && bK2_len[1] > 0 ? nrow_y_Xq_auc : 0] y1_z2q_auc[bK2_len[1]];
    vector[assoc_uses[3,2] == 1 && bK2_len[2] > 0 ? nrow_y_Xq_auc : 0] y2_z2q_auc[bK2_len[2]];
    vector[assoc_uses[3,3] == 1 && bK2_len[3] > 0 ? nrow_y_Xq_auc : 0] y3_z2q_auc[bK2_len[3]];
    int<lower=0> y1_z2q_id_auc[assoc_uses[3,1] == 1 && bK2_len[1] > 0 ? nrow_y_Xq_auc : 0];
    int<lower=0> y2_z2q_id_auc[assoc_uses[3,2] == 1 && bK2_len[2] > 0 ? nrow_y_Xq_auc : 0];
    int<lower=0> y3_z2q_id_auc[assoc_uses[3,3] == 1 && bK2_len[3] > 0 ? nrow_y_Xq_auc : 0];

  //---- data for calculating assoc*data interactions in GK quadrature

    // num assoc pars used in {ev/es/mv/ms}*data interactions
    int<lower=0,upper=a_K> a_K_data[M*4];

    // design matrix for interacting with ev/es/mv/ms at quadpoints
    matrix[sum(nrow_y_Xq[1:M]), sum(a_K_data)] y_Xq_data;

    // indexing specifying the rows of y_Xq_data that correspond to
    // each submodel
    int<lower=0> idx_q[3,2];

  //---- data for combining lower level units clustered within patients

    int<lower=0,upper=1> has_grp[M]; // 1 = has clustering below patient level
    int<lower=0,upper=4> grp_assoc; // 1=sum, 2=mean, 3=min, 4=max
    int<lower=0> grp_idx[nrow_e_Xq,2];
