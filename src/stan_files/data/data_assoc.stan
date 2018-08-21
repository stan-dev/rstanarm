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

    // fe design matrices

    matrix[assoc_uses[1,1] == 1 ? len_epts[1] : 0, yK[1]] y1_x_eta_epts;
    matrix[assoc_uses[1,2] == 1 ? len_epts[2] : 0, yK[2]] y2_x_eta_epts;
    matrix[assoc_uses[1,3] == 1 ? len_epts[3] : 0, yK[3]] y3_x_eta_epts;

    matrix[assoc_uses[1,1] == 1 ? len_qpts[1] : 0, yK[1]] y1_x_eta_qpts;
    matrix[assoc_uses[1,2] == 1 ? len_qpts[2] : 0, yK[2]] y2_x_eta_qpts;
    matrix[assoc_uses[1,3] == 1 ? len_qpts[3] : 0, yK[3]] y3_x_eta_qpts;

    matrix[assoc_uses[1,1] == 1 ? len_ipts[1] : 0, yK[1]] y1_x_eta_ipts;
    matrix[assoc_uses[1,2] == 1 ? len_ipts[2] : 0, yK[2]] y2_x_eta_ipts;
    matrix[assoc_uses[1,3] == 1 ? len_ipts[3] : 0, yK[3]] y3_x_eta_ipts;

    // re design matrices, group factor 1

    vector[assoc_uses[1,1] == 1 && bK1_len[1] > 0 ? len_epts[1] : 0] y1_z1_eta_epts[bK1_len[1]];
    vector[assoc_uses[1,2] == 1 && bK1_len[2] > 0 ? len_epts[2] : 0] y2_z1_eta_epts[bK1_len[2]];
    vector[assoc_uses[1,3] == 1 && bK1_len[3] > 0 ? len_epts[3] : 0] y3_z1_eta_epts[bK1_len[3]];

    vector[assoc_uses[1,1] == 1 && bK1_len[1] > 0 ? len_qpts[1] : 0] y1_z1_eta_qpts[bK1_len[1]];
    vector[assoc_uses[1,2] == 1 && bK1_len[2] > 0 ? len_qpts[2] : 0] y2_z1_eta_qpts[bK1_len[2]];
    vector[assoc_uses[1,3] == 1 && bK1_len[3] > 0 ? len_qpts[3] : 0] y3_z1_eta_qpts[bK1_len[3]];

    vector[assoc_uses[1,1] == 1 && bK1_len[1] > 0 ? len_ipts[1] : 0] y1_z1_eta_ipts[bK1_len[1]];
    vector[assoc_uses[1,2] == 1 && bK1_len[2] > 0 ? len_ipts[2] : 0] y2_z1_eta_ipts[bK1_len[2]];
    vector[assoc_uses[1,3] == 1 && bK1_len[3] > 0 ? len_ipts[3] : 0] y3_z1_eta_ipts[bK1_len[3]];

    // re design matrices, group factor 2

    vector[assoc_uses[1,1] == 1 && bK2_len[1] > 0 ? len_epts[1] : 0] y1_z2_eta_epts[bK2_len[1]];
    vector[assoc_uses[1,2] == 1 && bK2_len[2] > 0 ? len_epts[2] : 0] y2_z2_eta_epts[bK2_len[2]];
    vector[assoc_uses[1,3] == 1 && bK2_len[3] > 0 ? len_epts[3] : 0] y3_z2_eta_epts[bK2_len[3]];

    vector[assoc_uses[1,1] == 1 && bK2_len[1] > 0 ? len_qpts[1] : 0] y1_z2_eta_qpts[bK2_len[1]];
    vector[assoc_uses[1,2] == 1 && bK2_len[2] > 0 ? len_qpts[2] : 0] y2_z2_eta_qpts[bK2_len[2]];
    vector[assoc_uses[1,3] == 1 && bK2_len[3] > 0 ? len_qpts[3] : 0] y3_z2_eta_qpts[bK2_len[3]];

    vector[assoc_uses[1,1] == 1 && bK2_len[1] > 0 ? len_ipts[1] : 0] y1_z2_eta_ipts[bK2_len[1]];
    vector[assoc_uses[1,2] == 1 && bK2_len[2] > 0 ? len_ipts[2] : 0] y2_z2_eta_ipts[bK2_len[2]];
    vector[assoc_uses[1,3] == 1 && bK2_len[3] > 0 ? len_ipts[3] : 0] y3_z2_eta_ipts[bK2_len[3]];

    // ids for re design matrices, group factor 1

    int<lower=0> y1_z1_id_eta_epts[assoc_uses[1,1] == 1 && bK1_len[1] > 0 ? len_epts[1] : 0];
    int<lower=0> y2_z1_id_eta_epts[assoc_uses[1,2] == 1 && bK1_len[2] > 0 ? len_epts[2] : 0];
    int<lower=0> y3_z1_id_eta_epts[assoc_uses[1,3] == 1 && bK1_len[3] > 0 ? len_epts[3] : 0];

    int<lower=0> y1_z1_id_eta_qpts[assoc_uses[1,1] == 1 && bK1_len[1] > 0 ? len_qpts[1] : 0];
    int<lower=0> y2_z1_id_eta_qpts[assoc_uses[1,2] == 1 && bK1_len[2] > 0 ? len_qpts[2] : 0];
    int<lower=0> y3_z1_id_eta_qpts[assoc_uses[1,3] == 1 && bK1_len[3] > 0 ? len_qpts[3] : 0];

    int<lower=0> y1_z1_id_eta_ipts[assoc_uses[1,1] == 1 && bK1_len[1] > 0 ? len_ipts[1] : 0];
    int<lower=0> y2_z1_id_eta_ipts[assoc_uses[1,2] == 1 && bK1_len[2] > 0 ? len_ipts[2] : 0];
    int<lower=0> y3_z1_id_eta_ipts[assoc_uses[1,3] == 1 && bK1_len[3] > 0 ? len_ipts[3] : 0];

    // ids for re design matrices, group factor 1

    int<lower=0> y1_z2_id_eta_epts[assoc_uses[1,1] == 1 && bK2_len[1] > 0 ? len_epts[1] : 0];
    int<lower=0> y2_z2_id_eta_epts[assoc_uses[1,2] == 1 && bK2_len[2] > 0 ? len_epts[2] : 0];
    int<lower=0> y3_z2_id_eta_epts[assoc_uses[1,3] == 1 && bK2_len[3] > 0 ? len_epts[3] : 0];

    int<lower=0> y1_z2_id_eta_qpts[assoc_uses[1,1] == 1 && bK2_len[1] > 0 ? len_qpts[1] : 0];
    int<lower=0> y2_z2_id_eta_qpts[assoc_uses[1,2] == 1 && bK2_len[2] > 0 ? len_qpts[2] : 0];
    int<lower=0> y3_z2_id_eta_qpts[assoc_uses[1,3] == 1 && bK2_len[3] > 0 ? len_qpts[3] : 0];

    int<lower=0> y1_z2_id_eta_ipts[assoc_uses[1,1] == 1 && bK2_len[1] > 0 ? len_ipts[1] : 0];
    int<lower=0> y2_z2_id_eta_ipts[assoc_uses[1,2] == 1 && bK2_len[2] > 0 ? len_ipts[2] : 0];
    int<lower=0> y3_z2_id_eta_ipts[assoc_uses[1,3] == 1 && bK2_len[3] > 0 ? len_ipts[3] : 0];


  //---- data for calculating derivative of eta in GK quadrature

    // fe design matrices

    matrix[assoc_uses[2,1] == 1 ? len_epts[1] : 0, yK[1]] y1_x_eps_epts;
    matrix[assoc_uses[2,2] == 1 ? len_epts[2] : 0, yK[2]] y2_x_eps_epts;
    matrix[assoc_uses[2,3] == 1 ? len_epts[3] : 0, yK[3]] y3_x_eps_epts;

    matrix[assoc_uses[2,1] == 1 ? len_qpts[1] : 0, yK[1]] y1_x_eps_qpts;
    matrix[assoc_uses[2,2] == 1 ? len_qpts[2] : 0, yK[2]] y2_x_eps_qpts;
    matrix[assoc_uses[2,3] == 1 ? len_qpts[3] : 0, yK[3]] y3_x_eps_qpts;

    matrix[assoc_uses[2,1] == 1 ? len_ipts[1] : 0, yK[1]] y1_x_eps_ipts;
    matrix[assoc_uses[2,2] == 1 ? len_ipts[2] : 0, yK[2]] y2_x_eps_ipts;
    matrix[assoc_uses[2,3] == 1 ? len_ipts[3] : 0, yK[3]] y3_x_eps_ipts;

    // re design matrices, group factor 1

    vector[assoc_uses[2,1] == 1 && bK1_len[1] > 0 ? len_epts[1] : 0] y1_z1_eps_epts[bK1_len[1]];
    vector[assoc_uses[2,2] == 1 && bK1_len[2] > 0 ? len_epts[2] : 0] y2_z1_eps_epts[bK1_len[2]];
    vector[assoc_uses[2,3] == 1 && bK1_len[3] > 0 ? len_epts[3] : 0] y3_z1_eps_epts[bK1_len[3]];

    vector[assoc_uses[2,1] == 1 && bK1_len[1] > 0 ? len_qpts[1] : 0] y1_z1_eps_qpts[bK1_len[1]];
    vector[assoc_uses[2,2] == 1 && bK1_len[2] > 0 ? len_qpts[2] : 0] y2_z1_eps_qpts[bK1_len[2]];
    vector[assoc_uses[2,3] == 1 && bK1_len[3] > 0 ? len_qpts[3] : 0] y3_z1_eps_qpts[bK1_len[3]];

    vector[assoc_uses[2,1] == 1 && bK1_len[1] > 0 ? len_ipts[1] : 0] y1_z1_eps_ipts[bK1_len[1]];
    vector[assoc_uses[2,2] == 1 && bK1_len[2] > 0 ? len_ipts[2] : 0] y2_z1_eps_ipts[bK1_len[2]];
    vector[assoc_uses[2,3] == 1 && bK1_len[3] > 0 ? len_ipts[3] : 0] y3_z1_eps_ipts[bK1_len[3]];

    // re design matrices, group factor 2

    vector[assoc_uses[2,1] == 1 && bK2_len[1] > 0 ? len_epts[1] : 0] y1_z2_eps_epts[bK2_len[1]];
    vector[assoc_uses[2,2] == 1 && bK2_len[2] > 0 ? len_epts[2] : 0] y2_z2_eps_epts[bK2_len[2]];
    vector[assoc_uses[2,3] == 1 && bK2_len[3] > 0 ? len_epts[3] : 0] y3_z2_eps_epts[bK2_len[3]];

    vector[assoc_uses[2,1] == 1 && bK2_len[1] > 0 ? len_qpts[1] : 0] y1_z2_eps_qpts[bK2_len[1]];
    vector[assoc_uses[2,2] == 1 && bK2_len[2] > 0 ? len_qpts[2] : 0] y2_z2_eps_qpts[bK2_len[2]];
    vector[assoc_uses[2,3] == 1 && bK2_len[3] > 0 ? len_qpts[3] : 0] y3_z2_eps_qpts[bK2_len[3]];

    vector[assoc_uses[2,1] == 1 && bK2_len[1] > 0 ? len_ipts[1] : 0] y1_z2_eps_ipts[bK2_len[1]];
    vector[assoc_uses[2,2] == 1 && bK2_len[2] > 0 ? len_ipts[2] : 0] y2_z2_eps_ipts[bK2_len[2]];
    vector[assoc_uses[2,3] == 1 && bK2_len[3] > 0 ? len_ipts[3] : 0] y3_z2_eps_ipts[bK2_len[3]];

    // ids for re design matrices, group factor 1

    int<lower=0> y1_z1_id_eps_epts[assoc_uses[2,1] == 1 && bK1_len[1] > 0 ? len_epts[1] : 0];
    int<lower=0> y2_z1_id_eps_epts[assoc_uses[2,2] == 1 && bK1_len[2] > 0 ? len_epts[2] : 0];
    int<lower=0> y3_z1_id_eps_epts[assoc_uses[2,3] == 1 && bK1_len[3] > 0 ? len_epts[3] : 0];

    int<lower=0> y1_z1_id_eps_qpts[assoc_uses[2,1] == 1 && bK1_len[1] > 0 ? len_qpts[1] : 0];
    int<lower=0> y2_z1_id_eps_qpts[assoc_uses[2,2] == 1 && bK1_len[2] > 0 ? len_qpts[2] : 0];
    int<lower=0> y3_z1_id_eps_qpts[assoc_uses[2,3] == 1 && bK1_len[3] > 0 ? len_qpts[3] : 0];

    int<lower=0> y1_z1_id_eps_ipts[assoc_uses[2,1] == 1 && bK1_len[1] > 0 ? len_ipts[1] : 0];
    int<lower=0> y2_z1_id_eps_ipts[assoc_uses[2,2] == 1 && bK1_len[2] > 0 ? len_ipts[2] : 0];
    int<lower=0> y3_z1_id_eps_ipts[assoc_uses[2,3] == 1 && bK1_len[3] > 0 ? len_ipts[3] : 0];

    // ids for re design matrices, group factor 1

    int<lower=0> y1_z2_id_eps_epts[assoc_uses[2,1] == 1 && bK2_len[1] > 0 ? len_epts[1] : 0];
    int<lower=0> y2_z2_id_eps_epts[assoc_uses[2,2] == 1 && bK2_len[2] > 0 ? len_epts[2] : 0];
    int<lower=0> y3_z2_id_eps_epts[assoc_uses[2,3] == 1 && bK2_len[3] > 0 ? len_epts[3] : 0];

    int<lower=0> y1_z2_id_eps_qpts[assoc_uses[2,1] == 1 && bK2_len[1] > 0 ? len_qpts[1] : 0];
    int<lower=0> y2_z2_id_eps_qpts[assoc_uses[2,2] == 1 && bK2_len[2] > 0 ? len_qpts[2] : 0];
    int<lower=0> y3_z2_id_eps_qpts[assoc_uses[2,3] == 1 && bK2_len[3] > 0 ? len_qpts[3] : 0];

    int<lower=0> y1_z2_id_eps_ipts[assoc_uses[2,1] == 1 && bK2_len[1] > 0 ? len_ipts[1] : 0];
    int<lower=0> y2_z2_id_eps_ipts[assoc_uses[2,2] == 1 && bK2_len[2] > 0 ? len_ipts[2] : 0];
    int<lower=0> y3_z2_id_eps_ipts[assoc_uses[2,3] == 1 && bK2_len[3] > 0 ? len_ipts[3] : 0];


  //---- data for calculating integral of eta in GK quadrature

    // num. of nodes for GK quadrature for area under marker trajectory
    int<lower=0> auc_qnodes;
    int<lower=0> nrow_y_Xq_auc; // num. rows in long. predictor matrix at auc quadpoints
    vector[sum(assoc_uses[3,]) > 0 ? nrow_y_Xq_auc : 0] auc_qwts;

    // fe design matrices

    matrix[assoc_uses[3,1] == 1 ? len_epts[1] : 0, yK[1]] y1_x_auc_epts;
    matrix[assoc_uses[3,2] == 1 ? len_epts[2] : 0, yK[2]] y2_x_auc_epts;
    matrix[assoc_uses[3,3] == 1 ? len_epts[3] : 0, yK[3]] y3_x_auc_epts;

    matrix[assoc_uses[3,1] == 1 ? len_qpts[1] : 0, yK[1]] y1_x_auc_qpts;
    matrix[assoc_uses[3,2] == 1 ? len_qpts[2] : 0, yK[2]] y2_x_auc_qpts;
    matrix[assoc_uses[3,3] == 1 ? len_qpts[3] : 0, yK[3]] y3_x_auc_qpts;

    matrix[assoc_uses[3,1] == 1 ? len_ipts[1] : 0, yK[1]] y1_x_auc_ipts;
    matrix[assoc_uses[3,2] == 1 ? len_ipts[2] : 0, yK[2]] y2_x_auc_ipts;
    matrix[assoc_uses[3,3] == 1 ? len_ipts[3] : 0, yK[3]] y3_x_auc_ipts;

    // re design matrices, group factor 1

    vector[assoc_uses[3,1] == 1 && bK1_len[1] > 0 ? len_epts[1] : 0] y1_z1_auc_epts[bK1_len[1]];
    vector[assoc_uses[3,2] == 1 && bK1_len[2] > 0 ? len_epts[2] : 0] y2_z1_auc_epts[bK1_len[2]];
    vector[assoc_uses[3,3] == 1 && bK1_len[3] > 0 ? len_epts[3] : 0] y3_z1_auc_epts[bK1_len[3]];

    vector[assoc_uses[3,1] == 1 && bK1_len[1] > 0 ? len_qpts[1] : 0] y1_z1_auc_qpts[bK1_len[1]];
    vector[assoc_uses[3,2] == 1 && bK1_len[2] > 0 ? len_qpts[2] : 0] y2_z1_auc_qpts[bK1_len[2]];
    vector[assoc_uses[3,3] == 1 && bK1_len[3] > 0 ? len_qpts[3] : 0] y3_z1_auc_qpts[bK1_len[3]];

    vector[assoc_uses[3,1] == 1 && bK1_len[1] > 0 ? len_ipts[1] : 0] y1_z1_auc_ipts[bK1_len[1]];
    vector[assoc_uses[3,2] == 1 && bK1_len[2] > 0 ? len_ipts[2] : 0] y2_z1_auc_ipts[bK1_len[2]];
    vector[assoc_uses[3,3] == 1 && bK1_len[3] > 0 ? len_ipts[3] : 0] y3_z1_auc_ipts[bK1_len[3]];

    // re design matrices, group factor 2

    vector[assoc_uses[3,1] == 1 && bK2_len[1] > 0 ? len_epts[1] : 0] y1_z2_auc_epts[bK2_len[1]];
    vector[assoc_uses[3,2] == 1 && bK2_len[2] > 0 ? len_epts[2] : 0] y2_z2_auc_epts[bK2_len[2]];
    vector[assoc_uses[3,3] == 1 && bK2_len[3] > 0 ? len_epts[3] : 0] y3_z2_auc_epts[bK2_len[3]];

    vector[assoc_uses[3,1] == 1 && bK2_len[1] > 0 ? len_qpts[1] : 0] y1_z2_auc_qpts[bK2_len[1]];
    vector[assoc_uses[3,2] == 1 && bK2_len[2] > 0 ? len_qpts[2] : 0] y2_z2_auc_qpts[bK2_len[2]];
    vector[assoc_uses[3,3] == 1 && bK2_len[3] > 0 ? len_qpts[3] : 0] y3_z2_auc_qpts[bK2_len[3]];

    vector[assoc_uses[3,1] == 1 && bK2_len[1] > 0 ? len_ipts[1] : 0] y1_z2_auc_ipts[bK2_len[1]];
    vector[assoc_uses[3,2] == 1 && bK2_len[2] > 0 ? len_ipts[2] : 0] y2_z2_auc_ipts[bK2_len[2]];
    vector[assoc_uses[3,3] == 1 && bK2_len[3] > 0 ? len_ipts[3] : 0] y3_z2_auc_ipts[bK2_len[3]];

    // ids for re design matrices, group factor 1

    int<lower=0> y1_z1_id_auc_epts[assoc_uses[3,1] == 1 && bK1_len[1] > 0 ? len_epts[1] : 0];
    int<lower=0> y2_z1_id_auc_epts[assoc_uses[3,2] == 1 && bK1_len[2] > 0 ? len_epts[2] : 0];
    int<lower=0> y3_z1_id_auc_epts[assoc_uses[3,3] == 1 && bK1_len[3] > 0 ? len_epts[3] : 0];

    int<lower=0> y1_z1_id_auc_qpts[assoc_uses[3,1] == 1 && bK1_len[1] > 0 ? len_qpts[1] : 0];
    int<lower=0> y2_z1_id_auc_qpts[assoc_uses[3,2] == 1 && bK1_len[2] > 0 ? len_qpts[2] : 0];
    int<lower=0> y3_z1_id_auc_qpts[assoc_uses[3,3] == 1 && bK1_len[3] > 0 ? len_qpts[3] : 0];

    int<lower=0> y1_z1_id_auc_ipts[assoc_uses[3,1] == 1 && bK1_len[1] > 0 ? len_ipts[1] : 0];
    int<lower=0> y2_z1_id_auc_ipts[assoc_uses[3,2] == 1 && bK1_len[2] > 0 ? len_ipts[2] : 0];
    int<lower=0> y3_z1_id_auc_ipts[assoc_uses[3,3] == 1 && bK1_len[3] > 0 ? len_ipts[3] : 0];

    // ids for re design matrices, group factor 1

    int<lower=0> y1_z2_id_auc_epts[assoc_uses[3,1] == 1 && bK2_len[1] > 0 ? len_epts[1] : 0];
    int<lower=0> y2_z2_id_auc_epts[assoc_uses[3,2] == 1 && bK2_len[2] > 0 ? len_epts[2] : 0];
    int<lower=0> y3_z2_id_auc_epts[assoc_uses[3,3] == 1 && bK2_len[3] > 0 ? len_epts[3] : 0];

    int<lower=0> y1_z2_id_auc_qpts[assoc_uses[3,1] == 1 && bK2_len[1] > 0 ? len_qpts[1] : 0];
    int<lower=0> y2_z2_id_auc_qpts[assoc_uses[3,2] == 1 && bK2_len[2] > 0 ? len_qpts[2] : 0];
    int<lower=0> y3_z2_id_auc_qpts[assoc_uses[3,3] == 1 && bK2_len[3] > 0 ? len_qpts[3] : 0];

    int<lower=0> y1_z2_id_auc_ipts[assoc_uses[3,1] == 1 && bK2_len[1] > 0 ? len_ipts[1] : 0];
    int<lower=0> y2_z2_id_auc_ipts[assoc_uses[3,2] == 1 && bK2_len[2] > 0 ? len_ipts[2] : 0];
    int<lower=0> y3_z2_id_auc_ipts[assoc_uses[3,3] == 1 && bK2_len[3] > 0 ? len_ipts[3] : 0];


  //---- data for calculating assoc*data interactions in GK quadrature

    // num assoc pars used in {ev/es/mv/ms}*data interactions
    int<lower=0,upper=a_K> a_K_data[M*4];

    // design matrix for interacting with ev/es/mv/ms at quadpoints
    matrix[sum(len_epts[1:M]), sum(a_K_data)] y_x_data_qpts;
    matrix[sum(len_qpts[1:M]), sum(a_K_data)] y_x_data_qpts;
    matrix[sum(len_qpts[1:M]), sum(a_K_data)] y_x_data_qpts;
    matrix[sum(len_qpts[1:M]), sum(a_K_data)] y_x_data_qpts;

    // indexing specifying the rows of y_Xq_data that correspond to each submodel
    int<lower=0> idx_q[3,2];

  //---- data for combining lower level units clustered within patients

    int<lower=0,upper=1> has_grp[M]; // 1 = has clustering below patient level
    int<lower=0,upper=4> grp_assoc;  // 1 = sum, 2 = mean, 3 = min, 4 = max
    int<lower=0> grp_idx_epts      [len_epts,2];
    int<lower=0> grp_idx_qpts      [len_qpts ,2];
    int<lower=0> grp_idx_ipts[len_ipts,2];
