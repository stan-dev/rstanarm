    // !!! Be careful that indexing of has_assoc matches stan_jm.fit !!!

    // mark tracks indexing within a_beta vector, which is the
    // vector of association parameters
    int mark = 0;

    // mark2 tracks indexing within a_K_data vector, which is the
    // vector specifying the number of columns used for each possible
    // type of association term by data interaction
    int mark2 = 0;

    // mark3 tracks indexing within size_which_interactions vector
    int mark3 = 0;

    for (m in 1:M) {

      //----- etavalue and any interactions

      mark2 += 1;
      if (has_assoc[1,m]  == 1 || // etavalue
          has_assoc[9,m]  == 1 || // etavalue * data
          has_assoc[13,m] == 1 || // etavalue * etavalue
          has_assoc[14,m] == 1) { // etavalue * muvalue

        // declare and define eta at quadpoints for submodel m
#include /model/make_eta_tmp.stan

        // add etavalue and any interactions to event submodel eta
        if (has_assoc[1,m] == 1) { // etavalue
          vector[nrow_e_Xq] val;
          if (has_grp[m] == 0) { // no grouping factor clustered within patients
            val = eta_tmp;
          }
          else { // submodel has a grouping factor clustered within patients
            val = collapse_within_groups(eta_tmp, grp_idx, grp_assoc);
          }
          mark += 1;
          e_eta_q += a_beta[mark] * (val - a_xbar[mark]);
        }
        if (has_assoc[9,m] == 1) { // etavalue*data
          int J = a_K_data[mark2];
          int j_shift = (mark2 == 1) ? 0 : sum(a_K_data[1:(mark2-1)]);
          for (j in 1:J) {
            vector[nrow_e_Xq] val;
            int sel = j_shift + j;
            if (has_grp[m] == 0) {
              val = eta_tmp .* y_Xq_data[idx_q[m,1]:idx_q[m,2], sel];
            }
            else {
              val = collapse_within_groups(
                eta_tmp .* y_Xq_data[idx_q[m,1]:idx_q[m,2], sel],
                grp_idx, grp_assoc);
            }
            mark += 1;
            e_eta_q += a_beta[mark] * (val - a_xbar[mark]);
          }
        }
        mark3 += 1; // count even if assoc type isn't used
        if (has_assoc[13,m] == 1) { // etavalue*etavalue
          for (j in 1:size_which_interactions[mark3]) {
            int j_shift = (mark3 == 1) ? 0 : sum(size_which_interactions[1:(mark3-1)]);
            int sel = which_interactions[j+j_shift];
            vector[nrow_e_Xq] val;
#include /model/make_eta_tmp2.stan
            val = eta_tmp .* eta_tmp2;
            mark += 1;
            e_eta_q += a_beta[mark] * (val - a_xbar[mark]);
          }
        }
        mark3 += 1; // count even if assoc type isn't used
        if (has_assoc[14,m] == 1) { // etavalue*muvalue
          for (j in 1:size_which_interactions[mark3]) {
            int j_shift = (mark3 == 1) ? 0 : sum(size_which_interactions[1:(mark3-1)]);
            int sel = which_interactions[j+j_shift];
            vector[nrow_e_Xq] val;
            vector[nrow_y_Xq[sel]] mu_tmp2;
#include /model/make_eta_tmp2.stan
            mu_tmp2 = evaluate_mu(eta_tmp2, family[sel], link[sel]);
            val = eta_tmp .* mu_tmp2;
            mark += 1;
            e_eta_q += a_beta[mark] * (val - a_xbar[mark]);
          }
        }
      }
      else {
        mark3 += 2;
      }

      //----- etaslope and any interactions

      mark2 += 1;
      if ((has_assoc[2,m] == 1) || (has_assoc[10,m] == 1)) {

        // declare and define etaslope at quadpoints for submodel m
        vector[nrow_y_Xq[m]] dydt_eta_q;
        if (m == 1) {
          int bMat1_colshift = 0;
          int bMat2_colshift = 0;
          dydt_eta_q = evaluate_eta(y1_xq_eps, y1_z1q_eps, y1_z2q_eps,
                                    y1_z1q_id_eps, y1_z2q_id_eps,
                                    yGamma1, yBeta1, bMat1, bMat2,
                                    bMat1_colshift, bMat2_colshift, 0);
        }
        else if (m == 2) {
          int bMat1_colshift = bK1_len[1];
          int bMat2_colshift = bK2_len[1];
          dydt_eta_q = evaluate_eta(y2_xq_eps, y2_z1q_eps, y2_z2q_eps,
                                    y2_z1q_id_eps, y2_z2q_id_eps,
                                    yGamma2, yBeta2, bMat1, bMat2,
                                    bMat1_colshift, bMat2_colshift, 0);
        }
        else if (m == 3) {
          int bMat1_colshift = sum(bK1_len[1:2]);
          int bMat2_colshift = sum(bK2_len[1:2]);
          dydt_eta_q = evaluate_eta(y3_xq_eps, y3_z1q_eps, y3_z2q_eps,
                                    y3_z1q_id_eps, y3_z2q_id_eps,
                                    yGamma3, yBeta3, bMat1, bMat2,
                                    bMat1_colshift, bMat2_colshift, 0);
        }

        // add etaslope and any interactions to event submodel eta
        if (has_assoc[2,m] == 1) { // etaslope
          vector[nrow_e_Xq] val;
          if (has_grp[m] == 0) {
            val = dydt_eta_q;
          }
          else {
            val = collapse_within_groups(dydt_eta_q, grp_idx, grp_assoc);
          }
          mark += 1;
          e_eta_q += a_beta[mark] * (val - a_xbar[mark]);
        }
        if (has_assoc[10,m] == 1) { // etaslope*data
          int J = a_K_data[mark2];
          int j_shift = (mark2 == 1) ? 0 : sum(a_K_data[1:(mark2-1)]);
          for (j in 1:J) {
            vector[nrow_e_Xq] val;
            int sel = j_shift + j;
            if (has_grp[m] == 0) {
              val = dydt_eta_q .* y_Xq_data[idx_q[m,1]:idx_q[m,2], sel];
            }
            else {
              val = collapse_within_groups(
                      dydt_eta_q .* y_Xq_data[idx_q[m,1]:idx_q[m,2], sel],
                      grp_idx, grp_assoc);
            }
            mark += 1;
            e_eta_q += a_beta[mark] * (val - a_xbar[mark]);
          }
        }
      }

      //----- etaauc

      // add etaauc to event submodel eta
      if (has_assoc[3,m] == 1) { // etaauc
        vector[nrow_y_Xq_auc] eta_auc_tmp; // eta at all auc quadpoints (for submodel m)
        vector[nrow_y_Xq[m]] val; // eta following summation over auc quadpoints
        if (m == 1) {
          int bMat1_colshift = 0;
          int bMat2_colshift = 0;
          eta_auc_tmp = evaluate_eta(y1_xq_auc, y1_z1q_auc, y1_z2q_auc,
                                     y1_z1q_id_auc, y1_z2q_id_auc,
                                     yGamma1, yBeta1, bMat1, bMat2,
                                     bMat1_colshift, bMat2_colshift,
                                     intercept_type[1]);
        }
        else if (m == 2) {
          int bMat1_colshift = bK1_len[1];
          int bMat2_colshift = bK2_len[1];
          eta_auc_tmp = evaluate_eta(y2_xq_auc, y2_z1q_auc, y2_z2q_auc,
                                     y2_z1q_id_auc, y2_z2q_id_auc,
                                     yGamma2, yBeta2, bMat1, bMat2,
                                     bMat1_colshift, bMat2_colshift,
                                     intercept_type[2]);
        }
        else if (m == 3) {
          int bMat1_colshift = sum(bK1_len[1:2]);
          int bMat2_colshift = sum(bK2_len[1:2]);
          eta_auc_tmp = evaluate_eta(y3_xq_auc, y3_z1q_auc, y3_z2q_auc,
                                     y3_z1q_id_auc, y3_z2q_id_auc,
                                     yGamma3, yBeta3, bMat1, bMat2,
                                     bMat1_colshift, bMat2_colshift,
                                     intercept_type[3]);
        }
        mark += 1;
        for (r in 1:nrow_y_Xq[m]) {
          vector[auc_qnodes] val_tmp;
          vector[auc_qnodes] wgt_tmp;
          val_tmp = eta_auc_tmp[((r-1) * auc_qnodes + 1):(r * auc_qnodes)];
          wgt_tmp = auc_qwts[((r-1) * auc_qnodes + 1):(r * auc_qnodes)];
          val[r] = sum(wgt_tmp .* val_tmp);
        }
        e_eta_q += a_beta[mark] * (val - a_xbar[mark]);
      }

      //----- muvalue and any interactions

      mark2 += 1;
      if (has_assoc[4,m]  == 1 || // muvalue
          has_assoc[11,m] == 1 || // muvalue * data
          has_assoc[15,m] == 1 || // muvalue * etavalue
          has_assoc[16,m] == 1) { // muvalue * muvalue

        // declare and define mu for submodel m
        vector[nrow_y_Xq[m]] mu_tmp;
#include /model/make_eta_tmp.stan
        mu_tmp = evaluate_mu(eta_tmp, family[m], link[m]);

        // add muvalue and any interactions to event submodel eta
        if (has_assoc[4,m] == 1) { // muvalue
          vector[nrow_e_Xq] val;
          if (has_grp[m] == 0) {
            val = mu_tmp;
          }
          else {
            val = collapse_within_groups(mu_tmp, grp_idx, grp_assoc);
          }
          mark += 1;
          e_eta_q += a_beta[mark] * (val - a_xbar[mark]);
        }
        if (has_assoc[11,m] == 1) { // muvalue*data
          int tmp = a_K_data[mark2];
          int j_shift = (mark2 == 1) ? 0 : sum(a_K_data[1:(mark2-1)]);
          for (j in 1:tmp) {
            vector[nrow_e_Xq] val;
            int sel = j_shift + j;
            if (has_grp[m] == 0) {
              val = mu_tmp .* y_Xq_data[idx_q[m,1]:idx_q[m,2], sel];
            }
            else {
              val = collapse_within_groups(
                mu_tmp .* y_Xq_data[idx_q[m,1]:idx_q[m,2], sel],
                grp_idx, grp_assoc);
            }
            mark += 1;
            e_eta_q += a_beta[mark] * (val - a_xbar[mark]);
          }
        }
        mark3 += 1; // count even if assoc type isn't used
        if (has_assoc[15,m] == 1) { // muvalue*etavalue
          for (j in 1:size_which_interactions[mark3]) {
            int j_shift = (mark3 == 1) ? 0 : sum(size_which_interactions[1:(mark3-1)]);
            int sel = which_interactions[j+j_shift];
            vector[nrow_e_Xq] val;
#include /model/make_eta_tmp2.stan
            val = mu_tmp .* eta_tmp2;
            mark += 1;
            e_eta_q += a_beta[mark] * (val - a_xbar[mark]);
          }
        }
        mark3 += 1; // count even if assoc type isn't used
        if (has_assoc[16,m] == 1) { // muvalue*muvalue
          for (j in 1:size_which_interactions[mark3]) {
            int j_shift = (mark3 == 1) ? 0 : sum(size_which_interactions[1:(mark3-1)]);
            int sel = which_interactions[j+j_shift];
            vector[nrow_e_Xq] val;
            vector[nrow_y_Xq[sel]] mu_tmp2;
#include /model/make_eta_tmp2.stan
            mu_tmp2 = evaluate_mu(eta_tmp2, family[sel], link[sel]);
            val = mu_tmp .* mu_tmp2;
            mark += 1;
            e_eta_q += a_beta[mark] * (val - a_xbar[mark]);
          }
        }
      }
      else {
        mark3 += 2;
      }

      //----- muslope and any interactions

      mark2 += 1;
      if (has_assoc[5,m] == 1 || has_assoc[12,m] == 1) {
        reject("muslope association structure has been removed.");
      }

      //----- muauc

      // add muauc to event submodel eta
      if (has_assoc[6,m] == 1) { // muauc
        vector[nrow_y_Xq_auc] eta_auc_tmp; // eta at all auc quadpoints (for submodel m)
        vector[nrow_y_Xq_auc] mu_auc_tmp; // mu at all auc quadpoints (for submodel m)
        vector[nrow_y_Xq[m]] val; // mu following summation over auc quadpoints
        if (m == 1) {
          int bMat1_colshift = 0;
          int bMat2_colshift = 0;
          eta_auc_tmp = evaluate_eta(y1_xq_auc, y1_z1q_auc, y1_z2q_auc,
                                     y1_z1q_id_auc, y1_z2q_id_auc,
                                     yGamma1, yBeta1, bMat1, bMat2,
                                     bMat1_colshift, bMat2_colshift,
                                     intercept_type[1]);
        }
        else if (m == 2) {
          int bMat1_colshift = bK1_len[1];
          int bMat2_colshift = bK2_len[1];
          eta_auc_tmp = evaluate_eta(y2_xq_auc, y2_z1q_auc, y2_z2q_auc,
                                     y2_z1q_id_auc, y2_z2q_id_auc,
                                     yGamma2, yBeta2, bMat1, bMat2,
                                     bMat1_colshift, bMat2_colshift,
                                     intercept_type[2]);
        }
        else if (m == 3) {
          int bMat1_colshift = sum(bK1_len[1:2]);
          int bMat2_colshift = sum(bK2_len[1:2]);
          eta_auc_tmp = evaluate_eta(y3_xq_auc, y3_z1q_auc, y3_z2q_auc,
                                     y3_z1q_id_auc, y3_z2q_id_auc,
                                     yGamma3, yBeta3, bMat1, bMat2,
                                     bMat1_colshift, bMat2_colshift,
                                     intercept_type[3]);
        }
        mu_auc_tmp = evaluate_mu(eta_auc_tmp, family[m], link[m]);
        mark += 1;
        for (r in 1:nrow_y_Xq[m]) {
          vector[auc_qnodes] val_tmp;
          vector[auc_qnodes] wgt_tmp;
          val_tmp = mu_auc_tmp[((r-1) * auc_qnodes + 1):(r * auc_qnodes)];
          wgt_tmp = auc_qwts[((r-1) * auc_qnodes + 1):(r * auc_qnodes)];
          val[r] = sum(wgt_tmp .* val_tmp);
        }
        e_eta_q += a_beta[mark] * (val - a_xbar[mark]);
      }

    }

    //-----  shared random effects

    if (sum_size_which_b > 0) {
      reject("shared_b has been removed.");
    }
    if (sum_size_which_coef > 0) {
      reject("shared_coef has been removed.");
    }
