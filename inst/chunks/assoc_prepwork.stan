    // Linear predictor at quadrature time (NB does not include intercept)
    if (assoc_uses[1] == 1) {
      if (K > 0) y_eta_q = y_Xq_eta * beta;
      else y_eta_q = rep_vector(0.0, sum(nrow_y_Xq));
      y_eta_q = y_eta_q + csr_matrix_times_vector(sum(nrow_y_Xq), q, w_Zq_eta, v_Zq_eta, u_Zq_eta, b);
    }

    // Linear predictor at time plus epsilon (NB does not include intercept)
    if (assoc_uses[2] == 1) {
      if (K > 0) y_eta_q_eps = y_Xq_eps * beta;
      else y_eta_q_eps = rep_vector(0.0, sum(nrow_y_Xq));
      //if (has_offset == 1) y_eta_q_eps = y_eta_q_eps + y_offset; // how to handle offset?
      y_eta_q_eps = y_eta_q_eps + csr_matrix_times_vector(sum(nrow_y_Xq), q, w_Zq_eps, v_Zq_eps, u_Zq_eps, b);
    }

    // Linear predictor at auc quadpoints (NB does not include intercept)
    if (assoc_uses[3] == 1) {
      if (K > 0) y_eta_q_auc = y_Xq_auc * beta;
      else y_eta_q_auc = rep_vector(0.0, sum(nrow_y_Xq_auc));
      //if (has_offset == 1) y_eta_q_auc = y_eta_q_auc + y_offset; // how to handle offset?
      y_eta_q_auc = y_eta_q_auc + csr_matrix_times_vector(sum(nrow_y_Xq_auc), q, w_Zq_auc, v_Zq_auc, u_Zq_auc, b);
    }
