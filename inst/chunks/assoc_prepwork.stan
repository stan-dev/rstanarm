    # NB: the m loop take the linear predictors and/or expected
    # values evaluated in the previous step and, for submodel m, 
    # adds the intercept and stores the result in an element of an 
    # array rather than having all submodels in a single vector

    // Linear predictor at quadrature time
    if (assoc_uses[1] == 1) {
      if (K > 0) y_eta_q = y_Xq_eta * beta;
      else y_eta_q = rep_vector(0.0, sum(nrow_y_Xq));
      y_eta_q = y_eta_q + csr_matrix_times_vector(sum(nrow_y_Xq), q, w_Zq_eta, v_Zq_eta, u_Zq_eta, b);
      for (m in 1:M) {
        y_eta_q[idx_q[m,1]:idx_q[m,2]] = add_intercept(
          y_eta_q, m, idx_q, has_intercept, has_intercept_nob, 
          has_intercept_lob, has_intercept_upb, gamma_nob, 
          gamma_lob, gamma_upb, xbar, beta, KM);
      }
    }

    // Linear predictor at time plus epsilon
    if (assoc_uses[2] == 1) {
      if (K > 0) y_eta_q_eps = y_Xq_eps * beta;
      else y_eta_q_eps = rep_vector(0.0, sum(nrow_y_Xq));
      //if (has_offset == 1) y_eta_q_eps = y_eta_q_eps + y_offset; # how to handle offset?
      y_eta_q_eps = y_eta_q_eps + csr_matrix_times_vector(sum(nrow_y_Xq), q, w_Zq_eps, v_Zq_eps, u_Zq_eps, b);
      for (m in 1:M) {
        y_eta_q_eps[idx_q[m,1]:idx_q[m,2]] = add_intercept(
          y_eta_q_eps, m, idx_q, has_intercept, has_intercept_nob, 
          has_intercept_lob, has_intercept_upb, gamma_nob, 
          gamma_lob, gamma_upb, xbar, beta, KM);
      }
    }

    // Linear predictor at auc quadpoints
    if (assoc_uses[3] == 1) {
      if (K > 0) y_eta_q_auc = y_Xq_auc * beta;
      else y_eta_q_auc = rep_vector(0.0, sum(nrow_y_Xq_auc));
      //if (has_offset == 1) y_eta_q_auc = y_eta_q_auc + y_offset; # how to handle offset?
      y_eta_q_auc = y_eta_q_auc + csr_matrix_times_vector(sum(nrow_y_Xq_auc), q, w_Zq_auc, v_Zq_auc, u_Zq_auc, b);
      for (m in 1:M) {
        y_eta_q_auc[idx_qauc[m,1]:idx_qauc[m,2]] = add_intercept(
          y_eta_q_auc, m, idx_qauc, has_intercept, has_intercept_nob, 
          has_intercept_lob, has_intercept_upb, gamma_nob, 
          gamma_lob, gamma_upb, xbar, beta, KM);
      }
    }
