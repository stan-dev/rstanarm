    # NB: the m loop take the linear predictors and/or expected
    # values evaluated in the previous step and, for submodel m, 
    # adds the intercept and stores the result in an element of an 
    # array rather than having all submodels in a single vector

    // Linear predictor at quadrature time
    if (assoc_uses[1] == 1) {
      if (K > 0) y_eta_q = y_Xq_eta * beta;
      else y_eta_q = rep_vector(0.0, (M*nrow_y_Xq));
      y_eta_q = y_eta_q + csr_matrix_times_vector((M*nrow_y_Xq), q, w_Zq_eta, v_Zq_eta, u_Zq_eta, b);
      for (m in 1:M) {
        y_eta_qwide[m] = add_intercept(
          y_eta_q, m, nrow_y_Xq, has_intercept, has_intercept_nob, 
          has_intercept_lob, has_intercept_upb, gamma_nob, 
          gamma_lob, gamma_upb, xbar, beta, KM);
      }
    }

    // Linear predictor at time plus epsilon
    if (assoc_uses[2] == 1) {
      if (K > 0) y_eta_q_eps = y_Xq_eps * beta;
      else y_eta_q_eps = rep_vector(0.0, (M*nrow_y_Xq));
      //if (has_offset == 1) y_eta_q_eps = y_eta_q_eps + y_offset; # how to handle offset?
      y_eta_q_eps = y_eta_q_eps + csr_matrix_times_vector((M*nrow_y_Xq), q, w_Zq_eps, v_Zq_eps, u_Zq_eps, b);
      for (m in 1:M) {
        y_eta_qwide_eps[m] = add_intercept(
          y_eta_q_eps, m, nrow_y_Xq, has_intercept, has_intercept_nob, 
          has_intercept_lob, has_intercept_upb, gamma_nob, 
          gamma_lob, gamma_upb, xbar, beta, KM);
      }
    }

    // Linear predictor at auc quadpoints
    if (assoc_uses[3] == 1) {
      if (K > 0) y_eta_q_auc = y_Xq_auc * beta;
      else y_eta_q_auc = rep_vector(0.0, (M*nrow_y_Xq_auc));
      //if (has_offset == 1) y_eta_q_auc = y_eta_q_auc + y_offset; # how to handle offset?
      y_eta_q_auc = y_eta_q_auc + csr_matrix_times_vector((M*nrow_y_Xq_auc), q, w_Zq_auc, v_Zq_auc, u_Zq_auc, b);
      for (m in 1:M) {
        y_eta_qwide_auc[m] = add_intercept(
          y_eta_q_auc, m, nrow_y_Xq_auc, has_intercept, has_intercept_nob, 
          has_intercept_lob, has_intercept_upb, gamma_nob, 
          gamma_lob, gamma_upb, xbar, beta, KM);
      }
    }   

    // Expected value 
    if (assoc_uses[4] == 1) 
      for (m in 1:M) 
        y_qwide[m] = evaluate_mu(y_eta_qwide[m], family[m], link[m]);
      
    // Expected value at time plus epsilon
    if (assoc_uses[5] == 1)
      for (m in 1:M) 
        y_qwide_eps[m] = evaluate_mu(y_eta_qwide_eps[m], family[m], link[m]);

    // Expected value at auc quadpoints
    if (assoc_uses[6] == 1) 
      for (m in 1:M) 
        y_qwide_auc[m] = evaluate_mu(y_eta_qwide_auc[m], family[m], link[m]);
