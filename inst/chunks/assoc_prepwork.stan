    # NB: the m loop take the linear predictors and/or expected
    # values evaluated in the previous step and, for submodel m, 
    # adds the intercept and stores the result in an element of an 
    # array rather than having all submodels in a single vector

    // Linear predictor at quadrature time
    if (assoc_uses[1] == 1) {
      if (sum_y_K > 0) y_eta_q = y_Xq_eta * y_beta;
      else y_eta_q = rep_vector(0.0, (M*nrow_y_Xq));
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
