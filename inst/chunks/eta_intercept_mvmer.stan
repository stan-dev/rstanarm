    vector[NM[m]] eta_tmp;	           // eta for just one submodel 
    eta_tmp = eta[idx[m,1]:idx[m,2]];  // eta for just one submodel 
    if (has_intercept[m] == 1) {  // has intercept
      if (has_intercept_nob[m] == 1)
        eta_tmp = eta_tmp + gamma_nob[sum(has_intercept_nob[1:m])];
      else if (has_intercept_lob[m] == 1)
        eta_tmp = eta_tmp - min(eta_tmp) + gamma_lob[sum(has_intercept_lob[1:m])];
      else if (has_intercept_upb[m] == 1)
        eta_tmp = eta_tmp - max(eta_tmp) + gamma_upb[sum(has_intercept_upb[1:m])];					
    }
    else {  // no intercept, so model must have at least 1 predictor
      int K1 = idx_K[m,1];  // indexing for beta
      int K2 = idx_K[m,2];  // indexing for beta
      // correction to eta if model has no intercept (if X is centered)
      eta_tmp = eta_tmp + dot_product(xbar[K1:K2], beta[K1:K2]); 
    }
