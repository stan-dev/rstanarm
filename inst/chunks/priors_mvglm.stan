  // log-priors for coefficients
  if (prior_special_case == 1) {  // same prior type for all submodels (none, normal or student-t)
    if      (prior_dist[1] == 1) target += normal_lpdf(z_beta | 0, 1);
    else if (prior_dist[1] == 2) target += normal_lpdf(z_beta | 0, 1); // Student t  
  } 
  else {  // different prior types for each submodel
    for (m in 1:M) {
    int K1 = idx_K[m,1];  // indexing for beta vector
    int K2 = idx_K[m,2];  // indexing for beta vector
  if      (prior_dist[m] == 1) target += normal_lpdf(z_beta[K1:K2] | 0, 1);
  else if (prior_dist[m] == 2) target += normal_lpdf(z_beta[K1:K2] | 0, 1); // Student t
  else if (prior_dist[m] == 3) { // hs
	  int G1 = idx_global[m,1];  // indexing for global params
	  int G2 = idx_global[m,2];  // indexing for global params
	  int L1 = idx_local2[m,1];  // indexing for local params
	  int L2 = idx_local2[m,2];  // indexing for local params	  
    target += normal_lpdf(z_beta[K1:K2] | 0, 1);
    target += normal_lpdf(local2[1,L1:L2] | 0, 1);
    target += inv_gamma_lpdf(local2[2,L1:L2] | 0.5 * prior_df[K1:K2], 0.5 * prior_df[K1:K2]);
    target += normal_lpdf(global[G1] | 0, 1);
    target += inv_gamma_lpdf(global[G1+1] | 0.5 * global_prior_df[m], 0.5 * global_prior_df[m]);
  }
  else if (prior_dist[m] == 4) { // hs+
      int G1 = idx_global[m,1];  // indexing for global params
	  int G2 = idx_global[m,2];  // indexing for global params
	  int L1 = idx_local4[m,1];  // indexing for local params
	  int L2 = idx_local4[m,2];  // indexing for local params  
    target += normal_lpdf(z_beta[K1:K2] | 0, 1);
    target += normal_lpdf(local4[1,L1:L2] | 0, 1);
    target += inv_gamma_lpdf(local4[2,L1:L2] | 0.5 * prior_df[K1:K2], 0.5 * prior_df[K1:K2]);
    target += normal_lpdf(local4[3,L1:L2] | 0, 1);
    // unorthodox useage of prior_scale as another df hyperparameter
    target += inv_gamma_lpdf(local4[4,L1:L2] | 0.5 * prior_scale[K1:K2], 0.5 * prior_scale[K1:K2]);
    target += normal_lpdf(global[G1] | 0, 1);
    target += inv_gamma_lpdf(global[G1+1] | 0.5 * global_prior_df[m], 0.5 * global_prior_df[m]);
  }
  else if (prior_dist[m] == 5) { // laplace
    target += normal_lpdf(z_beta[K1:K2] | 0, 1);
    target += exponential_lpdf(mix[1][idx_mix[m,1]:idx_mix[m,2]] | 1);
  }
  else if (prior_dist[m] == 6) { // lasso
    target += normal_lpdf(z_beta[K1:K2] | 0, 1);
    target += exponential_lpdf(mix[1][idx_mix[m,1]:idx_mix[m,2]] | 1);
    target += chi_square_lpdf(ool[idx_ool[m]] | prior_df[K1]);
  }
  /* else prior_dist is 0 and nothing is added */	
  }
  }
  
  // log-priors for intercept parameters
  if (sum_has_intercept > 0) {
    int mark1 = 1; // indexing for unbounded intercepts
    int mark2 = 1; // indexing for lower bounded intercepts
    int mark3 = 1; // indexing for upper bounded intercepts
    for (m in 1:M) {
	if (has_intercept[m] == 1) {
      if (has_intercept_nob[m] == 1) {  // unbounded intercept
        gamma_lp(gamma_nob[mark1], prior_dist_for_intercept[m], 
                 prior_mean_for_intercept[m], prior_scale_for_intercept[m], 
                 prior_df_for_intercept[m]);
        mark1 = mark1 + 1;
      }
      // lower bounded intercept
      else if (has_intercept_lob[m] == 1) { // lower bounded intercept
        gamma_lp(gamma_lob[mark2], prior_dist_for_intercept[m], 
                 prior_mean_for_intercept[m], prior_scale_for_intercept[m], 
                 prior_df_for_intercept[m]);
        mark2 = mark2 + 1;
      }
      // upper bounded intercept
      else if (has_intercept_upb[m] == 1) { // upper bounded intercept
        gamma_lp(gamma_upb[mark3], prior_dist_for_intercept[m], 
                 prior_mean_for_intercept[m], prior_scale_for_intercept[m], 
                 prior_df_for_intercept[m]);
        mark3 = mark3 + 1;
      }
    }
	}
  }
