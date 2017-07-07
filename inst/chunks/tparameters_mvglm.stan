  vector[sum_has_aux] aux;  
  vector[K] beta;               // combined coefficients vector for all submodels
  vector[q] b_not_by_model;     // ordered by grouping factor (returned by make_b)
  vector[q] b;                  // ordered by submodel (structure of bdiag(Zmerge))
  vector[len_theta_L] theta_L; 
  
    // auxiliary parameters
  if (sum_has_aux > 0) {
    for (m in 1:M) {
      if (has_aux[m] == 1) {
        int mark = sum(has_aux[1:m]);
        if (prior_dist_for_aux[m] == 0) // none
          aux[mark] = aux_unscaled[mark];
        else {
          aux[mark] = prior_scale_for_aux[m] * aux_unscaled[mark];
          if (prior_dist_for_aux[m] <= 2) // normal or student_t
            aux[mark] = aux[mark] + prior_mean_for_aux[m];
        }
      }
    }
  }
  
  // coefficients
  if (prior_special_case == 1) { // same prior type for all submodels (none, normal or student-t)
  if      (prior_dist[1] == 0) beta = z_beta;
  else if (prior_dist[1] == 1) beta = z_beta .* prior_scale + prior_mean;
  else if (prior_dist[1] == 2) for (k in 1:K) {
    beta[k] = CFt(z_beta[k], prior_df[k]) * prior_scale[k] + prior_mean[k];
  }    
  }
  else {  // different prior types for each submodel
  for (m in 1:M) {
    int K1 = idx_K[m,1];  // indexing for beta vector
    int K2 = idx_K[m,2];  // indexing for beta vector
    if      (prior_dist[m] == 0) beta[K1:K2] = z_beta[K1:K2];
    else if (prior_dist[m] == 1) beta[K1:K2] = z_beta[K1:K2] .* prior_scale[K1:K2] + prior_mean[K1:K2];
    else if (prior_dist[m] == 2) for (k in K1:K2) {
      beta[k] = CFt(z_beta[k], prior_df[k]) * prior_scale[k] + prior_mean[k];
    }
    else if (prior_dist[m] == 3) {
	  int G1 = idx_global[m,1];  // indexing for global params
	  int G2 = idx_global[m,2];  // indexing for global params
	  int L1 = idx_local2[m,1];  // indexing for local params
	  int L2 = idx_local2[m,2];  // indexing for local params	  
      if (family[m] == 1) { // don't need is_continuous since family == 1 is gaussian in mvmer
		    int aux_mark = sum(has_aux[1:m]);
        beta[K1:K2] = hs_prior(z_beta[K1:K2], global[G1:G2], local2[,L1:L2], global_prior_scale[m], aux[aux_mark]);
	  }
      else beta[K1:K2] = hs_prior(z_beta[K1:K2], global[G1:G2], local2[,L1:L2], global_prior_scale[m], 1);
    }
    else if (prior_dist[m] == 4) {
      int G1 = idx_global[m,1];  // indexing for global params
	  int G2 = idx_global[m,2];  // indexing for global params
	  int L1 = idx_local4[m,1];  // indexing for local params
	  int L2 = idx_local4[m,2];  // indexing for local params
  	  if (family[m] == 1) { // don't need is_continuous since family == 1 is gaussian in mvmer
	  	int aux_mark = sum(has_aux[1:m]);
        beta[K1:K2] = hsplus_prior(z_beta[K1:K2], global[G1:G2], local4[,L1:L2], global_prior_scale[m], aux[aux_mark]);
	  }
      else beta[K1:K2] = hsplus_prior(z_beta[K1:K2], global[G1:G2], local4[,L1:L2], global_prior_scale[m], 1);
    }
    else if (prior_dist[m] == 5) // laplace
      beta[K1:K2] = prior_mean[K1:K2] + prior_scale[K1:K2] .* sqrt(2 * mix[1][idx_mix[m,1]:idx_mix[m,2]]) .* z_beta[K1:K2];
    else if (prior_dist[m] == 6) // lasso
      beta[K1:K2] = prior_mean[K1:K2] + ool[idx_ool[m]] * prior_scale[K1:K2] .* sqrt(2 * mix[1][idx_mix[m,1]:idx_mix[m,2]]) .* z_beta[K1:K2];  
  }
						 
  }
