  // this first block is same as tdata_glm.stan
  // but can't use that file because prior_dist throws error 
  len_z_T = 0;
  len_var_group = sum(p) * (t > 0);
  len_rho = sum(p) - t;
  pos = 1;
  for (i in 1:t) {
    if (p[i] > 1) {
      for (j in 1:p[i]) {
        delta[pos] = concentration[j];
        pos = pos + 1;
      }
    }
    for (j in 3:p[i]) len_z_T = len_z_T + p[i] - 1;
  }

  // indexing arrays
  for (m in 1:M) {
    // hs for each submodel
    hsM[m] = get_nvars_for_hs(prior_dist[m]);
  
    // indexing for hs global params
    if (prior_dist[m] == 3 || prior_dist[m] == 4) {
      idx_global[m,1] = len_global + 1;
      idx_global[m,2] = len_global + hsM[m]; 
      len_global      = len_global + hsM[m]; 
    } 
	else {
      idx_global[m,1] = 0;
	  idx_global[m,2] = 0;
    }
	
	// indexing for hs local params
	if (prior_dist[m] == 3) {
	  idx_local2[m,1] = len_local2 + 1;
	  idx_local2[m,2] = len_local2 + KM[m];
	  len_local2      = len_local2 + KM[m];
	} 
	else { 
	  idx_local2[m,1] = 0;
	  idx_local2[m,2] =	0;
	}
	
	// indexing for hs_plus local params
	if (prior_dist[m] == 4) {
	  idx_local4[m,1] = len_local4 + 1;
	  idx_local4[m,2] = len_local4 + KM[m];
	  len_local4      = len_local4 + KM[m];
	} 
	else { 
	  idx_local4[m,1] = 0;
	  idx_local4[m,2] =	0;
	}

    // indexing for shrinkage params
	if (prior_dist[m] == 5 || prior_dist[m] == 6) {
	  idx_mix[m,1] = len_mix + 1;
	  idx_mix[m,2] = len_mix + KM[m];
	  len_mix      = len_mix + KM[m];
	} 
	else {
	  idx_mix[m,1] = 0;
	  idx_mix[m,2] = 0;	
	}
	
    // indexing for one over lambda params
    if (prior_dist[m] == 6) {
	  idx_ool[m] = len_ool + 1;
	  len_ool = len_ool + 1;
	}
	else {
	  idx_ool[m] = 0;
	}
	
    // indexing for noise params
	if (family[m] == 8) {  // poisson-gamma mixture model
	  idx_noise[m,1] = len_noise + 1;
	  idx_noise[m,2] = len_noise + NM[m];
	  len_noise      = len_noise + NM[m];
	} 
	else {
	  idx_noise[m,1] = 0;
	  idx_noise[m,2] = 0;	
	}	
  }
  
  // transformations of outcome
  for (m in 1:M) {
    sum_log_y[m] = not_a_number();
    if (family[m] == 2 || family[m] == 3) {
      sum_log_y[m] = sum(log(y_real[idx_real[m,1]:idx_real[m,2]]));
    }
    if (family[m] == 3) {
      sqrt_y[idx_real[m,1]:idx_real[m,2]] = sqrt(y_real[idx_real[m,1]:idx_real[m,2]]);
      log_y[idx_real[m,1]:idx_real[m,2]] = log(y_real[idx_real[m,1]:idx_real[m,2]]);
    }
  }
