  real alpha[sum_has_intercept];
  real mean_PPD[M];
  int mark = 1;  // indexing for alphas within m loop
  for (m in 1:M) {
    if (has_intercept[m] == 1) {
	    if (dense_X) {  // dense X
        int K1 = idx_K[m,1];  // indexing for beta
        int K2 = idx_K[m,2];  // indexing for beta
        if (has_intercept_nob[m] == 1)
          alpha[mark] = gamma_nob[sum(has_intercept_nob[1:m])] - dot_product(xbar[K1:K2], beta[K1:K2]);
        else if (has_intercept_lob[m] == 1)
          alpha[mark] = gamma_lob[sum(has_intercept_lob[1:m])] - dot_product(xbar[K1:K2], beta[K1:K2]);
        else if (has_intercept_upb[m] == 1)
          alpha[mark] = gamma_upb[sum(has_intercept_upb[1:m])] - dot_product(xbar[K1:K2], beta[K1:K2]);		
	    }
	    else {  // sparse X
        if (has_intercept_nob[m] == 1)
          alpha[mark] = gamma_nob[sum(has_intercept_nob[1:m])];
        else if (has_intercept_lob[m] == 1)
          alpha[mark] = gamma_lob[sum(has_intercept_lob[1:m])];
        else if (has_intercept_upb[m] == 1)
          alpha[mark] = gamma_upb[sum(has_intercept_upb[1:m])];			  
	    }
	    mark = mark + 1;  
    }   
  }
  {
    int aux_mark = 1; // indexing for auxiliary parameters in m loop
    #include "make_eta.stan" // defines eta
    if (t > 0) {
      #include "eta_add_Zb.stan"
    } 
	
	  mark = 1;  // reset indexing for alphas before new m loop
    for (m in 1:M) {
      vector[NM[m]] eta_tmp;	           // eta for just one submodel 
      eta_tmp = eta[idx[m,1]:idx[m,2]];  // eta for just one submodel 
      if (has_intercept[m] == 1) {
        if      (has_intercept_lob[m] == 1) alpha[mark] = alpha[mark] - min(eta_tmp);
        else if (has_intercept_upb[m] == 1) alpha[mark] = alpha[mark] - max(eta_tmp);
	      mark = mark + 1;	  
	    }
	    #include "eta_intercept_mvmer.stan"	// adds intercept or shifts eta
      if (family[m] == 8) {  // poisson-gamma mixture
	      #include "eta_add_noise_mvmer.stan"
      }    
      eta_tmp = evaluate_mu(eta_tmp, family[m], link[m]);
	
      mean_PPD[m] = 0;	
      if (family[m] == 1) {  // gaussian
	      for (n in 1:NM[m]) 
	        mean_PPD[m] = mean_PPD[m] + normal_rng(eta_tmp[n], aux[aux_mark]);
      }
      else if (family[m] == 2) {  // gamma
	      for (n in 1:NM[m]) 
	        mean_PPD[m] = mean_PPD[m] + gamma_rng(aux[aux_mark], aux[aux_mark] / eta_tmp[n]);
      }
      else if (family[m] == 3) {  // inverse gaussian
        for (n in 1:NM[m]) 
          mean_PPD[m] = mean_PPD[m] + inv_gaussian_rng(eta_tmp[n], aux[aux_mark]);
      }
      else if (family[m] == 4) {  // bernoulli
        for (n in 1:NM[m]) 
          mean_PPD[m] = mean_PPD[m] + bernoulli_rng(eta_tmp[n]);
      } 
      else if (family[m] == 5) {  // binomial
	    // binomial with num trials > 1 has been removed	  
      }
      else if (family[m] == 6 || family[m] == 8) { 
        for (n in 1:NM[m]) {  // poisson or poisson-gamma
          if (eta_tmp[n] < poisson_max) mean_PPD[m] = mean_PPD[m] + poisson_rng(eta_tmp[n]);
	        else mean_PPD[m] = mean_PPD[m] + normal_rng(eta_tmp[n], sqrt(eta_tmp[n]));
        }
      }
      else if (family[m] == 7) for (n in 1:NM[m]) {  // negative binomial
        real gamma_temp;
        if (is_inf(aux[aux_mark])) gamma_temp = eta_tmp[n];
        else gamma_temp = gamma_rng(aux[aux_mark], aux[aux_mark] / eta_tmp[n]);
        if (gamma_temp < poisson_max) mean_PPD[m] = mean_PPD[m] + poisson_rng(gamma_temp);
        else mean_PPD[m] = mean_PPD[m] + normal_rng(gamma_temp, sqrt(gamma_temp));	
      }		  
      if (has_aux[m] == 1) aux_mark = aux_mark + 1;
	    mean_PPD[m] = mean_PPD[m] / NM[m];
    }
  }
