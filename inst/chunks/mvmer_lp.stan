	// accumulate log-likelihood for submodel m
	
	// unweighted log-likelihoods
    if (has_weights == 0 && prior_PD == 0) {  
	  int R1 = idx_real[m,1];  // indexing for outcome vector of reals
	  int R2 = idx_real[m,2];  // indexing for outcome vector of reals
	  int I1 = idx_int[m,1];   // indexing for outcome array of integers
	  int I2 = idx_int[m,2];   // indexing for outcome array of integers
      if (family[m] == 1) {  // gaussian
        if (link[m] == 1)      target += normal_lpdf(y_real[R1:R2] | eta_tmp, aux[aux_mark]);
        else if (link[m] == 2) target += lognormal_lpdf(y_real[R1:R2] | eta_tmp, aux[aux_mark]);
        else target += normal_lpdf(y_real[R1:R2] | divide_real_by_vector(1, eta_tmp), aux[aux_mark]);
      }
      else if (family[m] == 2) {  // gamma
        target += GammaReg(y_real[R1:R2], eta_tmp, aux[aux_mark], link[m], sum_log_y[m]);
      }
      else if (family[m] == 3) {  // inverse gaussian 
        target += inv_gaussian(y_real[R1:R2], linkinv_inv_gaussian(eta_tmp, link[m]), 
                               aux[aux_mark], sum_log_y[m], sqrt_y[R1:R2]);
      }
	    else if (family[m] == 4) {  // bernoulli
		    vector[N01[m,1]] eta0_tmp;
		    vector[N01[m,2]] eta1_tmp;
	      real dummy;  // irrelevant but useful for testing
		    eta0_tmp = segment(eta_tmp, 1, N01[m,1]);
		    eta1_tmp = segment(eta_tmp, (N01[m,1] + 1), N01[m,2]);
	      dummy = ll_bern_lp(eta0_tmp, eta1_tmp, link[m], N01[m,]);	  
	    }
	    else if (family[m] == 5) {  // binomial
	      // binomial with num trials > 1 has been removed	  
	    }
	    else if (family[m] == 6 || family[m] == 8) {  // poisson or poisson-gamma
          if (link[m] == 1) target += poisson_log_lpmf(y_int[I1:I2] | eta_tmp);
          else target += poisson_lpmf(y_int[I1:I2] | linkinv_count(eta_tmp, link[m]));
	    }
	    else if (family[m] == 7) {  // negative binomial
  	      if (link[m] == 1) target += neg_binomial_2_log_lpmf(y_int[I1:I2] | eta_tmp, aux[aux_mark]);
          else target += neg_binomial_2_lpmf(y_int[I1:I2] | linkinv_count(eta_tmp, link[m]), aux[aux_mark]);
	    }	    
    }

    // weighted log-likelihoods    
    else if (prior_PD == 0) { 
	  int R1 = idx_real[m,1];  // indexing for outcome vector of reals
	  int R2 = idx_real[m,2];  // indexing for outcome vector of reals
	  int I1 = idx_int[m,1];   // indexing for outcome array of integers
	  int I2 = idx_int[m,2];   // indexing for outcome array of integers
  	  vector[NM[m]] weights_tmp;	  
  	  vector[NM[m]] summands;
      weights_tmp = weights[idx[m,1]:idx[m,2]];	  
  	  if (family[m] == 1) {  // gaussian
  	    summands = pw_gauss(y_real[R1:R2], eta_tmp, aux[aux_mark], link[m]);
  	    target += dot_product(weights_tmp, summands);	  	    
  	  }
  	  else if (family[m] == 2) {  // gamma
  	    summands = pw_gamma(y_real[R1:R2], eta_tmp, aux[aux_mark], link[m]);
  	    target += dot_product(weights_tmp, summands);	  	    
  	  }
  	  else if (family[m] == 3) {  // inverse gaussian
  	    vector[NM[m]] log_y_tmp;	  
  	    vector[NM[m]] sqrt_y_tmp;  	    
    	log_y_tmp = log_y[idx[m,1]:idx[m,2]];
    	sqrt_y_tmp = sqrt_y[idx[m,1]:idx[m,2]];
  	    summands = pw_inv_gaussian(y_real[R1:R2], eta_tmp, aux[aux_mark], link[m], log_y_tmp, sqrt_y_tmp);
  	    target += dot_product(weights_tmp, summands);
  	  }
	  else if (family[m] == 4) {  // bernoulli
    	vector[N01[m,1]] weights0_tmp;
    	vector[N01[m,2]] weights1_tmp;
    	vector[N01[m,1]] eta0_tmp;
    	vector[N01[m,2]] eta1_tmp;
    	eta0_tmp = segment(eta_tmp, 1, N01[m,1]);
    	eta1_tmp = segment(eta_tmp, (N01[m,1] + 1), N01[m,2]);
    	weights0_tmp = segment(weights_tmp, 1, N01[m,1]);
    	weights1_tmp = segment(weights_tmp, (N01[m,1] + 1), N01[m,2]);		
        target += dot_product(weights0_tmp, pw_bern(0, eta0_tmp, link[m]));
        target += dot_product(weights1_tmp, pw_bern(1, eta1_tmp, link[m]));
  	  }
  	  else if (family[m] == 5) {  // binomial
	      // binomial with num trials > 1 has been removed
  	  }
  	  else if (family[m] == 6 || family[m] == 8) {  // poisson or poisson-gamma
        target += dot_product(weights_tmp, pw_pois(y_int[I1:I2], eta_tmp, link[m]));    		
  	  }
  	  else if (family[m] == 7) {  // negative binomial
        target += dot_product(weights_tmp, pw_nb(y_int[I1:I2], eta_tmp, aux[aux_mark], link[m]));
  	  }  	  
    }
