   # !!! Be careful that indexing of has_assoc matches stan_jm file
   for (m in 1:M) {

      // etavalue and any interactions
      mark2 = mark2 + 1; // count even if assoc type isn't used
      if (has_assoc[1,m] == 1) { # etavalue
        mark = mark + 1;
	      e_eta_q = e_eta_q + a_beta[mark] * y_eta_qwide[m];
      }	
      if (has_assoc[11,m] == 1) { # etavalue*data
  	    int tmp;
  	    int j_shift;
  	    if (mark2 == 1) j_shift = 0;
  	    else j_shift = sum(a_K_data[1:(mark2-1)]);
  	    tmp = a_K_data[mark2];  
        for (j in 1:tmp) {
          int sel;
          sel = j_shift + j;
          mark = mark + 1;
          e_eta_q = e_eta_q + (y_eta_qwide[m] .* y_Xq_data[sel]) * a_beta[mark];
        }
      }
      mark3 = mark3 + 1; // count even if assoc type isn't used
      if (has_assoc[15,m] == 1) { # etavalue*etavalue
        for (j in 1:size_which_interactions[mark3]) { 
          int sel;
    	    int j_shift;
     	    if (mark3 == 1) j_shift = 0;
    	    else j_shift = sum(size_which_interactions[1:(mark3-1)]);
    	    sel = which_interactions[j+j_shift];
  	      mark = mark + 1;
          e_eta_q = e_eta_q + (y_eta_qwide[m] .* y_eta_qwide[sel]) * a_beta[mark];  
       }
      }
      mark3 = mark3 + 1; // count even if assoc type isn't used
      if (has_assoc[16,m] == 1) { # etavalue*muvalue
        for (j in 1:size_which_interactions[mark3]) { 
  	      int sel;
    	    int j_shift;
    	    if (mark3 == 1) j_shift = 0;
    	    else j_shift = sum(size_which_interactions[1:(mark3-1)]);
    	    sel = which_interactions[j+j_shift];
  	      mark = mark + 1;
          e_eta_q = e_eta_q + (y_eta_qwide[m] .* y_qwide[sel]) * a_beta[mark];  
        }
      }
      
      // etaslope and any interactions
      mark2 = mark2 + 1;
      if ((has_assoc[2,m] == 1) || (has_assoc[12,m] == 1)) {
        vector[nrow_y_Xq] dydt_eta_q;
        dydt_eta_q = (y_eta_qwide_eps[m] - y_eta_qwide[m]) / eps;
        if (has_assoc[2,m] == 1) { # etaslope
          mark = mark + 1;
          e_eta_q = e_eta_q + a_beta[mark] * dydt_eta_q;
        }
        if (has_assoc[12,m] == 1) { # etaslope*data
    	    int tmp;
    	    int j_shift;
    	    if (mark2 == 1) j_shift = 0;
    	    else j_shift = sum(a_K_data[1:(mark2-1)]);
    	    tmp = a_K_data[mark2];  
          for (j in 1:tmp) {
            int sel;
            sel = j_shift + j;
            mark = mark + 1;
            e_eta_q = e_eta_q + (dydt_eta_q .* y_Xq_data[sel]) * a_beta[mark];
          }    
        }         
      }
      
      // etalag
      if (has_assoc[3,m] == 1) { # etalag
        mark = mark + 1;
        e_eta_q = e_eta_q + a_beta[mark] * y_eta_qwide_lag[m];          
      }  
      
      // etaauc
      if (has_assoc[4,m] == 1) { # etaauc
        vector[nrow_y_Xq] y_eta_q_auc_tmp;  
        mark = mark + 1;
        for (r in 1:nrow_y_Xq) {
          vector[auc_quadnodes] val_tmp;
          vector[auc_quadnodes] wgt_tmp;
          val_tmp = y_eta_qwide_auc[m,((r-1) * auc_quadnodes + 1):(r * auc_quadnodes)];
          wgt_tmp = auc_quadweights[((r-1) * auc_quadnodes + 1):(r * auc_quadnodes)];
          y_eta_q_auc_tmp[r] = sum(wgt_tmp .* val_tmp);
        }
        e_eta_q = e_eta_q + a_beta[mark] * y_eta_q_auc_tmp;          
      }       
      
      // muvalue and any interactions
      mark2 = mark2 + 1;
      if (has_assoc[5,m] == 1) { # muvalue
        mark = mark + 1;
        e_eta_q = e_eta_q + a_beta[mark] * y_qwide[m]; 
      }
      if (has_assoc[13,m] == 1) { # muvalue*data
  	    int tmp;
  	    int j_shift;
  	    if (mark2 == 1) j_shift = 0;
  	    else j_shift = sum(a_K_data[1:(mark2-1)]);
  	    tmp = a_K_data[mark2];  
        for (j in 1:tmp) {
          int sel;
          sel = j_shift + j;
          mark = mark + 1;
          e_eta_q = e_eta_q + (y_qwide[m] .* y_Xq_data[sel]) * a_beta[mark];
        }      
      } 
      mark3 = mark3 + 1; // count even if assoc type isn't used
      if (has_assoc[17,m] == 1) { # muvalue*etavalue
        for (j in 1:size_which_interactions[mark3]) { 
          int sel;
    	    int j_shift;
     	    if (mark3 == 1) j_shift = 0;
    	    else j_shift = sum(size_which_interactions[1:(mark3-1)]);
    	    sel = which_interactions[j+j_shift];
  	      mark = mark + 1;
          e_eta_q = e_eta_q + (y_qwide[m] .* y_eta_qwide[sel]) * a_beta[mark];  
       }
      }      
      mark3 = mark3 + 1; // count even if assoc type isn't used
      if (has_assoc[18,m] == 1) { # muvalue*muvalue
        for (j in 1:size_which_interactions[mark3]) { 
          int sel;
    	    int j_shift;
     	    if (mark3 == 1) j_shift = 0;
    	    else j_shift = sum(size_which_interactions[1:(mark3-1)]);
    	    sel = which_interactions[j+j_shift];
  	      mark = mark + 1;
          e_eta_q = e_eta_q + (y_qwide[m] .* y_qwide[sel]) * a_beta[mark];  
       }
      }      
      
      // muslope and any interactions
      mark2 = mark2 + 1;
      if ((has_assoc[6,m] == 1) || (has_assoc[14,m] == 1)) {
        vector[nrow_y_Xq] dydt_q;
        dydt_q = (y_qwide_eps[m] - y_qwide[m]) / eps;
        if (has_assoc[6,m] == 1) { # muslope
          mark = mark + 1;
          e_eta_q = e_eta_q + a_beta[mark] * dydt_q;          
        }
        if (has_assoc[14,m] == 1) { # muslope*data
    	    int tmp;
    	    int j_shift;
    	    if (mark2 == 1) j_shift = 0;
    	    else j_shift = sum(a_K_data[1:(mark2-1)]);
    	    tmp = a_K_data[mark2];  
          for (j in 1:tmp) {
            int sel;
            sel = j_shift + j;
            mark = mark + 1;
            e_eta_q = e_eta_q + (dydt_q .* y_Xq_data[sel]) * a_beta[mark];
          }          
        } 
      }
      
      // mulag
      if (has_assoc[7,m] == 1) { # mulag
        mark = mark + 1;
        e_eta_q = e_eta_q + a_beta[mark] * y_qwide_lag[m];          
      }

      // muauc
      if (has_assoc[8,m] == 1) { # muauc
        vector[nrow_y_Xq] y_qwide_auc_tmp;  
        mark = mark + 1;
        for (r in 1:nrow_y_Xq) {
          vector[auc_quadnodes] val_tmp;
          vector[auc_quadnodes] wgt_tmp;
          val_tmp = y_qwide_auc[m,((r-1) * auc_quadnodes + 1):(r * auc_quadnodes)];
          wgt_tmp = auc_quadweights[((r-1) * auc_quadnodes + 1):(r * auc_quadnodes)];
          y_qwide_auc_tmp[r] = sum(wgt_tmp .* val_tmp);
        }
        e_eta_q = e_eta_q + a_beta[mark] * y_qwide_auc_tmp;          
      }  

    }
    
    // shared random effects
  	if (sum_size_which_b > 0) {
  	  int mark_beg;  // used to define segment of a_beta
  	  int mark_end;
  	  matrix[nrow_e_Xq,sum_size_which_b] x_assoc_shared_b;	  
      mark_beg = mark + 1;	  
  	  mark_end = mark + sum_size_which_b;
  	  x_assoc_shared_b = make_x_assoc_shared_b(
  	    b, l, p, pmat, Npat, quadnodes, which_b_zindex,
  	    sum_size_which_b, size_which_b, t_i, M);
  	  e_eta_q = e_eta_q + x_assoc_shared_b * a_beta[mark_beg:mark_end];
  	  mark = mark + sum_size_which_b;
    }	
  	if (sum_size_which_coef > 0) {
  	  int mark_beg;  // used to define segment of a_beta
  	  int mark_end;
  	  matrix[nrow_e_Xq,sum_size_which_coef] x_assoc_shared_coef;	  
      mark_beg = mark + 1;	  
  	  mark_end = mark + sum_size_which_coef;
  	  x_assoc_shared_coef = make_x_assoc_shared_coef(
  	    b, y_beta, y_K, M, t_i, l, p, pmat, Npat, quadnodes,
  	    sum_size_which_coef, size_which_coef,
  	    which_coef_zindex, which_coef_xindex,
  	    y_has_intercept, y_has_intercept_unbound,
  	    y_has_intercept_lobound, y_has_intercept_upbound,
  	    y_gamma_unbound, y_gamma_lobound, y_gamma_upbound);
  	  e_eta_q = e_eta_q + x_assoc_shared_coef * a_beta[mark_beg:mark_end];
  	  mark = mark + sum_size_which_coef;
    }    
