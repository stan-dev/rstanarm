   # !!! Be careful that indexing of has_assoc matches stan_jm file
   for (m in 1:M) {

      // etavalue and any interactions
      mark2 = mark2 + 1; // count even if assoc type isn't used
      if (has_assoc[1,m] == 1) { # etavalue
        vector[nrow_e_Xq] val;    
        mark = mark + 1;
        if (has_clust == 1) val = clust_mat * y_eta_qwide[m];
        else val = y_eta_qwide[m];
	      e_eta_q = e_eta_q + a_beta[mark] * val;
      }	
      if (has_assoc[9,m] == 1) { # etavalue*data
  	    int tmp;
  	    int j_shift;
  	    if (mark2 == 1) j_shift = 0;
  	    else j_shift = sum(a_K_data[1:(mark2-1)]);
  	    tmp = a_K_data[mark2];  
        for (j in 1:tmp) {
          vector[nrow_e_Xq] val;    
          int sel;
          sel = j_shift + j;
          mark = mark + 1;
          if (has_clust == 1) val = clust_mat * (y_eta_qwide[m] .* y_Xq_data[sel]);
          else val = y_eta_qwide[m] .* y_Xq_data[sel];
          e_eta_q = e_eta_q + a_beta[mark] * val;
        }
      }
      mark3 = mark3 + 1; // count even if assoc type isn't used
      if (has_assoc[13,m] == 1) { # etavalue*etavalue
        for (j in 1:size_which_interactions[mark3]) { 
          vector[nrow_e_Xq] val;    
          int sel;
    	    int j_shift;
     	    if (mark3 == 1) j_shift = 0;
    	    else j_shift = sum(size_which_interactions[1:(mark3-1)]);
    	    sel = which_interactions[j+j_shift];
  	      mark = mark + 1;
          if (has_clust == 1) val = clust_mat * (y_eta_qwide[m] .* y_eta_qwide[sel]);
          else val = y_eta_qwide[m] .* y_eta_qwide[sel];
          e_eta_q = e_eta_q + a_beta[mark] * val;  
       }
      }
      mark3 = mark3 + 1; // count even if assoc type isn't used
      if (has_assoc[14,m] == 1) { # etavalue*muvalue
        for (j in 1:size_which_interactions[mark3]) { 
          vector[nrow_e_Xq] val;    
  	      int sel;
    	    int j_shift;
    	    if (mark3 == 1) j_shift = 0;
    	    else j_shift = sum(size_which_interactions[1:(mark3-1)]);
    	    sel = which_interactions[j+j_shift];
  	      mark = mark + 1;
          if (has_clust == 1) val = clust_mat * (y_eta_qwide[m] .* y_qwide[sel]);
          else val = y_eta_qwide[m] .* y_qwide[sel];  	      
          e_eta_q = e_eta_q + a_beta[mark] * val;  
        }
      }
      
      // etaslope and any interactions
      mark2 = mark2 + 1;
      if ((has_assoc[2,m] == 1) || (has_assoc[10,m] == 1)) {
        vector[nrow_y_Xq] dydt_eta_q;
        dydt_eta_q = (y_eta_qwide_eps[m] - y_eta_qwide[m]) / eps;
        if (has_assoc[2,m] == 1) { # etaslope
          vector[nrow_e_Xq] val;    
          mark = mark + 1;
          if (has_clust == 1) val = clust_mat * dydt_eta_q;
          else val = dydt_eta_q;          
          e_eta_q = e_eta_q + a_beta[mark] * val;
        }
        if (has_assoc[10,m] == 1) { # etaslope*data
    	    int tmp;
    	    int j_shift;
    	    if (mark2 == 1) j_shift = 0;
    	    else j_shift = sum(a_K_data[1:(mark2-1)]);
    	    tmp = a_K_data[mark2];  
          for (j in 1:tmp) {
            vector[nrow_e_Xq] val;    
            int sel;
            sel = j_shift + j;
            mark = mark + 1;
            if (has_clust == 1) val = clust_mat * (dydt_eta_q .* y_Xq_data[sel]);
            else val = dydt_eta_q .* y_Xq_data[sel];            
            e_eta_q = e_eta_q + a_beta[mark] * val;
          }    
        }         
      }
      
      // etaauc
      if (has_assoc[3,m] == 1) { # etaauc
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
      if (has_assoc[4,m] == 1) { # muvalue
        vector[nrow_e_Xq] val;    
        mark = mark + 1;
        if (has_clust == 1) val = clust_mat * y_qwide[m];
        else val = y_qwide[m];        
        e_eta_q = e_eta_q + a_beta[mark] * val; 
      }
      if (has_assoc[11,m] == 1) { # muvalue*data
  	    int tmp;
  	    int j_shift;
  	    if (mark2 == 1) j_shift = 0;
  	    else j_shift = sum(a_K_data[1:(mark2-1)]);
  	    tmp = a_K_data[mark2];  
        for (j in 1:tmp) {
          vector[nrow_e_Xq] val;    
          int sel;
          sel = j_shift + j;
          mark = mark + 1;
          if (has_clust == 1) val = clust_mat * (y_qwide[m] .* y_Xq_data[sel]);
          else val = y_qwide[m] .* y_Xq_data[sel];              
          e_eta_q = e_eta_q + a_beta[mark] * val;
        }      
      } 
      mark3 = mark3 + 1; // count even if assoc type isn't used
      if (has_assoc[15,m] == 1) { # muvalue*etavalue
        for (j in 1:size_which_interactions[mark3]) {
          vector[nrow_e_Xq] val;    
          int sel;
    	    int j_shift;
     	    if (mark3 == 1) j_shift = 0;
    	    else j_shift = sum(size_which_interactions[1:(mark3-1)]);
    	    sel = which_interactions[j+j_shift];
  	      mark = mark + 1;
          if (has_clust == 1) val = clust_mat * (y_qwide[m] .* y_eta_qwide[sel]);
          else val = y_qwide[m] .* y_eta_qwide[sel];        	      
          e_eta_q = e_eta_q + a_beta[mark] * val;  
       }
      }      
      mark3 = mark3 + 1; // count even if assoc type isn't used
      if (has_assoc[16,m] == 1) { # muvalue*muvalue
        for (j in 1:size_which_interactions[mark3]) { 
          vector[nrow_e_Xq] val;    
          int sel;
    	    int j_shift;
     	    if (mark3 == 1) j_shift = 0;
    	    else j_shift = sum(size_which_interactions[1:(mark3-1)]);
    	    sel = which_interactions[j+j_shift];
  	      mark = mark + 1;
          if (has_clust == 1) val = clust_mat * (y_qwide[m] .* y_qwide[sel]);
          else val = y_qwide[m] .* y_qwide[sel];   	      
          e_eta_q = e_eta_q + a_beta[mark] * val;  
       }
      }      
      
      // muslope and any interactions
      mark2 = mark2 + 1;
      if ((has_assoc[5,m] == 1) || (has_assoc[12,m] == 1)) {
        vector[nrow_y_Xq] dydt_q;
        dydt_q = (y_qwide_eps[m] - y_qwide[m]) / eps;
        if (has_assoc[5,m] == 1) { # muslope
          vector[nrow_e_Xq] val;    
          mark = mark + 1;
          if (has_clust == 1) val = clust_mat * dydt_q;
          else val = dydt_q;  
          e_eta_q = e_eta_q + a_beta[mark] * val;          
        }
        if (has_assoc[12,m] == 1) { # muslope*data
    	    int tmp;
    	    int j_shift;
    	    if (mark2 == 1) j_shift = 0;
    	    else j_shift = sum(a_K_data[1:(mark2-1)]);
    	    tmp = a_K_data[mark2];  
          for (j in 1:tmp) {
            vector[nrow_e_Xq] val;    
            int sel;
            sel = j_shift + j;
            mark = mark + 1;
            if (has_clust == 1) val = clust_mat * (dydt_q .* y_Xq_data[sel]);
            else val = dydt_q .* y_Xq_data[sel];             
            e_eta_q = e_eta_q + a_beta[mark] * val;
          }          
        } 
      }
      
      // muauc
      if (has_assoc[6,m] == 1) { # muauc
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
  	    b_not_by_model, l, p, pmat, Npat, quadnodes, which_b_zindex,
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
  	    b_not_by_model, beta, KM, M, t_i, l, p, pmat, Npat, quadnodes,
  	    sum_size_which_coef, size_which_coef,
  	    which_coef_zindex, which_coef_xindex,
  	    has_intercept, has_intercept_nob,
  	    has_intercept_lob, has_intercept_upb,
  	    gamma_nob, gamma_lob, gamma_upb);
  	  e_eta_q = e_eta_q + x_assoc_shared_coef * a_beta[mark_beg:mark_end];
  	  mark = mark + sum_size_which_coef;
    }    
