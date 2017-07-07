  /** 
  * Return the required number of local hs parameters
  *
  * @param prior_dist An integer indicating the prior distribution
  * @return An integer
  */
  int get_nvars_for_hs(int prior_dist) {
    int hs = 0;
    if (prior_dist == 3) hs = 2;
    else if (prior_dist == 4) hs = 4;
    return hs;
  }

  /** 
  * Generate betas using the primitive coefficients and prior information
  *
  * @param z_beta Vector of primitive coefficients
  * @param prior_dist An integer indicating the prior distribution
  * @param prior_{mean,scale,df} Vector of prior means, scales and dfs
  * @return A vector
  */
  vector generate_beta(vector z_beta, int prior_dist, vector prior_mean, 
                       vector prior_scale, vector prior_df, real[] global, 
                       vector[] local, real global_prior_scale, 
                       real[] one_over_lambda, vector[] mix) {
    vector[rows(z_beta)] beta;
    if      (prior_dist == 0) beta = z_beta;
    else if (prior_dist == 1) beta = z_beta .* prior_scale + prior_mean;
    else if (prior_dist == 2) for (k in 1:rows(z_beta)) {
      beta[k] = CFt(z_beta[k], prior_df[k]) * prior_scale[k] + prior_mean[k];
    }
    else if (prior_dist == 3) beta = hs_prior(z_beta, global, local, global_prior_scale, 1.0);
    else if (prior_dist == 4) beta = hsplus_prior(z_beta, global, local, global_prior_scale, 1.0);
    else if (prior_dist == 5) // laplace
      beta = prior_mean + prior_scale .* sqrt(2 * mix[1]) .* z_beta;
    else if (prior_dist == 6) // lasso
      beta = prior_mean + one_over_lambda[1] * prior_scale .* sqrt(2 * mix[1]) .* z_beta;
    return beta;
  }

  /** 
  * Generate auxiliary parameters using unscaled params and prior information
  *
  * @param aux_unscaled Vector of unscaled auxiliary params
  * @param prior_dist An integer indicating the prior distribution
  * @param prior_{mean,scale} Vector of prior means and scales
  * @return A vector
  */  
  vector generate_aux(vector aux_unscaled, int prior_dist, vector prior_mean, vector prior_scale) {
    vector[rows(aux_unscaled)] aux;
    if (prior_dist == 0) // none
      aux = aux_unscaled;
    else {
      aux = prior_scale .* aux_unscaled;
      if (prior_dist <= 2) // normal or student_t
        aux = aux + prior_mean;
    }
    return aux;	
  }
  
  /** 
  * Add intercept term to linear predictor for submodel m 
  *
  * @param eta Vector of linear predictors for all longitudinal submodels
  * @param m Integer specifying which submodel to return 'eta + intercept' for
  * @param len_eta_m Integer specifying the length of the returned eta vector
  *   (this is the same for each longitudinal submodel since it corresponds to 
  *    the number of quadrature points and individuals in the model, rather than
  *    the number of longitudinal observations)
  * @param has_intercept{_unbound,_lobound,_upbound} Integer arrays indicating 
  *   whether each submodel has (the specific type of) intercept
  * @param gamma{_unbound,_lobound,_upbound} Vector containing the parameters
  *   for each type of intercept
  * @param xbar Vector of predictor means
  * @param beta Vector of coefficients across all longitudinal submodels
  * @param K Integer array specifying the number of parameters in each
  *   longitudinal submodel
  * @return A vector of length len_eta_m
  */
  vector add_intercept(vector eta, int m, int len_eta_m, int[] has_intercept, 
                       int[] has_intercept_nob, int[] has_intercept_lob, 
                       int[] has_intercept_upb, real[] gamma_nob, 
                       real[] gamma_lob, real[] gamma_upb, 
                       vector xbar, vector beta, int[] K) {
    vector[len_eta_m] eta_m;
    eta_m = segment(eta, ((m-1) * len_eta_m) + 1, len_eta_m);
    if (has_intercept[m] == 1) {
      if (has_intercept_nob[m] == 1) 
        eta_m = eta_m + gamma_nob[sum(has_intercept_nob[1:m])];
      else if (has_intercept_lob[m] == 1)
        eta_m = eta_m - min(eta_m) + gamma_lob[sum(has_intercept_lob[1:m])];
	    else if (has_intercept_upb[m] == 1)
        eta_m = eta_m - max(eta_m) + gamma_upb[sum(has_intercept_upb[1:m])];
    } else {  // no intercept, so model must have at least 1 predictor
      int K1 = (m == 1) ? 1 : (sum(K[1:(m-1)]) + 1);
      int K2 = sum(K[1:m]);  
      // correction to eta if model has no intercept (and X is centered)
      eta_m = eta_m + dot_product(xbar[K1:K2], beta[K1:K2]);
    }
    return eta_m;
  } 

  /** 
  * Evaluate mu based on eta, family and link 
  *
  * @param eta Vector of linear predictors
  * @param family An integer indicating the family
  * @param link An integer indicating the link function (differs by family)
  * @return A vector
  */
  vector evaluate_mu(vector eta, int family, int link) {
    vector[rows(eta)] mu;
     if (family == 1) 
      mu = linkinv_gauss(eta, link);
    else if (family == 2) 
      mu = linkinv_gamma(eta, link);
    else if (family == 3)
      mu = linkinv_inv_gaussian(eta, link);
    else if (family == 4)
      mu = linkinv_bern(eta, link);	
    else if (family == 5)		  
      mu = linkinv_binom(eta, link);
    else if (family == 6 || family == 7 || family == 8)		  
      mu = linkinv_count(eta, link);      
    return mu;
  }

  /** 
  * Log-prior for coefficients
  *
  * @param z_beta Vector of primative coefficients
  * @param prior_dist Integer, the type of prior distribution
  * @param prior_scale Real, scale for the prior distribution
  * @param prior_df Real, df for the prior distribution
  * @param global_prior_df Real, df for the prior for the global hs parameter
  * @param local Vector of hs local parameters
  * @param global Real, the global parameter
  * @param mix Vector of shrinkage parameters
  * @param one_over_lambda Real
  * @return nothing
  */
  void beta_lp(vector z_beta, int prior_dist, vector prior_scale,
                 vector prior_df, real global_prior_df, vector[] local,
                 real[] global, vector[] mix, real[] one_over_lambda) {
    if      (prior_dist == 1) target += normal_lpdf(z_beta | 0, 1);
    else if (prior_dist == 2) target += normal_lpdf(z_beta | 0, 1); // Student t
    else if (prior_dist == 3) { // hs
      target += normal_lpdf(z_beta | 0, 1);
      target += normal_lpdf(local[1] | 0, 1);
      target += inv_gamma_lpdf(local[2] | 0.5 * prior_df, 0.5 * prior_df);
      target += normal_lpdf(global[1] | 0, 1);
      target += inv_gamma_lpdf(global[2] | 0.5 * global_prior_df, 0.5 * global_prior_df);
    }
    else if (prior_dist == 4) { // hs+
      target += normal_lpdf(z_beta | 0, 1);
      target += normal_lpdf(local[1] | 0, 1);
      target += inv_gamma_lpdf(local[2] | 0.5 * prior_df, 0.5 * prior_df);
      target += normal_lpdf(local[3] | 0, 1);
      // unorthodox useage of prior_scale as another df hyperparameter
      target += inv_gamma_lpdf(local[4] | 0.5 * prior_scale, 0.5 * prior_scale);
      target += normal_lpdf(global[1] | 0, 1);
      target += inv_gamma_lpdf(global[2] | 0.5 * global_prior_df, 0.5 * global_prior_df);
    }
    else if (prior_dist == 5) { // laplace
      target += normal_lpdf(z_beta | 0, 1);
      target += exponential_lpdf(mix[1] | 1);
    }
    else if (prior_dist == 6) { // lasso
      target += normal_lpdf(z_beta | 0, 1);
      target += exponential_lpdf(mix[1] | 1);
      target += chi_square_lpdf(one_over_lambda[1] | prior_df[1]);
    }
    else if (prior_dist == 7) { // product_normal
      target += normal_lpdf(z_beta | 0, 1);
    }
    /* else prior_dist is 0 and nothing is added */
  }  
  
  /** 
  * Log-prior for intercept parameters
  *
  * @param gamma Real, the intercept parameter
  * @param dist Integer, the type of prior distribution
  * @param mean Real, mean of prior distribution
  * @param scale Real, scale for the prior distribution
  * @param df Real, df for the prior distribution
  * @return nothing
  */  
  void gamma_lp(real gamma, int dist, real mean, real scale, real df) {
    if (dist == 1)  // normal
      target += normal_lpdf(gamma | mean, scale);
    else if (dist == 2)  // student_t
      target += student_t_lpdf(gamma | df, mean, scale);
    /* else dist is 0 and nothing is added */
  }

  /** 
  * Log-prior for auxiliary parameters
  *
  * @param aux_unscaled Vector (potentially of length 1) of unscaled 
  *   auxiliary parameter(s)
  * @param dist Integer specifying the type of prior distribution
  * @param scale Real specifying the scale for the prior distribution
  * @param df Real specifying the df for the prior distribution
  * @return nothing
  */
  void aux_lp(real aux_unscaled, int dist, real scale, real df) {
    if (dist > 0 && scale > 0) {
      if (dist == 1)
        target += normal_lpdf(aux_unscaled | 0, 1);
      else if (dist == 2)
        target += student_t_lpdf(aux_unscaled | df, 0, 1);
      else 
        target += exponential_lpdf(aux_unscaled | 1);
    }     
  }
  
  /** 
  * Log-prior for baseline hazard parameters
  *
  * @param aux_unscaled Vector (potentially of length 1) of unscaled 
  *   auxiliary parameter(s)
  * @param dist Integer specifying the type of prior distribution
  * @param scale Real specifying the scale for the prior distribution
  * @param df Real specifying the df for the prior distribution
  * @return nothing
  */
  void basehaz_lp(vector aux_unscaled, int dist, vector scale, vector df) {
    if (dist > 0) {
      if (dist == 1)
        target += normal_lpdf(aux_unscaled | 0, 1);
      else if (dist == 2)
        target += student_t_lpdf(aux_unscaled | df, 0, 1);
      else 
        target += exponential_lpdf(aux_unscaled | 1);
    }     
  }  

  /** 
  * Reorder the vector of group-specific coefficients
  *
  * @param b Vector whose elements are ordered in the following nested way:
  *   factor level (e.g. "subject ID" or "clinic ID") within grouping 
  *   factor (e.g. "subjects" or "clinics")
  * @param p An integer array with the number variables on the LHS of each |
  * @param pmat A matrix with the number variables on the LHS of each | in each
  *   longitudinal submodel. The rows correspond to each |, meaning the separate
  *   equations for each grouping variable, and the columns correspond to each
  *   longitudinal submodel. If subject ID is the only grouping variable then the
  *   matrix will have one row. If the joint model only has one longitudinal 
  *   submodel then the matrix will have one column.
  * @param p An integer array with the number of random coefficients for the 
  *   LHS of each |
  * @param qmat A matrix with the number of random coefficients for the LHS of each | 
  *   The rows correspond to each |, meaning the separate equations for 
  *   each grouping variable, and the columns correspond to each longitudinal
  *   submodel.
  * @param l An integer array with the number of levels for the factor(s) on 
  *   the RHS of each |
  * @param M An integer specifying the number of longitudinal submodels
  * @return A vector (containing the group-specific coefficients) whose whose 
  *   elements are ordered in the following nested way: factor level, within
  *   grouping factor, within longitudinal submodel
  */ 
  vector reorder_b(vector b, int[] p, int[,] pmat, int[] q1, int[] q2, 
                   int[,] qmat, int[] l, int M) {
    vector[rows(b)] b_new;
    # in the loops below:
    #   i = grouping factors, j = levels, m = models
    #   there are sum(q) random coefficients (ie, values within 
    #   vector b)
    
    # collection is ascending in terms of:
    #   sum q over 1:(i-1)
    #   sum q over 1:(j-1) within i 
    #   sum q over 1:(m-1) within j within i 
    for (i in 1:size(p)) { 
      int i_beg;
      int i_end;
      vector[q1[i]] b_i;
      if (i == 1) i_beg = 1;
      else i_beg = sum(q1[1:(i-1)]) + 1;
      i_end = sum(q1[1:i]);
      b_i = b[i_beg:i_end];
      
      for (j in 1:l[i]) {
        int j_beg;
        int j_end;
        vector[p[i]] b_ij;
        if (j == 1) j_beg = 1;
        else j_beg = (j-1) * p[i] + 1;
        j_end = j * p[i];
        b_ij = b_i[j_beg:j_end];
        
        for (m in 1:M) {
          int m_beg;
          int m_end;
          vector[pmat[i,m]] b_ijm;
          if (pmat[i,m] > 0) {
            int m_sum;
            int i_sum;
            int j_sum;
            int store_beg;
            int store_end;
            if (m == 1) m_beg = 1;
            else m_beg = sum(pmat[i, 1:(m-1)]) + 1;
            m_end = sum(pmat[i, 1:m]);
            b_ijm = b_ij[m_beg:m_end];
          
            # storage is ascending in terms of:
            #   sum q over 1:(m-1)
            #   sum q over 1:(i-1) within m 
            #   sum q over 1:(j-1) within i within m
            if (m == 1) m_sum = 0;
            else m_sum = sum(q2[1:(m-1)]);
            if (i == 1) i_sum = 0;
            else i_sum = sum(qmat[1:(i-1),m]);
            if (j == 1) j_sum = 0;
            else j_sum = (j-1) * pmat[i,m];
            store_beg = m_sum + i_sum + j_sum + 1;
            store_end = m_sum + i_sum + j_sum + pmat[i,m];
            b_new[store_beg:store_end] = b_ijm;          
          }
        }
      }  
    }
    return b_new;
  }  
  
  /** 
  * Create a design matrix for a shared random effects association
  * structure in the joint model
  *
  * @param b Vector of group-specific coefficients
  * @param l An integer array with the number of levels for the factor(s) on 
  *   the RHS of each |
  * @param p An integer array with the number of variables on the LHS of each |
  * @param pmat A matrix with the number variables on the LHS of each | in each
  *   longitudinal submodel. The rows correspond to each |, meaning the separate
  *   equations for each grouping variable, and the columns correspond to each
  *   longitudinal submodel. If subject ID is the only grouping variable then the
  *   matrix will have one row. If the joint model only has one longitudinal 
  *   submodel then the matrix will have one column.
  * @param Npat Integer specifying number of individuals represented 
  *   in vector b
  * @param quadnodes The number of quadrature nodes
  * @param which_b Integer array specifying the indices
  *   of the random effects to use in the association structure
  * @param sum_size_which_b Integer specifying total number of 
  *   random effects that are to be used in the association structure
  * @param size_which_b Integer array specifying number of random effects from
  *   each long submodel that are to be used in the association structure
  * @param t_i Integer specifying the index of the grouping factor that
  *   corresponds to the patient-level
  * @param M An integer specifying the number of longitudinal submodels
  * @return A matrix with the desired random effects represented
  *   in columns, and the individuals on the rows; the matrix is
  *   repeated (quadnodes + 1) times (bounded by rows)
  */  
  matrix make_x_assoc_shared_b(
    vector b, int[] l, int[] p, int[,] pmat, int Npat, int quadnodes,
    int[] which_b, int sum_size_which_b, int[] size_which_b, int t_i, int M) {
    int prior_shift;    // num. ranefs prior to subject-specific ranefs
    int start_store;
    int end_store;	
    matrix[Npat,sum_size_which_b] temp;
    matrix[(Npat*(quadnodes+1)),sum_size_which_b] x_assoc_shared_b;						  
    if (t_i == 1) prior_shift = 0;
    else prior_shift = sum(l[1:(t_i-1)]);
    for (i in 1:Npat) {
      int mark;
      int start_collect;  // index start of subject-specific ranefs for patient
      mark = 1;
      start_collect = prior_shift + (i - 1) * p[t_i];
      for (m in 1:M) {
        if (size_which_b[m] > 0) {
          int shift;  // num. subject-specific ranefs in prior submodels
          int j_shift; // shift in indexing of which_b vector
          if (m == 1) {
            shift = 0;
            j_shift = 0;
          }
          else {
            shift = sum(pmat[t_i, 1:(m-1)]);
            j_shift = sum(size_which_b[1:(m-1)]);
          }
          for (j in 1:size_which_b[m]) {
            int item_collect;   // subject-specific ranefs to select for current submodel
            item_collect = start_collect + shift + which_b[(j_shift + j)];
            temp[i,mark] = b[item_collect];
            mark = mark + 1;
          }
        }      
      }
    }
    for (i in 1:(quadnodes+1)) {
      start_store = (i - 1) * Npat + 1;
      end_store   = i * Npat;		
      x_assoc_shared_b[start_store:end_store,] = temp;
    }
  return x_assoc_shared_b;
  }
  
  /** 
  * Create a design matrix for a shared fixed + random effects association
  * structure in the joint model
  *
  * @param b Vector of group-specific coefficients
  * @param l An integer array with the number of levels for the factor(s) on 
  *   the RHS of each |
  * @param p An integer array with the number of variables on the LHS of each |
  * @param pmat A matrix with the number variables on the LHS of each | in each
  *   longitudinal submodel. The rows correspond to each |, meaning the separate
  *   equations for each grouping variable, and the columns correspond to each
  *   longitudinal submodel. If subject ID is the only grouping variable then the
  *   matrix will have one row. If the joint model only has one longitudinal 
  *   submodel then the matrix will have one column.
  * @param Npat Integer specifying number of individuals represented 
  *   in vector b
  * @param quadnodes The number of quadrature nodes
  * @param which_b Integer array specifying the indices
  *   of the random effects to use in the association structure
  * @param sum_size_which_b Integer specifying total number of 
  *   random effects that are to be used in the association structure
  * @param size_which_b Integer array specifying number of random effects from
  *   each long submodel that are to be used in the association structure
  * @param t_i Integer specifying the index of the grouping factor that
  *   corresponds to the patient-level
  * @param M An integer specifying the number of longitudinal submodels
  * @return A matrix with the desired random effects represented
  *   in columns, and the individuals on the rows; the matrix is
  *   repeated (quadnodes + 1) times (bounded by rows)
  */  
  matrix make_x_assoc_shared_coef(
    vector b, vector beta, int[] KM, int M, int t_i,
    int[] l, int[] p, int[,] pmat, int Npat, int quadnodes,
    int sum_size_which_coef, int[] size_which_coef,
    int[] which_coef_zindex, int[] which_coef_xindex,
    int[] has_intercept, int[] has_intercept_nob,
    int[] has_intercept_lob, int[] has_intercept_upb,
    real[] gamma_nob, real[] gamma_lob, real[] gamma_upb) {
      
    # in the loops below:
    #   t_i should only really ever equal 1 (since shared_coef association
    #       structure is not allowed if there is more than one clustering level)
    #   i = levels (ie, individuals)
    #   j = indices of the shared random effecs 
    #   m = models
    
    int t_shift;        // skip over group-level coefficients for earlier grouping factors
    int start_store;
    int end_store;	
    matrix[Npat,sum_size_which_coef] temp;
    matrix[(Npat*(quadnodes+1)),sum_size_which_coef] x_assoc_shared_coef;						  
    if (t_i == 1) t_shift = 0;
    else t_shift = sum(l[1:(t_i-1)]);
    for (i in 1:Npat) {
      int mark;        // counter for looping over shared coefficients
      int i_shift;     // skip over group-level coefficients for earlier levels
      mark = 1;
      i_shift = (i - 1) * p[t_i];
      for (m in 1:M) {
        if (size_which_coef[m] > 0) {  # if model has shared coefficients
          int j_shift;  // skip over elements of which_coef_zindex vector that are associated with earlier submodels
          int m_shift;  // skip over individual i's group-level coefficients for earlier submodels
          int shift_nb;
          int shift_lb;
          int shift_ub;
          int shift_beta;
          if (m == 1) {
            j_shift = 0; m_shift = 0; shift_nb = 0;
            shift_lb = 0; shift_ub = 0; shift_beta = 0;
          }
          else {
            j_shift = sum(size_which_coef[1:(m-1)]);
            m_shift = sum(pmat[t_i, 1:(m-1)]);
            shift_nb = sum(has_intercept_nob[1:(m-1)]); 
            shift_lb = sum(has_intercept_lob[1:(m-1)]); 
            shift_ub = sum(has_intercept_upb[1:(m-1)]); 
            shift_beta = sum(KM[1:(m-1)]);
          }
          for (j in 1:size_which_coef[m]) {
            int b_collect;      // group-level coefficients to extract for current i, j, m
            int beta_collect_m; // within-submodel index of fixed effect coefficient to extract
            int beta_collect;   // overall index of fixed effect coefficient to extract
            real coef;
            b_collect = t_shift + i_shift + m_shift + which_coef_zindex[(j_shift + j)];
            beta_collect_m = which_coef_xindex[(j_shift + j)];
            beta_collect = shift_beta + beta_collect_m;
            coef = b[b_collect];  // start with group-level coefficient
            if ((has_intercept[m] == 1) && (beta_collect == 1)) {
              # collect intercept
              if (has_intercept_nob[m] == 1)
                coef = coef + gamma_nob[sum(has_intercept_nob[1:m])];
              else if (has_intercept_lob[m] == 1)
                coef = coef + gamma_lob[sum(has_intercept_lob[1:m])];
              else if (has_intercept_upb[m] == 1)
                coef = coef + gamma_upb[sum(has_intercept_upb[1:m])];
            } 
            else if (has_intercept[m] == 1) {
              # collect fixed effect whilst recognising intercept term 
              # isn't in beta and correcting for that in the indexing
              coef = coef + beta[(beta_collect - 1)];
            }
            else
              coef = coef + beta[beta_collect];
              
            temp[i, mark] = coef;
            mark = mark + 1;  # move to next shared coefficient for individual i
          }
        }      
      }
    }
    # repeat the temp matrix quadnode times (ie, rbind)
    for (i in 1:(quadnodes+1)) {
      start_store = (i - 1) * Npat + 1;
      end_store   = i * Npat;		
      x_assoc_shared_coef[start_store:end_store, ] = temp;
    }
  return x_assoc_shared_coef;
  }  
