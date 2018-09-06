  /**
  * Scale a vector of auxiliary parameters based on prior information
  *
  * @param aux_unscaled A vector, the unscaled auxiliary parameters
  * @param prior_dist Integer, the type of prior distribution
  * @param prior_mean,prior_scale Vectors, the mean and scale 
  *   of the prior distribution
  * @return A vector, corresponding to the scaled auxiliary parameters
  */
  vector make_basehaz_coef(vector aux_unscaled, int prior_dist,
                           vector prior_mean, vector prior_scale) {
    vector[rows(aux_unscaled)] aux;
    if (prior_dist == 0) // none
      aux = aux_unscaled;
    else {
      aux = prior_scale .* aux_unscaled;
      if (prior_dist <= 2) // normal or student_t
        aux += prior_mean;
    }
    return aux;
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
  * Take the linear predictor and collapse across lower level
  * units of the grouping factor clustered within patients, using
  * the function specified by 'grp_assoc'
  *
  * @param eta The linear predictor evaluated for all the lower
  *   level units, having some length greater than N.
  * @param grp_idx An N-by-2 two dimensional array providing the
  *   beginning and ending index of the lower level units in eta that
  *   correspond to patient n (where n = 1,...,N).
  * @param grp_assoc The method for collapsing across the lower
  *   level units; 1=sum, 2=mean, 3=min, 4=max.
  * @return A vector
  */
  vector collapse_within_groups(vector eta, int[,] grp_idx,
                                int grp_assoc) {
    int N = size(grp_idx);
    vector[N] val;
    if (grp_assoc == 1) { // sum of lower level clusters
      for (n in 1:N)
        val[n] = sum(eta[grp_idx[n,1]:grp_idx[n,2]]);
    }
    else if (grp_assoc == 2) { // mean of lower level clusters
      for (n in 1:N)
        val[n] = mean(eta[grp_idx[n,1]:grp_idx[n,2]]);
    }
    else if (grp_assoc == 3) { // min of lower level clusters
      for (n in 1:N)
        val[n] = min(eta[grp_idx[n,1]:grp_idx[n,2]]);
    }
    else if (grp_assoc == 4) { // max of lower level clusters
      for (n in 1:N)
        val[n] = max(eta[grp_idx[n,1]:grp_idx[n,2]]);
    }
    return val;
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
  * @param qnodes The number of quadrature nodes
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
  *   repeated (qnodes + 1) times (bounded by rows)
  */
  matrix make_x_assoc_shared_b(
    vector b, int[] l, int[] p, int[,] pmat, int Npat, int qnodes,
    int[] which_b, int sum_size_which_b, int[] size_which_b, int t_i, int M) {
    int prior_shift; // num. ranefs prior to subject-specific ranefs
    int start_store;
    int end_store;
    matrix[Npat,sum_size_which_b] temp;
    matrix[(Npat*(qnodes+1)),sum_size_which_b] x_assoc_shared_b;
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
            mark += 1;
          }
        }
      }
    }
    for (i in 1:(qnodes+1)) {
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
  * @param qnodes The number of quadrature nodes
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
  *   repeated (qnodes + 1) times (bounded by rows)
  */
  matrix make_x_assoc_shared_coef(
    vector b, vector beta, int[] KM, int M, int t_i,
    int[] l, int[] p, int[,] pmat, int Npat, int qnodes,
    int sum_size_which_coef, int[] size_which_coef,
    int[] which_coef_zindex, int[] which_coef_xindex,
    int[] has_intercept, int[] has_intercept_nob,
    int[] has_intercept_lob, int[] has_intercept_upb,
    real[] gamma_nob, real[] gamma_lob, real[] gamma_upb) {

    // in the loops below:
    //   t_i should only really ever equal 1 (since shared_coef association
    //       structure is not allowed if there is more than one clustering level)
    //   i = levels (ie, individuals)
    //   j = indices of the shared random effecs
    //   m = models

    int t_shift;  // skip over group-level coefficients for earlier grouping factors
    int start_store;
    int end_store;
    matrix[Npat,sum_size_which_coef] temp;
    matrix[(Npat*(qnodes+1)),sum_size_which_coef] x_assoc_shared_coef;
    if (t_i == 1) t_shift = 0;
    else t_shift = sum(l[1:(t_i-1)]);
    for (i in 1:Npat) {
      int mark;    // counter for looping over shared coefficients
      int i_shift; // skip over group-level coefficients for earlier levels
      mark = 1;
      i_shift = (i - 1) * p[t_i];
      for (m in 1:M) {
        if (size_which_coef[m] > 0) {  // if model has shared coefficients
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
              // collect intercept
              if (has_intercept_nob[m] == 1)
                coef += gamma_nob[sum(has_intercept_nob[1:m])];
              else if (has_intercept_lob[m] == 1)
                coef += gamma_lob[sum(has_intercept_lob[1:m])];
              else if (has_intercept_upb[m] == 1)
                coef += gamma_upb[sum(has_intercept_upb[1:m])];
            }
            else if (has_intercept[m] == 1) {
              // collect fixed effect whilst recognising intercept term
              // isn't in beta and correcting for that in the indexing
              coef += beta[(beta_collect - 1)];
            }
            else
              coef += beta[beta_collect];

            temp[i, mark] = coef;
            mark += 1;  // move to next shared coefficient for individual i
          }
        }
      }
    }

    // repeat the temp matrix qnodes times (ie, rbind)
    for (i in 1:(qnodes+1)) {
      start_store = (i - 1) * Npat + 1;
      end_store   = i * Npat;
      x_assoc_shared_coef[start_store:end_store, ] = temp;
    }
  return x_assoc_shared_coef;
  }
