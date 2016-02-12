  /* for multiple .stan files */
  
  /** 
   * Create group-specific block-diagonal Cholesky factor, see section 2 of
   * https://cran.r-project.org/web/packages/lme4/vignettes/lmer.pdf
   * @param len_theta_L An integer indicating the length of returned vector, 
   *   which lme4 denotes as m
   * @param p An integer array with the number variables on the LHS of each |
   * @param dispersion Scalar standard deviation of the errors, calles sigma by lme4
   * @param tau Vector of scale parameters whose squares are proportional to the 
   *   traces of the relative covariance matrices of the group-specific terms
   * @param scale Vector of prior scales that are multiplied by elements of tau
   * @param zeta Vector of positive parameters that are normalized into simplexes
   *   and multiplied by the trace of the covariance matrix to produce variances
   * @param rho Vector of radii in the onion method for creating Cholesky factors
   * @param z_T Vector used in the onion method for creating Cholesky factors
   * @return A vector that corresponds to theta in lme4
   */
  vector make_theta_L(int len_theta_L, int[] p, real dispersion,
                      vector tau, vector scale, vector zeta,
                      vector rho, vector z_T) {
    vector[len_theta_L] theta_L;
    int zeta_mark;
    int rho_mark;
    int z_T_mark;
    int theta_L_mark;
    zeta_mark <- 1;
    rho_mark <- 1;
    z_T_mark <- 1;
    theta_L_mark <- 1;
    
    // each of these is a diagonal block of the implicit Cholesky factor
    for (i in 1:size(p)) { 
      int nc;
      nc <- p[i];
      if (nc == 1) { // "block" is just a standard deviation
        theta_L[theta_L_mark] <- tau[i] * scale[i] * dispersion;
        // unlike lme4, theta[theta_L_mark] includes the dispersion term in it
        theta_L_mark <- theta_L_mark + 1;
      }
      else { // block is lower-triangular               
        matrix[nc,nc] T_i; 
        real trace_T_i;
        vector[nc] pi; // variance = proportion of trace_T_i
        real std_dev;
        real T21;
        
        trace_T_i <- square(tau[i] * scale[i] * dispersion) * nc;
        // unlike lme4, T_i includes the dispersion term in it
        pi <- segment(zeta, zeta_mark, nc); // zeta ~ gamma(shape, 1)
        pi <- pi / sum(pi);                 // thus pi ~ dirichlet(shape)
        zeta_mark <- zeta_mark + nc;
        std_dev <- sqrt(pi[1] * trace_T_i);
        T_i[1,1] <- std_dev;
        
        // Put a correlation into T_i[2,1] and scale by std_dev
        std_dev <- sqrt(pi[2] * trace_T_i);
        T21 <- 2.0 * rho[rho_mark] - 1.0;
        rho_mark <- rho_mark + 1;
        T_i[2,2] <- std_dev * sqrt(1.0 - square(T21));
        T_i[2,1] <- std_dev * T21;
        
        for (r in 2:(nc - 1)) { // scaled onion method to fill T_i
          int rp1;
          vector[r] T_row;
          real scale_factor;
          T_row <- segment(z_T, z_T_mark, r);
          z_T_mark <- z_T_mark + r;
          rp1 <- r + 1;
          std_dev <- sqrt(pi[rp1] * trace_T_i);
          scale_factor <- sqrt(rho[rho_mark] / dot_self(T_row)) * std_dev;
          for(c in 1:r) T_i[rp1,c] <- T_row[c] * scale_factor;
          T_i[rp1,rp1] <- sqrt(1.0 - rho[rho_mark]) * std_dev;
          rho_mark <- rho_mark + 1;
        }
        
        // now vech T_i
        for (c in 1:nc) for (r in c:nc) {
          theta_L[theta_L_mark] <- T_i[r,c];
          theta_L_mark <- theta_L_mark + 1;
        }
      }
    }
    return theta_L;
  }
  
  /** 
  * Create group-specific coefficients, see section 2 of
  * https://cran.r-project.org/web/packages/lme4/vignettes/lmer.pdf
  *
  * @param z_b Vector whose elements are iid normal(0,sigma) a priori
  * @param theta Vector with covariance parameters as defined in lme4
  * @param p An integer array with the number variables on the LHS of each |
  * @param l An integer array with the number of levels for the factor(s) on 
  *   the RHS of each |
  * @return A vector of group-specific coefficients
  */
  vector make_b(vector z_b, vector theta_L, int[] p, int[] l) {
    vector[rows(z_b)] b;
    int b_mark;
    int theta_L_mark;
    b_mark <- 1;
    theta_L_mark <- 1;
    for (i in 1:size(p)) {
      int nc;
      nc <- p[i];
      if (nc == 1) {
        real theta_L_start;
        theta_L_start <- theta_L[theta_L_mark];
        for (s in b_mark:(b_mark + l[i] - 1)) 
          b[s] <- theta_L_start * z_b[s];
        b_mark <- b_mark + l[i];
        theta_L_mark <- theta_L_mark + 1;
      }
      else {
        matrix[nc,nc] T_i;
        T_i <- rep_matrix(0, nc, nc);
        for (c in 1:nc) {
          T_i[c,c] <- theta_L[theta_L_mark];
          theta_L_mark <- theta_L_mark + 1;
          for(r in (c+1):nc) {
            T_i[r,c] <- theta_L[theta_L_mark];
            theta_L_mark <- theta_L_mark + 1;
          }
        }
        for (j in 1:l[i]) {
          vector[nc] temp;
          temp <- T_i * segment(z_b, b_mark, nc);
          b_mark <- b_mark - 1;
          for (s in 1:nc) b[b_mark + s] <- temp[s];
          b_mark <- b_mark + nc + 1;
        }
      }
    }
    return b;
  }

  /** 
   * Prior on group-specific parameters
   *
   * @param z_b A vector of primitive coefficients
   * @param z_T A vector of primitives for the unit vectors in the onion method
   * @param rho A vector radii for the onion method
   * @param zeta A vector of primitives for the simplexes
   * @param tau A vector of scale parameters
   * @param regularization A real array of LKJ hyperparameters
   * @param delta A real array of concentration paramters
   * @param shape A vector of shape parameters
   * @param t An integer indicating the number of group-specific terms
   * @param p An integer array with the number variables on the LHS of each |
   * @return nothing
   */
  void decov_lp(vector z_b, vector z_T, vector rho, vector zeta, vector tau,
                real[] regularization, real[] delta, vector shape,
                int t, int[] p) {
    int pos_reg;
    int pos_rho;
    z_b ~ normal(0,1);
    z_T ~ normal(0,1);
    pos_reg <- 1;
    pos_rho <- 1;
    for (i in 1:t) if (p[i] > 1) {
      vector[p[i] - 1] shape1;
      vector[p[i] - 1] shape2;
      real nu;
      nu <- regularization[pos_reg] + 0.5 * (p[i] - 2);
      pos_reg <- pos_reg + 1;
      shape1[1] <- nu;
      shape2[1] <- nu;
      for (j in 2:(p[i]-1)) {
        nu <- nu - 0.5;
        shape1[j] <- 0.5 * j;
        shape2[j] <- nu;
      }
      rho[pos_rho:(pos_rho + p[i] - 2)] ~ beta(shape1,shape2);
      pos_rho <- pos_rho + p[i] - 1;
    }
    zeta ~ gamma(delta, 1);
    tau ~ gamma(shape, 1);
  }
  
  /** 
   * Elementwise square root
   *
   * @param y A vector of non-negative numbers
   * @return A vector of square roots
   */
  vector sqrt_vec(vector y) {
    vector[rows(y)] out;
    for (i in 1:rows(y)) out[i] <- sqrt(out[i]);
    return out;
  }

  /** 
   * Hierarchical shrinkage parameterization
   *
   * @param z_beta A vector of primitive coefficients
   * @param global A real array of positive numbers
   * @param local A vector array of positive numbers
   * @return A vector of coefficientes
   */
  vector hs_prior(vector z_beta, real[] global, vector[] local) {
    vector[rows(z_beta)] lambda;
    int K;
    K <- rows(z_beta);
    for (k in 1:K) lambda[k] <- local[1][k] * sqrt(local[2][k]);
    return z_beta .* lambda * global[1] * sqrt(global[2]);
  }

  /** 
   * Hierarchical shrinkage plus parameterization
   *
   * @param z_beta A vector of primitive coefficients
   * @param global A real array of positive numbers
   * @param local A vector array of positive numbers
   * @return A vector of coefficientes
   */
  vector hsplus_prior(vector z_beta, real[] global, vector[] local) {
    vector[rows(z_beta)] lambda;
    vector[rows(z_beta)] lambda_plus;
    int K;
    K <- rows(z_beta);
    for (k in 1:K) {
      lambda[k] <- local[1][k] * sqrt(local[2][k]);
      lambda_plus[k] <- local[3][k] * sqrt(local[4][k]);
    }
    return z_beta .* lambda .* lambda_plus * global[1] * sqrt(global[2]);
  }
  
  /** 
   * Divide a scalar by a vector
   *
   * @param x The scalar in the numerator
   * @param y The vector in the denominator
   * @return An elementwise vector
   */
  vector divide_real_by_vector(real x, vector y) {
    vector[rows(y)] ret;
    for (n in 1:rows(y)) ret[n] <- x / y[n];
    return ret;
  }

