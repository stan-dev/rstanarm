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
    int zeta_mark = 1;
    int rho_mark = 1;
    int z_T_mark = 1;
    int theta_L_mark = 1;

    // each of these is a diagonal block of the implicit Cholesky factor
    for (i in 1:size(p)) { 
      int nc = p[i];
      if (nc == 1) { // "block" is just a standard deviation
        theta_L[theta_L_mark] = tau[i] * scale[i] * dispersion;
        // unlike lme4, theta[theta_L_mark] includes the dispersion term in it
        theta_L_mark = theta_L_mark + 1;
      }
      else { // block is lower-triangular               
        matrix[nc,nc] T_i; 
        real std_dev;
        real T21;
        real trace_T_i = square(tau[i] * scale[i] * dispersion) * nc;
        vector[nc] pi = segment(zeta, zeta_mark, nc); // gamma(zeta | shape, 1)
        pi = pi / sum(pi);                            // thus dirichlet(pi | shape)
        
        // unlike lme4, T_i includes the dispersion term in it
        zeta_mark = zeta_mark + nc;
        std_dev = sqrt(pi[1] * trace_T_i);
        T_i[1,1] = std_dev;
        
        // Put a correlation into T_i[2,1] and scale by std_dev
        std_dev = sqrt(pi[2] * trace_T_i);
        T21 = 2.0 * rho[rho_mark] - 1.0;
        rho_mark = rho_mark + 1;
        T_i[2,2] = std_dev * sqrt(1.0 - square(T21));
        T_i[2,1] = std_dev * T21;
        
        for (r in 2:(nc - 1)) { // scaled onion method to fill T_i
          int rp1 = r + 1;
          vector[r] T_row = segment(z_T, z_T_mark, r);
          real scale_factor = sqrt(rho[rho_mark] / dot_self(T_row)) * std_dev;
          z_T_mark = z_T_mark + r;
          std_dev = sqrt(pi[rp1] * trace_T_i);
          for(c in 1:r) T_i[rp1,c] = T_row[c] * scale_factor;
          T_i[rp1,rp1] = sqrt(1.0 - rho[rho_mark]) * std_dev;
          rho_mark = rho_mark + 1;
        }
        
        // now vech T_i
        for (c in 1:nc) for (r in c:nc) {
          theta_L[theta_L_mark] = T_i[r,c];
          theta_L_mark = theta_L_mark + 1;
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
    int b_mark = 1;
    int theta_L_mark = 1;
    for (i in 1:size(p)) {
      int nc = p[i];
      if (nc == 1) {
        real theta_L_start = theta_L[theta_L_mark];
        for (s in b_mark:(b_mark + l[i] - 1)) 
          b[s] = theta_L_start * z_b[s];
        b_mark = b_mark + l[i];
        theta_L_mark = theta_L_mark + 1;
      }
      else {
        matrix[nc,nc] T_i = rep_matrix(0, nc, nc);
        for (c in 1:nc) {
          T_i[c,c] = theta_L[theta_L_mark];
          theta_L_mark = theta_L_mark + 1;
          for(r in (c+1):nc) {
            T_i[r,c] = theta_L[theta_L_mark];
            theta_L_mark = theta_L_mark + 1;
          }
        }
        for (j in 1:l[i]) {
          vector[nc] temp = T_i * segment(z_b, b_mark, nc);
          b_mark = b_mark - 1;
          for (s in 1:nc) b[b_mark + s] = temp[s];
          b_mark = b_mark + nc + 1;
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
    int pos_reg = 1;
    int pos_rho = 1;
    target += normal_lpdf(z_b | 0, 1);
    target += normal_lpdf(z_T | 0, 1);
    for (i in 1:t) if (p[i] > 1) {
      vector[p[i] - 1] shape1;
      vector[p[i] - 1] shape2;
      real nu = regularization[pos_reg] + 0.5 * (p[i] - 2);
      pos_reg = pos_reg + 1;
      shape1[1] = nu;
      shape2[1] = nu;
      for (j in 2:(p[i]-1)) {
        nu = nu - 0.5;
        shape1[j] = 0.5 * j;
        shape2[j] = nu;
      }
      target += beta_lpdf(rho[pos_rho:(pos_rho + p[i] - 2)] | shape1, shape2);
      pos_rho = pos_rho + p[i] - 1;
    }
    target += gamma_lpdf(zeta | delta, 1);
    target += gamma_lpdf(tau  | shape, 1);
  }
  
  /**
   * Hierarchical shrinkage parameterization
   *
   * @param z_beta A vector of primitive coefficients
   * @param global A real array of positive numbers
   * @param local A vector array of positive numbers
   * @param global_prior_scale A positive real number
   * @param error_scale 1 or sigma in the Gaussian case
   * @return A vector of coefficientes
   */
  vector hs_prior(vector z_beta, real[] global, vector[] local, 
                  real global_prior_scale, real error_scale) {
    vector[rows(z_beta)] lambda;
    int K;
    K = rows(z_beta);
    for (k in 1:K) lambda[k] = local[1][k] * sqrt(local[2][k]);
    return z_beta .* lambda * global[1] * sqrt(global[2]) * 
           global_prior_scale * error_scale;
  }

  /** 
   * Hierarchical shrinkage plus parameterization
   *
   * @param z_beta A vector of primitive coefficients
   * @param global A real array of positive numbers
   * @param local A vector array of positive numbers
   * @param global_prior_scale A positive real number
   * @param error_scale 1 or sigma in the Gaussian case
   * @return A vector of coefficientes
   */
  vector hsplus_prior(vector z_beta, real[] global, vector[] local, 
                      real global_prior_scale, real error_scale) {
    return z_beta .* (local[1] .* sqrt(local[2])) .* 
           (local[3] .* sqrt(local[4])) * global[1] * sqrt(global[2]) * 
           global_prior_scale * error_scale;
  }
  
  /** 
   * Divide a scalar by a vector
   *
   * @param x The scalar in the numerator
   * @param y The vector in the denominator
   * @return An elementwise vector
   */
  vector divide_real_by_vector(real x, vector y) {
    int K = rows(y); 
    vector[K] ret;
    for (k in 1:K) ret[k] = x / y[k];
    return ret;
  }

  /** 
   * Cornish-Fisher expansion for standard normal to Student t
   *
   * See result 26.7.5 of
   * http://people.math.sfu.ca/~cbm/aands/page_949.htm
   *
   * @param z A scalar distributed standard normal
   * @param df A scalar degrees of freedom
   * @return An (approximate) Student t variate with df degrees of freedom
   */
  real CFt(real z, real df) {
    real z2 = square(z);
    real z3 = z2 * z;
    real z5 = z2 * z3;
    real z7 = z2 * z5;
    real z9 = z2 * z7;
    real df2 = square(df);
    real df3 = df2 * df;
    real df4 = df2 * df2;
    return z + (z3 + z) / (4 * df) + (5 * z5 + 16 * z3 + 3 * z) / (96 * df2)
           + (3 * z7 + 19 * z5 + 17 * z3 - 15 * z) / (384 * df3)
           + (79 * z9 + 776 * z7 + 1482 * z5 - 1920 * z3 - 945 * z) / (92160 * df4);
  }

  /** 
   * Return two-dimensional array of group membership
   *
   * @param N An integer indicating the number of observations
   * @param t An integer indicating the number of grouping variables
   * @return An two-dimensional integer array of group membership
   */
  int[,] make_V(int N, int t, int[] v) {
    int V[t,N];
    int pos = 1;
    if (t > 0) for (j in 1:N) for (i in 1:t) {
      V[i,j] = v[pos];
      pos = pos + 1;
    }
    return V;
  }

