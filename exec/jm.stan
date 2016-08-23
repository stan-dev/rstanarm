/*
    This file is part of rstanarm.
    Copyright (C) 2015, 2016 Trustees of Columbia University
    
    rstanarm is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    rstanarm is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with rstanarm.  If not, see <http://www.gnu.org/licenses/>.
*/

# JM for a Gaussian, Gamma, or inverse Gaussian longitudinal outcome
functions {
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
  * Segment vector of group-specific coefficients into a vector
  * which contains only the group-specific coefficients for
  * the mth longitudinal submodel
  *
  * @param b Vector whose elements are to be separated
  * @param p An integer array with the number variables on the LHS of each |
  * @param p_y An integer array with the number of on the LHS of | in
  *   each of the longitudinal submodels
  * @param l An integer array with the number of levels for the factor(s) on 
  *   the RHS of each |
  * @param m An integer specifying which longitudinal submodel the 
  *   returned group-specific coefficients should relate to  
  * @param t An integer indicating the number of group-specific terms
  * @return A vectors containing the group-specific coefficients for 
  *   the mth longitudinal submodel
  */ 
  vector segment_b(vector b, int[] p, int[] p_y, int[] l, int m, int t) 
	vector[l[t] * p_y[m]] b_temp;
    int start_collect;
	int shift;
    int start_store;
    int end_store;
	if (m == 1) shift <- 0;
	else shift <- sum(p_y[1:(m-1)])
	for (i in 1:l[t]) {
	  start_collect <- (i - 1) * p[t] + 1 + shift;
      start_store <- (i - 1) * p_y[m] + 1;
      end_store <- (i - 1) * p_y[m] + p_y[m];
 	  b_temp[start_store:end_store] <- segment(b, start_collect, p_y[m]);
	}
	return b_temp;
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


  /** 
   * Apply inverse link function to linear predictor
   *
   * @param eta Linear predictor vector
   * @param link An integer indicating the link function
   * @return A vector, i.e. inverse-link(eta)
   */
  vector linkinv_gauss(vector eta, int link) {
    if (link < 1 || link > 3) reject("Invalid link");
    if (link < 3)  # link = identity or log 
      return eta; # return eta for log link too bc will use lognormal
    else {# link = inverse
      vector[rows(eta)] mu;
      for(n in 1:rows(eta)) mu[n] <- inv(eta[n]); 
      return mu;
    }
  }
  
  /** 
  * Apply inverse link function to linear predictor
  *
  * @param eta Linear predictor vector
  * @param link An integer indicating the link function
  * @return A vector, i.e. inverse-link(eta)
  */
  vector linkinv_gamma(vector eta, int link) {
    if (link < 1 || link > 3) reject("Invalid link");
    if (link == 1)  return eta;
    else if (link == 2) return exp(eta);
    else {
      vector[rows(eta)] mu;
      for(n in 1:rows(eta)) mu[n] <- inv(eta[n]); 
      return mu;
    }
  }
  
  /** 
  * Apply inverse link function to linear predictor
  *
  * @param eta Linear predictor vector
  * @param link An integer indicating the link function
  * @return A vector, i.e. inverse-link(eta)
  */
  vector linkinv_inv_gaussian(vector eta, int link) {
    if (link < 1 || link > 4) reject("Invalid link");
    if (link == 1)  return eta;
    else if (link == 2) return exp(eta);
    else {
      vector[rows(eta)] mu;
      if (link == 3) for( n in 1:rows(eta)) mu[n] <- inv(eta[n]);
      else for (n in 1:rows(eta)) mu[n] <- inv_sqrt(eta[n]);      
      return mu;
    }
  }
  
  /** 
  * Pointwise (pw) log-likelihood vector
  *
  * @param y The integer array corresponding to the outcome variable.
  * @param link An integer indicating the link function
  * @return A vector
  */
  vector pw_gauss(vector y, vector eta, real sigma, int link) {
    vector[rows(eta)] ll;
    if (link < 1 || link > 3) reject("Invalid link");
    if (link == 2) # link = log
      for (n in 1:rows(eta)) ll[n] <- lognormal_log(y[n], eta[n], sigma);
    else { # link = idenity or inverse
      vector[rows(eta)] mu;
      mu <- linkinv_gauss(eta, link);
      for (n in 1:rows(eta)) ll[n] <- normal_log(y[n], mu[n], sigma);
    }
    return ll;
  }
  
  real GammaReg_log(vector y, vector eta, real shape, 
                    int link, real sum_log_y) {
    real ret;
    if (link < 1 || link > 3) reject("Invalid link");
    ret <- rows(y) * (shape * log(shape) - lgamma(shape)) +
      (shape - 1) * sum_log_y;
    if (link == 2)      # link is log
      ret <- ret - shape * sum(eta) - shape * sum(y ./ exp(eta));
    else if (link == 1) # link is identity
      ret <- ret - shape * sum(log(eta)) - shape * sum(y ./ eta);
    else                # link is inverse
      ret <- ret + shape * sum(log(eta)) - shape * dot_product(eta, y);
    return ret;
  }
  
  /** 
  * Pointwise (pw) log-likelihood vector
  *
  * @param y The integer array corresponding to the outcome variable.
  * @param link An integer indicating the link function
  * @return A vector
  */
  vector pw_gamma(vector y, vector eta, real shape, int link) {
    vector[rows(eta)] ll;
    if (link < 1 || link > 3) reject("Invalid link");
    if (link == 3) { # link = inverse
      for (n in 1:rows(eta)) {
        ll[n] <- gamma_log(y[n], shape, shape * eta[n]);
      }
    }
    else if (link == 2) { # link = log
      for (n in 1:rows(eta)) {
        ll[n] <- gamma_log(y[n], shape, shape / exp(eta[n]));
      }
    }
    else { # link = identity
      for (n in 1:rows(eta)) {
        ll[n] <- gamma_log(y[n], shape, shape / eta[n]);
      }
    }
    return ll;
  }
  
  /** 
  * inverse Gaussian log-PDF (for data only, excludes constants)
  *
  * @param y The vector of outcomes
  * @param eta The vector of linear predictors
  * @param lambda A positive scalar nuisance parameter
  * @param link An integer indicating the link function
  * @return A scalar
  */
  real inv_gaussian_log(vector y, vector mu, real lambda, 
                        real sum_log_y, vector sqrt_y) {
    return 0.5 * rows(y) * log(lambda / (2 * pi())) - 
      1.5 * sum_log_y - 
      0.5 * lambda * dot_self( (y - mu) ./ (mu .* sqrt_y) );
  }
  
  /** 
  * Pointwise (pw) log-likelihood vector
  *
  * @param y The integer array corresponding to the outcome variable.
  * @param eta The linear predictors
  * @param lamba A positive scalar nuisance parameter
  * @param link An integer indicating the link function
  * @param log_y A precalculated vector of the log of y
  * @param sqrt_y A precalculated vector of the square root of y
  * @return A vector of log-likelihoods
  */
  vector pw_inv_gaussian(vector y, vector eta, real lambda, 
                         int link, vector log_y, vector sqrt_y) {
    vector[rows(y)] ll;
    vector[rows(y)] mu;
    if (link < 1 || link > 4) reject("Invalid link");
    mu <- linkinv_inv_gaussian(eta, link);
    for (n in 1:rows(y))
      ll[n] <- -0.5 * lambda * square( (y[n] - mu[n]) / (mu[n] * sqrt_y[n]) );
    ll <- ll + 0.5 * log(lambda / (2 * pi())) - 1.5 * log_y;
    return ll;
  }
  
  /** 
  * PRNG for the inverse Gaussian distribution
  *
  * Algorithm from wikipedia 
  *
  * @param mu The expectation
  * @param lambda The dispersion
  * @return A draw from the inverse Gaussian distribution
  */
  real inv_gaussian_rng(real mu, real lambda) {
    real z;
    real y;
    real x;
    real mu2;
    mu2 <- square(mu);
    y <- square(normal_rng(0,1));
    z <- uniform_rng(0,1);
    x <- mu + ( mu2 * y - mu * sqrt(4 * mu * lambda * y + mu2 * square(y)) )
      / (2 * lambda);
    if (z <= (mu / (mu + x))) return x;
    else return mu2 / x;
  }

  /** 
  * test function for csr_matrix_times_vector
  *
  * @param m Integer number of rows
  * @param n Integer number of columns
  * @param w Vector (see reference manual)
  * @param v Integer array (see reference manual)
  * @param u Integer array (see reference manual)
  * @param b Vector that is multiplied from the left by the CSR matrix
  * @return A vector that is the product of the CSR matrix and b
  */
  vector test_csr_matrix_times_vector(int m, int n, vector w, 
                                      int[] v, int[] u, vector b) {
    return csr_matrix_times_vector(m, n, w, v, u, b); 
  }
 
  /** 
  * Create a design matrix for a shared random effects association
  * structure in the joint model
  *
  * @param b Vector of group-specific coefficients
  * @param p An integer array with the number variables on the LHS of each |
  * @param Npat Integer specifying number of individuals represented 
  *   in vector b
  * @param quadnodes The number of quadrature nodes
  * @param which_b Integer array specifying the indices
  *   of the random effects to use in the association structure
  * @param size_which_b Integer specifying total number of 
  *   random effects which are to be used in the association structure
  * @return A matrix with the desired random effects represented
  *   in columns, and the individuals on the rows; the matrix is
  *   repeated (quadnodes + 1) times (bounded by rows)
  */
  matrix make_x_assoc_sh(vector b, int[] p, int Npat, int quadnodes,
	                  int[] which_b, int size_which_b) {
	int start_collect;
	int item_collect;
    matrix[Npat,size_which_b] temp;
	matrix[(Npat*(quadnodes+1)),size_which_b] x_assoc_sh;						  
	for (i in 1:Npat) {
	  start_collect <- (i - 1) * p[1];
	  for (j in 1:size_which_b) {
		item_collect <- start_collect + which_b[j];
		temp[i,j] <- b[item_collect];
      }
	}	
	x_assoc_sh <- temp;
	for (i in 1:quadnodes) 
	  x_assoc_sh <- append_row(x_assoc_sh, temp);
    return x_assoc_sh;
  }
  
}
data {
  // dimensions
  int<lower=1,upper=3> m;  // num. of long. submodels
  int<lower=0> Npat;  // num. individuals (equal to l[1] - 1)
  int<lower=0> y_N[3];  // num. of obs. in each long. submodel
  int<lower=0> y_K[3];  // num. of predictors in each long. submodel
  int<lower=0> sum_y_K; // total num. of predictors across all long submodels
  int<lower=0> e_K;   // num. of predictors in event submodel
  int<lower=0> a_K;   // num. of association parameters
  int<lower=0> quadnodes;  // num. of nodes for Gauss-Kronrod quadrature 
  int<lower=0> Npat_times_quadnodes;
  int<lower=0,upper=m> sum_y_has_intercept; // num. submodels w/ intercept
  
  // data for longitudinal submodel(s)
  int<lower=1> family[3];          // family
  int<lower=1> link[3];            // link function, varies by .stan file
  int<lower=0,upper=1> y_centre[3];         // 1 = yes for centred predictor matrix
  int<lower=0,upper=1> y_has_intercept[3];  // 1 = yes
  int<lower=0,upper=1> y_has_weights[3];    // 1 = Yes
  int<lower=0,upper=1> y_has_offset[3];     // 1 = Yes
  vector[y_N[1]*(m>0)] y1;           // outcome vector
  vector[y_N[2]*(m>1)] y2;               
  vector[y_N[3]*(m>2)] y3;                 
  vector[y_K[1]] y1_xbar[(y_centre[1]>0)*(m>0)];  // predictor means
  vector[y_K[2]] y2_xbar[(y_centre[2]>0)*(m>1)];
  vector[y_K[3]] y3_xbar[(y_centre[3]>0)*(m>2)];
  matrix[y_N[1],y_K[1]] y1_X;   // predictor matrix, possibly centred
  matrix[y_N[2],y_K[2]] y2_X;         
  matrix[y_N[3],y_K[3]] y3_X;           
  vector[y_N[1]*y_has_weights[1]] y1_weights;  
  vector[y_N[2]*y_has_weights[2]] y2_weights;  
  vector[y_N[3]*y_has_weights[3]] y3_weights;  
  vector[y_N[1]*y_has_offset[1]] y1_offset;  
  vector[y_N[2]*y_has_offset[2]] y2_offset;  
  vector[y_N[3]*y_has_offset[3]] y3_offset; 
  int<lower=0> num_non_zero[3];   // number of non-zero elements in the Z matrix
  vector[num_non_zero[1]*(m>0)] w1;  // non-zero elements in the implicit Z matrix
  vector[num_non_zero[2]*(m>1)] w2; 
  vector[num_non_zero[3]*(m>2)] w3; 
  int<lower=0> v1[num_non_zero[1]*(m>0)]; // column indices for w
  int<lower=0> v2[num_non_zero[2]*(m>1)]; 
  int<lower=0> v3[num_non_zero[3]*(m>2)];   
  int<lower=0> u1[(N[1]+1)*(m>0)];  // where the non-zeros start in each row 
  int<lower=0> u2[(N[2]+1)*(m>1)]; 
  int<lower=0> u3[(N[3]+1)*(m>2)]; 
  
  // data for event submodel
  int<lower=0,upper=1> basehaz_weibull;  // weibull baseline hazard
  int<lower=0,upper=1> e_centre;  // 1 = yes for centred predictor matrix
  int<lower=0,upper=1> e_has_intercept;  // 1 = yes
  int<lower=0> nrow_y_Xq;     // num. rows in long. predictor matrix at quad points
  int<lower=0> nrow_e_Xq;   // num. rows in event predictor matrix at quad points
  matrix[nrow_y_Xq,y_K[1]] y1_Xq;               // predictor matrix (long submodel) at quadpoints, possibly centred
  matrix[nrow_y_Xq,y_K[2]] y2_Xq;        
  matrix[nrow_y_Xq,y_K[3]] y3_Xq;               
  matrix[nrow_e_Xq,e_K] e_Xq;         // predictor matrix (event submodel) at quadpoints, possibly centred
  vector[nrow_e_Xq] e_times;          // event times and unstandardised quadrature points
  vector[nrow_e_Xq] e_d;              // event indicator, followed by dummy indicator for quadpoints
  vector[e_K] e_xbar[(e_centre>0)];   // predictor means (event submodel)
  int<lower=0> num_non_zero_Zq[3];    // number of non-zero elements in the Z matrix (at quadpoints)
  vector[num_non_zero_Zq[1]*(m>0)] w1_Zq;  // non-zero elements in the implicit Z matrix (at quadpoints)
  vector[num_non_zero_Zq[2]*(m>1)] w2_Zq; 
  vector[num_non_zero_Zq[3]*(m>2)] w3_Zq;  
  int<lower=0> v1_Zq[num_non_zero_Zq[1]*(m>0)]; // column indices for w (at quadpoints)
  int<lower=0> v2_Zq[num_non_zero_Zq[2]*(m>1)]; 
  int<lower=0> v3_Zq[num_non_zero_Zq[3]*(m>2)]; 
  int<lower=0> u1_Zq[(nrow_y_Xq+1)*(m>0)]; // where the non-zeros start in each row (at quadpoints)
  int<lower=0> u2_Zq[(nrow_y_Xq+1)*(m>1)];
  int<lower=0> u3_Zq[(nrow_y_Xq+1)*(m>2)];
  vector[Npat_times_quadnodes] quadweight_times_half_eventtime;
    
  // data for association structure
  int<lower=0,upper=1> assoc;              // 0 = no jm association structure, 1 = any jm association structure
  int<lower=0,upper=1> has_assoc_ev[m];    // eta value
  int<lower=0,upper=1> has_assoc_es[m];    // eta slope
  int<lower=0,upper=1> has_assoc_cv[m];    // mu value
  int<lower=0,upper=1> has_assoc_cs[m];    // mu slope
  int<lower=0,upper=m> sum_has_assoc_ev;   // num. long submodels linked via eta value
  int<lower=0,upper=m> sum_has_assoc_es;   // num. long submodels linked via eta slope
  int<lower=0,upper=m> sum_has_assoc_cv;   // num. long submodels linked via mu value
  int<lower=0,upper=m> sum_has_assoc_cs;   // num. long submodels linked via mu slope 
  int<lower=1> which_b_for_assoc[size_which_b_for_assoc];  // which random effects are shared 
  int<lower=0> size_which_b_for_assoc;     // num. of shared random effects
  
  // data for random effects model
  int<lower=0> t;     			// not used, set t = 1 for stan_jm
  int<lower=1> p[t];            // total num. random effects across all long. submodels
  int<lower=1> p_y[3];          // num. random effects in each long. submodel
  int<lower=1> l[t];            // num. individuals plus 1
  int<lower=0> q;               // = l * p
  int<lower=0> q_y1;            // = l * p_y1
  int<lower=0> q_y2;            // = l * p_y2
  int<lower=0> q_y3;            // = l * p_y3
  int<lower=0> len_theta_L;     // length of the theta_L vector

  // priors: 0 = none, 1 = normal, 2 = student_t, 3 = hs, 4 = hs_plus
  int<lower=0,upper=4> priorLong_dist;
  int<lower=0,upper=2> priorLong_dist_for_intercept;  
  int<lower=0,upper=4> priorEvent_dist;
  int<lower=0,upper=2> priorEvent_dist_for_intercept;
  int<lower=0,upper=4> priorAssoc_dist;
    
  // hyperparameters for priors, set to 0 if there is no prior
  vector[sum_K] priorLong_mean;
  vector[m]     priorLong_mean_for_intercept;
  vector[e_K]   priorEvent_mean;
  real          priorEvent_mean_for_intercept;
  vector[a_K]   priorAssoc_mean;
  vector<lower=0>[sum_K] priorLong_scale;
  vector<lower=0>[m]     priorLong_scale_for_intercept;
  vector<lower=0>[e_K]   priorEvent_scale;
  real<lower=0>          priorEvent_scale_for_intercept;
  vector<lower=0>[a_K]   priorAssoc_scale;
  vector<lower=0>[sum_K] priorLong_df;
  vector<lower=0>[m]     priorLong_df_for_intercept;
  vector<lower=0>[e_K]   priorEvent_df;
  real<lower=0>          priorEvent_df_for_intercept; 
  vector<lower=0>[a_K]   priorAssoc_df;
  vector<lower=0>[m]     priorLong_scale_for_dispersion;
  real<lower=0>          priorEvent_scale_for_weibull;
 
  // hyperparameters for random effects model
  vector<lower=0>[t] shape; 
  vector<lower=0>[t] scale;
  int<lower=0> len_concentration;
  real<lower=0> concentration[len_concentration];
  int<lower=0> len_regularization;
  real<lower=0> regularization[len_regularization];
   
  // flag indicating whether to draw from the prior
  int<lower=0,upper=1> prior_PD;  // 1 = yes

}
transformed data {
  vector[y_N[1] * (family[1] == 3)] sqrt_y1;
  vector[y_N[2] * (family[2] == 3)] sqrt_y2;
  vector[y_N[3] * (family[3] == 3)] sqrt_y3;
  vector[y_N[1] * (family[1] == 3)] log_y1;
  vector[y_N[2] * (family[2] == 3)] log_y2;
  vector[y_N[3] * (family[3] == 3)] log_y3;
  real sum_log_y1;
  real sum_log_y2;
  real sum_log_y3;
  vector[nrow_e_Xq] e_log_times;  // log of event times and unstandardised quadrature points
  int<lower=0,upper=1> y_t_any_124;
  int<lower=0,upper=1> y_t_all_124;
  int<lower=0,upper=1> e_t_any_124;  
  int<lower=0,upper=1> e_t_all_124;   
  int<lower=0,upper=1> a_t_any_124;  
  int<lower=0,upper=1> a_t_all_124; 
  int<lower=0> y_hs;
  int<lower=0> e_hs;                 
  int<lower=0> a_hs;                 
  int<lower=0> len_z_T;
  int<lower=0> len_var_group;
  int<lower=0> len_rho;
  real<lower=0> delta[len_concentration];
  int<lower=1> pos;
  
  // calculate transformations of outcome
  if      (family[1] == 1) sum_log_y1 <- not_a_number();
  else if (family[1] == 2) sum_log_y1 <- sum(log(y1));
  else {
    for (n in 1:y_N[1]) sqrt_y1[n] <- sqrt(y1[n]);
    log_y1 <- log(y1);
    sum_log_y1 <- sum(log_y1);
  }
  if      (family[2] == 1) sum_log_y2 <- not_a_number();
  else if (family[2] == 2) sum_log_y2 <- sum(log(y2));
  else {
    for (n in 1:y_N[2]) sqrt_y2[n] <- sqrt(y2[n]);
    log_y2 <- log(y2);
    sum_log_y2 <- sum(log_y2);
  }
  if      (family[3] == 1) sum_log_y3 <- not_a_number();
  else if (family[3] == 2) sum_log_y3 <- sum(log(y3));
  else {
    for (n in 1:y_N[3]) sqrt_y3[n] <- sqrt(y3[n]);
    log_y3 <- log(y3);
    sum_log_y3 <- sum(log_y3);
  }
  
  // calculate log of event times and unstandardised quadpoints
  e_log_times <- log(e_times);   
 
  // priors for longitudinal submodels  
  if (priorLong_dist <= 2) y_hs <- 0;
  else if (priorLong_dist == 3) y_hs <- 2;
  else if (priorLong_dist == 4) y_hs <- 4;
  if (priorLong_dist == 2) {
    y_t_any_124 <- 0;
    y_t_all_124 <- 1;
    for (k in 1:K) {
      if (priorLong_df[k] == 1 || priorLong_df[k] == 2 || priorLong_df[k] == 4)
        y_t_any_124 <- 1;
      else y_t_all_124 <- 0;
    }
  }
  else {
    t_any_124 <- 0;
    t_all_124 <- 0;
  }

  // priors for event submodel
  if (priorEvent_dist <= 2) e_hs <- 0;
  else if (priorEvent_dist == 3) e_hs <- 2;
  else if (priorEvent_dist == 4) e_hs <- 4;   
  if (priorEvent_dist == 2) {
    e_t_any_124 <- 0;
    e_t_all_124 <- 1;
    for (k in 1:e_K) {
      if (prior_df[k] == 1 || prior_df[k] == 2 || prior_df[k] == 4)
        e_t_any_124 <- 1;
      else e_t_all_124 <- 0;
    }
  }
  else {
    e_t_any_124 <- 0;
    e_t_all_124 <- 0;
  }
  
  // priors for association parameters
  if (priorAssoc_dist <= 2) a_hs <- 0;
  else if (priorAssoc_dist == 3) a_hs <- 2;
  else if (priorAssoc_dist == 4) a_hs <- 4;   
  if (priorAssoc_dist == 2) {
    a_t_any_124 <- 0;
    a_t_all_124 <- 1;
    for (k in 1:a_K) {
      if (prior_df[k] == 1 || prior_df[k] == 2 || prior_df[k] == 4)
        a_t_any_124 <- 1;
      else a_t_all_124 <- 0;
    }
  }
  else {
    a_t_any_124 <- 0;
    a_t_all_124 <- 0;
  }  
  
  // prior for covariance
  len_z_T <- 0;
  len_var_group <- sum(p) * (t > 0);
  len_rho <- sum(p) - t;
  pos <- 1;
  for (i in 1:t) {
    if (p[i] > 1) {
      for (j in 1:p[i]) {
        delta[pos] <- concentration[j];
        pos <- pos + 1;
      }
    }
    for (j in 3:p[i]) len_z_T <- len_z_T + p[i] - 1;
  }
  
}
parameters {
  // parameters for longitudinal submodel(s)
  real<lower=if_else(family[1] == 1 || link[1] == 2, 
        negative_infinity(), 0)> y1_gamma[y_has_intercept[1]];
  real<lower=if_else(family[2] == 1 || link[2] == 2, 
        negative_infinity(), 0)> y2_gamma[y_has_intercept[2]];  
  real<lower=if_else(family[3] == 1 || link[3] == 2, 
        negative_infinity(), 0)> y3_gamma[y_has_intercept[3]];
  vector[sum_y_K] y_z_beta;                // primative coefs (long submodels)
  vector<lower=0>[m] y_dispersion_unscaled; # interpretation depends on family!
  
  // parameters for event submodel
  real e_gamma[e_has_intercept];          // intercept (event model)
  vector[e_K] e_z_beta;                    // primative coefs (event submodel)
  real<lower=0> weibull_shape_unscaled[basehaz_weibull];  // unscaled weibull shape parameter 
  
  // parameters for association structure
  vector[a_K] a_z_beta;   // primative coefs
    
  // parameters for random effects model
  vector[q] z_b;
  vector[len_z_T] z_T;
  vector<lower=0,upper=1>[len_rho] rho;
  vector<lower=0>[len_concentration] zeta;
  vector<lower=0>[t] tau;  

  // parameters for priors
  real<lower=0> y_global[y_hs];
  vector<lower=0>[sum_y_K] y_local[y_hs];
  real<lower=0> e_global[e_hs];
  vector<lower=0>[e_K] e_local[e_hs];
  real<lower=0> a_global[a_hs];
  vector<lower=0>[e_K] a_local[a_hs];
  
}
transformed parameters {
  // parameters for longitudinal submodel(s)
  vector[sum_y_K] y_beta;
  vector[y_K[1]] y1_beta;
  vector[y_K[2]] y2_beta;
  vector[y_K[3]] y3_beta;
  vector[m] y_dispersion;
  
  // parameters for event submodel
  vector[e_K] e_beta; 
  real weibull_shape[basehaz_weibull]; 
  // parameters for GK quadrature  
  vector[nrow_y_Xq] y_eta_q[m];          // linear predictor (long submodel) evaluated at quadpoints
  vector[nrow_e_Xq] e_eta_q;      // linear predictor (event submodel) evaluated at quadpoints
  vector[nrow_e_Xq] ll_haz_q;     // log hazard contribution to the log likelihood for the event model at event time and quad points
  vector[Npat] ll_haz_eventtime;  // log hazard contribution to the log likelihood for the event model AT the event time only
  vector[Npat_times_quadnodes] ll_haz_quadtime;    // log hazard for the event model AT the quadrature points only
  vector[Npat_times_quadnodes] ll_surv_eventtime;  // log survival contribution to the log likelihood for the event model AT the event time
  real ll_event;                                   // log likelihood for the event model    
  
  // parameters for association structure  
  vector[a_K] a_beta;           

  // parameters for random effects model
  vector[q] b;
  vector[q_y1] b_y1;
  vector[q_y2] b_y2;
  vector[q_y3] b_y3;  
  vector[len_theta_L] theta_L; 

  // parameters for longitudinal submodel(s)
  if      (priorLong_dist == 0) y_beta <- y_z_beta;
  else if (priorLong_dist == 1) y_beta <- y_z_beta .* priorLong_scale + priorLong_mean;
  else if (priorLong_dist == 2) for (k in 1:sum_y_K) {
    real P;
    if (priorLong_df[k] == 1) {
      P <- Phi(y_z_beta[k]);
      y_beta[k] <- tan(pi() * (P - 0.5));
    }
    else if (priorLong_df[k] == 2) {
      P <- Phi(y_z_beta[k]);
      y_beta[k] <- 2 * (P - 0.5) / sqrt(2.0 * P * (1 - P));
    }
    else if (priorLong_df[k] == 4) {
      real q_a;
      P <- Phi(y_z_beta[k]);
      q_a <- sqrt(4.0 * P * (1 - P));
      q_a <- cos(acos(q_a) / 3) / q_a;
      y_beta[k] <- 2 * sqrt(q_a - 1);
      if (P < 0.5) y_beta[k] <- -y_beta[k];
    }
    else y_beta[k] <- y_z_beta[k];
    y_beta[k] <- y_beta[k] * priorLong_scale[k] + priorLong_mean[k];
  }
  else if (priorLong_dist == 3) y_beta <- hs_prior(y_z_beta, y_global, y_local);
  else if (priorLong_dist == 4) y_beta <- hsplus_prior(y_z_beta, y_global, y_local);
  
  if (m == 1) y1_beta <- y_beta
  else if (m > 1) {
    y1_beta <- segment(y_beta, 1, y_K[1])
	y2_beta <- segment(y_beta, (y_K[1] + 1), y_K[2])
	if (m > 2) 
	  y3_beta <- segment(y_beta, (y_K[1] + y_K[2] + 1), y_K[3])
  }
  
  for (j in 1:m) {
    if (priorLong_scale_for_dispersion[j] > 0)
	  y_dispersion[j] <-  priorLong_scale_for_dispersion[j] * y_dispersion_unscaled[j];
    else y_dispersion[j] <- y_dispersion_unscaled[j];
  }
  
  // parameters for event submodel
  if      (priorEvent_dist == 0) e_beta <- e_z_beta;
  else if (priorEvent_dist == 1) e_beta <- e_z_beta .* priorEvent_scale + priorEvent_mean;
  else if (priorEvent_dist == 2) for (k in 1:e_K) {
    real P;
    if (priorEvent_df[k] == 1) {
      P <- Phi(e_z_beta[k]);
      e_beta[k] <- tan(pi() * (P - 0.5));
    }
    else if (priorEvent_df[k] == 2) {
      P <- Phi(e_z_beta[k]);
      e_beta[k] <- 2 * (P - 0.5) / sqrt(2.0 * P * (1 - P));
    }
    else if (priorEvent_df[k] == 4) {
      real q_a;
      P <- Phi(e_z_beta[k]);
      q_a <- sqrt(4.0 * P * (1 - P));
      q_a <- cos(acos(q_a) / 3) / q_a;
      e_beta[k] <- 2 * sqrt(q_a - 1);
      if (P < 0.5) e_beta[k] <- -e_beta[k];
    }
    else e_beta[k] <- e_z_beta[k];
    e_beta[k] <- e_beta[k] * priorEvent_scale[k] + priorEvent_mean[k];
  }
  else if (priorEvent_dist == 3) e_beta <- hs_prior(e_z_beta, e_global, e_local);
  else if (priorEvent_dist == 4) e_beta <- hsplus_prior(e_z_beta, e_global, e_local);

  if (priorEvent_scale_for_weibull > 0)
    weibull_shape <-  priorEvent_scale_for_weibull * weibull_shape_unscaled;
  else weibull_shape <- weibull_shape_unscaled;   
  
   // parameters for association structure
  if      (priorAssoc_dist == 0) a_beta <- a_z_beta;
  else if (priorAssoc_dist == 1) a_beta <- a_z_beta .* priorAssoc_scale + priorAssoc_mean;
  else if (priorAssoc_dist == 2) for (k in 1:a_K) {
    real P;
    if (priorAssoc_df[k] == 1) {
      P <- Phi(a_z_beta[k]);
      a_beta[k] <- tan(pi() * (P - 0.5));
    }
    else if (priorAssoc_df[k] == 2) {
      P <- Phi(a_z_beta[k]);
      a_beta[k] <- 2 * (P - 0.5) / sqrt(2.0 * P * (1 - P));
    }
    else if (priorAssoc_df[k] == 4) {
      real q_a;
      P <- Phi(a_z_beta[k]);
      q_a <- sqrt(4.0 * P * (1 - P));
      q_a <- cos(acos(q_a) / 3) / q_a;
      a_beta[k] <- 2 * sqrt(q_a - 1);
      if (P < 0.5) a_beta[k] <- -a_beta[k];
    }
    else a_beta[k] <- a_z_beta[k];
    a_beta[k] <- a_beta[k] * priorAssoc_scale[k] + priorAssoc_mean[k];
  }
  else if (priorAssoc_dist == 3) a_beta <- hs_prior(a_z_beta, a_global, a_local);
  else if (priorAssoc_dist == 4) a_beta <- hsplus_prior(a_z_beta, a_global, a_local);

  // parameters for random effects model
  if (t > 0) {
    theta_L <- make_theta_L(len_theta_L, p, 
                            dispersion, tau, scale, zeta, rho, z_T);
    b <- make_b(z_b, theta_L, p, l);
	if (m>0) b_y1 <- segment_b(b, p, p_y, l, 1, t)
	if (m>1) b_y2 <- segment_b(b, p, p_y, l, 2, t)
	if (m>2) b_y3 <- segment_b(b, p, p_y, l, 3, t)
  }
  
  //===============
  // GK quadrature
  //===============
 
  // Longitudinal submodel(s): linear predictor at event and quad times
  if (y_K[1] > 0) y_eta_q[1] <- y1_Xq * y1_beta;
  else y_eta_q[1] <- rep_vector(0.0, nrow_Xq);
  if (y_has_offset[1] == 1) y_eta_q[1] <- y_eta_q[1] + y1_offset;
  y_eta_q[1] <- y_eta_q[1] + csr_matrix_times_vector(nrow_Xq, q_y1, w1_Zq, v1_Zq, u1_Zq, b_y1);
  if (y_has_intercept[1] == 1) {
    if (family[1] == 1 || link[1] == 2) y_eta_q[1] <- y_eta_q[1] + y1_gamma[1];
    else y_eta_q[1] <- y_eta_q[1] - min(y_eta_q[1]) + y1_gamma[1];
  }
  else if (y_centre[1] == 1) {
    // correction to eta if model has no intercept (if X is centered)
    y_eta_q[1] <- y_eta_q[1] + dot_product(y1_xbar, y1_beta); 
  }
  if (m > 1) {
    if (y_K[2] > 0) y_eta_q[2] <- y2_Xq * y2_beta;
    else y_eta_q[2] <- rep_vector(0.0, nrow_Xq);
    if (has_offset[2] == 1) y_eta_q[2] <- y_eta_q[2] + y2_offset;
    y_eta_q[2] <- y_eta_q[2] + csr_matrix_times_vector(nrow_Xq, q_y2, w2_Zq, v2_Zq, u2_Zq, b_y2);
	if (y_has_intercept[2] == 1) {
      if (family[2] == 1 || link[2] == 2) y_eta_q[2] <- y_eta_q[2] + y2_gamma[1];
      else y_eta_q[2] <- y_eta_q[2] - min(y_eta_q[2]) + y2_gamma[1];
    }
    else if (y_centre[2] == 1) {
      // correction to eta if model has no intercept (if X is centered)
      y_eta_q[2] <- y_eta_q[2] + dot_product(y2_xbar, y2_beta); 
    }
  }
  if (m > 2) {
    if (y_K[3] > 0) y_eta_q[3] <- y3_Xq * y3_beta;
    else y_eta_q[3] <- rep_vector(0.0, nrow_Xq);
    if (has_offset[3] == 1) y_eta_q[3] <- y_eta_q[3] + y3_offset;
    y_eta_q[3] <- y_eta_q[3] + csr_matrix_times_vector(nrow_Xq, q_y3, w3_Zq, v3_Zq, u3_Zq, b_y3);
 	if (y_has_intercept[3] == 1) {
      if (family[3] == 1 || link[3] == 2) y_eta_q[3] <- y_eta_q[3] + y3_gamma[1];
      else y_eta_q[3] <- y_eta_q[3] - min(y_eta_q[3]) + y3_gamma[1];
    }
    else if (y_centre[3] == 1) {
      // correction to eta if model has no intercept (if X is centered)
      y_eta_q[3] <- y_eta_q[3] + dot_product(y3_xbar, y3_beta); 
    }
  }

  // Event submodel: linear predictor at event and quad times
  if (e_K > 0) e_eta_q <- e_Xq * e_beta;
  else e_eta_q <- rep_vector(0.0, nrow_e_Xq);
  if (e_has_intercept == 1) {
    e_eta_q <- e_eta_q + e_gamma[1];
  }
  else if (e_centre == 1) {
    // correction to eta if model has no intercept (because X is centered)
    e_eta_q <- e_eta_q + dot_product(e_xbar, e_beta); 
  }
  if (assoc == 1) {
    int mark;
	mark <- 1;
    for (j in 1:m) {
      if (has_assoc_ev[j] == 1) {
	    e_eta_q <- e_eta_q + a_beta[mark] * y_eta_q[j];
        mark <- mark + 1;
      }		
      if (has_assoc_es[j] == 1) {
	    // NEED TO CALCULATE SLOPE
	    e_eta_q <- e_eta_q + a_beta[mark] * dydt_eta_q[j]; 
        mark <- mark + 1;
      }				
      if (has_assoc_cv[j] == 1) {
        vector[nrow_y_Xq] y_q;  // expected long. outcome at event and quad times   
	    if (link[j] == 1) y_q <- y_eta_q[j];
        else if (link[j] > 1) {
		  if (family[j] == 1) 
		    y_q <- linkinv_gauss(y_eta_q[j], link[j]);
          else if (family[j] == 2) 
		    y_q <- linkinv_gamma(y_eta_q[j], link[j]);
          else if (family[j] == 3) {
            y_q <- linkinv_inv_gaussian(y_eta_q[j], link[j]);
		}	  
	    e_eta_q <- e_eta_q + a_beta[mark] * y_q; 
        mark <- mark + 1;
      }				
      if (has_assoc_cs[j] == 1) {
	    // NEED TO CALCULATE SLOPE
	    e_eta_q <- e_eta_q + a_beta[mark] * dydt_q[j];
        mark <- mark + 1;
      }				
    }
	if (size_which_b_for_assoc > 0) {
	  int mark_end;  // used to define segment of a_beta
	  matrix[nrow_e_Xq,size_which_b_for_assoc] x_assoc_sh;	  
	  mark_end <- mark + size_which_b_for_assoc; 
	  x_assoc_sh <- make_x_assoc_sh(b, p, Npat, quadnodes, 
	                  which_b_for_assoc, size_which_b_for_assoc);
	  e_eta_q <- e_eta_q + x_assoc_sh * a_beta[mark:mark_end];
    }	
  }

 
  // Calculate log hazard at event times and unstandardised quadrature points 
  // NB assumes Weibull baseline hazard
  ll_haz_q <- e_d .* (log(weibull_shape) + (weibull_shape - 1) * e_log_times + e_eta_q);
					  
  // Partition event times and quad points
  ll_haz_eventtime <- segment(ll_haz_q, 1, Npat);
  ll_haz_quadtime  <- segment(ll_haz_q, (Npat + 1), Npat_times_quadnodes);
                
  // Log survival contribution to the likelihood (by summing over the 
  //   quadrature points to get the approximate integral) 
  ll_surv_eventtime <- quadweight_times_half_eventtime .* exp(ll_haz_quadtime);        

  // Log likelihood for event model 
  ll_event <- sum(ll_haz_eventtime) - sum(ll_surv_eventtime);				  

}
model {
  vector[y_N[1]] y1_eta;                                     
  vector[y_N[2]] y2_eta;                                     
  vector[y_N[3]] y3_eta;                                     
  
  // Longitudinal submodel(s): regression equations
  if (y_K[1] > 0) y1_eta <- y1_X * y1_beta;
  else y1_eta <- rep_vector(0.0, y_N[1]);
  if (y_has_offset[1] == 1) y1_eta <- y1_eta + y1_offset;
  y1_eta <- y1_eta + csr_matrix_times_vector(y_N[1], q_y1, w1, v1, u1, b_y1);
  if (y_has_intercept[1] == 1) {
    if (family[1] == 1 || link[1] == 2) y1_eta <- y1_eta + y1_gamma[1];
    else y1_eta <- y1_eta - min(y1_eta) + y1_gamma[1];
  }
  else if (y_centre[1] == 1) {
    // correction to eta if model has no intercept (if X is centered)
    y1_eta <- y1_eta + dot_product(y1_xbar, y1_beta); 
  }
  if (m > 1) {
    if (y_K[2] > 0) y2_eta <- y2_X * y2_beta;
    else y2_eta <- rep_vector(0.0, y_N[2]);
    if (has_offset[2] == 1) y2_eta <- y2_eta + y2_offset;
    y2_eta <- y2_eta + csr_matrix_times_vector(y_N[2], q_y2, w2, v2, u2, b_y2);
	if (y_has_intercept[2] == 1) {
      if (family[2] == 1 || link[2] == 2) y2_eta <- y2_eta + y2_gamma[1];
      else y2_eta <- y2_eta - min(y2_eta) + y2_gamma[1];
    }
    else if (y_centre[2] == 1) {
      // correction to eta if model has no intercept (if X is centered)
      y2_eta <- y2_eta + dot_product(y2_xbar, y2_beta); 
    }
  }
  if (m > 2) {
    if (y_K[3] > 0) y3_eta <- y3_X * y3_beta;
    else y3_eta <- rep_vector(0.0, y_N[3]);
    if (has_offset[3] == 1) y3_eta <- y3_eta + y3_offset;
    y3_eta <- y3_eta + csr_matrix_times_vector(y_N[3], q_y3, w3, v3, u3, b_y3);
 	if (y_has_intercept[3] == 1) {
      if (family[3] == 1 || link[3] == 2) y3_eta <- y3_eta + y3_gamma[1];
      else y3_eta <- y3_eta - min(y3_eta) + y3_gamma[1];
    }
    else if (y_centre[3] == 1) {
      // correction to eta if model has no intercept (if X is centered)
      y3_eta <- y3_eta + dot_product(y3_xbar, y3_beta); 
    }
  }
 
  // Log-likelihood for longitudinal submodel(s)
  if (has_weights == 0 && prior_PD == 0) { # unweighted log-likelihoods
    if (family[1] == 1) {
      if (link[1] == 1)      y1 ~ normal(y1_eta, y_dispersion[1]);
      else if (link[1] == 2) y1 ~ lognormal(y1_eta, y_dispersion[1]);
      else y1 ~ normal(divide_real_by_vector(1, y1_eta), y_dispersion[1]);
      // divide_real_by_vector() is defined in common_functions.stan
    }
    else if (family[1] == 2) {
      y1 ~ GammaReg(eta, y_dispersion[1], link[1], sum_log_y1);
    }
    else {
      y1 ~ inv_gaussian(linkinv_inv_gaussian(y1_eta, link[1]), 
                       y_dispersion[1], sum_log_y1, sqrt_y1);
    }
    if (m > 1) {	
      if (family[2] == 1) {
        if (link[2] == 1)      y2 ~ normal(y2_eta, y_dispersion[2]);
        else if (link[2] == 2) y2 ~ lognormal(y2_eta, y_dispersion[2]);
        else y2 ~ normal(divide_real_by_vector(1, y2_eta), y_dispersion[2]);
      }
      else if (family[2] == 2) {
        y2 ~ GammaReg(eta, y_dispersion[2], link[2], sum_log_y2);
      }
      else {
        y2 ~ inv_gaussian(linkinv_inv_gaussian(y2_eta, link[2]), 
                         y_dispersion[2], sum_log_y2, sqrt_y2);
      }
    }
	if (m > 2) {
      if (family[3] == 1) {
        if (link[3] == 1)      y3 ~ normal(y3_eta, y_dispersion[3]);
        else if (link[3] == 2) y3 ~ lognormal(y3_eta, y_dispersion[3]);
        else y3 ~ normal(divide_real_by_vector(1, y3_eta), y_dispersion[3]);
      }
      else if (family[3] == 2) {
        y3 ~ GammaReg(eta, y_dispersion[3], link[3], sum_log_y3);
      }
      else {
        y3 ~ inv_gaussian(linkinv_inv_gaussian(y3_eta, link[3]), 
                         y_dispersion[3], sum_log_y3, sqrt_y3);
      }
    }	
  }
  else if (prior_PD == 0) { # weighted log-likelihoods
    // WEIGHTS NOT YET HANDLED
  }
                           
  // Log-likelihood for event submodel   
  if (prior_PD == 0) increment_log_prob(ll_event); 

  // Log-priors for coefficients in longitudinal submodel(s)
  if (priorLong_dist == 1) y_z_beta ~ normal(0, 1);
  else if (priorLong_dist == 2) {
    if (t_all_124) y_z_beta ~ normal(0,1);
    else if (t_any_124) for (k in 1:sum_y_K) {
      if (priorLong_df[k] == 1 || priorLong_df[k] == 2 || priorLong_df[k] == 4)
        y_z_beta[k] ~ normal(0,1);
      else y_z_beta[k] ~ student_t(priorLong_df[k], 0, 1);
    }
    else y_z_beta ~ student_t(priorLong_df, 0, 1);
  }
  else if (priorLong_dist == 3) { // hs
    y_z_beta ~ normal(0,1);
    y_local[1] ~ normal(0,1);
    y_local[2] ~ inv_gamma(0.5 * priorLong_df, 0.5 * priorLong_df);
    y_global[1] ~ normal(0,1);
    y_global[2] ~ inv_gamma(0.5, 0.5);
  }
  else if (priorLong_dist == 4) { // hs+
    y_z_beta ~ normal(0,1);
    y_local[1] ~ normal(0,1);
    y_local[2] ~ inv_gamma(0.5 * priorLong_df, 0.5 * priorLong_df);
    y_local[3] ~ normal(0,1);
    // unorthodox useage of priorLong_scale as another df hyperparameter
    y_local[4] ~ inv_gamma(0.5 * priorLong_scale, 0.5 * priorLong_scale);
    y_global[1] ~ normal(0,1);
    y_global[2] ~ inv_gamma(0.5, 0.5);
  }
  // else priorLong_dist is 0 and nothing is added 
  
  // Log-prior for intercept in longitudinal submodel(s)
  if (y_has_intercept[1] == 1) {
    if (priorLong_dist_for_intercept == 1)  // normal
      y1_gamma ~ normal(priorLong_mean_for_intercept[1], priorLong_scale_for_intercept[1]);
    else if (priorLong_dist_for_intercept == 2)  // student_t
      y1_gamma ~ student_t(priorLong_df_for_intercept[1], priorLong_mean_for_intercept[1], 
                        priorLong_scale_for_intercept[1]);
    // else priorLong_dist is 0 and nothing is added 
  }
  if (y_has_intercept[2] == 1) {
    if (priorLong_dist_for_intercept == 1)  // normal
      y2_gamma ~ normal(priorLong_mean_for_intercept[2], priorLong_scale_for_intercept[2]);
    else if (priorLong_dist_for_intercept == 2)  // student_t
      y2_gamma ~ student_t(priorLong_df_for_intercept[2], priorLong_mean_for_intercept[2], 
                        priorLong_scale_for_intercept[2]);
  }
  if (y_has_intercept[3] == 1) {
    if (priorLong_dist_for_intercept == 1)  // normal
      y3_gamma ~ normal(priorLong_mean_for_intercept[3], priorLong_scale_for_intercept[3]);
    else if (priorLong_dist_for_intercept == 2)  // student_t
      y3_gamma ~ student_t(priorLong_df_for_intercept[3], priorLong_mean_for_intercept[3], 
                        priorLong_scale_for_intercept[3]);
  }  
      
  // Log-priors for coefficients in event submodel
  if (priorEvent_dist == 1) e_z_beta ~ normal(0, 1);
  else if (priorEvent_dist == 2) {
    if (e_t_all_124) e_z_beta ~ normal(0,1);
    else if (e_t_any_124) for (k in 1:e_K) {
      if (priorEvent_df[k] == 1 || priorEvent_df[k] == 2 || priorEvent_df[k] == 4)
        e_z_beta[k] ~ normal(0,1);
      else e_z_beta[k] ~ student_t(priorEvent_df[k], 0, 1);
    }
    else e_z_beta ~ student_t(priorEvent_df, 0, 1);
  }
  else if (priorEvent_dist == 3) { // hs
    e_z_beta ~ normal(0,1);
    e_local[1] ~ normal(0,1);
    e_local[2] ~ inv_gamma(0.5 * priorEvent_df, 0.5 * priorEvent_df);
    e_global[1] ~ normal(0,1);
    e_global[2] ~ inv_gamma(0.5, 0.5);
  }
  else if (priorEvent_dist == 4) { // hs+
    e_z_beta ~ normal(0,1);
    e_local[1] ~ normal(0,1);
    e_local[2] ~ inv_gamma(0.5 * priorEvent_df, 0.5 * priorEvent_df);
    e_local[3] ~ normal(0,1);
    // unorthodox useage of priorEvent_scale as another df hyperparameter
    e_local[4] ~ inv_gamma(0.5 * priorEvent_scale, 0.5 * priorEvent_scale);
    e_global[1] ~ normal(0,1);
    e_global[2] ~ inv_gamma(0.5, 0.5);
  }
  // else priorEvent_dist is 0 and nothing is added 
  
  // Log-prior for intercept in event submodel
  if (e_has_intercept == 1) {
    if (priorEvent_dist_for_intercept == 1)  // normal
      e_gamma ~ normal(priorEvent_mean_for_intercept, priorEvent_scale_for_intercept);
    else if (priorEvent_dist_for_intercept == 2)  // student_t
      e_gamma ~ student_t(priorEvent_df_for_intercept, priorEvent_mean_for_intercept, 
                        priorEvent_scale_for_intercept);
    // else priorEvent_dist is 0 and nothing is added 
  }    

  // Log-priors for association parameters
  if (priorAssoc_dist == 1) a_z_beta ~ normal(0, 1);
  else if (priorAssoc_dist == 2) {
    if (a_t_all_124) a_z_beta ~ normal(0,1);
    else if (a_t_any_124) for (k in 1:a_K) {
      if (priorAssoc_df[k] == 1 || priorAssoc_df[k] == 2 || priorAssoc_df[k] == 4)
        a_z_beta[k] ~ normal(0,1);
      else a_z_beta[k] ~ student_t(priorAssoc_df[k], 0, 1);
    }
    else a_z_beta ~ student_t(priorAssoc_df, 0, 1);
  }
  else if (priorAssoc_dist == 3) { // hs
    a_z_beta ~ normal(0,1);
    a_local[1] ~ normal(0,1);
    a_local[2] ~ inv_gamma(0.5 * priorAssoc_df, 0.5 * priorAssoc_df);
    a_global[1] ~ normal(0,1);
    a_global[2] ~ inv_gamma(0.5, 0.5);
  }
  else if (priorAssoc_dist == 4) { // hs+
    a_z_beta ~ normal(0,1);
    a_local[1] ~ normal(0,1);
    a_local[2] ~ inv_gamma(0.5 * priorAssoc_df, 0.5 * priorAssoc_df);
    a_local[3] ~ normal(0,1);
    // unorthodox useage of priorAssoc_scale as another df hyperparameter
    a_local[4] ~ inv_gamma(0.5 * priorAssoc_scale, 0.5 * priorAssoc_scale);
    a_global[1] ~ normal(0,1);
    a_global[2] ~ inv_gamma(0.5, 0.5);
  }
  // else priorAssoc_dist is 0 and nothing is added   
    
  // Log-prior for scale(s)
  for (j in 1:m)
	if (priorLong_scale_for_dispersion[j] > 0) y_dispersion_unscaled[j] ~ cauchy(0, 1);
    
  // Log-prior for Weibull shape
  if (basehaz_weibull == 1)
    weibull_shape_unscaled ~ cauchy(0,1);   
  
  // Prior for random effects model
  if (t > 0) decov_lp(z_b, z_T, rho, zeta, tau, 
                      regularization, delta, shape, t, p);
   
}
/*
generated quantities {
  real alpha[has_intercept];
  real mean_PPD;
  mean_PPD <- 0;
  if (has_intercept == 1)
    alpha[1] <- gamma[1] - dot_product(xbar, beta);
  {
  vector[N] eta;  // linear predictor
  if (K > 0) eta <- X * beta;
  else eta <- rep_vector(0.0, N);
  if (has_offset == 1) eta <- eta + offset;
    if (t > 0) eta <- eta + csr_matrix_times_vector(N, q, w, v, u, b);
    if (has_intercept == 1) {
      if (family == 1 || link == 2) eta <- eta + gamma[1];
      else {
        real min_eta;
        min_eta <- min(eta);
        alpha[1] <- alpha[1] - min_eta;
        eta <- eta - min_eta + gamma[1];
      }
    }
    else {
      // correction to eta if model has no intercept (because X is centered)
      eta <- eta + dot_product(xbar, beta); 
    }
    
    eta <- segment(eta, 1, N);
    
    if (family == 1) {
      if (link > 1) eta <- linkinv_gauss(eta, link);
      for (n in 1:N) mean_PPD <- mean_PPD + normal_rng(eta[n], dispersion);
    }
    else if (family == 2) {
      if (link > 1) eta <- linkinv_gamma(eta, link);
      for (n in 1:N) mean_PPD <- mean_PPD + gamma_rng(dispersion, dispersion / eta[n]);
    }
    else if (family == 3) {
      if (link > 1) eta <- linkinv_inv_gaussian(eta, link);
      for (n in 1:N) mean_PPD <- mean_PPD + inv_gaussian_rng(eta[n], dispersion);
    }
    mean_PPD <- mean_PPD / N;
  }
}
*/
