  // glmer stuff, see table 3 of
  // https://cran.r-project.org/web/packages/lme4/vignettes/lmer.pdf
  int<lower=0> t;               // num. terms (maybe 0) with a | in the glmer formula
  int<lower=1> p[t];            // num. variables on the LHS of each |
  int<lower=1> l[t];            // num. levels for the factor(s) on the RHS of each |
  int<lower=0> q;               // conceptually equals \sum_{i=1}^t p_i \times l_i
  int<lower=0> len_theta_L;     // length of the theta_L vector
  
  // interaction-specific inputs
  int<lower=1> n_one_way;       // num. one-way random intercepts
  int<lower=1> n_multi_way;     // num. multi-way interaction random intercepts
  int<lower=1> max_way;         // deepest multi-way interaction random intercepts
  int<lower=1> one_way_ix[n_one_way];      // positions in l[t] for one-way;
  int<lower=1> multi_way_ix[n_multi_way];      // positions in l[t] for multi-way;
  int<lower=1> multi_depth[n_multi_way];
  int main_multi_map[n_multi_way, max_way]; // zero-padded map from one to multi
  int<lower=1> len_multi_way_uniq; // number of unique multi-way terms;
  int<lower=1, upper = len_multi_way_uniq> depth_ind[n_multi_way];

  // hyperparameters for glmer stuff; if t > 0 priors are mandatory
  vector<lower=0>[t] shape; 
  vector<lower=0>[t] scale;
  int<lower=0> len_concentration;
  real<lower=0> concentration[len_concentration];
  int<lower=0> len_regularization;
  real<lower=0> regularization[len_regularization];
