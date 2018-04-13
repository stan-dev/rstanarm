  // population level dimensions
  int<lower=1,upper=3> M; // num submodels with data (limit of 3)
  int<lower=0,upper=1> has_aux[3]; // has auxiliary param
  int<lower=0,upper=1> has_weights; // has observation weights
  int<lower=0,upper=2> resp_type[3]; // 1=real,2=integer,0=none
  int<lower=0,upper=3> intercept_type[3]; // 1=unbounded,2=lob,3=upb,0=none
  int<lower=0> yNobs[3]; // num observations
  int<lower=0> yNeta[3]; // required length of eta
  int<lower=0> yK[3]; // num predictors

  // group level dimensions, for decov prior
  int<lower=0> t;    // num. terms (maybe 0) with a | in the glmer formula
  int<lower=1> p[t]; // num. variables on the LHS of each |
  int<lower=1> l[t]; // num. levels for the factor(s) on the RHS of each |
  int<lower=0> q;    // conceptually equals \sum_{i=1}^t p_i \times l_i
  int<lower=0> len_theta_L; // length of the theta_L vector

  // group level dimensions, for lkj prior

    // group factor 1
    int<lower=0> bN1; // num groups
    int<lower=0> bK1; // total num params
    int<lower=0> bK1_len[3]; // num params in each submodel
    int<lower=0> bK1_idx[3,2]; // beg/end index for group params

    // group factor 2
    int<lower=0> bN2; // num groups
    int<lower=0> bK2; // total num params
    int<lower=0> bK2_len[3]; // num params in each submodel
    int<lower=0> bK2_idx[3,2]; // beg/end index for group params
