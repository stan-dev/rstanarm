  // weights
  int<lower=0,upper=1> has_weights;  // 0 = No, 1 = Yes
  vector[N * has_weights] weights;
  
  // offset
  int<lower=0,upper=1> has_offset;  // 0 = No, 1 = Yes
  vector[N * has_offset] offset;
