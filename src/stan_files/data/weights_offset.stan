  // weights
  int<lower=0,upper=1> has_weights;  // 0 = No, 1 = Yes
  vector[has_weights ? N : 0] weights;
  
  // offset
  int<lower=0,upper=1> has_offset;  // 0 = No, 1 = Yes
  vector[has_offset ? N : 0] offset;
