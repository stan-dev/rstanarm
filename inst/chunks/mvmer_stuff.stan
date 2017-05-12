  int<lower=0> pmat[t,M]; // num. group-specific terms for each grouping factor (t) in each submodel (M)
  int<lower=0> qmat[t,M]; // = l * pmat --> num. group-specific coefs for each grouping factor in each submodel
  int<lower=0> q1[t];     // = rowsums(qmat) --> num. group-specific coefs for each grouping factor
  int<lower=0> q2[M];     // = colsums(qmat) --> num. group-specific coefs for each submodel
