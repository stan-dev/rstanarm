  // correction to eta if model has no intercept (because X is centered)
  eta = eta + dot_product(xbar, beta); 
