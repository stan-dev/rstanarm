  // correction to eta if model has no intercept (because X is centered)
  eta += dot_product(xbar, beta); 
