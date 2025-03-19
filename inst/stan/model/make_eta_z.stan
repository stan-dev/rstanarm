 if (family == 4 && z_dim > 0 && link_phi > 0) {
    eta_z = betareg_z * omega;
  }
  else if (family == 4 && z_dim == 0 && has_intercept_z == 1){
    eta_z = rep_vector(0.0, N); 
  }
