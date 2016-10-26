  if (link_phi > 1) {
    eta_z = eta_z - min(eta_z) + dot_product(zbar, omega);
  }
  else {
    eta_z = eta_z + dot_product(zbar, omega);
  }
