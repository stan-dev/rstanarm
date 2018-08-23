  if (link_phi > 1) {
    eta_z += dot_product(zbar, omega) - min(eta_z);
  }
  else {
    eta_z += dot_product(zbar, omega);
  }
