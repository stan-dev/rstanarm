	vector[len_epts] e_eta_epts = e_eta[idx_cpts[1,1]:idx_cpts[1,2]];
	vector[len_qpts] e_eta_qpts = e_eta[idx_cpts[2,1]:idx_cpts[2,2]];
	vector[len_ipts] e_eta_ipts = e_eta[idx_cpts[3,1]:idx_cpts[3,2]];

  real lhaz_epts = 0;  // summation of log hazard at event times
  real lsur_qpts = 0;  // summation of log surv based on qpts
  real lsur_ipts = 0;  // summation of log surv based on ipts

  // evaluate log hazard and log survival
  if (type == 5) { // exponential model
    if (len_epts > 0) {
      lhaz_epts = sum(e_eta_epts);
    }
    if (len_qpts > 0) {
      lsur_qpts = - dot_product(qwts, exp(e_eta_qpts));
    }
    if (len_ipts > 0) {
      lsur_ipts = - dot_product(iwts, exp(e_eta_ipts));
    }
  }
  else if (type == 1) { // weibull model
    real shape = e_aux[1];
    real log_shape = log(shape);
    if (len_epts > 0) {
      lhaz_epts = (len_epts * log_shape) + (shape - 1) * sum_log_t_events + sum(e_eta_epts);
    }
    if (len_qpts > 0) {
      vector[len_qpts] lhaz_qpts;
      lhaz_qpts = log_shape + (shape - 1) * log_qpts + e_eta_qpts;
      lsur_qpts = - dot_product(qwts, exp(lhaz_qpts));
    }
    if (len_ipts > 0) {
      vector[len_ipts] lhaz_ipts;
      lhaz_ipts = log_shape + (shape - 1) * log_ipts + e_eta_ipts;
      lsur_ipts = - dot_product(iwts, exp(lhaz_ipts));
    }
  }
  else if (type == 6) { // gompertz model
    real scale = e_aux[1];
    if (len_epts > 0) {
      lhaz_epts = scale * sum_t_events + sum(e_eta_epts);
    }
    if (len_qpts > 0) {
      vector[len_qpts] lhaz_qpts;
      lhaz_qpts = scale * qpts + e_eta_qpts;
      lsur_qpts = - dot_product(qwts, exp(lhaz_qpts));
    }
    if (len_ipts > 0) {
      vector[len_ipts] lhaz_ipts;
      lhaz_ipts = scale * ipts + e_eta_ipts;
      lsur_ipts = - dot_product(iwts, exp(lhaz_ipts));
    }
  }
  else if (type == 4) { // M-splines, on haz scale
    if (len_epts > 0)
      lhaz_epts = sum(log(basis_events * e_aux) + e_eta_epts);
    if (len_qpts > 0) {
      vector[len_qpts] lhaz_qpts;
      lhaz_qpts = log(basis_qpts * e_aux) + e_eta_qpts;
      lsur_qpts = - dot_product(qwts, exp(lhaz_qpts));
    }
    if (len_ipts > 0) {
      vector[len_ipts] lhaz_ipts;
      lhaz_ipts = log(basis_ipts * e_aux) + e_eta_ipts;
      lsur_ipts = - dot_product(iwts, exp(lhaz_ipts));
    }
  }
  else if (type == 2) { // B-splines, on log haz scale
    if (len_epts > 0) {
      lhaz_epts = sum(basis_events * e_aux + e_eta_epts);
    }
    if (len_qpts > 0) {
      vector[len_qpts] lhaz_qpts;
      lhaz_qpts = basis_qpts * e_aux + e_eta_qpts;
      lsur_qpts = - dot_product(qwts, exp(lhaz_qpts));
    }
    if (len_ipts > 0) {
      vector[len_ipts] lhaz_ipts;
      lhaz_ipts = basis_ipts * e_aux + e_eta_ipts;
      lsur_ipts = - dot_product(iwts, exp(lhaz_ipts));
    }
  }
  else {
    reject("Bug found: invalid baseline hazard.");
  }

  // increment target with log-lik for event submodel
  if (has_weights == 0 && prior_PD == 0) { // unweighted log likelihood
    target += lhaz_epts + lsur_qpts - lsur_ipts;
  }
  else if (prior_PD == 0) { // weighted log likelihood
    reject("Bug found: weights are not yet implemented.");
  }
