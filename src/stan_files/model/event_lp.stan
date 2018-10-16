  vector[len_cpts] lhaz; // log hazard at event and all quadrature times

  vector[nevent] lhaz_epts_event;
  vector[qevent] lhaz_qpts_event;
  vector[qlcens] lhaz_qpts_lcens;
  vector[qrcens] lhaz_qpts_rcens;
  vector[qicens] lhaz_qpts_icenl;
  vector[qicens] lhaz_qpts_icenu;
  vector[qdelay] lhaz_qpts_delay;

  // evaluate log hazard
  if (basehaz_type == 5) { // exponential model
    lhaz = exponential_log_haz(e_eta);
  }
  else if (basehaz_type == 1) { // weibull model
    real shape = e_aux[1];
    lhaz = weibull_log_haz(e_eta, cpts, shape);
  }
  else if (basehaz_type == 6) { // gompertz model
    real scale = e_aux[1];
    lhaz = gompertz_log_haz(e_eta, cpts, scale);
  }
  else if (basehaz_type == 4) { // M-splines, on haz scale
    lhaz = mspline_log_haz(e_eta, basis_cpts, e_aux);
  }
  else if (basehaz_type == 2) { // B-splines, on log haz scale
    lhaz = bspline_log_haz(e_eta, basis_cpts, e_aux);
  }
  else {
    reject("Bug found: invalid baseline hazard.");
  }

  // shift log hazard by a constant (used to center baseline hazard)
  lhaz = lhaz + norm_const;

  // split log hazard vector based on event types
  if (nevent > 0) lhaz_epts_event = lhaz[idx_cpts[1,1]:idx_cpts[1,2]];
  if (qevent > 0) lhaz_qpts_event = lhaz[idx_cpts[2,1]:idx_cpts[2,2]];
  if (qlcens > 0) lhaz_qpts_lcens = lhaz[idx_cpts[3,1]:idx_cpts[3,2]];
  if (qrcens > 0) lhaz_qpts_rcens = lhaz[idx_cpts[4,1]:idx_cpts[4,2]];
  if (qicens > 0) lhaz_qpts_icenl = lhaz[idx_cpts[5,1]:idx_cpts[5,2]];
  if (qicens > 0) lhaz_qpts_icenu = lhaz[idx_cpts[6,1]:idx_cpts[6,2]];
  if (qdelay > 0) lhaz_qpts_delay = lhaz[idx_cpts[7,1]:idx_cpts[7,2]];

  // increment target with log-lik contributions for event submodel
  if (has_weights == 0 && prior_PD == 0) { // unweighted log likelihood
    if (nevent > 0) target +=  lhaz_epts_event;
    if (qevent > 0) target +=  quadrature_log_surv(qwts_event, lhaz_qpts_event);
    if (qlcens > 0) target +=  quadrature_log_cdf (qwts_lcens, lhaz_qpts_lcens, qnodes);
    if (qrcens > 0) target +=  quadrature_log_surv(qwts_rcens, lhaz_qpts_rcens);
    if (qicens > 0) target +=  quadrature_log_cdf2(qwts_icenl, lhaz_qpts_icenl,
                                                   qwts_icenu, lhaz_qpts_icenu, qnodes);
    if (qdelay > 0) target += -quadrature_log_surv(qwts_delay, lhaz_qpts_delay);
  }
  else if (prior_PD == 0) { // weighted log likelihood
    reject("Bug found: weights are not yet implemented.");
  }
