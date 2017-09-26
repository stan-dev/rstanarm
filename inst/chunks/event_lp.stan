  vector[nrow_e_Xq] log_basehaz; // log baseline hazard AT event time and quadrature points
  vector[nrow_e_Xq] log_haz_q; // log hazard AT event time and quadrature points
  vector[Nevents] log_haz_etimes; // log hazard AT the event time only
  vector[Npat_times_qnodes] log_haz_qtimes; // log hazard AT the quadrature points
  
  // Log baseline hazard at event and quad times
  if (basehaz_type == 1) 
    log_basehaz = norm_const + log(e_aux[1]) + basehaz_X * (e_aux - 1) + e_gamma[1];
  else log_basehaz = norm_const + basehaz_X * e_aux;	
				  
  // Log hazard at event and quad times
  log_haz_q = log_basehaz + e_eta_q;
  log_haz_etimes = head(log_haz_q, Nevents);
  log_haz_qtimes = tail(log_haz_q, Npat_times_qnodes);

  // Log likelihood for event model
  if (has_weights == 0 && prior_PD == 0) { // unweighted log likelihood
    target += sum(log_haz_etimes) - 
	  dot_product(qwts, exp(log_haz_qtimes));
  } 
  else if (prior_PD == 0) { // weighted log likelihood
    target += dot_product(e_weights, log_haz_etimes) - 
      dot_product(e_weights_rep, qwts .* exp(log_haz_qtimes));
  }				    
