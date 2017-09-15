  vector[nrow_e_Xq] log_basehaz;  // baseline hazard evaluated at quadpoints
  vector[nrow_e_Xq] ll_haz_q;     // log hazard contribution to the log likelihood for the event model at event time and quad points
  vector[Npat] ll_haz_eventtime;  // log hazard contribution to the log likelihood for the event model AT the event time only
  vector[Npat_times_quadnodes] ll_surv_eventtime; // log survival contribution to the log likelihood for the event model AT the event time
  real ll_event;                  // log likelihood for the event model  
  
  // Log baseline hazard at event and quad times
  if (basehaz_type == 1) 
    log_basehaz = log(e_aux[1]) + basehaz_X * (e_aux - 1);
  else log_basehaz = basehaz_X * e_aux;	

  // Log hazard at event and quad times
  ll_haz_q = e_d .* (log_basehaz + e_eta_q);
					  
  // Log hazard contribution to the likelihood
  ll_haz_eventtime = segment(ll_haz_q, 1, Npat);

  // Log survival contribution to the likelihood (by summing over the 
  // quadrature points to get the approximate integral)
  // NB quadweight already incorporates the (b-a)/2 scaling such that the
  // integral is evaluated over limits (a,b) rather than (-1,+1)
  ll_surv_eventtime = quadweight .* 
    exp(segment(ll_haz_q, (Npat + 1), Npat_times_quadnodes));        

  // Log likelihood for event model
  if (has_weights == 0) { // unweighted log likelihood
    ll_event = sum(ll_haz_eventtime) - sum(ll_surv_eventtime);
  } 
  else {  // weighted log likelihood
    ll_event = sum(e_weights .* ll_haz_eventtime) - 
	  sum(e_weights_rep .* ll_surv_eventtime);
  }				    
                           
  // Log-likelihoods for event submodel  
  if (prior_PD == 0) target += ll_event;  
