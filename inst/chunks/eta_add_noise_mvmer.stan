	  int N1 = idx_noise[m,1]; // indexing for noise
	  int N2 = idx_noise[m,2]; // indexing for noise
      if      (link[m] == 1) eta_tmp = eta_tmp + log(aux[aux_mark]) + log(noise[1,N1:N2]);
      else if (link[m] == 2) eta_tmp = eta_tmp * aux[aux_mark] .* noise[1,N1:N2];
      else                   eta_tmp = eta_tmp + sqrt(aux[aux_mark]) + sqrt(noise[1,N1:N2]);
