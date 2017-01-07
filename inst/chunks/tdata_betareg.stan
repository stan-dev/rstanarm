  if (prior_dist_z <= 2) hs_z = 0;
  else if (prior_dist_z == 3) hs_z = 2;
  else if (prior_dist_z == 4) hs_z = 4;
  if (prior_dist_z == 2) {
    t_any_124_z = 0;
    t_all_124_z = 1;
    for (k in 1:z_dim) {
      if (prior_df_z[k] == 1 || prior_df_z[k] == 2 || prior_df_z[k] == 4)
        t_any_124_z = 1;
      else t_all_124_z = 0;
    }
  }
  else {
    t_any_124_z = 0;
    t_all_124_z = 0;
  }
