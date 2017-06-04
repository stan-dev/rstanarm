    vector[sum(nrow_y_Xq)*(assoc_uses[1]>0)]     y_eta_q;     // linear predictor (all long submodels) evaluated at quadpoints
    vector[sum(nrow_y_Xq)*(assoc_uses[2]>0)]     y_eta_q_eps; // linear predictor (all long submodels) evaluated at quadpoints plus time shift of epsilon
    vector[sum(nrow_y_Xq_auc)*(assoc_uses[3]>0)] y_eta_q_auc; // linear predictor (all long submodels) evaluated at auc quadpoints
    // mark tracks indexing within a_beta vector, which is the 
    // vector of association parameters
    int mark = 0;
    // mark2 tracks indexing within a_K_data vector, which is the 
    // vector specifying the number of columns used for each possible 
    // type of association term by data interaction
    int mark2 = 0;
    // mark3 tracks indexing within size_which_interactions vector
    int mark3 = 0;
