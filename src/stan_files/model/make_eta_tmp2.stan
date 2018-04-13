            vector[nrow_y_Xq[sel]] eta_tmp2;
            if (sel == 1) {
              int bMat1_colshift = 0;
              int bMat2_colshift = 0;
              eta_tmp2 = evaluate_eta(y1_xq_eta, y1_z1q_eta, y1_z2q_eta, 
                                      y1_z1q_id_eta, y1_z2q_id_eta, 
                                      yGamma1, yBeta1, bMat1, bMat2, 
                                      bMat1_colshift, bMat2_colshift, 
                                      intercept_type[1]);
            }
            else if (sel == 2) {
              int bMat1_colshift = bK1_len[1];
              int bMat2_colshift = bK2_len[1];
              eta_tmp2 = evaluate_eta(y2_xq_eta, y2_z1q_eta, y2_z2q_eta, 
                                      y2_z1q_id_eta, y2_z2q_id_eta, 
                                      yGamma2, yBeta2, bMat1, bMat2, 
                                      bMat1_colshift, bMat2_colshift, 
                                      intercept_type[2]);
            }
            else if (sel == 3) {
              int bMat1_colshift = sum(bK1_len[1:2]);
              int bMat2_colshift = sum(bK2_len[1:2]);
              eta_tmp2 = evaluate_eta(y3_xq_eta, y3_z1q_eta, y3_z2q_eta, 
                                      y3_z1q_id_eta, y3_z2q_id_eta, 
                                      yGamma3, yBeta3, bMat1, bMat2, 
                                      bMat1_colshift, bMat2_colshift, 
                                      intercept_type[3]);
            }                       
