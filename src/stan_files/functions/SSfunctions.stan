vector SS_asymp(vector input, matrix Phi);
vector SS_asympOff(vector input, matrix Phi);
vector SS_asympOrig(vector input, matrix Phi);
vector SS_biexp(vector input, matrix Phi);
vector SS_fol(vector Dose, vector input, matrix Phi);
vector SS_fpl(vector input, matrix Phi);
vector SS_gompertz(vector x, matrix Phi);
vector SS_logis(vector input, matrix Phi);
vector SS_micmen(vector input, matrix Phi);
vector SS_weibull(vector x, matrix Phi);

/* These functions (without the underscores) are all documented in R
   See also Appendix C of Pinheiro and Bates
https://books.google.com/books?id=3TVDAAAAQBAJ&lpg=PR3&dq=Pinheiro%20and%20Bates&pg=PA511#v=onepage&q&f=false
  These functions may be numerically unstable
*/

vector SS_asymp(vector input, matrix Phi) {
  // Phi[,1] = Asym, Phi[,2] = R0, Phi[,3] = lrc
  if (rows(Phi) > 1) {
    vector[rows(Phi)] Asym = Phi[,1];
    return Asym  + (Phi[,2] - Asym) .* exp(-exp(Phi[,3]) .* input);
  }
  else {
    real Asym = Phi[1,1];
    return Asym  + (Phi[1,2] - Asym) * exp(-exp(Phi[1,3]) * input);
  }
}

vector SS_asympOff(vector input, matrix Phi) {
  // Phi[,1] = Asym, Phi[,2] = lrc, Phi[,3] = c0
  if (rows(Phi) > 1) 
    return Phi[ ,1] .* (1 - exp(-exp(Phi[ ,2]) .* (input - Phi[ ,3])));
  else 
    return Phi[1,1]  * (1 - exp(-exp(Phi[1,2])  * (input - Phi[1,3])));
}

vector SS_asympOrig(vector input, matrix Phi) {
  // Phi[,1] = Asym, Phi[,2] = lrc
  if (rows(Phi) > 1)
    return Phi[ ,1] .* (1 - exp(-exp(Phi[ ,2]) .* input));
  else 
    return Phi[1,1]  * (1 - exp(-exp(Phi[1,2])  * input));
}

vector SS_biexp(vector input, matrix Phi) {
  // Phi[,1] = A1, Phi[,2] = lrc1, Phi[,3] = A2, Phi[,4] = lrc2
  if (rows(Phi) > 1)
    return Phi[ ,1] .* exp(-exp(Phi[ ,2]) .* input) + 
           Phi[ ,3] .* exp(-exp(Phi[ ,4]) .* input);
  else          
    return Phi[1,1]  * exp(-exp(Phi[1,2])  * input) + 
           Phi[1,3]  * exp(-exp(Phi[1,4])  * input);
}

vector SS_fol(vector Dose, vector input, matrix Phi) {
  // Phi[,1] = lKe, Phi[,2] = lKa, Phi[,3] = lCl
  int Phi_rows = rows(Phi);
  if (Phi_rows > 1) {
    vector[Phi_rows] lKe = Phi[,1];
    vector[Phi_rows] lKa = Phi[,2];
    vector[Phi_rows] exp_lKe = exp(lKe);
    vector[Phi_rows] exp_lKa = exp(lKa);
    return Dose .* exp(lKe + lKa - Phi[,3]) .* 
      (exp(-exp_lKe .* input) - exp(-exp_lKa .* input)) ./ (exp_lKa - exp_lKe);
  }
  else {
    real lKe = Phi[1,1];
    real lKa = Phi[1,2];
    real exp_lKe = exp(lKe);
    real exp_lKa = exp(lKa);
    return Dose * exp(lKe + lKa - Phi[1,3]) .* 
      (exp(-exp_lKe * input) - exp(-exp_lKa * input)) / (exp_lKa - exp_lKe);
  }
}
  
vector SS_fpl(vector input, matrix Phi) {
  // Phi[,1] = A, Phi[,2] = B, Phi[,3] = xmid, Phi[,4] = scal
  // input is generally data so cannot switch signs
  if (rows(Phi) > 1) {
    vector[rows(Phi)] A = Phi[,1];
    return A + (Phi[,2] - A) ./ (1 + exp((Phi[,3] - input) ./ exp(Phi[,4])));
  }
  else {
    real A = Phi[1,1];
    return A + rep_vector(Phi[1,2] - A, rows(input)) 
      ./ (1 + exp((Phi[1,3] - input) / exp(Phi[1,4])));
  }
}

vector SS_gompertz(vector x, matrix Phi) {
  // Phi[,1] = Asym, Phi[,2] = b2, Phi[,3] = b3
  vector[rows(x)] out;
  if (rows(Phi) > 1) for (i in 1:rows(x))
    out[i] = Phi[i,1] * exp(-Phi[i,2] * Phi[i,3] ^ x[i]);
  else {
    real Asym = Phi[1,1];
    real b2 = Phi[1,2];
    real b3 = Phi[1,3];
    for (i in 1:rows(x)) out[i] = Asym * exp(-b2 * b3 ^ x[i]);
  }
  return out;
}

vector SS_logis(vector input, matrix Phi) {
  // Phi[,1] = Asym, Phi[,2] = xmid, Phi[,3] = scal
  // input is typically data so cannot switch signs of everything
  if (rows(Phi) > 1)
    return Phi[,1] ./ (1 + exp( (Phi[,2] - input) ./ exp(Phi[,3])));
  else 
    return rep_vector(Phi[1,1], rows(input)) ./ 
      (1 + exp( (Phi[1,2] - input) / exp(Phi[1,3])));
}

vector SS_micmen(vector input, matrix Phi) {
  // Phi[,1] = Vm, Phi[,2] = K
  if (rows(Phi) > 1)
    return Phi[ ,1] .* input ./ (Phi[ ,2] + input);
  else
    return Phi[1,1]  * input ./ (Phi[1,2] + input);
}

vector SS_weibull(vector x, matrix Phi) {
  // Phi[,1] = Asym, Phi[,2] = Drop, Phi[,3] = lrc, Phi[,4] = pwr
  vector[rows(x)] out;
  if (rows(Phi) > 1) for (i in 1:rows(x))
    out[i] = Phi[i,1] - Phi[i,2] * exp(-exp(Phi[i,3]) * x[i] ^ Phi[i,4]);
  else {
    real Asym = Phi[1,1];
    real Drop = Phi[1,2];
    real lrc = Phi[1,3];
    real pwr = Phi[1,4];
    for (i in 1:rows(x)) 
      out[i] = Asym - Drop * exp(-exp(lrc) * x[i] ^ pwr);
  }
  return out;
}

matrix reshape_vec(vector x, int Rows, int Cols) {
  matrix[Rows, Cols] out;
  int pos = 1;
  if (rows(x) != Rows * Cols) reject("x is the wrong length");
  for (c in 1:Cols) for (r in 1:Rows) {
    out[r,c] = x[pos];
    pos += 1;
  }
  return out;
}
