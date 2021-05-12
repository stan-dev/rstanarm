vector SS_asymp(vector input, matrix Phi_);
vector SS_asympOff(vector input, matrix Phi_);
vector SS_asympOrig(vector input, matrix Phi_);
vector SS_biexp(vector input, matrix Phi_);
vector SS_fol(vector Dose, vector input, matrix Phi_);
vector SS_fpl(vector input, matrix Phi_);
vector SS_gompertz(vector x, matrix Phi_);
vector SS_logis(vector input, matrix Phi_);
vector SS_micmen(vector input, matrix Phi_);
vector SS_weibull(vector x, matrix Phi_);

/* These functions (without the underscores) are all documented in R
   See also Appendix C of Pinheiro and Bates
https://books.google.com/books?id=3TVDAAAAQBAJ&lpg=PR3&dq=Pinheiro%20and%20Bates&pg=PA511#v=onepage&q&f=false
  These functions may be numerically unstable
*/

vector SS_asymp(vector input, matrix Phi_) {
  // Phi_[,1] = Asym, Phi_[,2] = R0, Phi_[,3] = lrc
  if (rows(Phi_) > 1) {
    vector[rows(Phi_)] Asym = Phi_[,1];
    return Asym  + (Phi_[,2] - Asym) .* exp(-exp(Phi_[,3]) .* input);
  }
  else {
    real Asym = Phi_[1,1];
    return Asym  + (Phi_[1,2] - Asym) * exp(-exp(Phi_[1,3]) * input);
  }
}

vector SS_asympOff(vector input, matrix Phi_) {
  // Phi_[,1] = Asym, Phi_[,2] = lrc, Phi_[,3] = c0
  if (rows(Phi_) > 1) 
    return Phi_[ ,1] .* (1 - exp(-exp(Phi_[ ,2]) .* (input - Phi_[ ,3])));
  else 
    return Phi_[1,1]  * (1 - exp(-exp(Phi_[1,2])  * (input - Phi_[1,3])));
}

vector SS_asympOrig(vector input, matrix Phi_) {
  // Phi_[,1] = Asym, Phi_[,2] = lrc
  if (rows(Phi_) > 1)
    return Phi_[ ,1] .* (1 - exp(-exp(Phi_[ ,2]) .* input));
  else 
    return Phi_[1,1]  * (1 - exp(-exp(Phi_[1,2])  * input));
}

vector SS_biexp(vector input, matrix Phi_) {
  // Phi_[,1] = A1, Phi_[,2] = lrc1, Phi_[,3] = A2, Phi_[,4] = lrc2
  if (rows(Phi_) > 1)
    return Phi_[ ,1] .* exp(-exp(Phi_[ ,2]) .* input) + 
           Phi_[ ,3] .* exp(-exp(Phi_[ ,4]) .* input);
  else          
    return Phi_[1,1]  * exp(-exp(Phi_[1,2])  * input) + 
           Phi_[1,3]  * exp(-exp(Phi_[1,4])  * input);
}

vector SS_fol(vector Dose, vector input, matrix Phi_) {
  // Phi_[,1] = lKe, Phi_[,2] = lKa, Phi_[,3] = lCl
  int Phi__rows = rows(Phi_);
  if (Phi__rows > 1) {
    vector[Phi__rows] lKe = Phi_[,1];
    vector[Phi__rows] lKa = Phi_[,2];
    vector[Phi__rows] exp_lKe = exp(lKe);
    vector[Phi__rows] exp_lKa = exp(lKa);
    return Dose .* exp(lKe + lKa - Phi_[,3]) .* 
      (exp(-exp_lKe .* input) - exp(-exp_lKa .* input)) ./ (exp_lKa - exp_lKe);
  }
  else {
    real lKe = Phi_[1,1];
    real lKa = Phi_[1,2];
    real exp_lKe = exp(lKe);
    real exp_lKa = exp(lKa);
    return Dose * exp(lKe + lKa - Phi_[1,3]) .* 
      (exp(-exp_lKe * input) - exp(-exp_lKa * input)) / (exp_lKa - exp_lKe);
  }
}
  
vector SS_fpl(vector input, matrix Phi_) {
  // Phi_[,1] = A, Phi_[,2] = B, Phi_[,3] = xmid, Phi_[,4] = scal
  // input is generally data so cannot switch signs
  if (rows(Phi_) > 1) {
    vector[rows(Phi_)] A = Phi_[,1];
    return A + (Phi_[,2] - A) ./ (1 + exp((Phi_[,3] - input) ./ exp(Phi_[,4])));
  }
  else {
    real A = Phi_[1,1];
    return A + rep_vector(Phi_[1,2] - A, rows(input)) 
      ./ (1 + exp((Phi_[1,3] - input) / exp(Phi_[1,4])));
  }
}

vector SS_gompertz(vector x, matrix Phi_) {
  // Phi_[,1] = Asym, Phi_[,2] = b2, Phi_[,3] = b3
  vector[rows(x)] out;
  if (rows(Phi_) > 1) for (i in 1:rows(x))
    out[i] = Phi_[i,1] * exp(-Phi_[i,2] * Phi_[i,3] ^ x[i]);
  else {
    real Asym = Phi_[1,1];
    real b2 = Phi_[1,2];
    real b3 = Phi_[1,3];
    for (i in 1:rows(x)) out[i] = Asym * exp(-b2 * b3 ^ x[i]);
  }
  return out;
}

vector SS_logis(vector input, matrix Phi_) {
  // Phi_[,1] = Asym, Phi_[,2] = xmid, Phi_[,3] = scal
  // input is typically data so cannot switch signs of everything
  if (rows(Phi_) > 1)
    return Phi_[,1] ./ (1 + exp( (Phi_[,2] - input) ./ exp(Phi_[,3])));
  else 
    return rep_vector(Phi_[1,1], rows(input)) ./ 
      (1 + exp( (Phi_[1,2] - input) / exp(Phi_[1,3])));
}

vector SS_micmen(vector input, matrix Phi_) {
  // Phi_[,1] = Vm, Phi_[,2] = K
  if (rows(Phi_) > 1)
    return Phi_[ ,1] .* input ./ (Phi_[ ,2] + input);
  else
    return Phi_[1,1]  * input ./ (Phi_[1,2] + input);
}

vector SS_weibull(vector x, matrix Phi_) {
  // Phi_[,1] = Asym, Phi_[,2] = Drop, Phi_[,3] = lrc, Phi_[,4] = pwr
  vector[rows(x)] out;
  if (rows(Phi_) > 1) for (i in 1:rows(x))
    out[i] = Phi_[i,1] - Phi_[i,2] * exp(-exp(Phi_[i,3]) * x[i] ^ Phi_[i,4]);
  else {
    real Asym = Phi_[1,1];
    real Drop = Phi_[1,2];
    real lrc = Phi_[1,3];
    real pwr = Phi_[1,4];
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
