  int<lower = 0, upper = 1> is_vae;
  int<lower = 0> depth_vae;
  int<lower = 0> p_vae;
  int<lower = 0> ncoefs_vae; // excludes padding
  int<lower = 0> K_vae;      // includes padding
  matrix[K_vae, K_vae] W_vae[depth_vae];
  vector[K_vae] B_vae[depth_vae];