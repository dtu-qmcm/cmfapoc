data {
  int<lower=1> N_meas;
  int<lower=1> N_mets;
  int<lower=1> N_samples;
  int<lower=1> N_met_samples;
  array[N_meas] int<lower=1, upper=N_mets> mets;
  array[N_meas] int<lower=1, upper=N_met_samples> met_samples;
  array[N_meas] real p_theoretical;
  vector[N_meas] y;
}

transformed data {
  vector[N_meas] rescaled_y = y./mean(y);
}

parameters {
  array[N_mets] real ln_offset;
  array[N_mets] real ln_sigma;
  array[N_met_samples] real<lower=0> height_T;
}

transformed parameters {
  array[N_meas] real sigma;
  array[N_mets] real ldl;
  vector[N_meas] estimated_height;
  ldl = exp(ln_offset);
  for (i in 1:N_meas){
    sigma[i] = exp(ln_sigma[mets[i]]);
    estimated_height[i] = p_theoretical[i]*height_T[met_samples[i]] + ldl[mets[i]];
  }
}

model {
  rescaled_y ~ lognormal(log(estimated_height), sigma);
  ln_sigma ~ normal(-2, 2);
  ln_offset ~ normal(0, 2);
  height_T ~ lognormal(0, 3);
}

generated quantities {
  array[N_mets] real true_offset;
  vector[N_meas] rescaled_height_offset = estimated_height .* mean(y) - y;
  for (i in 1:N_mets){
    true_offset[i] = ldl[i]*mean(y);
  }
}
