data {
  
  int<lower=0> N_ED;                 // number of eds
  int<lower=0> N;                    // number patients
  
  matrix[N, 8] x_abc;
  
  array[N] int<lower=1,upper=N_ED> ed;     // ED for individual i
  array[N] int<lower=0,upper=1> y;         // outcome for individual i
  
}

parameters {
  
  vector[8] tau;
  
  vector[N_ED] ed_effect;
  real<lower=0> sigma_ed;

}

model {
  
  ed_effect ~ normal(0, sigma_ed);
  sigma_ed ~ student_t(3, 0, 2.5);

  tau ~ normal(0, 1); // was 2.5

  y ~ bernoulli_logit(ed_effect[ed] + x_abc * tau);
  
}

generated quantities {
  
  array[7] real lOR;
  
  lOR[1] = tau[2];                                            //  a=1, b=0, c=0
  lOR[2] = tau[3];                                            //  a=0, b=1, c=0
  lOR[3] = tau[4];                                            //  a=0, b=0, c=1
  lOR[4] = tau[2] + tau[3] + tau[5];                          //  a=1, b=1, c=0
  lOR[5] = tau[2] + tau[4] + tau[6];                          //  a=1, b=0, c=1
  lOR[6] = tau[3] + tau[4] + tau[7];                          //  a=0, b=1, c=1
  lOR[7] = tau[2]+tau[3]+tau[4]+tau[5]+tau[6]+tau[7]+tau[8];  //  a=1, b=1, c=1
  
}
