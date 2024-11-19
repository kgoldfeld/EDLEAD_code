data {
  
  int<lower=0> N_ED;                 // number of eds
  int<lower=0> N;                    // number patients
  
  matrix[N, 8] x_abc;
  
  array[N] int<lower=1,upper=N_ED> ed;     // ED for individual i
  array[N] int<lower=0,upper=1> y;         // outcome for individual i
  
}

parameters {
  
  vector[8] z;

  real delta_m;
  real<lower = 0> sigma_m;
  
  real delta_x;
  real<lower=0> sigma_x;
  
  vector[N_ED] ed_effect;
  real<lower=0> sigma_ed;
  
}

transformed parameters {
  
  vector[8] tau;
  
  tau[1] = z[1];
  
  for (i in 2:4){
    tau[i] = sigma_m * z[i] + delta_m;
  }
  
  for (i in 5:7){
    tau[i] = sigma_x * z[i] + delta_x;
  }
  
  tau[8] = z[8];
  
}

model {
  
  sigma_m ~ student_t(3, 0, 2.5); 
  sigma_x ~ student_t(3, 0, 2.5); 

  // delta_m ~ normal(0, 1); 
  // delta_x ~ normal(0, 1); 
  
  delta_m ~ normal(0, 5); 
  delta_x ~ normal(0, 5); 

  ed_effect ~ normal(0, sigma_ed);
  sigma_ed ~ student_t(3, 0, 2.5);

  z ~ std_normal();

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
