functions {
  real gc_modeled(real beta_1, real plane){
    return beta_1 * plane;
  }
}

data {
  // Total number of data points 
  int N;
  
  // Number of entries in each level of the hierarchy
  int J_1;
  
  // Index arrays to keep track of hierarchical structure
  int index_1[N];

  // Plane gas measurements
  real<lower=0> plane[N];
  
  // GC gas data
  real<lower=0> gc_sim[N];
  
  // Likelihood sigma estimate using plane tropopause height
  real sigma[N];
  
  // Info for priors

  // Tau mean comes from delta t between the replicate and plane sample
  // It is now a vector! 
  real<lower=0> tau_[N];
  
  // other
  //
  real mu_beta;
  real sigma_beta;
 
}

parameters {
  real<lower=0> beta_;
  //vector[J_1] beta_1;
  vector[J_1] beta_1_tilde;

}

transformed parameters {  

  //vector[J_1] beta_1 = beta_ + tau_[J_1] * beta_1_tilde;


  vector[J_1] beta_1;
  real<lower=0> mu_gc[N];
  real epsilon = 1.0e-6;
  for (i in 1:N){
    beta_1[index_1[i]] = fmax((beta_ + tau_[i] * beta_1_tilde[index_1[i]]), epsilon);
    //beta_1[index_1[i]] = beta_ + tau_[i] * beta_1_tilde[index_1[i]];
    mu_gc[i] = gc_modeled(beta_1[index_1[i]], plane[i]);
  }

  
  vector[J_1] alpha_1;
  for (j in 1:J_1){
    alpha_1[j] = 1.0 / beta_1[j];  
  }
  real alpha = 1.0/beta_;
  
}

model { 

  beta_ ~ normal(mu_beta, sigma_beta);
  //beta_1 ~ normal(beta_, tau_);
  beta_1_tilde ~ normal(0, 1);
  
  gc_sim ~ lognormal(log(mu_gc), sigma); 
    
}

generated quantities {
  real gc_sim_ppc[N];
  for (i in 1:N) {
    gc_sim_ppc[i] = lognormal_rng(log(gc_modeled(beta_1[index_1[i]], plane[i])), sigma[i]);
  }
}