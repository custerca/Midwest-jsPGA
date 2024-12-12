functions{
  // Thermal response curve function
  real thermal_response_curve(real temp, real CTmax, real Topt, real sigma){
    real TRC;
          if (temp <= Topt) 
              TRC = exp(-((temp - Topt)/(2*sigma))^2);
          else if (CTmax >= temp && temp > Topt)
              TRC = 1 - ((temp - Topt)/(Topt - CTmax))^2;
          else
              TRC = 0.0001;
      return TRC;
  }
  
} 

// begin data section
data{
  // dimension parameters
  int<lower=0> N; // number of rows
  int<lower=0> J; // number of gears
  int<lower=0> K; // number of species
  int<lower=0> P; // number of covariates
  // vector<lower=1>[J] alpha; // Dirichlet prior parameter
  int M; // number of basis functions
  // data
  array[K] int J_num; // Number of gears per species
  array[N] int<lower=0> Y; // total caught
  // array[N*K] int<lower=0> y_vec; // vector of Y
  array[N] int species;
  array[N] int<lower=0> E; // effort
  array[N] int gears;

  matrix[N,P] X; // design matrix abundance
  
  matrix[N,M] Psi; // basis function matrix
  
  // priors for TRC...not actually priors anymore...fixed effect simulated through numerical integration
  vector[N] temp; // observed temperatures
  vector[N] CTmax_prior;
  vector[N] Topt_prior;
  vector[N] sigma; // scaling parameter for thermal response curves
    
}

// begin transformed data section
transformed data{
  vector<lower=0,upper=1>[N] TRCvec; 

  for(n in 1:N){
    TRCvec[n] = thermal_response_curve(temp[n], CTmax_prior[n], Topt_prior[n], sigma[n]);;
  }

}

// begin parameters section
parameters {
//   row_vector[K] beta_0; // intercept
  matrix[P,K] beta; // abundance parameters
  // array[K] theta; // catchability parameters
  real recip_phi; // inverse of negative binomial overdispersion parameter
  
  // optimization of code to estimate spatial basis coefficient's covariance structure
  // https://mc-stan.org/docs/2_19/stan-users-guide/multivariate-hierarchical-priors-section.html
  matrix[K,M] z_A; // standard normal
  cholesky_factor_corr[K] L_Sigma; // Cholesky decomposition of sigma
  vector<lower=0, upper=pi() / 2>[K] tau_unif;

  simplex[J_num[1]] theta1;
  simplex[J_num[2]] theta2;
  simplex[J_num[3]] theta3;
  simplex[J_num[4]] theta4;
  simplex[J_num[5]] theta5;
  simplex[J_num[6]] theta6;
  simplex[J_num[7]] theta7;
  simplex[J_num[8]] theta8;

}

// begin transformed parameters section
transformed parameters{
  vector<lower=0>[K] tau = 2.5 * tan(tau_unif);
  real<lower=0> phi; // negative binomial overdispersion parameter
  matrix[M,K] A; // spatial basis coefficients


  A = (diag_pre_multiply(tau,L_Sigma) * z_A)';
  
  phi = 1 / recip_phi; 

}

// begin model section
model {
  vector[N] Etilde;
  vector[N] lambda;

  // priors 
  to_vector(beta) ~ normal(0, 10); // beta prior
  to_vector(z_A) ~ std_normal(); 
  L_Sigma ~ lkj_corr_cholesky(2); //https://distribution-explorer.github.io/multivariate_continuous/lkj.html
  recip_phi ~ cauchy(0.,5); // inverse of overdispersion parameter prior

  // for(k in 1:K){
  //   theta[k] ~ dirichlet(rep_vector(1,J_num[k])); // catchability prior
  // }

  theta1 ~ dirichlet(rep_vector(1,J_num[1]));
  theta2 ~ dirichlet(rep_vector(1,J_num[2]));
  theta3 ~ dirichlet(rep_vector(1,J_num[3]));
  theta4 ~ dirichlet(rep_vector(1,J_num[4]));
  theta5 ~ dirichlet(rep_vector(1,J_num[5]));
  theta6 ~ dirichlet(rep_vector(1,J_num[6]));
  theta7 ~ dirichlet(rep_vector(1,J_num[7]));
  theta8 ~ dirichlet(rep_vector(1,J_num[8]));

  for(n in 1:N){
    if(species[n]==1)
      Etilde[n] = E[n] * theta1[gears[n]] * J;
    else if(species[n]==2)
      Etilde[n] = E[n] * theta2[gears[n]] * J;
    else if(species[n]==3)
      Etilde[n] = E[n] * theta3[gears[n]] * J;
    else if(species[n]==4)
      Etilde[n] = E[n] * theta4[gears[n]] * J;
    else if(species[n]==5)
      Etilde[n] = E[n] * theta5[gears[n]] * J;
    else if(species[n]==6)
      Etilde[n] = E[n] * theta6[gears[n]] * J;
    else if(species[n]==7)
      Etilde[n] = E[n] * theta7[gears[n]] * J;
    else if(species[n]==8)
      Etilde[n] = E[n] * theta8[gears[n]] * J;
  }

  for(n in 1:N){
    lambda[n] = (X[n] * beta[,species[n]]) + (Psi[n] * A[,species[n]]) + log(TRCvec[n]);
  }

  Y ~ neg_binomial_2_log(log(Etilde + 1e-10) + lambda, phi); // modeling catch data
  
}

generated quantities{
  corr_matrix[K] Sigma; // Correlation matrix Sigma for spatial basis coefficients covariance structure
  cov_matrix[K] Pi; // Full covariance matrix T*A*T
  Pi = diag_pre_multiply(tau, L_Sigma) * diag_pre_multiply(tau, L_Sigma)';
  Sigma = multiply_lower_tri_self_transpose(L_Sigma); 
}