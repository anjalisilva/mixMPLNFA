data{
  int<lower=1> d;// Dimension of theta
  int<lower=1> q;// Number of latent factor
  int<lower=0> N;//Sample size
  int y[N,d];//array of Y
  vector[d] mu;
  matrix[d,d] Psi;
  matrix [d,q] Lambda;
  vector[d] normfactors;
  vector[q] mu_u;
  matrix[q,q] I_u;
  
}
parameters{
  matrix[N,d] theta;
  matrix[N,q] u;
}
model{ //Change values for the priors as appropriate
  for (n in 1:N){
    matrix[q,1] inter;
    vector[d] inter2;
    u[n,]~multi_normal(mu_u,I_u);
    inter=to_matrix(u[n,],q,1);
    inter2=to_vector(Lambda*inter);
    theta[n,]~multi_normal(mu+inter2,Psi);
  }
    for (k in 1:d){
      vector[N] z;
      z=exp(normfactors[k]+theta[,k]);
      y[,k]~poisson(z);
    }
  }
