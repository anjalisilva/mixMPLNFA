rm(list=ls())
#### Libraries needed
library(mvtnorm)

set.seed(100)
####
p<-2 ### Dimensionality of latent variable
d<-5 ### Dimensionality of observed variable
G<-2 ### Number of Groups
N<-1000 ### Sample size
tot_data<-5 ### Total number of datasets

###Storing output
dat<-list() ###Stores information about all run
par<-list() ###Stores information about parameters at every run



lambda_1<-matrix(runif(p*d),nrow=d)
lambda_2<-matrix(runif(p*d),nrow=d)

psi_1<-diag(d)*runif(d)
psi_2<-diag(d)*runif(1)

sigma_1<-lambda_1%*%t(lambda_1)+psi_1
sigma_2<-lambda_2%*%t(lambda_2)+psi_2

mu_1<-c(6, 3, 3, 6 ,3)
mu_2<-c(5, 3, 5, 3, 5)

pi_g<-c(0.32,0.68)



par[[1]]<-NULL # for pi
par[[2]]<-list() # for mu
par[[3]]<-list() # for Lambda
par[[4]]<-list() # for Psi
par[[5]]<-list() # for Sigma
names(par)<-c("pi","mu","Lambda","Psi","Sigma")

par[[1]]<-pi_g
par[[2]]<-list(mu_1,mu_2)
par[[3]]<-list(lambda_1,lambda_2)
par[[4]]<-list(psi_1,psi_2)
par[[5]]<-list(sigma_1,sigma_2)


for (run in 1:tot_data) {

  # set.seed(run)
  Y2 <- Y <- matrix(0, ncol = d, nrow = N)
  X <- matrix(0, ncol = d, nrow = N)
  U <- rmvnorm(N, mean = rep(0, p), sigma = diag(p))
  
  
  z <- t(rmultinom(n = N, size = 1, prob = pi_g))
  
  for (i in 1:N) {
    grp <- which(z[i, ] == 1)
    X[i, ] <- rmvnorm(n = 1, 
                      mean = par[[2]][[grp]] + par[[3]][[grp]] %*% U[i, ],
                      sigma = par[[4]][[grp]])
    for (j in 1:d) {
      Y[i, j] <- rpois(1, exp(X[i, j]))
    }
  }
  
  norms <- log(edgeR::calcNormFactors(Y))
  
  for (i in 1:N){
    for (j in 1:d){
      Y2[i,j] <- rpois(1, exp(X[i, j] + norms[j])) 
    }
  }
  
  
  dat[[run]] <- list()
  dat[[run]][[1]] <- par
  dat[[run]][[2]] <- Y
  dat[[run]][[3]] <- X
  dat[[run]][[4]] <- z
  dat[[run]][[5]] <- U
  
}


save.image("sim_data_ASilva.Rdata")