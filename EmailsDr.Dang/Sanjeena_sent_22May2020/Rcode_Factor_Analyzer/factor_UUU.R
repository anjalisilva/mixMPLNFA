rm(list=ls())
load("sim_data.Rdata")

parallel_FA<-function(run,G=2,pmax=2){

###### Parameter Updates ####

Y<-dat[[run]][[2]]
true<-dat[[run]][[4]]
true_par<-dat[[run]][[1]]

true_par


N<-nrow(Y)
d<-ncol(Y)

lib_mat<-rep(1,d)



#### Initialization ###
mu<-list()
psi<-list()
lambda<-list()
sigma<-list()
isigma<-list()
m<-list()
S<-list()
P<-list()
Q<-list()

###Other intermediate items initialized
start<-list()
Sk<-array(0, c(d,d,G) )
GX<-list()
dGX<-list()
z_S<-list()


k_means<-kmeans(log(Y+1),centers=G,nstart=100)$cluster  
z<-mclust::unmap(k_means)
pi_g<-colSums(z)/N

for (g in 1:G){
  obs<-which(z[,g]==1)
  mu[[g]]<-colMeans(log(Y[obs,]+1/6)) 
  sigma[[g]]<-var(log(Y[obs,]+1/6))
  isigma[[g]]<-solve(sigma[[g]])
  temp<-eigen(sigma[[g]])
  lambda[[g]]<-matrix(NA,ncol=pmax,nrow=d)
  for (q in 1:pmax){
    lambda[[g]][,q]<-temp$vectors[,q]*sqrt(temp$values[q]) 
  }
  psi[[g]]<-diag(sigma[[g]]-lambda[[g]]%*%t(lambda[[g]]))*diag(d)
}


for (g in 1:G){
  start[[g]]<-log(Y+1/6) ###Starting value for M
  m[[g]]<-log(Y+1/6)
  S[[g]]<-list()
  for (i in 1:N){
  S[[g]][[i]]<-diag(d)*0.000000001
  }
}

checks<-0
it<-1
aloglik<-NULL
loglik<-NULL
aloglik[c(1,2,3)]<-0
it_max<-1000


while (checks==0){


  for (g in 1:G){
    GX[[g]]<-list()
    dGX[[g]]<-list()
    z_S[[g]]<-list()
    for (i in 1:N){
      dGX[[g]][[i]]<-diag(exp(log(lib_mat)+start[[g]][i,])+0.5*diag(S[[g]][[i]]),d)+isigma[[g]]
      S[[g]][[i]]<-solve(dGX[[g]][[i]])
      z_S[[g]][[i]]<-z[i,g]*S[[g]][[i]]
      GX[[g]][[i]]<-Y[i,]-exp(start[[g]][i,]+log(lib_mat)+0.5*diag(S[[g]][[i]]))-(isigma[[g]])%*%(start[[g]][i,]-mu[[g]])
      m[[g]][i, ]<-start[[g]][i,]+S[[g]][[i]]%*% GX[[g]][[i]]
    }
    start[[g]]<-m[[g]]
    ####Updating mu
    mu[[g]]<-colSums(z[,g]*m[[g]])/sum(z[,g])
    
    ####Updating Sample covariance
    mu_mat<-matrix(rep(mu[[g]],N),nrow=N,byrow=TRUE)
    res<-m[[g]]-mu_mat
    temp <- cov.wt(res, wt=z[,g], center=FALSE, method="ML")
    Sk[,,g]<- temp$cov
  
    ###Updating Beta
    beta<-t(lambda[[g]])%*%isigma[[g]]
    Q[[g]]<-diag(pmax)-beta%*%lambda[[g]]
    temp3<-solve(psi[[g]])
   
    P[[g]]<-(m[[g]]-mu_mat)%*%t(Q[[g]]%*%t(lambda[[g]])%*%temp3)
    
    diff<-1
    min_diff<-10^{-6}
    while(diff==1){
    
    lambda_old<-lambda[[g]]
    psi_old<-psi[[g]]
    theta_old<-diag(pmax)-beta%*%lambda[[g]]+beta%*%Sk[,,g]%*%t(beta)
    ####Updating Lambda
    lambda[[g]]<-Sk[,,g]%*%t(beta)%*%solve(theta_old)
    
    ###Updating Beta with new lambda
    beta<-t(lambda[[g]])%*%isigma[[g]]
  
    ###Updating Psi
    psi[[g]]<-diag(Sk[,,g]-lambda[[g]]%*%beta%*%Sk[,,g]+Reduce("+",z_S[[g]])/sum(z[,g]))*diag(d)
  
    diff_lam<-norm(lambda[[g]]-lambda_old,type="F")
    diff_psi<-norm(psi[[g]]-psi_old,type="F")
    if (abs(diff_lam)<min_diff&abs(diff_psi)<min_diff) diff<-0
    
    sigma[[g]]<-lambda[[g]]%*%t(lambda[[g]])+psi[[g]]
    isigma[[g]]<-solve(sigma[[g]]) ###Ideally for high dimensional data, we should be using Woodbury identity so it avoids inverting a p by p matrix
  
    }
    
  }

  pi_g<-colSums(z)/N
  lib_mat_full<-matrix(1,ncol=d,nrow=N) ###Matrix containing normaization factor
  ### Some useful functions
  fun_five<-function(x,y=isigma[[g]]){
    sum(diag(x%*%y))
  }
  
   F<-matrix(NA,ncol=G,nrow=N)
  
  for (g in 1:G){
    two<-rowSums(exp(m[[g]]+log(lib_mat_full)+0.5*matrix(unlist(lapply(S[[g]],diag)),ncol=d,byrow=TRUE)))
    five<-0.5*unlist(lapply(S[[g]],fun_five))
    six<-0.5*log(unlist(lapply(S[[g]],det)))
    F[,g]<-pi_g[g]*exp(rowSums(m[[g]]*Y)-two-rowSums(lfactorial(Y))+rowSums(log(lib_mat_full)*Y)-0.5*mahalanobis(m[[g]],center=mu[[g]],cov=sigma[[g]])-five+six-0.5*log(det(sigma[[g]]))-d/2)
  }  
  
  
loglik[it]<-sum(log(rowSums(F)))
z<-F/rowSums(F)


if (it>2){
  #Aitkaine's stopping criterion
  if ((loglik[it-1]-loglik[it-2])==0) checks<-1 else{
    a<-(loglik[it]-loglik[it-1])/(loglik[it-1]-loglik[it-2])
    add_to<-(1/(1-a)*(loglik[it]-loglik[it-1]))
    # }
    aloglik[it]<-loglik[it-1]+add_to
    if (abs(aloglik[it]-loglik[it-1])<0.01) checks<-1 else checks<-checks
  }
}	
print(it)
it<-it+1
if (it==it_max) checks<-1   

  
}
return(list(pi_g=pi_g,mu=mu,sigma=sigma,lambda=lambda,psi=psi,z=z,loglik=loglik,kmeans=k_means,true=true))
}

if (.Platform$OS.type == "unix") {
  library("doMC")
  registerDoMC(6) #parallel::detectCores() gives the total number of cores available
}
library("plyr")

total_run<-200
ptm<-proc.time()
output_all<-list()
output_all <- foreach(run = 1:total_run, .errorhandling = "pass") %dopar% {
  parallel_FA(run)
}
proc.time()-ptm



ARI<-NULL
for (i in 1:total_run){
ARI[i]<- mclust::adjustedRandIndex(mclust::map(output_all[[i]]$true),mclust::map(output_all[[i]]$z))
}

pi_all<-matrix(0,ncol=2,nrow=total_run)
mean_all<-matrix(0,ncol=d*2,nrow=total_run)
var_all<-matrix(0,ncol=d*d*2,nrow=total_run)


for (i in 1:total_run){
  or<-order(output_all[[i]]$pi_g)
  pi_all[i,]<-c(output_all[[i]]$pi_g[or[[1]]],output_all[[i]]$pi_g[or[[2]]])
  mean_all[i,]<-c(output_all[[i]]$mu[[or[1]]],output_all[[i]]$mu[[or[2]]])
  var_all[i,]<-c(c(output_all[[i]]$sigma[[or[1]]]),c(output_all[[i]]$sigma[[or[2]]]))
}

av_mean<-round(apply(mean_all,2,mean),2)
av_mean[1:d]
av_mean[-c(1:d)]
sd_mean<-round(apply(mean_all,2,sd),2)
sd_mean[1:d]
sd_mean[-c(1:d)]

av_sigma<-round(apply(var_all,2,mean),2)
round(dat[[1]][[1]]$Sigma[[1]],2)
round(dat[[1]][[1]]$Sigma[[2]],2)
matrix(av_sigma[1:(d*d)],d,d)
matrix(av_sigma[-c(1:(d*d))],d,d)

sd_sigma<-round(apply(var_all,2,sd),2)
matrix(sd_sigma[1:(d*d)],d,d)
matrix(sd_sigma[-c(1:(d*d))],d,d)

save.image("Routput.Rdata")

