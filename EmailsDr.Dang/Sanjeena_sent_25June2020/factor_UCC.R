rm(list=ls())
load("sim_data.Rdata")


G<-2
pmax<-2

fun_lambda_g<-function(g){Sk[,,g]%*%t(beta[[g]])%*%solve(bigTheta[[g]])}
fun_psi_ucc<-function(g){1/d*ng[g]/N*sum(diag(Sk[,,g]-lambdanew[[g]]%*%beta[[g]]%*%Sk[,,g]))}

#parallel_FA<-function(run,G=2,pmax=2){

###### Parameter Updates ####

Y<-dat[[run]][[2]]
true<-dat[[run]][[4]]
true_par<-dat[[run]][[1]]
theta_old<-list()
lambdanew<-list()
psinew<-list()

true_par
beta<-list()
bigTheta<-list()


N<-nrow(Y)
d<-ncol(Y)

lib_mat<-rep(1,d)
modelName<-"UCC"



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
ng<-colSums(z)

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
}

  
  diff<-1
  min_diff<-10^{-6}
  while(diff==1){
    
    ##########Model UCU###########################
    if (substr(modelName,1,3)=="UCC"){
      #Updating Lambda and beta
      lambda_old<-lambda
      psi_old<-psi
      for (g in 1:G){
        beta[[g]]<-t(lambda[[g]])%*%solve(lambda[[g]]%*%t(lambda[[g]])+psi[[g]])
        bigTheta[[g]]<-diag(pmax)-beta[[g]]%*%lambda[[g]]+beta[[g]]%*%Sk[,,g]%*%t(beta[[g]])
        lambdanew[[g]]<-fun_lambda_g(g)
      }		
      psinew[[1]]<-sum(sapply(1:G,fun_psi_ucc))*diag(d)
      for (g in 1:G){		
        #Complete Sigma		
        sigma[[g]]<-(lambdanew[[g]]%*%t(lambdanew[[g]])+psinew[[1]])
        lambda[[g]]<-lambdanew[[g]]
        psi[[g]]<-psinew[[1]]	
      }
      		
      par<-G*(d*pmax-0.5*pmax*(pmax-1))+1
    }
  
 
    
    
    diff_lam<-NULL
    diff_psi<-NULL
for (g in 1:G){
  diff_lam[g]<-norm(lambda[[g]]-lambda_old[[g]],type="F")
  diff_psi[g]<-norm(psi[[g]]-psi_old[[g]],type="F")
}
  if (all(abs(diff_lam)<min_diff)&all(abs(diff_psi)<min_diff)) diff<-0
  
  }
  
for (g in 1:G){
  isigma[[g]]<-solve(sigma[[g]])
}
  pi_g<-colSums(z)/N
  ng<-colSums(z)
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
    if (abs(aloglik[it]-loglik[it-1])<0.0001) checks<-1 else checks<-checks
  }
}	
print(it)
it<-it+1
if (it==it_max) checks<-1   
}


