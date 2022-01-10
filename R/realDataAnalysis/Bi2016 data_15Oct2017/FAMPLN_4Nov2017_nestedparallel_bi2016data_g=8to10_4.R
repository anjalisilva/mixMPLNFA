# 4 November 2017
# testing for G=8:10, q=1:3, models=(2)

#### Packages ####
#rm(list=ls())
#install.packages("rstan", type = "source")
#devtools::install_github("stan-dev/rstan", ref = "develop", subdir = "rstan/rstan")
#options(max.print=999999)
library(rstan)
#install.packages("Rcpp", type = "source")
library(Rcpp)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(mvtnorm)
library(mclust)
library(edgeR)
#install.packages("coda")
library("coda")

#### Reading in data ####

bi2016<-as.matrix(read.csv("median_counts_DEgenes_Bietal2016data.csv", header = TRUE, row.names=1))
dim(bi2016) # 570  12

which(rowSums(bi2016)==0) #none = no rows are all zero


#### functions ####

zvalue_calculation<-function(theta_Stan,y,G,mu_g,Sig_g,PI, normalizefactors){
  d<-ncol(y)
  n<-nrow(y)
  forz<-matrix(NA,ncol=G,nrow=n)
  for (g in 1:G){
    for (i in 1:n){
      x<-theta_Stan[[g]][i,]
      # for zig calculation (the numerator part)
      forz[i,g]<-PI[g]*exp(t(y[i,])%*%(x+normalizefactors)-sum(exp(x+normalizefactors))-sum(lfactorial(y[i,]))-
                             d/2*log(2*pi)-1/2*log(det(Sig_g[((g-1)*d+1):(g*d),]))-0.5*t(x-mu_g[g,])%*%solve(Sig_g[((g-1)*d+1):(g*d),])%*%(x-mu_g[g,]))
      
    }
    # check which forz == 0 and rowSums(forz)==0 and which of these
    # have both equalling to 0 (because 0/0 =NaN)
    if (G==1){
      errorpossible<-Reduce(intersect, list(which(forz==0),which(rowSums(forz)==0)))
      zvalue<-forz/rowSums(forz)
      zvalue[errorpossible,]<-1
    }else {zvalue<-forz/rowSums(forz)}
    
  }
  
  # check which forz == 0 and rowSums(forz)==0 and which of these
  # have both equalling to 0 (because 0/0 =NaN)
  if (G==1){
    errorpossible<-Reduce(intersect, list(which(forz==0),which(rowSums(forz)==0)))
    zvalue<-forz/rowSums(forz)
    zvalue[errorpossible,]<-1
  }else {zvalue<-forz/rowSums(forz)}
  
  
  return(zvalue)
}

calc_likelihood<-function(z, PI, y, mu_g, G, Sig_g, theta_Stan, normalizefactors){ 
  n<-nrow(y)
  like<-matrix(NA, nrow=n, ncol=G)
  for (g in 1:G){
    for (i in 1:n){
      x<-theta_Stan[[g]][i,]
      d<-ncol(y)
      like[i,g]<-(z[i,g] *(log(PI[g]) +
                             t(y[i,])%*%(x+normalizefactors)-sum(exp(x+normalizefactors))-sum(lfactorial(y[i,]))-
                             d/2*log(2*pi)-1/2*log(det(Sig_g[((g-1)*d+1):(g*d),]))-0.5*t(x-mu_g[g,])%*%solve(Sig_g[((g-1)*d+1):(g*d),])%*%(x-mu_g[g,])))
    }
  }
  loglike<-sum(rowSums(like))
  return(loglike)
}

stanrun<-function(model, gmin,gmax,y,mu_all_outer, it_outer, sigma_all_outer,numb_iterations, n_chain, normalizefacs){
  fitrstan<-list()
  d<-ncol(y)
  
  
  for (g in gmin:gmax){
    data1=list(d=ncol(y),N=nrow(y),y=y,mu=mu_all_outer[[it_outer-1]][g,],Sigma=sigma_all_outer[[it_outer-1]][[g]], normfactors=as.vector(normalizefacs))
    stanproceed<-0
    try=1
    
    while (!stanproceed){
      
      #cat("\nRstan generating sample at outer iteration", it_outer, "for g: ",g , "try: ", try)
      #cat("\nNumber of iterations is", numb_iterations, "\n")
      fitrstan[[g]]<-sampling(object=model,
                              data=data1,
                              iter=numb_iterations, chains = n_chain, verbose=FALSE, refresh=-1)
      
      if (all(summary(fitrstan[[g]])$summary[,"Rhat"] < 1.1) == TRUE && all(summary(fitrstan[[g]])$summary[,"n_eff"]>100) == TRUE){
        stanproceed<-1
      } else if(all(summary(fitrstan[[g]])$summary[,"Rhat"] < 1.1) != TRUE || all(summary(fitrstan[[g]])$summary[,"n_eff"]>100) != TRUE){
        if(try == 10){ # stop after 10 tries
          stanproceed = 1
        }
        numb_iterations = numb_iterations+100
        try=try+1
      }
    }
  } # close g loop
  
  
  results <- list(fitrstan = fitrstan,
                  numb_iterations = numb_iterations)
  class(results) <- "RStan"
  return(results)
  
  return(results)
}

initializationrun<-function(gmodel, y, init_method, init_iterations, n_chain, numb_iterations, initialization=NA, normalizefactors, modelname, q){
  z<-init_runs<-list()
  logL_init<-vector()
  n<-nrow(y)
  d<-ncol(y)
  
  for(iterations in 1:init_iterations){
    if (init_method=="kmeans" | is.na(init_method)){
      if (!require(mclust)) suppressWarnings(install.packages('mclust')) # loading needed packages
      suppressWarnings(library(mclust))
      z[[iterations]]<-unmap(kmeans(log(y+1/3),gmodel)$cluster)
    }else if (init_method=="random"){
      if(gmodel==1){ # generating z if g=1
        z[[iterations]] <- as.matrix(rep.int(1, times=n), ncol=gmodel, nrow=n)
      } else { # generating z if g>1
        z_conv=0
        while(!z_conv){ # ensure that dimension of z is same as G (i.e. 
          # if one column contains all 0s, then generate z again)
          z[[iterations]] <- t(rmultinom(n, size = 1, prob=rep(1/gmodel,gmodel))) 
          if(length(which(colSums(z[[iterations]])>0)) ==gmodel){
            z_conv=1
          }
        }
      }
    }else if (init_method=="medoids"){
      if (!require(cluster)) install.packages('cluster') 
      library(cluster)
      
      if (!require(mclust)) suppressWarnings(install.packages('mclust')) # loading needed packages
      suppressWarnings(library(mclust))
      
      z[[iterations]]<-unmap(pam(log(y+1/3),k=gmodel)$cluster)
    }else if (init_method=="clara"){
      if (!require(cluster)) install.packages('cluster') 
      library(cluster)
      
      z[[iterations]]<-unmap(clara(log(y+1/3),k=gmodel)$cluster)
    }else if (init_method=="fanny"){
      if (!require(cluster)) install.packages('cluster') 
      library(cluster)
      
      z[[iterations]]<-unmap(fanny(log(y+1/3),k=gmodel)$cluster)
    }
    
    init_runs[[iterations]]=cluster_mpln(y=y,z=z[[iterations]],G=gmodel,n_chain=n_chain,numb_iterations=numb_iterations, initialization="init", normalizefactors=normalizefactors, modelname=modelname, q=q)
    logL_init[iterations] <- unlist(tail((init_runs[[iterations]]$loglikelihood), n=1)) 
  }
  
  initialization<-init_runs[[which(logL_init==max(logL_init, na.rm = TRUE))[1]]]
  return(initialization)
}  

# Model selection Functions 
# BIC function 
BIC_function<-function(ll, k, n, run, qmin, qmax, gmin, gmax, modelnames, models){
  
  BICvalues = -2*ll+ (k* log(n))
  colnames(BICvalues) = c(paste0("q=", qmin:qmax))
  rownames(BICvalues) = c(paste0("g=", gmin:gmax))
  
  BICmodel_g = seq(gmin, gmax, 1)[which(min(BICvalues) == BICvalues, arr.ind = T)[1]] # obtaining row (which is g)
  BICmodel_q = seq(qmin, qmax, 1)[which(min(BICvalues) == BICvalues, arr.ind = T)[2]] # obtaining column (which is q)
  BICmodel_m = modelnames[models[which(min(BICvalues) == BICvalues, arr.ind = T)[3]]] # obtaining model name
  BICmodel_labels = run[[which(min(BICvalues) == BICvalues, arr.ind = T)[1]]][[which(min(BICvalues) == BICvalues, arr.ind = T)[2]]][[which(min(BICvalues) == BICvalues, arr.ind = T)[3]]]$allresults$clusterlabels # obtaining model labels
  BICMessage = NA
  
  if (max(BICmodel_labels)!=BICmodel_g){
    BICmodel_g<-max(BICmodel_labels)
    BICMessage<-"Spurious or empty cluster resulted."
  }
  
  dimnames(BICvalues)[[3]]=modelnames[models]
  
  BICresults<-list(allBICvalues = BICvalues,
                   BICmodelselected_g=BICmodel_g,
                   BICmodelselected_q=BICmodel_q,
                   BICmodelselected_modelname=BICmodel_m,
                   BICmodelselected_labels=BICmodel_labels,
                   BICMessage=BICMessage)
  class(BICresults) <- "BIC"
  return(BICresults)
}  

# ICL function
ICL_function<-function(bIc, run, qmin, qmax, gmin, gmax, modelnames, models){
  ICLvalues<-array(dim=c(gmax-gmin+1,qmax-qmin+1,length(models)))
  
  for (g in 1:(gmax-gmin+1)){
    for (qvalue in 1:(qmax-qmin+1)){
      for (numbmodel in 1:length(models)){
        z<-run[[g]][[qvalue]][[numbmodel]]$allresults$probaPost
        mapz<-mclust::unmap(run[[g]][[qvalue]][[numbmodel]]$allresults$clusterlabels)
        forICL<-function(g){sum(log(z[which(mapz[,g]==1),g]))}
        ICLvalues[g,qvalue,numbmodel] <- bIc$allBICvalues[g,qvalue,numbmodel] + sum(sapply(1:ncol(mapz),forICL))
      }
    }
  }
  
  colnames(ICLvalues) = c(paste0("q=", qmin:qmax))
  rownames(ICLvalues) = c(paste0("g=", gmin:gmax))
  ICLvalues[which(ICLvalues== -Inf)] = NA
  
  ICLmodel_g = seq(gmin, gmax, 1)[which(min(ICLvalues, na.rm=TRUE) == ICLvalues, arr.ind = T)[1]] # obtaining row (which is g)
  ICLmodel_q = seq(qmin, qmax, 1)[which(min(ICLvalues, na.rm=TRUE) == ICLvalues, arr.ind = T)[2]] # obtaining column (which is q)
  ICLmodel_m = modelnames[models[which(min(ICLvalues, na.rm=TRUE) == ICLvalues, arr.ind = T)[3]]] # obtaining model name
  ICLmodel_labels = run[[which(min(ICLvalues, na.rm=TRUE) == ICLvalues, arr.ind = T)[1]]][[which(min(ICLvalues, na.rm=TRUE) == ICLvalues, arr.ind = T)[2]]][[which(min(ICLvalues, na.rm=TRUE) == ICLvalues, arr.ind = T)[3]]]$allresults$clusterlabels
  ICLMessage = NA
  
  if (max(ICLmodel_labels)!=ICLmodel_g){
    ICLmodel_g<-max(ICLmodel_labels)
    ICLMessage<-"Spurious or empty cluster resulted."
  }
  
  dimnames(ICLvalues)[[3]]=modelnames[models]
  
  ICLresults<-list(allICLvalues=ICLvalues,
                   ICLmodelselected_g=ICLmodel_g,
                   ICLmodelselected_q=ICLmodel_q,
                   ICLmodelselected_modelname=ICLmodel_m,
                   ICLmodelselected_labels=ICLmodel_labels,
                   ICLMessage=ICLMessage)
  class(ICLresults) <- "ICL"
  return(ICLresults)
}

# AIC function
AIC_function<-function(ll, k, run, qmin, qmax, gmin, gmax, modelnames, models){
  
  AICvalues = -2*ll+ 2*k
  colnames(AICvalues) = c(paste0("q=", qmin:qmax))
  rownames(AICvalues) = c(paste0("g=", gmin:gmax))
  
  AICmodel_g = seq(gmin, gmax, 1)[which(min(AICvalues) == AICvalues, arr.ind = T)[1]] # obtaining row (which is g)
  AICmodel_q = seq(qmin, qmax, 1)[which(min(AICvalues) == AICvalues, arr.ind = T)[2]] # obtaining column (which is q)
  AICmodel_m = modelnames[models[which(min(AICvalues) == AICvalues, arr.ind = T)[3]]] # obtaining model name
  AICmodel_labels = run[[which(min(AICvalues) == AICvalues, arr.ind = T)[1]]][[which(min(AICvalues) == AICvalues, arr.ind = T)[2]]][[which(min(AICvalues) == AICvalues, arr.ind = T)[3]]]$allresults$clusterlabels
  AICMessage = NA
  
  if (max(AICmodel_labels)!=AICmodel_g){
    AICmodel_g = max(AICmodel_labels)
    AICMessage = "Spurious or empty cluster resulted."
  }
  
  dimnames(AICvalues)[[3]] = modelnames[models]
  
  AICresults<-list(allAICvalues = AICvalues,
                   AICmodelselected_g=AICmodel_g,
                   AICmodelselected_q=AICmodel_q,
                   AICmodelselected_modelname=AICmodel_m,
                   AICmodelselected_labels=AICmodel_labels,
                   AICMessage=AICMessage)
  
  class(AICresults) <- "AIC"
  return(AICresults)
}  

# AIC3 function
AIC3_function<-function(ll, k, run, qmin, qmax, gmin, gmax, modelnames, models){
  
  AIC3values <- -2*ll+ 3*k
  colnames(AIC3values)<-c(paste0("q=", qmin:qmax))
  rownames(AIC3values)<-c(paste0("g=", gmin:gmax))
  
  AIC3model_g<-seq(gmin, gmax, 1)[which(min(AIC3values) == AIC3values, arr.ind = T)[1]] # obtaining row (which is g)
  AIC3model_q<-seq(qmin, qmax, 1)[which(min(AIC3values) == AIC3values, arr.ind = T)[2]] # obtaining column (which is q)
  AIC3model_m<-modelnames[models[which(min(AIC3values) == AIC3values, arr.ind = T)[3]]] # obtaining model name
  AIC3model_labels<-run[[which(min(AIC3values) == AIC3values, arr.ind = T)[1]]][[which(min(AIC3values) == AIC3values, arr.ind = T)[2]]][[which(min(AIC3values) == AIC3values, arr.ind = T)[3]]]$allresults$clusterlabels
  AIC3Message<-NA
  
  if (max(AIC3model_labels)!=AIC3model_g){
    AIC3model_g<-max(AIC3model_labels)
    AIC3Message<-"Spurious or empty cluster resulted."
  }
  
  dimnames(AIC3values)[[3]] = modelnames[models]
  
  AIC3results<-list(allAIC3values = AIC3values,
                    AIC3modelselected_g=AIC3model_g,
                    AIC3modelselected_q=AIC3model_q,
                    AIC3modelselected_modelname=AIC3model_m,
                    AIC3modelselected_labels=AIC3model_labels,
                    AIC3Message=AIC3Message)
  
  class(AIC3results) <- "AIC3"
  return(AIC3results)
}

calculate_parameters<-function(Gmin, Gmax, qmin, qmax, modelnames, models, allruns){
  
  obs = ll = k = mu_para = sigma_para = pi_para = gchoice = array(dim=c(Gmax-Gmin+1,qmax-qmin+1,length(models)))
  colnames(obs) = colnames(ll) = colnames(k) = colnames(mu_para) = colnames(gchoice) = colnames(sigma_para) = colnames(pi_para) = c(paste('q =',qmin:qmax))
  rownames(obs) = rownames(ll) = rownames(k) = rownames(mu_para) = rownames(gchoice) = rownames(sigma_para) = rownames(pi_para) = c(paste('g =',Gmin:Gmax))
  dimnames(obs)[[3]] = dimnames(ll)[[3]] = dimnames(gchoice)[[3]] = dimnames(k)[[3]] = modelnames[models]
  totaltime = list()
  
  
  for(g in 1:(Gmax-Gmin+1)) {
    for(qvalue in 1:(qmax-qmin+1)){
      for (numbmodel in 1:length(models)){
        inputmodel = modelnames[models[numbmodel]]
        
        totaltime[[numbmodel]] = allruns[[g]][[qvalue]][[numbmodel]]$totaltime
        
        obs[g,qvalue,numbmodel] = nrow(allruns[[g]][[qvalue]][[numbmodel]]$dataset)
        p = ncol(allruns[[g]][[qvalue]][[numbmodel]]$dataset)
        
        # calculating log-likelihood
        ll[g,qvalue,numbmodel]<-unlist(tail(allruns[[g]][[qvalue]][[numbmodel]]$allresults$loglikelihood, n=1)) # save the final log-likelihood
        
        if ((Gmax-Gmin+1) == 1){ # if running for only one specific G
          G=Gmin
        }else{  # if running for only multiple G
          G=g
        }
        
        if ((qmax-qmin+1) == 1){# if running for only one specific G
          Q=qmin
        }else{ # if running for only multiple G
          Q=qvalue
        }
        
        
        #calculating number of parameters of mu, sigma, and pi
        mu_para[g,qvalue,numbmodel]<-p*G
        
        # because if you have g-1 parameters, you can do 1-these to get the last one
        pi_para[g,qvalue,numbmodel]<-G-1
        
        if(inputmodel == "CCC"){
          sigma_para[g,qvalue,numbmodel] = ((p*Q)-(Q*(Q-1))/2)+1
        } else if (inputmodel == "CCU"){
          sigma_para[g,qvalue,numbmodel] = ((p*Q)-(Q*(Q-1))/2)+p
        } else if (inputmodel == "CUC") {
          sigma_para[g,qvalue,numbmodel] = ((p*Q)-(Q*(Q-1))/2)+G
        } else if (inputmodel == "CUU") {
          sigma_para[g,qvalue,numbmodel] = ((p*Q)-(Q*(Q-1))/2)+(G*p)
        } else if (inputmodel == "UCC") {
          sigma_para[g,qvalue,numbmodel] = (G*((p*Q)-(Q*(Q-1))/2))+1
        } else if (inputmodel == "UCU") {
          sigma_para[g,qvalue,numbmodel] = (G*((p*Q)-(Q*(Q-1))/2))+p
        } else if (inputmodel == "UUC") {
          sigma_para[g,qvalue,numbmodel] = (G*((p*Q)-(Q*(Q-1))/2))+G
        } else if (inputmodel == "UUU") {
          sigma_para[g,qvalue,numbmodel] = (G*((p*Q)-(Q*(Q-1))/2))+(G*p)
        }
        
        # number of parameters
        k[g,qvalue,numbmodel]<-mu_para[g,qvalue,numbmodel]+sigma_para[g,qvalue,numbmodel]+pi_para[g,qvalue,numbmodel] # total parameters are
        
        # keeping track of g for slope heuristics
        gchoice[g,qvalue,numbmodel]<-g 
        
      } # end of numbmodel loop
    } # end of qvalue loop
  } # end of g loop    
  
  results_parameters<-list(parameters= k,
                           gchoice = gchoice,
                           loglikelihood = ll,
                           totaltime = totaltime,
                           n = obs)
  
  class(results_parameters) <- "MPLN"
  return(results_parameters)
  
}

fun_FA<-function(G,Sg,psi,lambda,ng,p,q,modelName){
  n<-sum(ng)
  d<-p
  
  ###Specifying empty list needed at some point  
  gamma<-list()
  bigTheta<-list()
  lambdanew<-list()
  psinew<-list()
  Sigma<-list() 
  
  fun_sgav<-function(g){ng[g]/n*Sg[[g]]}
  Sg_av<-matrix(rowSums(sapply(1:G,fun_sgav)),d,d)
  
  ###################Functions for decomposition of Covariance matrix#
  ###Lambda for UUU,UUC,UCU,UCC#
  fun_lambda_g<-function(g){Sg[[g]]%*%t(gamma[[g]])%*%solve(bigTheta[[g]])}
  #####Lambda for CUC #
  fun_lambda_cuc_1<-function(g){ng[g]/(diag(psi[[g]])[1])*(Sg[[g]]%*%t(gamma[[g]]))}	
  fun_lambda_cuc_2<-function(g){ng[g]/(diag(psi[[g]])[1])*bigTheta[[g]]}
  fun_lambda_ccu_ccc<-function(){Sg_av%*%t(gamma[[1]])%*%solve(bigTheta_av)}
  fun_lambda_cuu_1<-function(g){ng[g]*solve(psi[[g]])%*%Sg[[g]]%*%t(gamma[[g]])}
  fun_lambda_cuu_2<-function(g){ng[g]*diag(solve(psi[[g]]))[k]*bigTheta[[g]]}
  fun_sgav<-function(g){ng[g]/n*Sg[[g]]}
  
  
  
  #############Updates for Psi matrix / vector #
  fun_psi_uuu<-function(g){diag(Sg[[g]]-lambdanew[[g]]%*%gamma[[g]]%*%Sg[[g]])*diag(d)}
  fun_psi_uuc<-function(g){mean(diag(Sg[[g]]-lambdanew[[g]]%*%gamma[[g]]%*%Sg[[g]]))}
  fun_psi_ucu<-function(g){ng[g]/n*diag(Sg[[g]]-lambdanew[[g]]%*%gamma[[g]]%*%Sg[[g]])}
  fun_psi_ucc<-function(g){1/d*ng[g]/n*sum(diag(Sg[[g]]-lambdanew[[g]]%*%gamma[[g]]%*%Sg[[g]]))}
  fun_psi_cuc<-function(g){mean(diag(Sg[[g]]-2*lambdanew[[1]]%*%gamma[[g]]%*%Sg[[g]]+lambdanew[[1]]%*%bigTheta[[g]]%*%t(lambdanew[[1]])))}
  fun_psi_cuu<-function(g){diag(Sg[[g]]-2*lambdanew[[1]]%*%gamma[[g]]%*%Sg[[g]]+lambdanew[[1]]%*%bigTheta[[g]]%*%t(lambdanew[[1]]))}
  fun_psi_ccu<-function(){diag(Sg_av-lambdanew[[1]]%*%gamma[[1]]%*%Sg_av)}
  fun_psi_ccc<-function(){mean(diag(Sg_av-lambdanew[[1]]%*%gamma[[1]]%*%Sg_av))}	
  
  
  ##########Model CCC#
  if (modelName=="CCC"){	
    #Updating Lambda and Gamma
    gamma[[1]]<-t(lambda[[1]])%*%solve(lambda[[1]]%*%t(lambda[[1]])+psi[[1]])
    bigTheta[[1]]<-diag(q)-gamma[[1]]%*%lambda[[1]]+gamma[[1]]%*%Sg_av%*%t(gamma[[1]])
    bigTheta_av<-bigTheta[[1]]
    ####Lambda for CCU, CCC #
    lambdanew[[1]]<-fun_lambda_ccu_ccc()
    psinew[[1]]<-fun_psi_ccc()*diag(d)
    #Complete Sigma	
    for (g in 1:G){	
      Sigma[[g]]<-(lambdanew[[1]]%*%t(lambdanew[[1]])+psinew[[1]])
    }
    
    
    psi[[1]]<-psinew[[1]]	
    lambda[[1]]<-lambdanew[[1]]		
    par<-(d*q-0.5*q*(q-1))+1
  }
  ###
  if (modelName=="UUU"){
    #Updating Lambda and Gamma
    for (g in 1:G){
      gamma[[g]]<-t(lambda[[g]])%*%solve(lambda[[g]]%*%t(lambda[[g]])+psi[[g]])
      bigTheta[[g]]<-diag(q)-gamma[[g]]%*%lambda[[g]]+gamma[[g]]%*%Sg[[g]]%*%t(gamma[[g]])
      lambdanew[[g]]<-fun_lambda_g(g)
      psinew[[g]]<-fun_psi_uuu(g)
      #Complete Sigma		
      Sigma[[g]]<-(lambdanew[[g]]%*%t(lambdanew[[g]])+psinew[[g]])
      lambda[[g]]<-lambdanew[[g]]
      psi[[g]]<-psinew[[g]]			
    }
    par<-G*(d*q-0.5*q*(q-1))+G*d
  }
  ###
  
  ##########Model UUC#
  if (modelName=="UUC"){
    #Updating Lambda and Gamma
    for (g in 1:G){
      gamma[[g]]<-t(lambda[[g]])%*%solve(lambda[[g]]%*%t(lambda[[g]])+psi[[g]])
      bigTheta[[g]]<-diag(q)-gamma[[g]]%*%lambda[[g]]+gamma[[g]]%*%Sg[[g]]%*%t(gamma[[g]])
      lambdanew[[g]]<-fun_lambda_g(g)
      psinew[[g]]<-fun_psi_uuc(g)*diag(d)
      #Complete Sigma		
      Sigma[[g]]<-(lambdanew[[g]]%*%t(lambdanew[[g]])+psinew[[g]])
      lambda[[g]]<-lambdanew[[g]]
      psi[[g]]<-psinew[[g]]			
    }
    par<-G*(d*q-0.5*q*(q-1))+G
  }
  ###
  
  ##########Model UCU#
  if (modelName=="UCU"){
    #Updating Lambda and Gamma
    for (g in 1:G){
      gamma[[g]]<-t(lambda[[g]])%*%solve(lambda[[g]]%*%t(lambda[[g]])+psi[[1]])
      bigTheta[[g]]<-diag(q)-gamma[[g]]%*%lambda[[g]]+gamma[[g]]%*%Sg[[g]]%*%t(gamma[[g]])
      lambdanew[[g]]<-fun_lambda_g(g)
    }		
    psinew[[1]]<-rowSums(sapply(1:G,fun_psi_ucu))*diag(d)
    for (g in 1:G){		
      #Complete Sigma		
      Sigma[[g]]<-(lambdanew[[g]]%*%t(lambdanew[[g]])+psinew[[1]])
      lambda[[g]]<-lambdanew[[g]]
    }
    psi[[1]]<-psinew[[1]]			
    par<-G*(d*q-0.5*q*(q-1))+d
  }
  
  ###
  
  ##########Model UCU#
  if (modelName=="UCC"){
    #Updating Lambda and Gamma
    for (g in 1:G){
      gamma[[g]]<-t(lambda[[g]])%*%solve(lambda[[g]]%*%t(lambda[[g]])+psi[[1]])
      bigTheta[[g]]<-diag(q)-gamma[[g]]%*%lambda[[g]]+gamma[[g]]%*%Sg[[g]]%*%t(gamma[[g]])
      lambdanew[[g]]<-fun_lambda_g(g)
    }		
    psinew[[1]]<-sum(sapply(1:G,fun_psi_ucc))*diag(d)
    for (g in 1:G){		
      #Complete Sigma		
      Sigma[[g]]<-(lambdanew[[g]]%*%t(lambdanew[[g]])+psinew[[1]])
      lambda[[g]]<-lambdanew[[g]]
    }
    psi[[1]]<-psinew[[1]]			
    par<-G*(d*q-0.5*q*(q-1))+1
  }
  #
  
  
  ##########Model CUU#
  if (modelName=="CUU"){
    #Updating Lambda and Gamma
    for (g in 1:G){
      gamma[[g]]<-t(lambda[[1]])%*%solve(lambda[[1]]%*%t(lambda[[1]])+psi[[g]])
      bigTheta[[g]]<-diag(q)-gamma[[g]]%*%lambda[[1]]+gamma[[g]]%*%Sg[[g]]%*%t(gamma[[g]])
    }
    for_lam_1<-matrix(rowSums(sapply(1:G,fun_lambda_cuu_1)),d,q)
    for_lam_2<-list()
    for (k in 1:d){
      if (q>1) for_lam_2[[k]]<-solve(matrix(rowSums(sapply(1:G,fun_lambda_cuu_2)),q,q)) else for_lam_2[[k]]<-solve(matrix(sum(sapply(1:G,fun_lambda_cuu_2)),q,q))
    }
    lambdanew[[1]]<-matrix(NA,d,q)
    for (k in 1:d){
      lambdanew[[1]][k,]<-for_lam_1[k,]%*%for_lam_2[[k]]
    }
    for (g in 1:G){
      psinew[[g]]<-fun_psi_cuu(g)*diag(d)
      #Complete Sigma		
      Sigma[[g]]<-(lambdanew[[1]]%*%t(lambdanew[[1]])+psinew[[g]])
      lambda[[1]]<-lambdanew[[1]]
      psi[[g]]<-psinew[[g]]			
    }
    par<-(d*q-0.5*q*(q-1))+G*d
  }
  #
  ##########Model CUC#
  if (modelName=="CUC"){
    #Updating Lambda and Gamma
    for (g in 1:G){
      gamma[[g]]<-t(lambda[[1]])%*%solve(lambda[[1]]%*%t(lambda[[1]])+psi[[g]])
      bigTheta[[g]]<-diag(q)-gamma[[g]]%*%lambda[[1]]+gamma[[g]]%*%Sg[[g]]%*%t(gamma[[g]])
    }	
    # print("lambda before update")
    # print(lambda)	
    if (q>1){	
      lambdanew[[1]]<-matrix(rowSums(sapply(1:G,fun_lambda_cuc_1)),d,q)%*%solve(matrix(rowSums(sapply(1:G,fun_lambda_cuc_2)),q,q)	)
    } else lambdanew[[1]]<-matrix(rowSums(sapply(1:G,fun_lambda_cuc_1)),d,q)%*%solve(matrix(sum(sapply(1:G,fun_lambda_cuc_2)),q,q)	)
    lambda[[1]]<-lambdanew[[1]]	
    
    
    ###
    for (g in 1:G){	
      psinew[[g]]<-fun_psi_cuc(g)*diag(d)
      #Complete Sigma		
      Sigma[[g]]<-(lambdanew[[1]]%*%t(lambdanew[[1]])+psinew[[g]])
      psi[[g]]<-psinew[[g]]	
    }
    
    par<-(d*q-0.5*q*(q-1))+G
  }
  ###
  
  
  ##########Model CCU#
  if (modelName=="CCU"){
    #Updating Lambda and Gamma
    gamma[[1]]<-t(lambda[[1]])%*%solve(lambda[[1]]%*%t(lambda[[1]])+psi[[1]])
    fun_sgav<-function(g){ng[g]/n*Sg[[g]]}
    Sg_av<-matrix(rowSums(sapply(1:G,fun_sgav)),d,d)
    bigTheta[[1]]<-diag(q)-gamma[[1]]%*%lambda[[1]]+gamma[[1]]%*%Sg_av%*%t(gamma[[1]])
    bigTheta_av<-bigTheta[[1]]
    ####Lambda for CCU, CCC #
    lambdanew[[1]]<-fun_lambda_ccu_ccc()
    psinew[[1]]<-fun_psi_ccu()*diag(d)
    #Complete Sigma	
    for (g in 1:G){	
      Sigma[[g]]<-(lambdanew[[1]]%*%t(lambdanew[[1]])+psinew[[1]])
    }
    psi[[1]]<-psinew[[1]]	
    lambda[[1]]<-lambdanew[[1]]		
    par<-(d*q-0.5*q*(q-1))+d
  }
  ###
  FAresults<-list(Psi=psi,Lambda=lambda,Sigma=Sigma)
  class(FAresults) <- "FactorAnalyzer"
  return(FAresults)
}

cluster_mpln<-function(y,z,G,n_chain,numb_iterations, initialization, normalizefactors, modelname, q){
  
  d<-ncol(y)
  n<-nrow(y)
  
  obs = PI = norm_mu_outer = logL = vector()
  median_mu_outer = median_Sigma_outer = norm_Sigma_outer = list()
  mu_all_outer = Sg_FA = Sigma = Sigma_all_outer = Psi_all_outer = Lambda_all_outer = list()
  it_outer<-2 # the starting value of interation for outer loop
  conv_outer<-0
  
  # for initializing loading matrix (Lambda) = p*q matrix
  if (substr(modelname, 1, 1) == "U"){ 
    mat<-matrix(ncol=q, nrow=d)
    for (model in 1:q){
      mat[,model]= runif(d, min=model-1, max=model)}
    lambda = rep(list(mat), G)
  } else if(substr(modelname, 1, 1) == "C"){
    mat<-matrix(ncol=q, nrow=d)
    for (model in 1:q){
      mat[,model]= runif(d, min=model-1, max=model)}
    lambda = rep(list(mat), G=1)}
  
  if (substr(modelname, 2, 2) == "U"){ # for initializing psi error variance
    psi <- rep(list(diag(d)), G) 
    if (substr(modelname, 3, 3) == "U"){ # for initializing psi isotropic
      for (k in 1:G){
        diag(psi[[k]]) <- 1:d*k
      }
      # if substr(modelnames, 3, 3) == "C", no need to do anything as
      # diagonal elements will be the same anyway
    }
  } else if(substr(modelname, 2, 2) == "C"){ # for initializing psi error variance
    psi <- rep(list(diag(d)), G=1) 
    if (substr(modelname, 3, 3) == "U"){ # for initializing psi isotropic
      diag(psi[[1]]) <- 1:d}
    # if substr(modelnames, 3, 3) == "C", no need to do anything as
    # diagonal elements will be the same anyway
  }
  
  # for initialization of Sigma
  for (gcurrent in 1:G){ 
    if (substr(modelname, 1, 1) == "U"){ # indexing lambda
      lam <- gcurrent
    } else if(substr(modelname, 1, 1) == "C"){
      lam <- 1 }
    
    if (substr(modelname, 2, 2) == "U"){ # indexing psi 
      psii<-gcurrent
    } else if(substr(modelname, 2, 2) == "C"){ 
      psii<-1 }
    
    Sigma[[gcurrent]]<-(lambda[[lam]]%*%t(lambda[[lam]]))+psi[[psii]] # initializing Sg
  }
  
  # if initialization hasn't been done
  if (all(is.na(initialization))==TRUE || all(initialization =="init")){
    mu_all_outer[[1]]<-mu_g<- matrix(log(mean(y)), ncol=d, nrow=G) 
    Sigma_all_outer[[1]]<-Sigma
    Psi_all_outer[[1]]<-psi
    Lambda_all_outer[[1]]<-lambda
  } else{  
    # if initialization has been done
    logL[1]<-calc_likelihood(z=initialization$probaPost, PI=initialization$proportion, y=y, mu_g=initialization$finalmu, G=1, Sig_g=do.call("rbind",initialization$finalsigma), theta_Stan=initialization$thetaRStan, normalizefactors=normalizefactors)
    mu_all_outer[[1]]<-mu_g<-initialization$finalmu
    Sigma_all_outer[[1]]<-initialization$finalsigma
    Psi_all_outer[[1]]<-initialization$finalpsi
    Lambda_all_outer[[1]]<-initialization$finallambda 
    z=initialization$probaPost
    numb_iterations = initialization$FinalRstan_iterations
  }
  
  mod = stan_model("MPLN.stan")
  
  while(!conv_outer){ 
    for(g in 1:G){
      obs[g]<-sum(z[,g]) # number of observations in each group
      PI[g]<-obs[g]/n  # obtain probability of each group
    }
    
    theta_Stan<-E_theta2<-list()
    rstan_results<-stanrun(model=mod, gmin=1,gmax=G,y=y,mu_all_outer=mu_all_outer, it_outer=it_outer, sigma_all_outer=Sigma_all_outer, numb_iterations=numb_iterations, n_chain=n_chain, normalizefacs=normalizefactors)
    
    fit = rstan_results$fitrstan
    numb_iterations = rstan_results$numb_iterations
    
    
    for (g in 1:G){
      tt<-as.matrix(fit[[g]])
      theta_Stan[[g]]<-matrix(NA,nrow=n,ncol=d)
      E_theta2[[g]]<-list()
      
      for (i in 1:n){
        zz<-c(1:(d-1))*n+i
        theta_mat<-tt[,c(i,zz)]
        theta_Stan[[g]][i,]<-colMeans(theta_mat-normalizefactors)
        E_theta2[[g]][[i]]<-z[i,g]*t(tt[,c(i,zz)])%*%tt[,c(i,zz)]/((0.5*numb_iterations)*n_chain)
      } 
      
      
      mu_g[g,]<-colSums(z[,g]*theta_Stan[[g]])/sum(z[,g])
      Sg_FA[[g]]<-Reduce("+",E_theta2[[g]])/sum(z[,g])-mu_g[g,]%*%t(mu_g[g,])
      colnames(Sg_FA[[g]])<-paste(1:d)
      rownames(Sg_FA[[g]])<-paste(1:d)
    }
    
    # calculating Sigma, Lambda and Psi
    calcfact<-fun_FA(G=G,Sg=Sg_FA,ng=colSums(z),psi=Psi_all_outer[[it_outer-1]],lambda=Lambda_all_outer[[it_outer-1]],p=ncol(y),q=q,modelName=modelname)
    
    
    mu_all_outer[[it_outer]]<-mu_g # save the final converged mu 
    Sigma_all_outer[[it_outer]]<-calcfact$Sigma # save the final converged sigma
    Psi_all_outer[[it_outer]]<-calcfact$Psi
    Lambda_all_outer[[it_outer]]<-calcfact$Lambda
    
    logL[it_outer]<-calc_likelihood(z=z, PI=PI, y=y, mu_g=mu_all_outer[[it_outer]], G=G, Sig_g=do.call("rbind", Sigma_all_outer[[it_outer]]), theta_Stan=theta_Stan, normalizefactors=normalizefactors)
    
    
    # convergence of outer loop
    norm_mu_outer[it_outer]<-norm((mu_all_outer[[it_outer]]-mu_all_outer[[it_outer-1]]),type="F")
    norm_Sigma_outer[[it_outer]]<-list()
    for (gcurrent in 1:G){ 
      norm_Sigma_outer[[it_outer]][[gcurrent]]<-norm(matrix(unlist(mapply('-',Sigma_all_outer[[it_outer]],Sigma_all_outer[[it_outer-1]],SIMPLIFY=FALSE)),ncol=ncol(y), nrow=ncol(y)),type="F")
    }
    
    median_mu_outer[[it_outer]]<-median(norm_mu_outer[2:it_outer])
    median_Sigma_outer[[it_outer]]<-list()
    for (gcurrent in 1:g){ 
      median_Sigma_outer[[it_outer]][[gcurrent]]<-median(norm_Sigma_outer[[2]][[gcurrent]]:norm_Sigma_outer[[it_outer]][[gcurrent]])
    }
    
    threshold_outer<-10
    
    # check for convergence
    if(it_outer>(threshold_outer+1)){
      
      check2<-vector()
      for (gcurrent in 1:G){ 
        check2[gcurrent]<-(abs(median_Sigma_outer[[it_outer-threshold_outer]][[gcurrent]]-median_Sigma_outer[[it_outer]][[gcurrent]]))
      }
      if((abs(median_mu_outer[[it_outer-threshold_outer]]-median_mu_outer[[it_outer]])<5 && (all(check2<5)==TRUE) ) || it_outer>100){
        programclust<-vector()
        programclust<-map(z)
        
        # checking for spurious clusters and getting rid of them
        #keep<-as.numeric(names(which(table(programclust)>5))) 
        #if (length(keep) !=length(unique(programclust))){
        #  z<-as.matrix(z[,keep])
        #  z<-z/rowSums(z)
        #  programclust<-map(z)
        #}
        
        # checking for empty clusters
        J <- 1:ncol(z)
        K <- as.logical(match(J, sort(unique(programclust)), nomatch = 0))
        if(length(J[!K])>0){ # J[!K] tells which are empty clusters
          z<-z[,-J[!K]]
          programclust<-map(z)
        }
        
        results <- list(finalmu = mu_all_outer[[it_outer]]+ matrix(rep(normalizefactors,nrow(mu_all_outer[[it_outer]])),byrow=TRUE,ncol=ncol(mu_all_outer[[it_outer]])), 
                        finalsigma = Sigma_all_outer[[it_outer]], 
                        finallambda = Lambda_all_outer[[it_outer]],
                        finalpsi = Psi_all_outer[[it_outer]],
                        allmu = lapply(mu_all_outer, function(x) (x+matrix(rep(normalizefactors,nrow(mu_all_outer[[it_outer]])),byrow=TRUE,ncol=ncol(mu_all_outer[[it_outer]])))), 
                        allsigma = Sigma_all_outer, 
                        alllambda = Lambda_all_outer,
                        allpsi = Psi_all_outer,
                        clusterlabels = programclust,
                        iterations = it_outer, 
                        FinalRstan_iterations = numb_iterations,
                        proportion = PI, 
                        loglikelihood = logL,
                        probaPost = z,
                        thetaRStan = theta_Stan,
                        q = q)
        
        conv_outer<-1
      } # end of conv_outer == 1 loop
    } # end of it_outer>20 loop
    
    # if running for initialization, need to stop after 10 iterations
    if(it_outer==10 && all(is.na(initialization) !=TRUE)){
      if(all(initialization == "init")){
        programclust<-vector()
        programclust<-map(z)
        
        # will not check for spurious and empty clusters 
        results <- list(finalmu=mu_all_outer[[it_outer]]+ matrix(rep(normalizefactors,nrow(mu_all_outer[[it_outer]])),byrow=TRUE,ncol=ncol(mu_all_outer[[it_outer]])), 
                        finalsigma=Sigma_all_outer[[it_outer]], 
                        finallambda = Lambda_all_outer[[it_outer]],
                        finalpsi = Psi_all_outer[[it_outer]],
                        allmu = lapply(mu_all_outer, function(x) (x+matrix(rep(normalizefactors,nrow(mu_all_outer[[it_outer]])),byrow=TRUE,ncol=ncol(mu_all_outer[[it_outer]])))), 
                        allsigma = Sigma_all_outer, 
                        alllambda = Lambda_all_outer,
                        allpsi = Psi_all_outer,
                        clusterlabels = programclust,
                        iterations = it_outer, 
                        FinalRstan_iterations = numb_iterations,
                        proportion = PI, 
                        loglikelihood = logL,
                        probaPost = z,
                        thetaRStan = theta_Stan,
                        thetaRStan = theta_Stan,
                        q=q)
        
        conv_outer<-1
      }
    }
    
    if(conv_outer!=1){ # only update until outer convergence is not reached, not after
      z<-zvalue_calculation(theta_Stan=theta_Stan,y=y,G=G,mu_g=mu_g,Sig_g=do.call("rbind", Sigma_all_outer[[it_outer]]),PI=PI, normalizefactors=normalizefactors)
      it_outer<-it_outer+1 # updating outer loop iteration
      numb_iterations = numb_iterations+10
    }
  }# ending while loop
  
  if(all(is.na(initialization))==TRUE){
    class(results) <- "cluster_initMPLN"
  }else{
    class(results) <- "cluster_MPLN"
  }
  return(results)
}    

main_mpln<-function(y, G, n_chain, numb_iterations=NA, membership=NA, init_method=NA, init_iterations=NA, normalize=NA, modelname=NA, q=NA){
  ptm<-proc.time()  
  
  if (typeof(y) != "double" & typeof(y) != "integer"){
    stop("Dataset type needs to be integer");}
  
  d<-ncol(y)
  n<-nrow(y)
  
  if(length(which(apply(y, 1, function(x) all(x==0))==TRUE))!=0){
    cat("\nDataset row(s)", c(which(apply(y, 1, function(x) all(x==0))==TRUE)), "will be removed as this/these contain(s) all zeros")
    if(all(is.na(membership)==FALSE)){membership<-membership[-c(which(apply(y, 1, function(x) all(x==0))==TRUE))]}
    y<-y[-c(which(apply(y, 1, function(x) all(x==0))==TRUE)),]
    n<-nrow(y)
  }
  
  if(all(is.na(membership)==TRUE)){
    membership<-"Not provided"}
  
  # generating normalization factors
  if(is.na(normalize) == FALSE) {
    if (!require(edgeR)) install.packages('edgeR') # loading needed package
    library(edgeR)
    norm_factors<-log(as.vector(calcNormFactors(as.matrix(y), method = "TMM")))
  } else {norm_factors<-rep(0,d)}
  #cat("\nNormalize factors in main_mpln are: ",norm_factors)
  
  
  # clustering
  if(init_iterations!=0){
    initializeruns=initializationrun(gmodel=G, y=y, init_method=init_method, init_iterations=init_iterations, n_chain=n_chain, numb_iterations=numb_iterations, initialization=NA, normalizefactors=norm_factors, modelname=modelname, q=q)
    allruns=cluster_mpln(y=y,z=NA,G=G,n_chain=n_chain,numb_iterations=NA, initialization=initializeruns,normalizefactors=norm_factors, modelname=modelname, q=q)
  }else if(init_iterations == 0){
    allruns=cluster_mpln(y=y,z=unmap(kmeans(log(y+1/3),G)$cluster),G=G, n_chain=n_chain,numb_iterations=numb_iterations, initialization=NA, normalizefactors=norm_factors, modelname=modelname, q=q)
  }
  
  final = proc.time()-ptm
  RESULTS = list(dataset = y,
                 dimensionality = d,
                 normalization_factors = norm_factors,
                 G = G,
                 initalization_method = init_method,
                 allresults = allruns,
                 truelabels = membership,
                 totaltime = final)
  
  class(RESULTS) = "MPLN"
  return(RESULTS)
}

downstream_analysis = function(dataset, G=G, modelnames, models, q){  
  ptm = proc.time() 
  Gmax = max(G)
  Gmin = min(G)
  modelnames = modelnames
  qmax = max(q)
  qmin = min(q)
  
  calpara<-calculate_parameters(Gmin=Gmin, Gmax=Gmax, qmin=qmin, qmax=qmax, modelnames=modelnames, models=models, allruns=dataset)
  k<-calpara$parameters
  gchoice<-calpara$gchoice
  ll = calpara$loglikelihood
  totaltime_eachg = calpara$totaltime
  n = calpara$n
  
  # starting model selection
  bic<-BIC_function(ll=ll, k=k, n=n, run=dataset, qmin=qmin, qmax=qmax, gmin=Gmin, gmax=Gmax, modelnames=modelnames, models=models)
  icl<-ICL_function(bIc=bic, run=dataset, qmin=qmin, qmax=qmax, gmin=Gmin, gmax=Gmax, modelnames=modelnames, models=models)
  aic<-AIC_function(ll=ll,k=k, run=dataset, qmin=qmin, qmax=qmax, gmin=Gmin, gmax=Gmax, modelnames=modelnames, models=models)
  aic3<-AIC3_function(ll=ll,k=k, run=dataset, qmin=qmin, qmax=qmax, gmin=Gmin, gmax=Gmax, modelnames=modelnames, models=models)
  
  # for Djump and DDSE
  if((Gmax-Gmin+1) > 10 ) {
    if (!require(capushe)) install.packages('capushe') # loading needed package
    library(capushe)
    
    # adapted based on HTSCluster package 2.0.8 (25 Oct 2016)
    
    message("Note: diagnostic plots for results corresponding to model selection via slope heuristics (Djump and DDSE) should be examined to ensure that sufficiently complex models have been considered.")
    
    for(g in 1:(Gmax-Gmin+1)) {
      for(qvalue in 1:(qmax-qmin+1)){
        for (numbmodel in 1:length(models)){
          mat <- cbind(cbind(gchoice[g,qvalue,numbmodel]), cbind(k[g,qvalue,numbmodel])/n[g,qvalue,numbmodel], cbind(k[g,qvalue,numbmodel])/n[g,qvalue,numbmodel], cbind(ll[g,qvalue,numbmodel]))
          ResCapushe <- capushe(mat, n)
          DDSEmodel<- ResCapushe@DDSE@model
          Djumpmodel<- ResCapushe@Djump@model
        }      
      }
    }
    
    final = (proc.time()-ptm) + Reduce("+",totaltime_eachg)
    
    RESULTS <- list(dataset = dataset[[1]][[1]][[1]]$dataset, 
                    normalization_factors = dataset[[1]][[1]][[1]]$normalization_factors,
                    dimensionality = ncol(dataset[[1]][[1]][[1]]$dataset), 
                    gmin = Gmin,
                    gmax = Gmax,
                    qmin = qmin,
                    qmax = qmax,
                    initalization = dataset[[1]][[1]][[1]]$initalization_method, 
                    allresults = dataset,
                    loglikelihood = ll, 
                    numbofparameters = k,
                    truelabels = dataset[[1]][[1]][[1]]$truelabels,
                    models = c(modelnames[models]),
                    ICL.all = icl,
                    BIC.all = bic,
                    AIC.all = aic,
                    AIC3.all = aic3,
                    SlopeHeuristics = ResCapushe,
                    Djumpmodelselected = ResCapushe@Djump@model,
                    DDSEmodelselected = ResCapushe@DDSE@model,
                    totaltime = final)
    
  }else if((Gmax-Gmin+1) < 10 ){# end of Djump and DDSE
    final = (proc.time()-ptm) + Reduce("+",totaltime_eachg)
    RESULTS <- list(dataset = dataset[[1]][[1]][[1]]$dataset,
                    normalization_factors = dataset[[1]][[1]][[1]]$normalization_factors,
                    dimensionality = ncol(dataset[[1]][[1]][[1]]$dataset),
                    gmin = Gmin,
                    gmax = Gmax,
                    qmin = qmin,
                    qmax = qmax,
                    initalization = dataset[[1]][[1]][[1]]$initalization_method,
                    allresults = dataset,
                    loglikelihood = ll, 
                    numbofparameters = k,
                    truelabels = dataset[[1]][[1]][[1]]$truelabels, 
                    models = c(modelnames[models]),
                    ICL.all = icl,
                    BIC.all = bic,
                    AIC.all = aic,
                    AIC3.all = aic3,
                    SlopeHeuristics = "Not used",
                    Djumpmodelselected = "Not used",
                    DDSEmodelselected = "Not used",
                    totaltime = final)
  }
  
  class(RESULTS) <- "MPLN_FA"
  return(RESULTS)
} # end of MPLN function

final_results = function(allresults, G, modelnames, q, models){
  
  eachdataset <- vector("list", length = length(allresults)) # obtaining number of datasets
  names(eachdataset) <- paste("dataset", 1:length(allresults), sep = "")
  
  for (number in 1:length(allresults)){
    eachdataset[[number]] = downstream_analysis(dataset=allresults[[number]], G=G, modelnames=modelnames, models=models, q=q)
  }
  
  class(eachdataset) <- "MPLN_finalresults"
  return(eachdataset)
}


#### Running the model ####

G_testing = c(8:10) # total number of cluster sizes to test 
q_testing = c(1:3) # total number of latent factors
models_testing = c(4) # total number of models
modelnames_consider=c("CCC","CCU","CUC","CUU","UCC","UCU","UUC","UUU") # identitiy of model


library(foreach)
opts <- list(chunkSize=12)
results2 <-
  foreach(a=G_testing, .options.nws=opts) %:% # for the cluster size G
  foreach(c=q_testing, .options.nws=opts) %:% # for each q
  foreach(m=models_testing, .inorder=FALSE) %dopar% { # for each model
    set.seed(1)
    obj =  main_mpln(y=bi2016, G=a, n_chain=3, numb_iterations=1000, membership=NA, init_method="kmeans", init_iterations=5, normalize="TMM", modelname=modelnames_consider[m], q=c)
  }
save.image("FAMPLN_4Nov2017_nestedparallel_bi2016data_g=8to10_4.RData")

#results
#finalresults = final_results(allresults=results2, G=G_testing, modelnames=modelnames_consider, models=models_testing, q=q_testing)
#save.image("FAMPLN_26Oct2017_nestedparallel_bi2016data_g=4to7_2.RData")



#run_FA_g2q2CCC<-list()
#for (i in 1:1){
#  set.seed(i)
#  run_FA_g2q2CCC[[i]]<-main_mpln(y=datasets_FA[[1]]$dataset, G=2, n_chain=1, numb_iterations=70, membership=datasets_FA[[1]]$truemembership, init_method="kmeans", init_iterations=1, normalize="TMM", modelnames=c("CCC","CCU","CUC","CUU","UCC","UCU","UUC","UUU"), qmin=2, qmax=2)
#  save.image(file="FAMPLN_27May2017_g2q2CCC.RData")
#} 

#for (i in 1:10){
#  cat("\nRun", i)
#  cat("\nAdjusted Rand Index is", adjustedRandIndex(try_nointernal[[i]]$truelabels, try_nointernal[[i]]$BIC.all$BICmodelselected_labels), "\n")
#  print("Mean is")
#  print(try_nointernal[[i]]$allresults$`g=2`$`q=1`$`model=UUU`$finalmu)
#  print("Sigma is")
#  print(try_nointernal[[i]]$allresults$`g=2`$`q=1`$`model=UUU`$finalsigma)
#  cat("Normalization factors are: ", try_nointernal[[i]]$normalization_factors, "\n")
#  print(try_nointernal[[i]]$totaltime)
#}

# checking clustering of labels
#table(run_FA[[2]]$truelabels,datasets_FA[[1]]$truemembership)

