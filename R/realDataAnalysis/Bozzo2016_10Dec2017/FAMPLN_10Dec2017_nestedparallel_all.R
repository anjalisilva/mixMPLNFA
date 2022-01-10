# 10 Dec November 2017
# testing for G=2, q=1:3, models=(2)

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

bozzo2016<-as.matrix(read.csv("bozzo_mediancounts.csv", header = TRUE, row.names=1))
dim(bozzo2016) # 1336 6

which(rowSums(bozzo2016)==0) #none = no rows are all zero


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

initialization<-function(mod, gmodel, y, init_method, init_iterations, n_chain, numb_iterations, normalizefactors, modelname, q){
  z<-init_runs<-list()
  logL_init<-vector()
  n<-nrow(y)
  d<-ncol(y)
  
  for(iterations in 1:init_iterations){
    set.seed(iterations)
    
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
    
    init_runs[[iterations]]=initialization_run(y=y,z=z[[iterations]],G=gmodel,n_chain=n_chain,numb_iterations=numb_iterations, initialization="init", normalizefactors=normalizefactors, modelname=modelname, q=q)
    logL_init[iterations] <- unlist(tail((init_runs[[iterations]]$loglikelihood), n=1)) 
  }
  
  initialization<-init_runs[[which(logL_init==max(logL_init, na.rm = TRUE))[1]]]
  return(initialization)
}  

initialization_run<-function(mod, y,z,G,n_chain,numb_iterations, initialization, normalizefactors, modelname, q){
  
  d<-ncol(y)
  n<-nrow(y)
  mod = stan_model("MPLN.stan")
  
  obs = PI = norm_mu_outer = logL = vector()
  median_mu_outer = median_Sigma_outer = norm_Sigma_outer = list()
  mu_all_outer = Sg_FA = Sigma = Sigma_all_outer = Psi_all_outer = Lambda_all_outer = list()
  it_outer<-2 # the starting value of interation for outer loop
  conv_outer_init<-0
  
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
  mu_all_outer[[1]]<-mu_g<- matrix(log(mean(y)), ncol=d, nrow=G) 
  Sigma_all_outer[[1]]<-Sigma
  Psi_all_outer[[1]]<-psi
  Lambda_all_outer[[1]]<-lambda
  

  # do only one iteration

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
    
    
    z<-zvalue_calculation(theta_Stan=theta_Stan,y=y,G=G,mu_g=mu_g,Sig_g=do.call("rbind", Sigma_all_outer[[it_outer]]),PI=PI, normalizefactors=normalizefactors)
    
    
     programclust<-vector()
     programclust<-map(z)
        
        
        # will not check for spurious and empty clusters 
        results <- list(finalmu=mu_all_outer[[it_outer]]+ matrix(rep(normalizefactors,nrow(mu_all_outer[[it_outer]])),byrow=TRUE,ncol=ncol(mu_all_outer[[it_outer]])), 
                        finalsigma=Sigma_all_outer[[it_outer]], 
                        finallambda = Lambda_all_outer[[it_outer]],
                        finalpsi = Psi_all_outer[[it_outer]],
                        clusterlabels = programclust,
                        iterations = it_outer, 
                        FinalRstan_iterations = numb_iterations,
                        proportion = PI, 
                        loglikelihood = logL,
                        probaPost = z,
                        thetaRStan = theta_Stan,
                        q=q)
        
  
  class(results) <- "cluster_initMPLN"
  return(results)
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
  BICmodel_labels = run[[which(min(BICvalues) == BICvalues, arr.ind = T)[1]]][[which(min(BICvalues) == BICvalues, arr.ind = T)[3]]][[1]][[which(min(BICvalues) == BICvalues, arr.ind = T)[2]]][[1]]$allresults$clusterlabels # obtaining model labels
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
        
        if(g==1){
          z = matrix(rep(1, 1336),1336,1)
          mapz = mclust::unmap(z)
          forICL<-function(g){sum(log(z[which(mapz[,g]==1),g]))}
          ICLvalues[g,qvalue,numbmodel] <- bIc$allBICvalues[g,qvalue,numbmodel] + sum(sapply(1:ncol(mapz),forICL))
        }else{
          z<-run[[g]][[numbmodel]][[1]][[qvalue]][[1]]$allresults$probaPost
          mapz<-mclust::unmap(run[[g]][[numbmodel]][[1]][[qvalue]][[1]]$allresults$clusterlabels)
          forICL<-function(g){sum(log(z[which(mapz[,g]==1),g]))}
          ICLvalues[g,qvalue,numbmodel] <- bIc$allBICvalues[g,qvalue,numbmodel] + sum(sapply(1:ncol(mapz),forICL))
          
        }
        
      }
    }
  }
  
  colnames(ICLvalues) = c(paste0("q=", qmin:qmax))
  rownames(ICLvalues) = c(paste0("g=", gmin:gmax))
  ICLvalues[which(ICLvalues== -Inf)] = NA
  
  ICLmodel_g = seq(gmin, gmax, 1)[which(min(ICLvalues, na.rm=TRUE) == ICLvalues, arr.ind = T)[1]] # obtaining row (which is g)
  ICLmodel_q = seq(qmin, qmax, 1)[which(min(ICLvalues, na.rm=TRUE) == ICLvalues, arr.ind = T)[2]] # obtaining column (which is q)
  ICLmodel_m = modelnames[models[which(min(ICLvalues, na.rm=TRUE) == ICLvalues, arr.ind = T)[3]]] # obtaining model name
  ICLmodel_labels = run[[which(min(ICLvalues, na.rm=TRUE) == ICLvalues, arr.ind = T)[1]]][[which(min(ICLvalues, na.rm=TRUE) == ICLvalues, arr.ind = T)[3]]][[1]][[which(min(ICLvalues, na.rm=TRUE) == ICLvalues, arr.ind = T)[2]]][[1]]$allresults$clusterlabels
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
  AICmodel_labels = run[[which(min(AICvalues) == AICvalues, arr.ind = T)[1]]][[which(min(AICvalues) == AICvalues, arr.ind = T)[3]]][[1]][[which(min(AICvalues) == AICvalues, arr.ind = T)[2]]][[1]]$allresults$clusterlabels
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
  AIC3model_labels<-run[[which(min(AIC3values) == AIC3values, arr.ind = T)[1]]][[which(min(AIC3values) == AIC3values, arr.ind = T)[3]]][[1]][[which(min(AIC3values) == AIC3values, arr.ind = T)[2]]][[1]]$allresults$clusterlabels
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
  totaltime = vector("list", length = (Gmax-Gmin+1) )
  

    
  for(g in 1:(Gmax-Gmin+1)) {
    totaltime[[g]] = list()
    for(qvalue in 1:(qmax-qmin+1)){
      totaltime[[g]][[qvalue]] = list()
      for (numbmodel in 1:length(models)){
        totaltime[[g]][[qvalue]][[numbmodel]] = list()
        inputmodel = modelnames[models[numbmodel]]

        
        # calculating log-likelihood
        if (g==1){   # from another run which was not run in this format, for G=1 only
          
          totaltime[[g]][[qvalue]][[numbmodel]] = NA
          
          obs[g,qvalue,numbmodel] = 1336
          p = 6
          
          ll[1,1,1] = -42359.63
          ll[1,1,2] = -42359.36   
          ll[1,1,3] = -42359.63
          ll[1,1,4] = -42359.63
          
          ll[1,2,1] = -41176.03
          ll[1,2,2] = -41176.20   
          ll[1,2,3] = -41176.03
          ll[1,2,4] = -41176.03
          
          ll[1,3,1] = -40537.04
          ll[1,3,2] = -40536.76  
          ll[1,3,3] = -40537.04
          ll[1,3,4] = -40537.04
        }else{
          totaltime[[g]][[qvalue]][[numbmodel]] = allruns[[g]][[numbmodel]][[1]][[qvalue]][[1]]$totaltime
          
          obs[g,qvalue,numbmodel] = nrow(allruns[[g]][[numbmodel]][[1]][[qvalue]][[1]]$dataset)
          p = ncol(allruns[[g]][[numbmodel]][[1]][[qvalue]][[1]]$dataset)
          
          ll[g,qvalue,numbmodel]<-unlist(tail(allruns[[g]][[numbmodel]][[1]][[qvalue]][[1]]$allresults$loglikelihood, n=1)) # save the final log-likelihood
        }
        
        if ((Gmax-Gmin+1) == 1){ # if running for only one specific G
          G=Gmin
        }else{  # if running for only multiple G
          G=c(Gmin:Gmax)[g]
        }
        
        if ((qmax-qmin+1) == 1){# if running for only one specific G
          Q=qmin
        }else{ # if running for only multiple G
          Q=c(qmin:qmax)[qvalue]
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
        gchoice[g,qvalue,numbmodel]<-G
        
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

cluster_mpln<-function(mod, y,z,G,n_chain,numb_iterations, initialization, normalizefactors, modelname, q){
  
  d<-ncol(y)
  n<-nrow(y)
  
  mod = stan_model("MPLN.stan")
  
  obs = PI = norm_mu_outer = logL = vector()
  median_mu_outer = median_Sigma_outer = norm_Sigma_outer = list()
  mu_all_outer = Sg_FA = Sigma = Sigma_all_outer = Psi_all_outer = Lambda_all_outer = list()
  it_outer<-2 # the starting value of interation for outer loop
  conv_outer<-0
  

  
  # if initialization hasn't been done
  if (all(is.na(initialization))==TRUE || all(initialization =="init")){
    
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
    
    threshold_outer<-20
    
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
                        q = q)
        
        conv_outer<-1
      } # end of conv_outer == 1 loop
    } # end of it_outer>20 loop
    

    
    if(conv_outer!=1){ # only update until outer convergence is not reached, not after
      z<-zvalue_calculation(theta_Stan=theta_Stan,y=y,G=G,mu_g=mu_g,Sig_g=do.call("rbind", Sigma_all_outer[[it_outer]]),PI=PI, normalizefactors=normalizefactors)
      it_outer<-it_outer+1 # updating outer loop iteration
      numb_iterations = numb_iterations+10
    }
  }# ending while loop
  
  class(results) <- "cluster_MPLN"

  return(results)
}    

main_mpln<-function(mod, y, G, n_chain, numb_iterations=NA, membership=NA, init_method=NA, init_iterations=NA, normalize=NA, modelname=NA, q=NA){
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
    initializeruns=initialization(mod=mod,gmodel=G, y=y, init_method=init_method, init_iterations=init_iterations, n_chain=n_chain, numb_iterations=numb_iterations, normalizefactors=norm_factors, modelname=modelname, q=q)
    allruns=cluster_mpln(mod=mod,y=y,z=NA,G=G,n_chain=n_chain,numb_iterations=NA, initialization=initializeruns,normalizefactors=norm_factors, modelname=modelname, q=q)
  }else if(init_iterations == 0){
    allruns=cluster_mpln(mod=mod,y=y,z=unmap(kmeans(log(y+1/3),G)$cluster),G=G, n_chain=n_chain,numb_iterations=numb_iterations, initialization=NA, normalizefactors=norm_factors, modelname=modelname, q=q)
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
  Gmax = max(G)
  Gmin = min(G)
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
    
    final = totaltime_eachg
    
    RESULTS <- list(dataset = dataset[[Gmax]][[1]][[1]][[1]][[1]]$dataset, 
                    normalization_factors = dataset[[Gmax]][[1]][[1]][[1]][[1]]$normalization_factors,
                    dimensionality = ncol(dataset[[Gmax]][[1]][[1]][[1]][[1]]$dataset), 
                    gmin = Gmin,
                    gmax = Gmax,
                    qmin = qmin,
                    qmax = qmax,
                    initalization = dataset[[Gmax]][[1]][[1]][[1]][[1]]$initalization_method, 
                    allresults = dataset,
                    loglikelihood = ll, 
                    numbofparameters = k,
                    truelabels = dataset[[Gmax]][[1]][[1]][[1]][[1]]$truelabels,
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
    final = totaltime_eachg 
    RESULTS <- list(dataset = dataset[[Gmax]][[1]][[1]][[1]][[1]]$dataset,
                    normalization_factors = dataset[[Gmax]][[1]][[1]][[1]][[1]]$normalization_factors,
                    dimensionality = ncol(dataset[[Gmax]][[1]][[1]][[1]][[1]]$dataset),
                    gmin = Gmin,
                    gmax = Gmax,
                    qmin = qmin,
                    qmax = qmax,
                    initalization = dataset[[Gmax]][[1]][[1]][[1]][[1]]$initalization_method,
                    allresults = dataset,
                    loglikelihood = ll, 
                    numbofparameters = k,
                    truelabels = dataset[[Gmax]][[1]][[1]][[1]][[1]]$truelabels, 
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


#### Running the model ####

mod = stan_model("MPLN.stan")


G_testing = c(1:9) # total number of cluster sizes to test 
q_testing = c(1:3) # total number of latent factors
models = c(2,4,6,8) # total number of models
modelnames=c("CCC","CCU","CUC","CUU","UCC","UCU","UUC","UUU") # identitiy of model


#allresults_FA=list()
allresults_FA[[1]]=list()
allresults_FA[[1]][[1]]=results2 # g=1
allresults_FA[[1]][[2]]=results2
allresults_FA[[1]][[3]]=results2
allresults_FA[[1]][[4]]=results2
allresults_FA[[2]]=list()
allresults_FA[[2]][[1]]=results2 # g=2
allresults_FA[[2]][[2]]=results2
allresults_FA[[2]][[3]]=results2
allresults_FA[[2]][[4]]=results2
allresults_FA[[3]]=list()
allresults_FA[[3]][[1]]=results2 # g= 3
allresults_FA[[3]][[2]]=results2
allresults_FA[[3]][[3]]=results2
allresults_FA[[3]][[4]]=results2
allresults_FA[[4]]=list()
allresults_FA[[4]][[1]]=results2 # g= 4
allresults_FA[[4]][[2]]=results2
allresults_FA[[4]][[3]]=results2
allresults_FA[[4]][[4]]=results2
allresults_FA[[5]]=list()
allresults_FA[[5]][[1]]=results2 # g= 5
allresults_FA[[5]][[2]]=results2
allresults_FA[[5]][[3]]=results2
allresults_FA[[5]][[4]]=results2
allresults_FA[[6]]=list()
allresults_FA[[6]][[1]]=results2 # g=6
allresults_FA[[6]][[2]]=results2
allresults_FA[[6]][[3]]=results2
allresults_FA[[6]][[4]]=results2
allresults_FA[[7]]=list()
allresults_FA[[7]][[1]]=results2 # g=7
allresults_FA[[7]][[2]]=results2
allresults_FA[[7]][[3]]=results2
allresults_FA[[7]][[4]]=results2
allresults_FA[[8]]=list()
allresults_FA[[8]][[1]]=results2 # g=8
allresults_FA[[8]][[2]]=results2
allresults_FA[[8]][[3]]=results2
allresults_FA[[8]][[4]]=results2
allresults_FA[[9]]=list()
allresults_FA[[9]][[1]]=results2 # g=9
allresults_FA[[9]][[2]]=results2
allresults_FA[[9]][[3]]=results2
allresults_FA[[9]][[4]]=results2

# format allruns[[g]][[numbmodel]][[1]][[qvalue]][[1]]
finalresults = downstream_analysis(dataset=allresults_FA, G=c(1:9), modelnames=modelnames, models=models, q=q_testing)
save.image("FAMPLN_10Dec2017_nestedparallel_all.RData")

pairs(log(bozzo2016), col=finalresults$BIC.all$BICmodelselected_labels)
# the log-likelihood values
# format allruns[[g]][[numbmodel]][[1]][[qvalue]][[1]]

# g = 2
plot(finalresults$allresults[[2]][[1]][[1]][[1]][[1]]$allresults$loglikelihood[-1], type="l")
plot(finalresults$allresults[[2]][[1]][[1]][[2]][[1]]$allresults$loglikelihood[-1], type="l")
plot(finalresults$allresults[[2]][[1]][[1]][[3]][[1]]$allresults$loglikelihood[-1], type="l")

plot(finalresults$allresults[[2]][[2]][[1]][[1]][[1]]$allresults$loglikelihood[-1], type="l")
plot(finalresults$allresults[[2]][[2]][[1]][[2]][[1]]$allresults$loglikelihood[-1], type="l")
plot(finalresults$allresults[[2]][[2]][[1]][[3]][[1]]$allresults$loglikelihood[-1], type="l")

plot(finalresults$allresults[[2]][[3]][[1]][[1]][[1]]$allresults$loglikelihood[-1], type="l")
plot(finalresults$allresults[[2]][[3]][[1]][[2]][[1]]$allresults$loglikelihood[-1], type="l")
plot(finalresults$allresults[[2]][[3]][[1]][[3]][[1]]$allresults$loglikelihood[-1], type="l")
  # UUU (4)TH MODEL 
plot(finalresults$allresults[[2]][[4]][[1]][[1]][[1]]$allresults$loglikelihood[-1], type="l")
plot(finalresults$allresults[[2]][[4]][[1]][[2]][[1]]$allresults$loglikelihood[-1], type="l")
plot(finalresults$allresults[[2]][[4]][[1]][[3]][[1]]$allresults$loglikelihood[-1], type="l")


# g = 3
plot(finalresults$allresults[[3]][[1]][[1]][[1]][[1]]$allresults$loglikelihood[-1], type="l")
plot(finalresults$allresults[[3]][[1]][[1]][[2]][[1]]$allresults$loglikelihood[-1], type="l")
plot(finalresults$allresults[[3]][[1]][[1]][[3]][[1]]$allresults$loglikelihood[-1], type="l")

plot(finalresults$allresults[[3]][[2]][[1]][[1]][[1]]$allresults$loglikelihood[-1], type="l")
plot(finalresults$allresults[[3]][[2]][[1]][[2]][[1]]$allresults$loglikelihood[-1], type="l")
plot(finalresults$allresults[[3]][[2]][[1]][[3]][[1]]$allresults$loglikelihood[-1], type="l")

plot(finalresults$allresults[[3]][[3]][[1]][[1]][[1]]$allresults$loglikelihood[-1], type="l")
plot(finalresults$allresults[[3]][[3]][[1]][[2]][[1]]$allresults$loglikelihood[-1], type="l")
plot(finalresults$allresults[[3]][[3]][[1]][[3]][[1]]$allresults$loglikelihood[-1], type="l")

plot(finalresults$allresults[[3]][[4]][[1]][[1]][[1]]$allresults$loglikelihood[-1], type="l")
plot(finalresults$allresults[[3]][[4]][[1]][[2]][[1]]$allresults$loglikelihood[-1], type="l")
plot(finalresults$allresults[[3]][[4]][[1]][[3]][[1]]$allresults$loglikelihood[-1], type="l")

# g = 4
plot(finalresults$allresults[[4]][[1]][[1]][[1]][[1]]$allresults$loglikelihood[-1], type="l")
plot(finalresults$allresults[[4]][[1]][[1]][[2]][[1]]$allresults$loglikelihood[-1], type="l")
plot(finalresults$allresults[[4]][[1]][[1]][[3]][[1]]$allresults$loglikelihood[-1], type="l")

plot(finalresults$allresults[[4]][[2]][[1]][[1]][[1]]$allresults$loglikelihood[-1], type="l")
plot(finalresults$allresults[[4]][[2]][[1]][[2]][[1]]$allresults$loglikelihood[-1], type="l")
plot(finalresults$allresults[[4]][[2]][[1]][[3]][[1]]$allresults$loglikelihood[-1], type="l")

plot(finalresults$allresults[[4]][[3]][[1]][[1]][[1]]$allresults$loglikelihood[-1], type="l")
plot(finalresults$allresults[[4]][[3]][[1]][[2]][[1]]$allresults$loglikelihood[-1], type="l")
plot(finalresults$allresults[[4]][[3]][[1]][[3]][[1]]$allresults$loglikelihood[-1], type="l")

plot(finalresults$allresults[[4]][[4]][[1]][[1]][[1]]$allresults$loglikelihood[-1], type="l")
plot(finalresults$allresults[[4]][[4]][[1]][[2]][[1]]$allresults$loglikelihood[-1], type="l")
plot(finalresults$allresults[[4]][[4]][[1]][[3]][[1]]$allresults$loglikelihood[-1], type="l")


# g = 5
plot(finalresults$allresults[[5]][[1]][[1]][[1]][[1]]$allresults$loglikelihood[-1], type="l")
plot(finalresults$allresults[[5]][[1]][[1]][[2]][[1]]$allresults$loglikelihood[-1], type="l")
plot(finalresults$allresults[[5]][[1]][[1]][[3]][[1]]$allresults$loglikelihood[-1], type="l")

plot(finalresults$allresults[[5]][[2]][[1]][[1]][[1]]$allresults$loglikelihood[-1], type="l")
plot(finalresults$allresults[[5]][[2]][[1]][[2]][[1]]$allresults$loglikelihood[-1], type="l")
plot(finalresults$allresults[[5]][[2]][[1]][[3]][[1]]$allresults$loglikelihood[-1], type="l")

plot(finalresults$allresults[[5]][[3]][[1]][[1]][[1]]$allresults$loglikelihood[-1], type="l")
plot(finalresults$allresults[[5]][[3]][[1]][[2]][[1]]$allresults$loglikelihood[-1], type="l")
plot(finalresults$allresults[[5]][[3]][[1]][[3]][[1]]$allresults$loglikelihood[-1], type="l")

plot(finalresults$allresults[[5]][[4]][[1]][[1]][[1]]$allresults$loglikelihood[-1], type="l")
plot(finalresults$allresults[[5]][[4]][[1]][[2]][[1]]$allresults$loglikelihood[-1], type="l")
plot(finalresults$allresults[[5]][[4]][[1]][[3]][[1]]$allresults$loglikelihood[-1], type="l")

# g = 6
plot(finalresults$allresults[[6]][[1]][[1]][[1]][[1]]$allresults$loglikelihood[-1], type="l")
plot(finalresults$allresults[[6]][[1]][[1]][[2]][[1]]$allresults$loglikelihood[-1], type="l")
plot(finalresults$allresults[[6]][[1]][[1]][[3]][[1]]$allresults$loglikelihood[-1], type="l")

plot(finalresults$allresults[[6]][[2]][[1]][[1]][[1]]$allresults$loglikelihood[-1], type="l")
plot(finalresults$allresults[[6]][[2]][[1]][[2]][[1]]$allresults$loglikelihood[-1], type="l")
plot(finalresults$allresults[[6]][[2]][[1]][[3]][[1]]$allresults$loglikelihood[-1], type="l")

plot(finalresults$allresults[[6]][[3]][[1]][[1]][[1]]$allresults$loglikelihood[-1], type="l")
plot(finalresults$allresults[[6]][[3]][[1]][[2]][[1]]$allresults$loglikelihood[-1], type="l")
plot(finalresults$allresults[[6]][[3]][[1]][[3]][[1]]$allresults$loglikelihood[-1], type="l")

plot(finalresults$allresults[[6]][[4]][[1]][[1]][[1]]$allresults$loglikelihood[-1], type="l")
plot(finalresults$allresults[[6]][[4]][[1]][[2]][[1]]$allresults$loglikelihood[-1], type="l")
plot(finalresults$allresults[[6]][[4]][[1]][[3]][[1]]$allresults$loglikelihood[-1], type="l")

# g = 7
plot(finalresults$allresults[[7]][[1]][[1]][[1]][[1]]$allresults$loglikelihood[-1], type="l")
plot(finalresults$allresults[[7]][[1]][[1]][[2]][[1]]$allresults$loglikelihood[-1], type="l")
plot(finalresults$allresults[[7]][[1]][[1]][[3]][[1]]$allresults$loglikelihood[-1], type="l")

plot(finalresults$allresults[[7]][[2]][[1]][[1]][[1]]$allresults$loglikelihood[-1], type="l")
plot(finalresults$allresults[[7]][[2]][[1]][[2]][[1]]$allresults$loglikelihood[-1], type="l")
plot(finalresults$allresults[[7]][[2]][[1]][[3]][[1]]$allresults$loglikelihood[-1], type="l")

plot(finalresults$allresults[[7]][[3]][[1]][[1]][[1]]$allresults$loglikelihood[-1], type="l")
plot(finalresults$allresults[[7]][[3]][[1]][[2]][[1]]$allresults$loglikelihood[-1], type="l")
plot(finalresults$allresults[[7]][[3]][[1]][[3]][[1]]$allresults$loglikelihood[-1], type="l")

plot(finalresults$allresults[[7]][[4]][[1]][[1]][[1]]$allresults$loglikelihood[-1], type="l")
plot(finalresults$allresults[[7]][[4]][[1]][[2]][[1]]$allresults$loglikelihood[-1], type="l")
plot(finalresults$allresults[[7]][[4]][[1]][[3]][[1]]$allresults$loglikelihood[-1], type="l")

# g = 8
plot(finalresults$allresults[[8]][[1]][[1]][[1]][[1]]$allresults$loglikelihood[-1], type="l")
plot(finalresults$allresults[[8]][[1]][[1]][[2]][[1]]$allresults$loglikelihood[-1], type="l")
plot(finalresults$allresults[[8]][[1]][[1]][[3]][[1]]$allresults$loglikelihood[-1], type="l")

plot(finalresults$allresults[[8]][[2]][[1]][[1]][[1]]$allresults$loglikelihood[-1], type="l")
plot(finalresults$allresults[[8]][[2]][[1]][[2]][[1]]$allresults$loglikelihood[-1], type="l")
plot(finalresults$allresults[[8]][[2]][[1]][[3]][[1]]$allresults$loglikelihood[-1], type="l")

plot(finalresults$allresults[[8]][[3]][[1]][[1]][[1]]$allresults$loglikelihood[-1], type="l")
plot(finalresults$allresults[[8]][[3]][[1]][[2]][[1]]$allresults$loglikelihood[-1], type="l")
plot(finalresults$allresults[[8]][[3]][[1]][[3]][[1]]$allresults$loglikelihood[-1], type="l")

plot(finalresults$allresults[[8]][[4]][[1]][[1]][[1]]$allresults$loglikelihood[-1], type="l")
plot(finalresults$allresults[[8]][[4]][[1]][[2]][[1]]$allresults$loglikelihood[-1], type="l")
plot(finalresults$allresults[[8]][[4]][[1]][[3]][[1]]$allresults$loglikelihood[-1], type="l")

# g = 9
plot(finalresults$allresults[[9]][[1]][[1]][[1]][[1]]$allresults$loglikelihood[-1], type="l")
plot(finalresults$allresults[[9]][[1]][[1]][[2]][[1]]$allresults$loglikelihood[-1], type="l")
plot(finalresults$allresults[[9]][[1]][[1]][[3]][[1]]$allresults$loglikelihood[-1], type="l")

plot(finalresults$allresults[[9]][[2]][[1]][[1]][[1]]$allresults$loglikelihood[-1], type="l")
plot(finalresults$allresults[[9]][[2]][[1]][[2]][[1]]$allresults$loglikelihood[-1], type="l")
plot(finalresults$allresults[[9]][[2]][[1]][[3]][[1]]$allresults$loglikelihood[-1], type="l")

plot(finalresults$allresults[[9]][[3]][[1]][[1]][[1]]$allresults$loglikelihood[-1], type="l")
plot(finalresults$allresults[[9]][[3]][[1]][[2]][[1]]$allresults$loglikelihood[-1], type="l")
plot(finalresults$allresults[[9]][[3]][[1]][[3]][[1]]$allresults$loglikelihood[-1], type="l")

plot(finalresults$allresults[[9]][[4]][[1]][[1]][[1]]$allresults$loglikelihood[-1], type="l")
plot(finalresults$allresults[[9]][[4]][[1]][[2]][[1]]$allresults$loglikelihood[-1], type="l")
plot(finalresults$allresults[[9]][[4]][[1]][[3]][[1]]$allresults$loglikelihood[-1], type="l")

# Checking convergence
library(coda)
# format allruns[[g]][[numbmodel]][[1]][[qvalue]][[1]]

# g = 1
all(heidel.diag(finalresults$allresults[[1]][[1]][[1]][[1]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)
all(heidel.diag(finalresults$allresults[[1]][[1]][[1]][[2]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)
all(heidel.diag(finalresults$allresults[[1]][[1]][[1]][[3]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)

all(heidel.diag(finalresults$allresults[[1]][[2]][[1]][[1]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)
all(heidel.diag(finalresults$allresults[[1]][[2]][[1]][[2]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)
all(heidel.diag(finalresults$allresults[[1]][[2]][[1]][[3]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)

all(heidel.diag(finalresults$allresults[[1]][[3]][[1]][[1]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)
all(heidel.diag(finalresults$allresults[[1]][[3]][[1]][[2]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)
all(heidel.diag(finalresults$allresults[[1]][[3]][[1]][[3]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)

all(heidel.diag(finalresults$allresults[[1]][[4]][[1]][[1]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)
all(heidel.diag(finalresults$allresults[[1]][[4]][[1]][[2]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)
all(heidel.diag(finalresults$allresults[[1]][[4]][[1]][[3]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)


# g = 2
all(heidel.diag(finalresults$allresults[[2]][[1]][[1]][[1]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)
all(heidel.diag(finalresults$allresults[[2]][[1]][[1]][[2]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)
all(heidel.diag(finalresults$allresults[[2]][[1]][[1]][[3]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)

all(heidel.diag(finalresults$allresults[[2]][[2]][[1]][[1]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)
all(heidel.diag(finalresults$allresults[[2]][[2]][[1]][[2]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)
all(heidel.diag(finalresults$allresults[[2]][[2]][[1]][[3]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)

all(heidel.diag(finalresults$allresults[[2]][[3]][[1]][[1]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)
all(heidel.diag(finalresults$allresults[[2]][[3]][[1]][[2]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)
all(heidel.diag(finalresults$allresults[[2]][[3]][[1]][[3]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)

all(heidel.diag(finalresults$allresults[[2]][[4]][[1]][[1]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)
all(heidel.diag(finalresults$allresults[[2]][[4]][[1]][[2]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)
all(heidel.diag(finalresults$allresults[[2]][[4]][[1]][[3]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)

# g = 3
all(heidel.diag(finalresults$allresults[[3]][[1]][[1]][[1]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)
all(heidel.diag(finalresults$allresults[[3]][[1]][[1]][[2]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)
all(heidel.diag(finalresults$allresults[[3]][[1]][[1]][[3]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)

all(heidel.diag(finalresults$allresults[[3]][[2]][[1]][[1]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)
all(heidel.diag(finalresults$allresults[[3]][[2]][[1]][[2]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)
all(heidel.diag(finalresults$allresults[[3]][[2]][[1]][[3]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)

all(heidel.diag(finalresults$allresults[[3]][[3]][[1]][[1]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)
all(heidel.diag(finalresults$allresults[[3]][[3]][[1]][[2]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)
all(heidel.diag(finalresults$allresults[[3]][[3]][[1]][[3]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)

all(heidel.diag(finalresults$allresults[[3]][[4]][[1]][[1]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)
all(heidel.diag(finalresults$allresults[[3]][[4]][[1]][[2]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)
all(heidel.diag(finalresults$allresults[[3]][[4]][[1]][[3]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)

# g = 4
all(heidel.diag(finalresults$allresults[[4]][[1]][[1]][[1]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)
all(heidel.diag(finalresults$allresults[[4]][[1]][[1]][[2]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)
all(heidel.diag(finalresults$allresults[[4]][[1]][[1]][[3]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)

all(heidel.diag(finalresults$allresults[[4]][[2]][[1]][[1]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)
all(heidel.diag(finalresults$allresults[[4]][[2]][[1]][[2]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)
all(heidel.diag(finalresults$allresults[[4]][[2]][[1]][[3]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)

all(heidel.diag(finalresults$allresults[[4]][[3]][[1]][[1]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)
all(heidel.diag(finalresults$allresults[[4]][[3]][[1]][[2]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)
all(heidel.diag(finalresults$allresults[[4]][[3]][[1]][[3]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)

all(heidel.diag(finalresults$allresults[[4]][[4]][[1]][[1]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)
all(heidel.diag(finalresults$allresults[[4]][[4]][[1]][[2]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)
all(heidel.diag(finalresults$allresults[[4]][[4]][[1]][[3]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)


# g = 5
all(heidel.diag(finalresults$allresults[[5]][[1]][[1]][[1]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)
all(heidel.diag(finalresults$allresults[[5]][[1]][[1]][[2]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)
all(heidel.diag(finalresults$allresults[[5]][[1]][[1]][[3]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)

all(heidel.diag(finalresults$allresults[[5]][[2]][[1]][[1]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)
all(heidel.diag(finalresults$allresults[[5]][[2]][[1]][[2]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)
all(heidel.diag(finalresults$allresults[[5]][[2]][[1]][[3]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)

all(heidel.diag(finalresults$allresults[[5]][[3]][[1]][[1]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)
all(heidel.diag(finalresults$allresults[[5]][[3]][[1]][[2]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)
all(heidel.diag(finalresults$allresults[[5]][[3]][[1]][[3]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)

all(heidel.diag(finalresults$allresults[[5]][[4]][[1]][[1]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)
all(heidel.diag(finalresults$allresults[[5]][[4]][[1]][[2]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)
all(heidel.diag(finalresults$allresults[[5]][[4]][[1]][[3]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)

# g = 6
all(heidel.diag(finalresults$allresults[[6]][[1]][[1]][[1]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)
all(heidel.diag(finalresults$allresults[[6]][[1]][[1]][[2]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)
all(heidel.diag(finalresults$allresults[[6]][[1]][[1]][[3]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)

all(heidel.diag(finalresults$allresults[[6]][[2]][[1]][[1]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)
all(heidel.diag(finalresults$allresults[[6]][[2]][[1]][[2]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)
all(heidel.diag(finalresults$allresults[[6]][[2]][[1]][[3]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)

all(heidel.diag(finalresults$allresults[[6]][[3]][[1]][[1]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)
all(heidel.diag(finalresults$allresults[[6]][[3]][[1]][[2]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)
all(heidel.diag(finalresults$allresults[[6]][[3]][[1]][[3]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)

all(heidel.diag(finalresults$allresults[[6]][[4]][[1]][[1]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)
all(heidel.diag(finalresults$allresults[[6]][[4]][[1]][[2]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)
all(heidel.diag(finalresults$allresults[[6]][[4]][[1]][[3]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)

# g = 7
all(heidel.diag(finalresults$allresults[[7]][[1]][[1]][[1]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)
all(heidel.diag(finalresults$allresults[[7]][[1]][[1]][[2]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)
all(heidel.diag(finalresults$allresults[[7]][[1]][[1]][[3]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)

all(heidel.diag(finalresults$allresults[[7]][[2]][[1]][[1]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)
all(heidel.diag(finalresults$allresults[[7]][[2]][[1]][[2]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)
all(heidel.diag(finalresults$allresults[[7]][[2]][[1]][[3]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)

all(heidel.diag(finalresults$allresults[[7]][[3]][[1]][[1]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)
all(heidel.diag(finalresults$allresults[[7]][[3]][[1]][[2]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)
all(heidel.diag(finalresults$allresults[[7]][[3]][[1]][[3]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)

all(heidel.diag(finalresults$allresults[[7]][[4]][[1]][[1]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)
all(heidel.diag(finalresults$allresults[[7]][[4]][[1]][[2]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)
all(heidel.diag(finalresults$allresults[[7]][[4]][[1]][[3]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)

# g = 8
all(heidel.diag(finalresults$allresults[[8]][[1]][[1]][[1]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)
all(heidel.diag(finalresults$allresults[[8]][[1]][[1]][[2]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)
all(heidel.diag(finalresults$allresults[[8]][[1]][[1]][[3]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)

all(heidel.diag(finalresults$allresults[[8]][[2]][[1]][[1]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)
all(heidel.diag(finalresults$allresults[[8]][[2]][[1]][[2]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)
all(heidel.diag(finalresults$allresults[[8]][[2]][[1]][[3]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)

all(heidel.diag(finalresults$allresults[[8]][[3]][[1]][[1]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)
all(heidel.diag(finalresults$allresults[[8]][[3]][[1]][[2]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)
all(heidel.diag(finalresults$allresults[[8]][[3]][[1]][[3]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)

all(heidel.diag(finalresults$allresults[[8]][[4]][[1]][[1]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)
all(heidel.diag(finalresults$allresults[[8]][[4]][[1]][[2]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)
all(heidel.diag(finalresults$allresults[[8]][[4]][[1]][[3]][[1]]$allresults$loglikelihood[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1)


# Model selection
finalresults$BIC.all$BICmodelselected_g
min(finalresults$BIC.all$allBICvalues, na.rm = T) # g=5, q=3, UUU
finalresults$ICL.all$ICLmodelselected_g
min(finalresults$ICL.all$allICLvalues, na.rm = T) # g=5, q=3, UUU


# plotting model selection
colvector=c(rep(1,27), rep(2,27), rep(3,27), rep(4,27))
labelvector=c(rep(c(1:27),4))
plot(finalresults$BIC.all$allBICvalues, 
     xlab= "The model",
     ylab= "BIC value",
     col= colvector+1, pch = colvector, cex = 0.5, lty = "solid", lwd = 5,ylim=c(min(finalresults$BIC.all$allBICvalues)-100,max(finalresults$BIC.all$allBICvalues)+200),  xaxt="n" )

text(finalresults$BIC.all$allBICvalues, labels=labelvector, cex= 0.8, pos=3)
par(xpd=TRUE)
legend(25,86200,c("M = CC","M = CU","M = UC","M = UU"), horiz = T, lty = "solid", pch=c(1:4), col = c(2:5),  bty = "n", cex=0.75, lwd = 2.5)
abline(77828.00,0,lty=2, col=5)
# Note 1 is G=1, q=1, 
# Note 2 is G=2, q=1,
# Note 3 is G=3, q=1

editedres = array(dim=c(10,3,4))
editedres[1:9,1:3,1:4] =finalresults$BIC.all$allBICvalues
editedres[10,1:3,1] = c(81474.77, 80112.72, 78890.82)
editedres[10,1:3,2] = c(81814.32, 80259.27, 79079.83)
editedres[10,1:3,3] = c(82480.68, 79765.68, 78651.19)
editedres[10,1:3,4] = c(81751.39, 79838.72, 78887.44)

colvector=c(rep(1,30), rep(2,30), rep(3,30), rep(4,30))
labelvector=c(rep(c(1:30),4))
plot(editedres, 
     xlab= "The model",
     ylab= "BIC value",
     col= colvector+1, pch = colvector, cex = 0.5, lty = "solid", lwd = 5, ylim=c(min(editedres)-300,max(editedres)+300),  xaxt="n" )
text(editedres, labels=labelvector, cex= 0.8, pos=3)
par(xpd=TRUE)
legend(25,86200,c("M = CC","M = CU","M = UC","M = UU"), horiz = T, lty = "solid", pch=c(1:4), col = c(2:5),  bty = "n", cex=0.75, lwd = 2.5)
abline(77828.00,0,lty=2, col=5)



# plotting loglikelihood
colvector=c(rep(1,24), rep(2,24), rep(3,24), rep(4,24))
#labelvector=c("1,1,1", "1,2,1","1,3,1","2,1,1","2,2,1","2,3,1","3,1,1","3,2,1","3,3,1", "1,1,2", "1,2,2","1,3,2","2,1,2","2,2,2","2,3,2","3,1,2","3,2,2","3,3,2",   "1,1,3", "1,2,3","1,3,3","2,1,3","2,2,3","2,3,3","3,1,3","3,2,3","3,3,3",   "1,1,4", "1,2,4","1,3,4","2,1,4","2,2,4","2,3,4","3,1,4","3,2,4","3,3,4")
labelvector=c(rep(1:24,4))
plot(finalresults$loglikelihood, 
     xlab= "The model",
     ylab= "LogL value",
     col= colvector+1, pch = colvector, cex = 0.5, lty = "solid", lwd = 5,ylim=c(min(finalresults$loglikelihood)-100,max(finalresults$loglikelihood)+100) )

text(finalresults$loglikelihood, labels=labelvector, cex= 0.8, pos=3)
abline(-38413.78,0,lty=2, col=5)


# save labels
mat_labels=cbind(bozzo2016,finalresults$ICL.all$ICLmodelselected_labels)
ordervector_nb = list()
anothervector_nb = list()
for (i in 1:5){
  ordervector_nb[[i]]=which(mat_labels[,7]==i)
  anothervector_nb[[i]]=rep(i,length(which(mat_labels[,7]==i)))
}


tableoflabels = matrix(ncol=5, nrow=723)
colnames(tableoflabels) = c("cluster_1", "cluster_2","cluster_3","cluster_4","cluster_5")
tableoflabels[1:length(unlist(ordervector_nb[[1]])),1] = unlist(names(ordervector_nb[[1]]))
tableoflabels[1:length(unlist(ordervector_nb[[2]])),2] = unlist(names(ordervector_nb[[2]]))
tableoflabels[1:length(unlist(ordervector_nb[[3]])),3] = unlist(names(ordervector_nb[[3]]))
tableoflabels[1:length(unlist(ordervector_nb[[4]])),4] = unlist(names(ordervector_nb[[4]]))
tableoflabels[1:length(unlist(ordervector_nb[[5]])),5] = unlist(names(ordervector_nb[[5]]))
write.csv(tableoflabels, "FAMPLN_25Dec2017_bozzo_all_g=5.RData.csv")

# heatmap
vec<-unlist(ordervector_nb)
colorsvector<-unlist(anothervector_nb)

library(gplots)
library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

heatmap.2(as.matrix(bozzo2016[vec,]),dendrogram="none",trace="none",scale="row",
          Rowv=FALSE,  col = rev(redgreen(75)), RowSideColors=col_vector[colorsvector+1])
par(xpd=TRUE)
legend(-0.34, 0.7,      
       legend = paste0("Cluster ", unique(colorsvector)),
       col = unique(col_vector[colorsvector+1]), 
       lty= 1,             
       lwd = 5,           
       cex=.5,  xpd = TRUE, horiz = FALSE
)
legend('left','groups',paste0("Cluster ", unique(colorsvector)), lty = c(1,2,3),
       col=unique(col_vector[colorsvector]),ncol=3,bty ="n")

# matplot of expression
mat_labels=cbind(bozzo2016,finalresults$ICL.all$ICLmodelselected_labels)
matplot(t(log(bozzo2016+1)), type="l" ,pch=1,col = 1, xlab="sample", ylab="Counts", cex=1, lty=2, xaxt="n", xlim=c(1,ncol(bozzo2016)) )
matplot(t(bozzo2016), type="l" ,pch=1,col = 1, xlab="sample", ylab="Counts", cex=1, lty=2, xaxt="n", xlim=c(1,ncol(bozzo2016)) ) #plot

#### Line plots ####
# arranged by dark vs non-dark
mat_labels=cbind(bozzo2016,finalresults$ICL.all$ICLmodelselected_labels)

par(mfrow=c(2,3))
toplot_1=mat_labels[which(mat_labels[,7]==1),c(1:6)]
toplot1=rbind(log(toplot_1+1), colMeans(log(toplot_1+1)))
toplot1_space=cbind(toplot1[,c(1:3)],rep(NA,723),toplot1[,c(4:6)])
matplot(t(toplot1_space), type="l" ,pch=1,col=c(rep(1,722),7), xlab="Sample", ylab="Expression", cex=1, lty=c(rep(2,722),1),lwd=c(rep(1,722),3), xaxt="n", xlim=c(1,ncol(toplot1_space)), main="Cluster 1")
axis(1,at=c(1:3, 5:7),labels=c("DE", "DI", "DM","NDE", "NDI", "NDM"))
#matplot(t((toplot_1+1)), type="l" ,pch=1,col = 1, xlab="Sample", ylab="Expression", cex=1, lty=2, xaxt="n", xlim=c(1,ncol(toplot_1)))

toplot_1=mat_labels[which(mat_labels[,7]==2),c(1:6)]
toplot1=rbind(log(toplot_1+1), colMeans(log(toplot_1+1)))
toplot1_space=cbind(toplot1[,c(1:3)],rep(NA,60),toplot1[,c(4:6)])
matplot(t(toplot1_space), type="l" ,pch=1,col=c(rep(2,59),7), xlab="Sample", ylab="Expression", cex=1, lty=c(rep(2,59),1),lwd=c(rep(1,59),3), xaxt="n", xlim=c(1,ncol(toplot1_space)), main="Cluster 2")
axis(1,at=c(1:3, 5:7),labels=c("DE", "DI", "DM","NDE", "NDI", "NDM"))
#matplot(t((toplot_2+1)), type="l" ,pch=1,col = 2, xlab="Sample", ylab="Expression", cex=1, lty=2, xaxt="n", xlim=c(1,ncol(toplot_2)))

toplot_1=mat_labels[which(mat_labels[,7]==3),c(1:6)]
toplot1=rbind(log(toplot_1+1), colMeans(log(toplot_1+1)))
toplot1_space=cbind(toplot1[,c(1:3)],rep(NA,181),toplot1[,c(4:6)])
matplot(t(toplot1_space), type="l" ,pch=1,col=c(rep(3,180),7), xlab="Sample", ylab="Expression", cex=1, lty=c(rep(2,180),1),lwd=c(rep(1,180),3), xaxt="n", xlim=c(1,ncol(toplot1_space)), main="Cluster 3")
axis(1,at=c(1:3, 5:7),labels=c("DE", "DI", "DM","NDE", "NDI", "NDM"))
#matplot(t((toplot_3+1)), type="l" ,pch=1,col = 3, xlab="Sample", ylab="Expression", cex=1, lty=2, xaxt="n", xlim=c(1,ncol(toplot_3)))

toplot_1=mat_labels[which(mat_labels[,7]==4),c(1:6)]
toplot1=rbind(log(toplot_1+1), colMeans(log(toplot_1+1)))
toplot1_space=cbind(toplot1[,c(1:3)],rep(NA,112),toplot1[,c(4:6)])
matplot(t(toplot1_space), type="l" ,pch=1,col=c(rep(4,111),7), xlab="Sample", ylab="Expression", cex=1, lty=c(rep(2,111),1),lwd=c(rep(1,111),3), xaxt="n", xlim=c(1,ncol(toplot1_space)), main="Cluster 4")
axis(1,at=c(1:3, 5:7),labels=c("DE", "DI", "DM","NDE", "NDI", "NDM"))
#matplot(t((toplot_4+1)), type="l" ,pch=1,col = 4, xlab="Sample", ylab="Expression", cex=1, lty=2, xaxt="n", xlim=c(1,ncol(toplot_4)))

toplot_1=mat_labels[which(mat_labels[,7]==5),c(1:6)]
toplot1=rbind(log(toplot_1+1), colMeans(log(toplot_1+1)))
toplot1_space=cbind(toplot1[,c(1:3)],rep(NA,265),toplot1[,c(4:6)])
matplot(t(toplot1_space), type="l" ,pch=1,col=c(rep(6,264),7), xlab="Sample", ylab="Expression", cex=1, lty=c(rep(2,264),1),lwd=c(rep(1,264),3), xaxt="n", xlim=c(1,ncol(toplot1_space)), main="Cluster 5")
axis(1,at=c(1:3, 5:7),labels=c("DE", "DI", "DM","NDE", "NDI", "NDM"))
#matplot(t((toplot_5+1)), type="l" ,pch=1,col = 6, xlab="Sample", ylab="Expression", cex=1, lty=2, xaxt="n", xlim=c(1,ncol(toplot_5)))


#par(xpd=NA)
axis(1,at=c(1:6),labels=c("DE", "DI", "DM", "NDE", "NDI", "NDM"))
par(xpd=TRUE)
legend(1.3,80500,c("BIC","ICL","AIC","AIC3"), horiz = T, pch = c(1), lty = c(2), col=1:4,  bty = "n", cex=1)
legend(1.3,80500,c("BIC","ICL","AIC","AIC3"), horiz = T, pch = c(1), lty = c(2), col=1:4,  bty = "y")


# arranged by growth stage
mat_labels2=cbind(bozzo2016[,c(1,4,2,5,3,6)],finalresults$ICL.all$ICLmodelselected_labels)

par(mfrow=c(2,3))
toplot_2=mat_labels2[which(mat_labels2[,7]==1),c(1:6)]
toplot_2=rbind(log(toplot_2+1), colMeans(log(toplot_2+1)))
toplot_2_space=cbind(toplot_2[,c(1:2)],rep(NA,723),toplot_2[,c(3:4)],rep(NA,723),toplot_2[,c(5:6)])
matplot(t(toplot_2_space), type="l" ,pch=1,col=c(rep(1,722),7), xlab="Sample", ylab="Expression", cex=1, lty=c(rep(2,722),1),lwd=c(rep(1,722),3), xaxt="n", xlim=c(1,ncol(toplot_2_space)), main="Cluster 1")
axis(1,at=c(1:2, 4:5, 7:8),labels=c("DE", "NDE", "DI", "NDI", "DM", "NDM"))


toplot_2=mat_labels[which(mat_labels2[,7]==2),c(1:6)]
toplot_2=rbind(log(toplot_2+1), colMeans(log(toplot_2+1)))
toplot_2_space=cbind(toplot_2[,c(1:2)],rep(NA,60),toplot_2[,c(3:4)],rep(NA,60),toplot_2[,c(5:6)])
matplot(t(toplot_2_space), type="l" ,pch=1,col=c(rep(2,59),7), xlab="Sample", ylab="Expression", cex=1, lty=c(rep(2,59),1),lwd=c(rep(1,59),3), xaxt="n", xlim=c(1,ncol(toplot_2_space)), main="Cluster 2")
axis(1,at=c(1:2, 4:5, 7:8),labels=c("DE", "NDE", "DI", "NDI", "DM", "NDM"))

toplot_2=mat_labels[which(mat_labels2[,7]==3),c(1:6)]
toplot_2=rbind(log(toplot_2+1), colMeans(log(toplot_2+1)))
toplot_2_space=cbind(toplot_2[,c(1:2)],rep(NA,181),toplot_2[,c(3:4)],rep(NA,181),toplot_2[,c(5:6)])
matplot(t(toplot_2_space), type="l" ,pch=1,col=c(rep(3,180),7), xlab="Sample", ylab="Expression", cex=1, lty=c(rep(2,180),1),lwd=c(rep(1,180),3), xaxt="n", xlim=c(1,ncol(toplot_2_space)), main="Cluster 3")
axis(1,at=c(1:2, 4:5, 7:8),labels=c("DE", "NDE", "DI", "NDI", "DM", "NDM"))

toplot_2=mat_labels[which(mat_labels2[,7]==4),c(1:6)]
toplot_2=rbind(log(toplot_2+1), colMeans(log(toplot_2+1)))
toplot_2_space=cbind(toplot_2[,c(1:2)],rep(NA,112),toplot_2[,c(3:4)],rep(NA,112),toplot_2[,c(5:6)])
matplot(t(toplot_2_space), type="l" ,pch=1,col=c(rep(4,111),7), xlab="Sample", ylab="Expression", cex=1, lty=c(rep(2,111),1),lwd=c(rep(1,111),3), xaxt="n", xlim=c(1,ncol(toplot_2_space)), main="Cluster 4")
axis(1,at=c(1:2, 4:5, 7:8),labels=c("DE", "NDE", "DI", "NDI", "DM", "NDM"))

toplot_2=mat_labels[which(mat_labels2[,7]==5),c(1:6)]
toplot_2=rbind(log(toplot_2+1), colMeans(log(toplot_2+1)))
toplot_2_space=cbind(toplot_2[,c(1:2)],rep(NA,265),toplot_2[,c(3:4)],rep(NA,265),toplot_2[,c(5:6)])
matplot(t(toplot_2_space), type="l" ,pch=1,col=c(rep(6,264),7), xlab="Sample", ylab="Expression", cex=1, lty=c(rep(2,264),1),lwd=c(rep(1,264),3), xaxt="n", xlim=c(1,ncol(toplot_2_space)), main="Cluster 5")
axis(1,at=c(1:2, 4:5, 7:8),labels=c("DE", "NDE", "DI", "NDI", "DM", "NDM"))

