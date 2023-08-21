
## functions for updating lambda
fun_lambda_g<-function(g){Sk[,,g]%*%t(beta.var[[g]])%*%solve(bigTheta[[g]])}
fun_lambda_cuu_1<-function(g){ng[g]*solve(psi[[g]])%*%Sk[,,g]%*%t(beta.var[[g]])}
fun_lambda_cuu_2<-function(g){ng[g]*diag(solve(psi[[g]]))[k]*bigTheta[[g]]}
fun_lambda_cuc_1<-function(g){ng[g]/(diag(psi[[g]])[1])*(Sk[,,g]%*%t(beta.var[[g]]))}
fun_lambda_cuc_2<-function(g){ng[g]/(diag(psi[[g]])[1])*bigTheta[[g]]}
fun_lambda_ccu_ccc<-function(){Sg_av%*%t(beta.var[[1]])%*%solve(bigTheta_av)}

## functions for updating psi
fun_psi_uuu<-function(g){diag(Sk[,,g]-lambdanew[[g]]%*%beta.var[[g]]%*%Sk[,,g]+Reduce("+",zS[[g]])/sum(z[,g]))*diag(d)}
fun_psi_uuc<-function(g){mean(diag(Sk[,,g]-lambdanew[[g]]%*%beta.var[[g]]%*%Sk[,,g]+Reduce("+",zS[[g]])/ng[g]))}
fun_psi_ucu<-function(g){ng[g]/N*diag(Sk[,,g]-lambdanew[[g]]%*%beta.var[[g]]%*%Sk[,,g])+diag(Reduce("+",zS[[g]])/N)}
fun_psi_ucc<-function(g){1/d*ng[g]/N*sum(diag(Sk[,,g]-lambdanew[[g]]%*%beta.var[[g]]%*%Sk[,,g]))+sum(diag(Reduce("+",zS[[g]])))/(N*d)}
fun_psi_cuu<-function(g){diag(Sk[,,g]-2*lambdanew[[1]]%*%beta.var[[g]]%*%Sk[,,g]+lambdanew[[1]]%*%bigTheta[[g]]%*%t(lambdanew[[1]])+Reduce("+",zS[[g]])/sum(z[,g]))*diag(d)}
fun_psi_cuc<-function(g){mean(diag(Sk[,,g]-2*lambdanew[[1]]%*%beta.var[[g]]%*%Sk[,,g]+lambdanew[[1]]%*%bigTheta[[g]]%*%t(lambdanew[[1]])+Reduce("+",zS[[g]])/ng[g]))}
fun_psi_ccu<-function(){diag(Sg_av-lambdanew[[1]]%*%beta.var[[1]]%*%Sg_av)+diag(Reduce("+",zS[[1]])/N)}
fun_psi_ccc<-function(){mean(diag(Sg_av-lambdanew[[1]]%*%beta.var[[1]]%*%Sg_av))+sum(diag(Reduce("+",zS[[1]])))/(N*d)}

fun_sgav<-function(g){ng[g]/N*Sk[,,g]}



## model updates
modelUpdates<-function(modelName,
                        zS,
                        ng,
                        z,
                        lambda,
                        isigma,
                        clustersize,
                        pmax.var,
                        Sk,
                        psi){
  beta.var <- list()
  bigTheta <- list()
  sigma.var<-list()

  Sk<<-Sk
  zS<<-zS
  ng<<-ng
  z<<-z
  psi<<-psi
  lambda<<-lambda



  ##### MODEL UUU #####
  if (substr(modelName,1,3)=="UUU"){
    #Updating Lambda and beta.var
    for (g in 1:clustersize){
      beta.var[[g]]<-t(lambda[[g]])%*%isigma[[g]]
      bigTheta[[g]]<-diag(pmax.var)-beta.var[[g]]%*%lambda[[g]]+beta.var[[g]]%*%Sk[,,g]%*%t(beta.var[[g]])
    }
    beta.var<<-beta.var
    bigTheta<<-bigTheta
    for(g in 1:clustersize){
      lambdanew[[g]]<-fun_lambda_g(g)
    }
    lambdanew<<-lambdanew
    for(g in 1:clustersize){
      psinew[[g]]<-fun_psi_uuu(g)
    }
    for (g in 1:clustersize){
      #Complete sigma.var
      sigma.var[[g]]<-(lambdanew[[g]]%*%t(lambdanew[[g]])+psinew[[g]])
      lambda[[g]]<-lambdanew[[g]]
      psi[[g]]<-psinew[[g]]
    }
    #sigma.var<<-sigma.var
    par<-clustersize*(d*pmax.var-0.5*pmax.var*(pmax.var-1))+clustersize*d
  }

  ##### MODEL UUC #####
  if (substr(modelName,1,3)=="UUC"){
    #Updating Lambda and beta.var
    for (g in 1:clustersize){
      beta.var[[g]]<-t(lambda[[g]])%*%isigma[[g]]
      bigTheta[[g]]<-diag(pmax.var)-beta.var[[g]]%*%lambda[[g]]+beta.var[[g]]%*%Sk[,,g]%*%t(beta.var[[g]])
    }
    beta.var<<-beta.var
    bigTheta<<-bigTheta
    for (g in 1:clustersize){
      lambdanew[[g]]<-fun_lambda_g(g)
    }
    lambdanew<<-lambdanew
    for (g in 1:clustersize){
      psinew[[g]]<-fun_psi_uuc(g)*diag(d)
    }
    psi<-psinew ###This line was missing. So, psi was never updated
    for (g in 1:clustersize){
      #Complete sigma.var
      sigma.var[[g]]<-(lambdanew[[g]]%*%t(lambdanew[[g]])+psinew[[g]])
      lambda[[g]]<-lambdanew[[g]]

    }
    par<-clustersize*(d*pmax.var-0.5*pmax.var*(pmax.var-1))+clustersize
  }

  ##### MODEL UCU #####
  if (substr(modelName,1,3)=="UCU"){
    #Updating Lambda and beta
    for (g in 1:clustersize){
      beta.var[[g]]<-t(lambda[[g]])%*%isigma[[g]]
      bigTheta[[g]]<-diag(pmax.var)-beta.var[[g]]%*%lambda[[g]]+beta.var[[g]]%*%Sk[,,g]%*%t(beta.var[[g]])
    }
    beta.var<<-beta.var
    bigTheta<<-bigTheta
    for(g in 1:clustersize){
      lambdanew[[g]]<-fun_lambda_g(g)
    }
    lambdanew<<-lambdanew
    psinew[[1]]<-rowSums(sapply(1:clustersize,fun_psi_ucu))*diag(d)
    for (g in 1:clustersize){
      #Complete sigma.var
      sigma.var[[g]]<-(lambdanew[[g]]%*%t(lambdanew[[g]])+psinew[[1]])
      lambda[[g]]<-lambdanew[[g]]
      psi[[g]]<-psinew[[1]]
    }
    par<-clustersize*(d*pmax.var-0.5*pmax.var*(pmax.var-1))+d
  }

  ##### MODEL UCC #####
  if (substr(modelName,1,3)=="UCC"){
    #Updating Lambda and beta
    for (g in 1:clustersize){
      beta.var[[g]]<-t(lambda[[g]])%*%isigma[[g]]
      bigTheta[[g]]<-diag(pmax.var)-beta.var[[g]]%*%lambda[[g]]+beta.var[[g]]%*%Sk[,,g]%*%t(beta.var[[g]])
    }
    beta.var<<-beta.var
    bigTheta<<-bigTheta
    for(g in 1:clustersize){
      lambdanew[[g]]<-fun_lambda_g(g)
    }
    lambdanew<<-lambdanew
    psinew[[1]]<-sum(sapply(1:clustersize,fun_psi_ucc))*diag(d)
    for (g in 1:clustersize){
      #Complete sigma.var
      sigma.var[[g]]<-(lambdanew[[g]]%*%t(lambdanew[[g]])+psinew[[1]])
      lambda[[g]]<-lambdanew[[g]]
      psi[[g]]<-psinew[[1]]
    }
    par<-clustersize*(d*pmax.var-0.5*pmax.var*(pmax.var-1))+1
  }

  ##### MODEL CUU #####
  if (substr(modelName,1,3)=="CUU"){

    #Updating Lambda and beta.var
    for (g in 1:clustersize){
      beta.var[[g]]<-t(lambda[[1]])%*%isigma[[g]]
      bigTheta[[g]]<-diag(pmax.var)-beta.var[[g]]%*%lambda[[1]]+beta.var[[g]]%*%Sk[,,g]%*%t(beta.var[[g]])
    }
    beta.var<<-beta.var
    bigTheta<<-bigTheta
    for_lam_1<-matrix(rowSums(sapply(1:clustersize,fun_lambda_cuu_1)),d,pmax.var)
    for_lam_2<-list()
    for (k in 1:d){
      k<<-k
      if (pmax.var>1) for_lam_2[[k]]<-solve(matrix(rowSums(sapply(1:clustersize,fun_lambda_cuu_2)),pmax.var,pmax.var)) else for_lam_2[[k]]<-solve(matrix(sum(sapply(1:clustersize,fun_lambda_cuu_2)),pmax.var,pmax.var))
    }
    lambdanew[[1]]<-matrix(NA,d,pmax.var)
    for (k in 1:d){
      lambdanew[[1]][k,]<-for_lam_1[k,]%*%for_lam_2[[k]]
    }
    lambdanew<<-lambdanew
    for (g in 1:clustersize){
      psinew[[g]]<-fun_psi_cuu(g)*diag(d)
      #Complete sigma.var
      sigma.var[[g]]<-(lambdanew[[1]]%*%t(lambdanew[[1]])+psinew[[g]])
      lambda[[g]]<-lambdanew[[1]]
    }
    psi<-psinew ###This line was missing. So, psi was never updated
    par<-(d*pmax.var-0.5*pmax.var*(pmax.var-1))+clustersize*d
  }

  ##### MODEL CUC #####
  if (substr(modelName,1,3)=="CUC"){
    #Updating Lambda and beta.var
    for (g in 1:clustersize){
      beta.var[[g]]<-t(lambda[[1]])%*%isigma[[g]]
      bigTheta[[g]]<-diag(pmax.var)-beta.var[[g]]%*%lambda[[1]]+beta.var[[g]]%*%Sk[,,g]%*%t(beta.var[[g]])
    }
    beta.var<<-beta.var
    bigTheta<<-bigTheta
    if (pmax.var>1){
      lambdanew[[1]]<-matrix(rowSums(sapply(1:clustersize,fun_lambda_cuc_1)),d,pmax.var)%*%solve(matrix(rowSums(sapply(1:clustersize,fun_lambda_cuc_2)),pmax.var,pmax.var)	)
    } else lambdanew[[1]]<-matrix(rowSums(sapply(1:clustersize,fun_lambda_cuc_1)),d,pmax.var)%*%solve(matrix(sum(sapply(1:clustersize,fun_lambda_cuc_2)),pmax.var,pmax.var)	)
    lambda[[1]]<-lambdanew[[1]]
    lambdanew<<-lambdanew

    ################################
    for (g in 1:clustersize){
      psinew[[g]]<-fun_psi_cuc(g)*diag(d)
      #Complete sigma.var
      sigma.var[[g]]<-(lambdanew[[1]]%*%t(lambdanew[[1]])+psinew[[g]])
      lambda[[g]]<-lambdanew[[1]]
      psi[[g]]<-psinew[[g]]
    }
    par<-(d*pmax.var-0.5*pmax.var*(pmax.var-1))+clustersize
  }

  ##### MODEL CCU #####
  if (substr(modelName,1,3)=="CCU"){
    #Updating Lambda and beta.var
    beta.var[[1]]<-t(lambda[[1]])%*%isigma[[1]]

    Sg_av<-matrix(rowSums(sapply(1:clustersize,fun_sgav)),d,d)
    bigTheta[[1]]<-diag(pmax.var)-beta.var[[1]]%*%lambda[[1]]+beta.var[[1]]%*%Sg_av%*%t(beta.var[[1]])
    bigTheta_av<-bigTheta[[1]]


    beta.var<<-beta.var
    bigTheta_av<<-bigTheta_av
    Sg_av<<-Sg_av
    ####Lambda for CCU, CCC ############
    lambdanew[[1]]<-fun_lambda_ccu_ccc()
    lambdanew<<-lambdanew
    psinew[[1]]<-fun_psi_ccu()*diag(d)
    bigtheta_old<-bigTheta[[1]]
    #Complete sigma.var
    for (g in 1:clustersize){
      sigma.var[[g]]<-(lambdanew[[1]]%*%t(lambdanew[[1]])+psinew[[1]])
      lambda[[g]]<-lambdanew[[1]]
      psi[[g]]<-psinew[[1]]
      bigTheta[[g]]<-bigtheta_old
    }
    bigTheta<<-bigTheta
    par<-(d*pmax.var-0.5*pmax.var*(pmax.var-1))+d
  }

  ##### MODEL CCC #####
  if (substr(modelName,1,3) == "CCC") {

    #Updating Lambda and beta.var
    beta.var[[1]] <- t(lambda[[1]]) %*% isigma[[1]]

    Sg_av <- matrix(rowSums(sapply(1:clustersize, fun_sgav)), d, d)
    bigTheta[[1]] <- diag(pmax.var) - beta.var[[1]] %*% lambda[[1]] +
                     beta.var[[1]] %*% Sg_av %*% t(beta.var[[1]])
    bigTheta_av <- bigTheta[[1]]

    beta.var <<- beta.var
    bigTheta_av <<- bigTheta_av
    Sg_av <<- Sg_av
    # Lambda for CCU, CCC ############
    lambdanew[[1]] <- fun_lambda_ccu_ccc()
    lambdanew <<- lambdanew
    psinew[[1]] <- fun_psi_ccc()*diag(d)
    bigtheta_old <- bigTheta[[1]]
    # Complete sigma.var
    for (g in 1:clustersize) {
      sigma.var[[g]]<-(lambdanew[[1]]%*%t(lambdanew[[1]])+psinew[[1]])
      lambda[[g]]<-lambdanew[[1]]
      psi[[g]]<-psinew[[1]]
      bigTheta[[g]]<-bigtheta_old
    }
    bigTheta<<-bigTheta
    par<-(d*pmax.var-0.5*pmax.var*(pmax.var-1))+1
  }
  return(list(sigma.var = sigma.var, psi = psi, lambda = lambda, par = par,isigma=isigma))
}
