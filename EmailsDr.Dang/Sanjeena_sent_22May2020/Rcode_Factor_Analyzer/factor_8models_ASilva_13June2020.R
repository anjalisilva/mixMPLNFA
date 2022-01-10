rm(list=ls())
library(edgeR)
load("sim_data_ASilva.Rdata")


parallel_FA <- function(run = 1, 
                        G = 2, 
                        pmax = 2, 
                        modelname = "CUU") {

  ####Parameter updates #
  set.seed(run)
  Y <- dat[[run]][[2]]
  true <- dat[[run]][[4]]
  true_par <- dat[[run]][[1]]
  
  # true_par
  
  
  N <- nrow(Y)
  d <- ncol(Y)
  
  lib_mat <- rep(1, d)
  # lib_mat <- log(as.vector(edgeR::calcNormFactors(as.matrix(Y), method = "TMM")))
  
  
  #### Initialization #
  mu <- psi <- lambda <- sigma <- isigma <- list()
  m <- S <- P <- Q <- list()
  n_g <- vector()
  diff_psi <- diff_lam <- list()
  
  #### Other intermediate items initialized #
  Sk <- array(0, c(d, d, G))
  start <- GX <- dGX <- z_S <- list()

  
  
  k_means <- stats::kmeans(log(Y + 1/6), centers = G, nstart = 100)$cluster  
  z <- mclust::unmap(k_means)
  pi_g <- colSums(z) / N

  
  for (g in 1:G) { # initalization
    obs <- which(z[, g] == 1)
    n_g[g] <- length(obs)
    mu[[g]] <- colMeans(log(Y[obs, ] + 1/6)) 
    sigma[[g]] <- var(log(Y[obs, ] + 1/6))
    isigma[[g]] <- solve(sigma[[g]])
  }

    
  # Initializing loading matrix (Lambda) = d * p matrix
  if (substr(modelname, 1, 1) == "U") { 
    for (g in 1:G) {
      lambda[[g]] <- matrix(NA, ncol = pmax, nrow = d)
      temp <- eigen(sigma[[g]])
      for (model in 1:pmax) {
        lambda[[g]][, model] <- temp$vectors[, model] * sqrt(temp$values[model]) 
      }
    }
  } else if(substr(modelname, 1, 1) == "C") {
      lambda <- matrix(ncol = pmax, nrow = d)
      temp <- eigen(sigma[[1]])
      for (model in 1:pmax) {
        lambda[, model] <- temp$vectors[, model] * sqrt(temp$values[model])
      } 
      lambda <- rep(list(lambda), times = 1) 
    }

    
  
  
  # Initializing psi error variance = d * d matrix
  if (substr(modelname, 1, 1) == "U") { # lambda unconstrained
   if (substr(modelname, 2, 2) == "U") { # psi unconstrained
    for (g in 1:G) {
      psi[[g]] <- diag(sigma[[g]]-lambda[[g]]%*%t(lambda[[g]]))*diag(d)
    }
    if (substr(modelname, 3, 3) == "C") {
      cat("\n UUC")
      for (g in 1:G) {
        diag(psi[[g]]) <- 1
      }
    } else {
        cat("\n UUU")
        for (g in 1:G) {
         diag(psi[[g]]) <- g
        }
      }
   } else if(substr(modelname, 2, 2) == "C") { # for initializing psi error variance
     cat("\n UCC or UCU")
     psi <- list(diag(d))
     # if substr(modelnames, 3, 3) == "C", no need to do anything as
     # diagonal elements will be the same anyway
     
     # if substr(modelnames, 3, 3) == "U", no need to do anything as
     # only one psi is present
   }
  } else if (substr(modelname, 1, 1) == "C") { 
    if (substr(modelname, 2, 2) == "U") {
      for (g in 1:G) {
        psi[[g]] <- diag(sigma[[g]] - lambda[[1]] %*% t(lambda[[1]]))*diag(d)
      }
      if (substr(modelname, 3, 3) == "C") {
        cat("\n CUC")
        for (g in 1:G) {
          diag(psi[[g]]) <- 1
        }
      } else {
        cat("\n CUU")
      }
    } else if(substr(modelname, 2, 2) == "C") { # for initializing psi error variance
      cat("\n CCU or CCC")
      psi <- list(diag(d))
      # if substr(modelnames, 3, 3) == "C", no need to do anything as
      # diagonal elements will be the same anyway
      
      # if substr(modelnames, 3, 3) == "U", no need to do anything as
      # only one psi is present
    }
  }
 
  
  
  
  # Initializing sigma = d * d matrix
  for (g in 1:G) {
    
    if (substr(modelname, 1, 1) == "U") {    
      if (substr(modelname, 2, 2) == "U") {
        cat("\nUU_")
        sigma[[g]] <- (lambda[[g]] %*% t(lambda[[g]])) + psi[[g]]
        isigma[[g]] <- solve(sigma[[g]])
      } else if (substr(modelname, 2, 2) == "C") {
        cat("\nUC_")
        sigma[[g]] <- (lambda[[g]] %*% t(lambda[[g]])) + psi[[1]]
        isigma[[g]] <- solve(sigma[[g]])
      }
    } else if (substr(modelname, 1, 1) == "C") {
      if (substr(modelname, 2, 2) == "U") {
        cat("\nCU_")
        sigma[[g]] <- (lambda[[1]] %*% t(lambda[[1]])) + psi[[g]]
        isigma[[g]] <- solve(sigma[[g]])
      } else if (substr(modelname, 2, 2) == "C") {
        cat("\nCC_")
        sigma[[g]] <- (lambda[[1]] %*% t(lambda[[1]])) + psi[[1]]
        isigma[[g]] <- solve(sigma[[g]])
      }
    }
    
    
    start[[g]] <- log(Y + 1/6) ###Starting value for M
    m[[g]] <- log(Y + 1/6)
    S[[g]] <- list()
    for (i in 1:N) {
      S[[g]][[i]] <- diag(d) * 0.000000001
    }
  }
  
  
  
  checks <- 0
  it <- 1
  aloglik <- loglik <- NULL
  aloglik[c(1, 2, 3)] <- 0
  it_max <- 1000
  
  # start clustering 
  while (checks == 0) {
  cat("\n Iteration:", it)
  
    for (g in 1:G) {
      GX[[g]] <- list()
      dGX[[g]] <- list()
      z_S[[g]] <- list()
      for (i in 1:N) {
        dGX[[g]][[i]] <- diag(exp(log(lib_mat) + 
                              start[[g]][i, ]) + 
                              0.5 * diag(S[[g]][[i]]), d) +
                              isigma[[g]]
        S[[g]][[i]] <- solve(dGX[[g]][[i]])
        z_S[[g]][[i]] <- z[i, g] * S[[g]][[i]]
        GX[[g]][[i]] <- Y[i, ] - 
                        exp(start[[g]][i, ] + log(lib_mat) + 
                        0.5 * diag(S[[g]][[i]])) - (isigma[[g]]) %*% 
                        (start[[g]][i, ] - mu[[g]])
        m[[g]][i, ] <- start[[g]][i, ] + S[[g]][[i]] %*% GX[[g]][[i]]
      }
      start[[g]] <- m[[g]]
    }
    cat(" Done line 187")
      
    for (g in 1:G) {  
      ####Updating mu
      mu[[g]] <- colSums(z[, g] * m[[g]]) / sum(z[, g])
      
      ####Updating Sample covariance
      mu_mat <- matrix(rep(mu[[g]], N), nrow = N, byrow = TRUE)
      res <- m[[g]] - mu_mat
      temp <- cov.wt(res, wt = z[, g], center = FALSE, method = "ML")
      Sk[,,g] <- temp$cov
    
      ### Updating Beta
      if (substr(modelname, 1, 1) == "U") {
       beta <- t(lambda[[g]]) %*% isigma[[g]]     # ****
       Q[[g]] <- diag(pmax) - beta %*% lambda[[g]] # ****
      } else if (substr(modelname, 1, 1) == "C") {
        beta <- t(lambda[[1]]) %*% isigma[[g]]     # ****
        Q[[g]] <- diag(pmax) - beta %*% lambda[[1]] # ****
      }
      
      if (substr(modelname, 2, 2) == "U") { # if psi is unconstrained
        temp3 <- solve(psi[[g]]) # ****
      } else if (substr(modelname, 2, 2) == "C") { # if psi is constrained
        temp3 <- solve(psi[[1]]) # ****
      }
     
      if (substr(modelname, 1, 1) == "U") { # if lambda is unconstrained
        P[[g]] <- (m[[g]] - mu_mat) %*% t(Q[[g]] %*% t(lambda[[g]]) %*% temp3)  # ****
      } else if (substr(modelname, 1, 1) == "C") { # if lambda is constrained
        P[[g]] <- (m[[g]] - mu_mat) %*% t(Q[[g]] %*% t(lambda[[1]]) %*% temp3)  # ****  
      }
      
      
      diff <- 1
      min_diff <- 10^{-6}
      diffcount <- 1
      while(diff == 1) {
      cat("\n Internal loop, it: ", diffcount)
      
        if (substr(modelname, 1, 1) == "U") { # if lambda is unconstrained
          lambda_old <- lambda[[g]]  # ****
          theta_old <- diag(pmax) - beta %*% lambda[[g]] + beta %*% Sk[,,g] %*% t(beta)  # ****
        } else if (substr(modelname, 1, 1) == "C") { # if lambda is constrained
          lambda_old <- lambda[[1]]  # ****
          theta_old <- diag(pmax) - beta %*% lambda[[1]] + beta %*% Sk[,,g] %*% t(beta)  # ****
        }
        
        
        if (substr(modelname, 2, 2) == "U") { # if psi is unconstrained
          psi_old <- psi[[g]]  # ****
        } else if (substr(modelname, 2, 2) == "C") { # if psi is constrained
          psi_old <- psi[[1]]  # ****
        }
        
        
        
        
        
        
        ####Updating Lambda
        if (modelname == "UUU" || 
            modelname == "UUC" || 
            modelname == "UCU" ||
            modelname == "UCC") { 
          lambda[[g]] <- Sk[,,g] %*% t(beta) %*% solve(theta_old)  # ****
        } else if (modelname == "CUU" ||
                   modelname == "CUC" || 
                   modelname == "CCU" || 
                   modelname == "CCC") {
         
           # equation 5, calculating r 
          rMatrix <- list()
          if (substr(modelname, 2, 2) == "U") { # if psi is unconstrained
            for (g in 1:G) {
              rMatrix[[g]] <- n_g[g] * solve(psi[[g]]) %*% 
                (stats::cov.wt(res, wt = z[, g], center = FALSE, method = "ML")$cov) %*% t(beta)
            }
            rMatrix <- Reduce('+', rMatrix) 
          } else if (substr(modelname, 2, 2) == "C") { # if psi is constrained
            for (g in 1:G) {
              rMatrix[[g]] <- n_g[g] * solve(psi[[1]]) %*% 
                (stats::cov.wt(res, wt = z[, g], center = FALSE, method = "ML")$cov) %*% t(beta)
            }
            rMatrix <- Reduce('+', rMatrix) 
          }
          
         # equation 5, calculating lambda_i
          for (k in 1:d) {
            internalMat <- list()
            for (g in 1:G) {
              if (substr(modelname, 1, 1) == "U") { # if lambda is unconstrained
                theta_old <- diag(pmax) - beta %*% lambda[[g]] + beta %*% Sk[,,g] %*% t(beta)
              } else if (substr(modelname, 1, 1) == "C") { # if lambda is constrained
                theta_old <- diag(pmax) - beta %*% lambda[[1]] + beta %*% Sk[,,g] %*% t(beta)
              }
              
              if (substr(modelname, 2, 2) == "U") { # if psi is unconstrained
               internalMat[[g]] <- (n_g[g]/diag(psi[[g]])[k]) * theta_old
              } else  if (substr(modelname, 2, 2) == "C") { # if psi is constrained
                internalMat[[g]] <- (n_g[g]/diag(psi[[1]])[k]) * theta_old
              }
              
            }
            internalMat <- Reduce('+', internalMat) 
            lambda[[1]][k, ] <- rMatrix[k, ] %*% solve(internalMat)
            # non-conformable arguments ISSUE **************
          }
        }
         
        
        
        
        ###Updating Beta with new lambda
        if (substr(modelname, 1, 1) == "U") {  # if lambda is unconstrained
          beta <- t(lambda[[g]]) %*% isigma[[g]]  # ****
        } else if (substr(modelname, 1, 1) == "C") {  # if lambda is constrained
          beta <- t(lambda[[1]]) %*% isigma[[g]]  # ****
        }
        
        
      
        
        
        
        ### Updating Psi
        if (modelname == "UUU" || 
            modelname == "CUU") { 
          
          # psi unconstrained
          if (substr(modelname, 1, 1) == "U") {  # if lambda is unconstrained
            psi[[g]] <- diag(Sk[,,g] - lambda[[g]] %*% beta %*% Sk[,,g] + 
                        Reduce("+", z_S[[g]])/sum(z[, g])) * diag(d)  # ****
          } else if (substr(modelname, 1, 1) == "C") {  # if lambda is constrained
            psi[[g]] <- diag(Sk[,,g] - lambda[[1]] %*% beta %*% Sk[,,g] + 
                        Reduce("+", z_S[[g]])/sum(z[, g])) * diag(d)  # ****
          }
            
        } else if (modelname == "UUC" ||
                   modelname == "CUC") {
          # psi unconstrained
          if (substr(modelname, 1, 1) == "U") {  # if lambda is unconstrained
            psi[[g]] <- (1 / d * sum(diag(Sk[,,g] - lambda[[g]] %*% beta %*% Sk[,,g] + 
                        Reduce("+", z_S[[g]])/sum(z[, g])))) * diag(d)  # ****
          } else if (substr(modelname, 1, 1) == "C") {  # if lambda is constrained
            psi[[g]] <- (1 / d * sum(diag(Sk[,,g] - lambda[[1]] %*% beta %*% Sk[,,g] + 
                        Reduce("+", z_S[[g]])/sum(z[, g])))) * diag(d)  # ****
          }
            
        } else {
          # UCU; UCC; CCU; CCC
          
          firstTerm <- secondTerm <- thirdTerm <- list()
          for (g in 1:G) {
           firstTerm[[g]] <- n_g[g]/N * Sk[,,g]
           
           if (substr(modelname, 1, 1) == "U") {  # if lambda is unconstrained
            secondTerm[[g]] <- n_g[g]/N * lambda[[g]] %*% beta %*% Sk[,,g]
           } else if (substr(modelname, 1, 1) == "C") {  # if lambda is constrained
             secondTerm[[g]] <- n_g[g]/N * lambda[[1]] %*% beta %*% Sk[,,g]
           }
           
           thirdTerm[[g]] <- Reduce("+", z_S[[g]])
          }
          
          if (modelname == "UCU" ||
              modelname == "CCU") {
                # psi is constrained
                psi[[1]] <- diag(Reduce("+", firstTerm) - Reduce("+", secondTerm) + 
                            (1 / N * Reduce("+", thirdTerm))) * diag(d) # ****
          }
          if (modelname == "UCC" ||
              modelname == "CCC") {
              # psi is constrained
              psi[[1]] <- (1 / d * sum(diag(Reduce("+", firstTerm) - Reduce("+", secondTerm) + 
                          (1/N * Reduce("+", thirdTerm))))) * diag(d) # ****
              
          }
        }
        
        
        
        if (substr(modelname, 1, 1) == "U") {  # if lambda is unconstrained
          diff_lam[[diffcount]] <- norm(lambda[[g]] - lambda_old, type = "F")  # ****
        } else if (substr(modelname, 1, 1) == "C") {  # if lambda is constrained
          diff_lam[[diffcount]] <- norm(lambda[[1]] - lambda_old, type = "F")  # ****
        }
        
        if (substr(modelname, 2, 2) == "U") { # if psi is unconstrained
          diff_psi[[diffcount]] <- norm(psi[[g]] - psi_old, type = "F")  # ****
        } else if (substr(modelname, 2, 2) == "C") { # if psi is constrained
          diff_psi[[diffcount]] <- norm(psi[[1]] - psi_old, type = "F")  # ****
        }
        if (abs(diff_lam[[diffcount]]) < min_diff & abs(diff_psi[[diffcount]]) < min_diff) {
          diff <- 0
        } 
        # cat("\n Value of diff_psi:", abs(diff_psi[[diffcount]]))
        # cat("\n Value of diff_lam:", abs(diff_lam[[diffcount]]))
        # if(it > 2) {
        # par(mfrow = c(1, 2))
          # plot(unlist(diff_psi), type ="l", main = paste("diff psi; it", it, "; internal it:", diffcount))
          # plot(unlist(diff_lam), type ="l", main = paste("diff lam; it", it, "; internal it:", diffcount))
        # }
        diffcount <- diffcount + 1 

        
        ### Updating Sigma
        if (modelname == "UUU" || modelname == "UUC") {
          sigma[[g]] <- lambda[[g]] %*% t(lambda[[g]]) + psi[[g]]  # ****
        } else if (modelname == "UCU" || modelname == "UCC") {
          sigma[[g]] <- lambda[[g]] %*% t(lambda[[g]]) + psi[[1]]  # ****
        } else if (modelname == "CUU" || modelname == "CUC") {
          sigma[[g]] <- lambda[[1]] %*% t(lambda[[1]]) + psi[[g]]  # ****
        } else if (modelname == "CCU" || modelname == "CCC") {
          sigma[[g]] <- lambda[[1]] %*% t(lambda[[1]]) + psi[[1]]  # ****
        }
          
        isigma[[g]] <- solve(sigma[[g]]) ### Ideally for high dimensional data, 
        # we should be using Woodbury identity so it avoids inverting a p by p matrix
      }
        
      }
    cat(" Done line 395")
    
  
    pi_g <- colSums(z) / N
    lib_mat_full <- matrix(1, ncol = d, nrow = N) ###Matrix containing normaization factor
    ### Some useful functions
    fun_five <- function(x, y = isigma[[g]]) {
      sum(diag(x %*% y))
    }
    
    F <- matrix(NA, ncol = G, nrow = N)
    
    for (g in 1:G) {
      two <- rowSums(exp(m[[g]] + log(lib_mat_full) +
             0.5 * matrix(unlist(lapply(S[[g]], diag)), ncol=d, byrow=TRUE)))
      five <- 0.5 * unlist(lapply(S[[g]], fun_five))
      six <- 0.5 * log(unlist(lapply(S[[g]], det)))
      F[, g] <- pi_g[g] * exp(rowSums(m[[g]] * Y) - two - rowSums(lfactorial(Y)) + 
                rowSums(log(lib_mat_full) * Y) - 0.5 * 
                stats::mahalanobis(m[[g]], center = mu[[g]], cov = sigma[[g]]) - 
                five + six - 0.5 * log(det(sigma[[g]])) - d / 2)
    }  
    
    
  loglik[it] <- sum(log(rowSums(F)))
  par(mfrow = c(1, 1))
  plot(loglik, type = "l", ylab = "logL", main = paste("logL; it", it))
  z <- F / rowSums(F)
  
  
  if (it > 2) {
    #Aitkaine's stopping criterion
    if ((loglik[it - 1] - loglik[it - 2]) == 0) checks <- 1 else {
      a <- (loglik[it] - loglik[it - 1]) / (loglik[it - 1] - loglik[it - 2])
      add_to <- (1 / (1 - a) * (loglik[it] - loglik[it - 1]))
      # }
      aloglik[it] <- loglik[it - 1] + add_to
      if (abs(aloglik[it] - loglik[it - 1]) < 0.01) checks <- 1 else checks <- checks
    }
  }	

  it <- it + 1
  if (it == it_max) checks <- 1   
  
    
  }
  
  
  return(list(pi_g = pi_g,
              mu = mu,
              sigma = sigma,
              lambda = lambda,
              psi = psi,
              z = z,
              loglik = loglik,
              kmeans = k_means,
              true = true))
}

if (.Platform$OS.type == "unix") {
  library("doMC")
  registerDoMC(6) #parallel::detectCores() gives the total number of cores available
}
library("plyr")

total_run <- 5
ptm<-proc.time()
output_all<-list()
output_all <- foreach(run = 1:total_run, .errorhandling = "pass") %dopar% {
  parallel_FA(run)
}
proc.time()-ptm



ARI <- NULL
for (i in 1:total_run) {
  ARI[i] <- mclust::adjustedRandIndex(mclust::map(output_all[[i]]$true), mclust::map(output_all[[i]]$z))
}

pi_all <- matrix(0, ncol = 2, nrow = total_run)
mean_all <- matrix(0, ncol = d * 2, nrow = total_run)
var_all <- matrix(0, ncol = d * d * 2, nrow = total_run)


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

save.image("Routput_ASilva.Rdata")

