
## functions for updating lambda
funLambdag <- function(g, Sk, betaVar, bigTheta) {
  calcValue <- Sk[,,g] %*%
    t(betaVar[[g]]) %*%
    solve(bigTheta[[g]])
  return(calcValue)
}

funLambdaCUU1 <- function(g, ng, psi, Sk, betaVar) {
  calcValue <- ng[g] *
    solve(psi[[g]]) %*% Sk[,,g] %*%
    t(betaVar[[g]])
  return(calcValue)
}

funLambdaCUU2 <- function(g, ng, psi, k, bigTheta) {
  calcValue <- ng[g] *
    diag(solve(psi[[g]]))[k] *
    bigTheta[[g]]
  return(calcValue)
}

funLambdaCUC1 <- function(g, ng, psi, Sk, betaVar) {
  calcValue <- ng[g] / (diag(psi[[g]])[1]) *
    (Sk[,,g] %*% t(betaVar[[g]]))
  return(calcValue)
}


funLambdaCUC2 <- function(g, ng, psi, bigTheta) {
  calcValue <- ng[g] / (diag(psi[[g]])[1]) * bigTheta[[g]]
  return(calcValue)
}


funLambdaCCUnCCC <- function(Sgav, betaVar, bigThetaav) {
  calcValue <- Sgav %*% t(betaVar[[1]]) %*% solve(bigThetaav)
  return(calcValue)
}

## functions for updating psi
funpsiUUU <- function(g, Sk, lambdanew, betaVar, z, dimensionality, zS) {
  calcValue <- diag(Sk[,,g] - lambdanew[[g]] %*%
  betaVar[[g]] %*% Sk[,,g] +
  Reduce("+", zS[[g]]) / sum(z[, g])) * diag(dimensionality)
  return(calcValue)
}

funpsiUUC <- function(g, Sk, lambdanew, betaVar, zS, ng) {
  calcValue <- mean(diag(Sk[,,g] - lambdanew[[g]] %*%
  betaVar[[g]]%*%Sk[,,g] + Reduce("+", zS[[g]]) / ng[g]))
  return(calcValue)
}

funpsiUCU <- function(g, nObservations, ng, Sk, lambdanew, betaVar, zS) {
  calcValue <- ng[g] / nObservations * diag(Sk[,,g] - lambdanew[[g]] %*%
  betaVar[[g]] %*% Sk[,,g]) +
  diag(Reduce("+", zS[[g]]) / nObservations)
  return(calcValue)
}

funPsiUCC <- function(g, ng, dimensionality, nObservations, Sk, lambdanew, zS, betaVar) {
  calcValue <- (1 / dimensionality) * (ng[g] / nObservations) *
    sum(diag(Sk[,,g] - lambdanew[[g]] %*%
    betaVar[[g]] %*% Sk[,,g])) +
    sum(diag(Reduce("+", zS[[g]]))) / (nObservations * dimensionality)
  return(calcValue)
}

funpsiCUU <- function(g, Sk, zS, z, dimensionality, lambdanew, betaVar, bigTheta) {
  calcValue <- diag(Sk[,,g] - 2 * lambdanew[[1]] %*%
    betaVar[[g]] %*% Sk[, , g] +
    lambdanew[[1]] %*% bigTheta[[g]] %*%
    t(lambdanew[[1]]) +
    Reduce("+", zS[[g]]) / sum(z[, g])) * diag(dimensionality)
  return(calcValue)
}

funpsiCUC <- function(g, Sk, lambdanew, betaVar, bigTheta, zS, ng) {
  calcValue <- mean(diag(Sk[,,g] - 2 * lambdanew[[1]] %*%
    betaVar[[g]] %*% Sk[,,g] +
    lambdanew[[1]] %*% bigTheta[[g]] %*%
    t(lambdanew[[1]]) + Reduce("+", zS[[g]]) / ng[g]))
  return(calcValue)
}

funpsiCCU <- function(Sgav, lambdanew, betaVar, zS, nObservations) {
  calcValue <- diag(Sgav - lambdanew[[1]] %*%
    betaVar[[1]] %*% Sgav) +
    diag(Reduce("+", zS[[1]]) / nObservations)
  return(calcValue)
}

funpsiCCC <- function(Sgav, lambdanew, betaVar, zS, nObservations, dimensionality) {
  calcValue <- mean(diag(Sgav - lambdanew[[1]] %*%
    betaVar[[1]] %*% Sgav)) +
    sum(diag(Reduce("+", zS[[1]])))/(nObservations * dimensionality)
  return(calcValue)
}

funSgav <- function(g, ng, nObservations, Sk) {
  calcValue <- ng[g] / nObservations * Sk[,,g]
  return(calcValue)
}

modelUpdates <- function(modelName,
                         zS,
                         ng,
                         z,
                         lambda,
                         isigma,
                         clustersize,
                         pmaxVar,
                         Sk,
                         psi,
                         dimensionality,
                         nObservations) {

  betaVar <- bigTheta <- sigmaVar <-
    lambdanew <- psinew <- list()
  # Sk <<- Sk
  # zS <<- zS
  # ng <<- ng
  # z <<- z
  # psi <<- psi
  # lambda <<- lambda

  if (substr(modelName, 1, 3) == "UUU") {
    for (g in seq_along(1:clustersize)) {
      betaVar[[g]] <- t(lambda[[g]]) %*% isigma[[g]]
      bigTheta[[g]] <- diag(pmaxVar) - betaVar[[g]] %*%
        lambda[[g]] + betaVar[[g]] %*% Sk[, , g] %*%
        t(betaVar[[g]])
    }
    # betaVar <<- betaVar
    bigTheta <<- bigTheta
    for (g in seq_along(1:clustersize)) {
      lambdanew[[g]] <- funLambdag(g = g,
                                   Sk = Sk,
                                   betaVar = betaVar,
                                   bigTheta = bigTheta)
    }

    # lambdanew <<- lambdanew
    for (g in seq_along(1:clustersize)) {
      psinew[[g]] <- funpsiUUU(g = g,
                               Sk = Sk,
                               lambdanew = lambdanew,
                               betaVar = betaVar,
                               z = z,
                               dimensionality = dimensionality,
                               zS = zS)
    }
    for (g in seq_along(1:clustersize)) {
      sigmaVar[[g]] <- (lambdanew[[g]] %*% t(lambdanew[[g]]) +
                          psinew[[g]])
      lambda[[g]] <- lambdanew[[g]]
      psi[[g]] <- psinew[[g]]
    }
    par <- clustersize *
      (dimensionality * pmaxVar - 0.5 * pmaxVar * (pmaxVar - 1)) + clustersize * dimensionality
  }
  if (substr(modelName, 1, 3) == "UUC") {
    for (g in seq_along(1:clustersize)) {
      betaVar[[g]] <- t(lambda[[g]]) %*% isigma[[g]]
      bigTheta[[g]] <- diag(pmaxVar) - betaVar[[g]] %*%
        lambda[[g]] + betaVar[[g]] %*% Sk[, , g] %*%
        t(betaVar[[g]])
    }

    # betaVar <<- betaVar
    # bigTheta <<- bigTheta
    for (g in seq_along(1:clustersize)) {
      lambdanew[[g]] <- funLambdag(g = g,
                                   Sk = Sk,
                                   betaVar = betaVar,
                                   bigTheta = bigTheta)
    }

    # lambdanew <<- lambdanew
    for (g in seq_along(1:clustersize)) {
      psinew[[g]] <- funpsiUUC(g = g,
                               Sk = Sk,
                               lambdanew = lambdanew,
                               betaVar = betaVar,
                               zS = zS,
                               ng = ng) * diag(dimensionality)
    }
    psi <- psinew
    for (g in seq_along(1:clustersize)) {
      sigmaVar[[g]] <- (lambdanew[[g]] %*% t(lambdanew[[g]]) +
                          psinew[[g]])
      lambda[[g]] <- lambdanew[[g]]
    }
    par <- clustersize *
      (dimensionality * pmaxVar - 0.5 * pmaxVar * (pmaxVar - 1)) + clustersize
  }
  if (substr(modelName, 1, 3) == "UCU") {
    for (g in seq_along(1:clustersize)) {
      betaVar[[g]] <- t(lambda[[g]]) %*% isigma[[g]]
      bigTheta[[g]] <- diag(pmaxVar) - betaVar[[g]] %*%
        lambda[[g]] + betaVar[[g]] %*% Sk[, , g] %*%
        t(betaVar[[g]])
    }
    # betaVar <<- betaVar
    # bigTheta <<- bigTheta
    for (g in seq_along(1:clustersize)) {
      lambdanew[[g]] <- funLambdag(g = g,
                                   Sk = Sk,
                                   betaVar = betaVar,
                                   bigTheta = bigTheta)
    }

    # lambdanew <<- lambdanew
    psinew[[1]] <- rowSums(sapply(1:clustersize, funpsiUCU,
                                  nObservations = nObservations, ng = ng, Sk = Sk,
                                  lambdanew = lambdanew,
                                  betaVar = betaVar,
                                  zS = zS)) * diag(dimensionality)
    for (g in seq_along(1:clustersize)) {
      sigmaVar[[g]] <- (lambdanew[[g]] %*% t(lambdanew[[g]]) +
                          psinew[[1]])
      lambda[[g]] <- lambdanew[[g]]
      psi[[g]] <- psinew[[1]]
    }
    par <- clustersize *
      (dimensionality * pmaxVar - 0.5 * pmaxVar * (pmaxVar - 1)) + dimensionality
  }
  if (substr(modelName, 1, 3) == "UCC") {
    for (g in seq_along(1:clustersize)) {
      betaVar[[g]] <- t(lambda[[g]]) %*% isigma[[g]]
      bigTheta[[g]] <- diag(pmaxVar) - betaVar[[g]] %*%
        lambda[[g]] + betaVar[[g]] %*% Sk[, , g] %*%
        t(betaVar[[g]])
    }

    # betaVar <<- betaVar
    # bigTheta <<- bigTheta
    for (g in seq_along(1:clustersize)) {
      lambdanew[[g]] <- funLambdag(g = g,
                                   Sk = Sk,
                                   betaVar = betaVar,
                                   bigTheta = bigTheta)
    }

    # lambdanew <<- lambdanew
    psinew[[1]] <- sum(sapply(1:clustersize, funPsiUCC,
                              ng = ng, dimensionality = dimensionality, nObservations = nObservations, Sk = Sk,
                              lambdanew = lambdanew, zS = zS, betaVar = betaVar),
                       na.rm = TRUE) * diag(dimensionality)

    for (g in seq_along(1:clustersize)) {
      sigmaVar[[g]] <- (lambdanew[[g]] %*% t(lambdanew[[g]]) +
                          psinew[[1]])
      lambda[[g]] <- lambdanew[[g]]
      psi[[g]] <- psinew[[1]]
    }
    par <- clustersize *
      (dimensionality * pmaxVar - 0.5 * pmaxVar * (pmaxVar - 1)) + 1
  }

  if (substr(modelName, 1, 3) == "CUU") {
    for (g in seq_along(1:clustersize)) {
      betaVar[[g]] <- t(lambda[[1]]) %*% isigma[[g]]
      bigTheta[[g]] <- diag(pmaxVar) - betaVar[[g]] %*%
        lambda[[1]] + betaVar[[g]] %*% Sk[, , g] %*%
        t(betaVar[[g]])
    }

    # betaVar <<- betaVar
    # bigTheta <<- bigTheta
    forLam1 <- matrix(data =
                        rowSums(sapply(1:clustersize, funLambdaCUU1,
                                       ng = ng, psi = psi, Sk = Sk, betaVar = betaVar)),
                      nrow = dimensionality, ncol = pmaxVar)
    forLam2 <- list()
    for (k in seq_along(1:dimensionality)) {
      # cat("k = ", k, "\n")
      # k <<- k
      if (pmaxVar > 1) {

        # solve issue anticipation
        # forLam2[[k]] <- solve(matrix(rowSums(sapply(1:clustersize,
        #                                              funLambdaCUU2)), pmaxVar, pmaxVar))
        matLambdaCUU2 <- matrix(data =
                                  rowSums(sapply(1:clustersize,
                                                 funLambdaCUU2, ng = ng, psi = psi,
                                                 k = k, bigTheta = bigTheta)),
                                nrow = pmaxVar,
                                ncol = pmaxVar)
        forLam2[[k]] <- tryCatch(solve(matLambdaCUU2), error = function(err) NA)
        if(all(is.na(forLam2[[k]]))) {
          # cat("\n solve issue matLambdaCUU2 line 593")
          forLam2[[k]] <- diag(ncol(matLambdaCUU2)) # if error with inverse
        }
      } else {
        # solve issue anticipation
        # forLam2[[k]] <- solve(matrix(sum(sapply(1:clustersize,
        #                                              funLambdaCUU2)), pmaxVar, pmaxVar))
        matLambdaCUU2 <- matrix(data = sum(sapply(1:clustersize,
                                                  funLambdaCUU2, ng = ng, psi = psi,
                                                  k = k, bigTheta = bigTheta)),
                                nrow = pmaxVar,
                                ncol = pmaxVar)
        forLam2[[k]] <- tryCatch(solve(matLambdaCUU2), error = function(err) NA)
        if(all(is.na(forLam2[[k]]))) {
          # cat("\n solve issue matLambdaCUU2 in line 607")
          forLam2[[k]] <- diag(ncol(matLambdaCUU2)) # if error with inverse
        }
      }
    }

    lambdanew[[1]] <- matrix(NA, dimensionality, pmaxVar)
    for (k in seq_along(1:dimensionality)) {
      lambdanew[[1]][k, ] <- forLam1[k, ] %*% forLam2[[k]]
    }

    # lambdanew <<- lambdanew
    for (g in seq_along(1:clustersize)) {
      psinew[[g]] <- funpsiCUU(g = g,
                               Sk = Sk,
                               zS = zS,
                               z = z,
                               dimensionality = dimensionality,
                               lambdanew = lambdanew,
                               betaVar = betaVar,
                               bigTheta = bigTheta) * diag(dimensionality)
      sigmaVar[[g]] <- (lambdanew[[1]] %*% t(lambdanew[[1]]) +
                          psinew[[g]])
      lambda[[g]] <- lambdanew[[1]]
    }

    psi <- psinew
    par <- (dimensionality * pmaxVar - 0.5 * pmaxVar * (pmaxVar - 1)) +
      clustersize * dimensionality
  }
  if (substr(modelName, 1, 3) == "CUC") {
    for (g in seq_along(1:clustersize)) {
      betaVar[[g]] <- t(lambda[[1]]) %*% isigma[[g]]
      bigTheta[[g]] <- diag(pmaxVar) - betaVar[[g]] %*%
        lambda[[1]] + betaVar[[g]] %*% Sk[, , g] %*%
        t(betaVar[[g]])
    }

    # betaVar <<- betaVar
    # bigTheta <<- bigTheta
    if (pmaxVar > 1) {
      # solve issue anticipation
      # lambdanew[[1]] <- matrix(rowSums(sapply(1:clustersize, funLambdaCUC1)), dimensionality, pmaxVar) %*%
      #                          solve(matrix(rowSums(sapply(1:clustersize,
      #                                funLambdaCUC2)), pmaxVar, pmaxVar))

      matLambdaCUC2 <- matrix(data = rowSums(sapply(1:clustersize,
                                                    funLambdaCUC2, ng = ng,
                                                    psi = psi, bigTheta = bigTheta)),
                              nrow = pmaxVar,
                              ncol = pmaxVar)

      imatLambdaCUC2 <-  tryCatch(solve(matLambdaCUC2), error = function(err) NA)
      if(all(is.na(imatLambdaCUC2))) {
        imatLambdaCUC2 <- diag(ncol(matLambdaCUC2)) # if error with inverse
      }

      lambdanew[[1]] <- matrix(data = rowSums(sapply(1:clustersize, funLambdaCUC1,
                                                     ng = ng, psi = psi, Sk = Sk, betaVar = betaVar)),
                               nrow = dimensionality,
                               ncol = pmaxVar) %*%
        imatLambdaCUC2


    } else {
      # solve issue anticipation
      # lambdanew[[1]] <- matrix(rowSums(sapply(1:clustersize, funLambdaCUC1)), dimensionality, pmaxVar) %*%
      #                          solve(matrix(sum(sapply(1:clustersize, funLambdaCUC2)),
      #                          pmaxVar, pmaxVar))

      matLambdaCUC2 <- matrix(data = sum(sapply(1:clustersize, funLambdaCUC2,
                                                ng = ng, psi = psi, bigTheta = bigTheta)),
                              nrow = pmaxVar,
                              ncol = pmaxVar)
      imatLambdaCUC2 <-  tryCatch(solve(matLambdaCUC2), error = function(err) NA)
      if(all(is.na(imatLambdaCUC2))) {
        imatLambdaCUC2 <- diag(ncol(matLambdaCUC2)) # if error with inverse
      }

      lambdanew[[1]] <- matrix(data = rowSums(sapply(1:clustersize, funLambdaCUC1,
                                                     ng = ng, psi = psi, Sk = Sk, betaVar = betaVar)),
                               nrow = dimensionality,
                               ncol = pmaxVar) %*%
        imatLambdaCUC2
    }
    lambda[[1]] <- lambdanew[[1]]
    # lambdanew <<- lambdanew
    for (g in 1:clustersize) {
      psinew[[g]] <- funpsiCUC(g = g,
                               Sk = Sk,
                               lambdanew = lambdanew,
                               betaVar = betaVar,
                               bigTheta = bigTheta,
                               zS = zS,
                               ng = ng) * diag(dimensionality)
      sigmaVar[[g]] <- (lambdanew[[1]] %*% t(lambdanew[[1]]) +
                          psinew[[g]])
      lambda[[g]] <- lambdanew[[1]]
      psi[[g]] <- psinew[[g]]
    }
    par <- (dimensionality * pmaxVar - 0.5 * pmaxVar * (pmaxVar - 1)) +
      clustersize
  }
  if (substr(modelName, 1, 3) == "CCU") {
    betaVar[[1]] <- t(lambda[[1]]) %*% isigma[[1]]
    Sgav <- matrix(data = rowSums(sapply(1:clustersize, funSgav,
                                         ng = ng, nObservations = nObservations, Sk = Sk)),
                   nrow = dimensionality,
                   ncol = dimensionality)
    bigTheta[[1]] <- diag(pmaxVar) - betaVar[[1]] %*% lambda[[1]] +
      betaVar[[1]] %*% Sgav %*% t(betaVar[[1]])
    bigThetaav <- bigTheta[[1]]
    # betaVar <<- betaVar
    # bigThetaav <<- bigThetaav
    # Sgav <<- Sgav
    lambdanew[[1]] <- funLambdaCCUnCCC(Sgav = Sgav,
                                       betaVar = betaVar,
                                       bigThetaav = bigThetaav)
    # lambdanew <<- lambdanew
    psinew[[1]] <- funpsiCCU(Sgav = Sgav,
                             lambdanew = lambdanew,
                             betaVar = betaVar,
                             zS = zS,
                             nObservations = nObservations) * diag(dimensionality)
    bigthetaOld <- bigTheta[[1]]
    for (g in seq_along(1:clustersize)) {
      sigmaVar[[g]] <- (lambdanew[[1]] %*% t(lambdanew[[1]]) +
                          psinew[[1]])
      lambda[[g]] <- lambdanew[[1]]
      psi[[g]] <- psinew[[1]]
      bigTheta[[g]] <- bigthetaOld
    }
    # bigTheta <<- bigTheta
    par <- (dimensionality * pmaxVar - 0.5 * pmaxVar * (pmaxVar - 1)) + dimensionality
  }
  if (substr(modelName, 1, 3) == "CCC") {
    betaVar[[1]] <- t(lambda[[1]]) %*% isigma[[1]]
    Sgav <- matrix(data = rowSums(sapply(1:clustersize, funSgav,
                                         ng = ng, nObservations = nObservations, Sk = Sk)),
                   nrow = dimensionality,
                   ncol = dimensionality)
    bigTheta[[1]] <- diag(pmaxVar) - betaVar[[1]] %*% lambda[[1]] +
      betaVar[[1]] %*% Sgav %*% t(betaVar[[1]])
    bigThetaav <- bigTheta[[1]]
    # betaVar <<- betaVar
    # bigThetaav <<- bigThetaav
    # Sgav <<- Sgav
    lambdanew[[1]] <- funLambdaCCUnCCC(Sgav = Sgav,
                                       betaVar = betaVar,
                                       bigThetaav = bigThetaav)
    # lambdanew <<- lambdanew
    psinew[[1]] <- funpsiCCC(Sgav = Sgav,
                             lambdanew = lambdanew,
                             betaVar = betaVar,
                             zS = zS,
                             nObservations = nObservations,
                             dimensionality = dimensionality) * diag(dimensionality)
    bigthetaOld <- bigTheta[[1]]
    for (g in seq_along(1:clustersize)) {
      sigmaVar[[g]] <- (lambdanew[[1]] %*% t(lambdanew[[1]]) +
                          psinew[[1]])
      lambda[[g]] <- lambdanew[[1]]
      psi[[g]] <- psinew[[1]]
      bigTheta[[g]] <- bigthetaOld
    }
    # bigTheta <<- bigTheta
    par <- (dimensionality * pmaxVar - 0.5 * pmaxVar * (pmaxVar - 1)) +
      1
  }

  returnVals <- list(sigmaVar = sigmaVar,
                     psi = psi,
                     lambda = lambda,
                     par = par,
                     isigma = isigma,
                     bigTheta = bigTheta)
  class(returnVals) <- "modelUpdates"
  return(returnVals)
}



# [END]
