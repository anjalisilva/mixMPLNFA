
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
funpsiUUU <- function(g, Sk, lambdanew, betaVar, z, d, zS) {
  calcValue <- diag(Sk[,,g] - lambdanew[[g]] %*%
  betaVar[[g]] %*% Sk[,,g] +
  Reduce("+", zS[[g]]) / sum(z[, g])) * diag(d)
  return(calcValue)
}

funpsiUUC <- function(g, Sk, lambdanew, betaVar, zS, ng) {
  calcValue <- mean(diag(Sk[,,g] - lambdanew[[g]] %*%
  betaVar[[g]]%*%Sk[,,g] + Reduce("+", zS[[g]]) / ng[g]))
  return(calcValue)
}

funpsiUCU <- function(g, N, ng, Sk, lambdanew, betaVar, zS) {
  calcValue <- ng[g] / N * diag(Sk[,,g] - lambdanew[[g]] %*%
  betaVar[[g]] %*% Sk[,,g]) +
  diag(Reduce("+", zS[[g]]) / N)
  return(calcValue)
}

funPsiUCC <- function(g, ng, d, N, Sk, lambdanew, zS, betaVar) {
  calcValue <- (1 / d) * (ng[g] / N) *
    sum(diag(Sk[,,g] - lambdanew[[g]] %*%
    betaVar[[g]] %*% Sk[,,g])) +
    sum(diag(Reduce("+", zS[[g]]))) / (N * d)
  return(calcValue)
}

funpsiCUU <- function(g, Sk, zS, z, d, lambdanew, betaVar, bigTheta) {
  calcValue <- diag(Sk[,,g] - 2 * lambdanew[[1]] %*%
    betaVar[[g]] %*% Sk[, , g] +
    lambdanew[[1]] %*% bigTheta[[g]] %*%
    t(lambdanew[[1]]) +
    Reduce("+", zS[[g]]) / sum(z[, g])) * diag(d)
  return(calcValue)
}

funpsiCUC <- function(g, Sk, lambdanew, betaVar, bigTheta, zS, ng) {
  calcValue <- mean(diag(Sk[,,g] - 2 * lambdanew[[1]] %*%
    betaVar[[g]] %*% Sk[,,g] +
    lambdanew[[1]] %*% bigTheta[[g]] %*%
    t(lambdanew[[1]]) + Reduce("+", zS[[g]]) / ng[g]))
  return(calcValue)
}

funpsiCCU <- function(Sgav, lambdanew, betaVar, zS, N) {
  calcValue <- diag(Sgav - lambdanew[[1]] %*%
    betaVar[[1]] %*% Sgav) +
    diag(Reduce("+", zS[[1]]) / N)
  return(calcValue)
}

funpsiCCC <- function(Sgav, lambdanew, betaVar, zS, N, d) {
  calcValue <- mean(diag(Sgav - lambdanew[[1]] %*%
    betaVar[[1]] %*% Sgav)) +
    sum(diag(Reduce("+", zS[[1]])))/(N * d)
  return(calcValue)
}

funSgav <- function(g, ng, N, Sk) {
  calcValue <- ng[g] / N * Sk[,,g]
  return(calcValue)
}

# [END]
