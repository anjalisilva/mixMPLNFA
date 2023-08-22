
## functions for updating lambda
funLambdag <- function(g, Sk, betaVar, bigTheta) {
  Sk[,,g] %*% t(betaVar[[g]]) %*% solve(bigTheta[[g]])
}

funLambdaCUU1 <- function(g){ng[g] * solve(psi[[g]]) %*% Sk[,,g] %*% t(betaVar[[g]])}
funLambdaCUC1 <- function(g){ng[g] / (diag(psi[[g]])[1]) * (Sk[,,g] %*% t(betaVar[[g]]))}
funLambdaCUC2 <- function(g){ng[g] / (diag(psi[[g]])[1]) * bigTheta[[g]]}
funLambdaCCUnCCC <- function(){Sgav %*% t(betaVar[[1]]) %*% solve(bigThetaav)}

## functions for updating psi
funpsiUUU <- function(g){diag(Sk[,,g] - lambdanew[[g]] %*%
                                betaVar[[g]] %*% Sk[,,g] +
                                Reduce("+", zS[[g]]) / sum(z[, g])) * diag(d)}
funpsiUUC <- function(g){mean(diag(Sk[,,g] - lambdanew[[g]] %*%
                                     betaVar[[g]]%*%Sk[,,g] +
                                     Reduce("+", zS[[g]]) / ng[g]))}
funpsiUCU <- function(g){ng[g] / N*diag(Sk[,,g] - lambdanew[[g]] %*%
                                          betaVar[[g]] %*% Sk[,,g]) +
                                          diag(Reduce("+", zS[[g]]) / N)}
funPsiUCC <- function(g){(1 / d) * (ng[g] / N) *
                         sum(diag(Sk[,,g] - lambdanew[[g]] %*%
                         betaVar[[g]] %*% Sk[,,g])) +
                         sum(diag(Reduce("+", zS[[g]]))) / (N * d)}
fun_psi_cuu <- function(g){diag(Sk[,,g] - 2*lambdanew[[1]] %*%
                                betaVar[[g]]%*%Sk[,,g] +
                                lambdanew[[1]] %*% bigTheta[[g]] %*%
                                t(lambdanew[[1]]) +
                                Reduce("+", zS[[g]]) / sum(z[,g])) * diag(d)}
funpsiCUC <- function(g){mean(diag(Sk[,,g] - 2 * lambdanew[[1]] %*%
                                     betaVar[[g]] %*% Sk[,,g] +
                                     lambdanew[[1]] %*% bigTheta[[g]] %*%
                                     t(lambdanew[[1]]) + Reduce("+", zS[[g]]) / ng[g]))}
funpsiCCU <- function(){diag(Sgav - lambdanew[[1]] %*% betaVar[[1]] %*% Sgav) +
                          diag(Reduce("+", zS[[1]]) / N)}
funpsiCCC <- function(){mean(diag(Sgav - lambdanew[[1]] %*%
                          betaVar[[1]] %*% Sgav)) +
                          sum(diag(Reduce("+", zS[[1]])))/(N * d)}

funSgav <- function(g){ng[g] / N * Sk[,,g]}
