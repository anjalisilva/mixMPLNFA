
## functions for updating lambda
funLambdag <- function(g){Sk[,,g] %*% t(beta.var[[g]]) %*% solve(bigTheta[[g]])}
funLambdaCUU1 <- function(g){ng[g] * solve(psi[[g]]) %*% Sk[,,g] %*% t(beta.var[[g]])}
funLambdaCUC1 <- function(g){ng[g] / (diag(psi[[g]])[1]) * (Sk[,,g] %*% t(beta.var[[g]]))}
funLambdaCUC2 <- function(g){ng[g] / (diag(psi[[g]])[1]) * bigTheta[[g]]}
funLambdaCCUnCCC <- function(){Sgav %*% t(beta.var[[1]]) %*% solve(bigThetaav)}

## functions for updating psi
fun_psi_uuu <- function(g){diag(Sk[,,g] - lambdanew[[g]] %*%
                                beta.var[[g]] %*% Sk[,,g] +
                                Reduce("+", zS[[g]]) / sum(z[, g])) * diag(d)}
fun_psi_uuc <- function(g){mean(diag(Sk[,,g] - lambdanew[[g]] %*%
                                     beta.var[[g]]%*%Sk[,,g] +
                                     Reduce("+", zS[[g]]) / ng[g]))}
fun_psi_ucu <- function(g){ng[g] / N*diag(Sk[,,g] - lambdanew[[g]] %*%
                                          beta.var[[g]] %*% Sk[,,g]) +
                                          diag(Reduce("+", zS[[g]]) / N)}
funPsiUCC <- function(g){(1 / d) * (ng[g] / N) *
                         sum(diag(Sk[,,g] - lambdanew[[g]] %*%
                         beta.var[[g]] %*% Sk[,,g])) +
                         sum(diag(Reduce("+", zS[[g]]))) / (N * d)}
fun_psi_cuu <- function(g){diag(Sk[,,g] - 2*lambdanew[[1]] %*%
                                beta.var[[g]]%*%Sk[,,g] +
                                lambdanew[[1]] %*% bigTheta[[g]] %*%
                                t(lambdanew[[1]]) +
                                Reduce("+", zS[[g]]) / sum(z[,g])) * diag(d)}
fun_psi_cuc <- function(g){mean(diag(Sk[,,g] - 2 * lambdanew[[1]] %*%
                                     beta.var[[g]] %*% Sk[,,g] +
                                     lambdanew[[1]] %*% bigTheta[[g]] %*%
                                     t(lambdanew[[1]]) + Reduce("+", zS[[g]]) / ng[g]))}
fun_psi_ccu <- function(){diag(Sgav - lambdanew[[1]] %*% beta.var[[1]] %*% Sgav) +
                          diag(Reduce("+", zS[[1]]) / N)}
fun_psi_ccc <- function(){mean(diag(Sgav - lambdanew[[1]] %*%
                          beta.var[[1]] %*% Sgav)) +
                          sum(diag(Reduce("+", zS[[1]])))/(N * d)}

funSgav <- function(g){ng[g] / N * Sk[,,g]}
