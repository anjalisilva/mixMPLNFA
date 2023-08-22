#' Clustering via Parsimonious Mixtures of MPLN Factor Analyzers Family
#'
#' Performs clustering using parsimonious mixtures of multivariate
#' Poisson-log normal factor analyzers family (PMPLNFA) via
#' variational Gaussian approximations. Model selection can
#' be done using AIC, AIC3, BIC and ICL.
#'
#' @details Starting values (argument: initMethod) play an
#'     important role in the successful operation of this algorithm.
#'
#' @param dataset A dataset of class matrix and type integer such that
#'    rows correspond to observations and columns correspond to variables.
#'    The dataset have dimensions n x d, where n is the total number of
#'    observations and d is the dimensionality. If rowSums are zero, these
#'    rows will be removed prior to cluster analysis.
#' @param membership A numeric vector of length nrow(dataset) containing the
#'    cluster membership of each observation. If not available,
#'    leave as "none".
#' @param gmin A positive integer specifying the minimum number of components
#'    to be considered in the clustering run.
#' @param gmax A positive integer, >= gmin, specifying the maximum number of
#'    components to be considered in the clustering run.
#' @param initMethod An algorithm for initialization. Current options are
#'    "kmeans", "random", "medoids", "clara", or "fanny". Default is "kmeans".
#' @param nInitIterations A positive integer or zero, specifying the number
#'    of initialization runs to be performed. This many runs, each with 10
#'    iterations, will be performed via MPLNClust and values from the run with
#'    highest log-likelihood will be used as initialization values. Default is 2.
#' @param normalize A string with options "Yes" or "No" specifying
#'     if normalization should be performed. Currently, normalization factors
#'     are calculated using TMM method of edgeR package. Default is "Yes".
#'
#' @return Returns an S3 object of class mixMPLNFA with results.
#' \itemize{
#'   \item dataset - The input dataset on which clustering is performed.
#'   \item dimensionality - Dimensionality of the input dataset.
#'   \item normalizationFactors - A vector of normalization factors used
#'      for input dataset.
#'   \item gmin - Minimum number of components/clusters considered in the clustering
#'      run.
#'   \item gmax - Maximum number of components/clusters considered in the clustering
#'      run.
#'   \item initalizationMethod - Method used for initialization.
#'   \item allResults - A list with all results.
#'   \item logLikelihood - A vector with value of final log-likelihoods for
#'      each component/cluster size.
#'   \item numbParameters - A vector with number of parameters for each
#'      component/cluster size.
#'   \item trueLabels - The vector of true labels, if provided by user.
#'   \item ICLresults - A list with all ICL model selection results.
#'   \item BICresults - A list with all BIC model selection results.
#'   \item AICresults - A list with all AIC model selection results.
#'   \item AIC3results - A list with all AIC3 model selection results.
#'   \item slopeHeuristics - If more than 10 models are considered, slope heuristic
#'      results as obtained via capushe::capushe().
#'   \item DjumpModelSelected - If more than 10 models are considered, slope heuristic
#'      results as obtained via capushe::capushe().
#'   \item DDSEModelSelected - If more than 10 models are considered, slope heuristic
#'      results as obtained via capushe::capushe().
#'   \item totalTime - Total time used for clustering and model selection.
#' }
#'
#'
#' @examples
#' \dontrun{
#'   modelChoice <- c("UUU", "UUC", "UCU", "UCC", "CUU", "CUC", "CCU", "CCC")
#' }
#'
#' @author {Anjali Silva, \email{anjali@alumni.uoguelph.ca}, Andrea Payne, \email{andreapayne@cmail.carleton.ca},
#' Sanjeena Dang, \email{sanjeenadang@cunet.carleton.ca}. }
#'
#' @references
#' Aitchison, J. and C. H. Ho (1989). The multivariate Poisson-log normal distribution.
#' \emph{Biometrika} 76.
#'
#' Akaike, H. (1973). Information theory and an extension of the maximum likelihood
#' principle. In \emph{Second International Symposium on Information Theory}, New York, NY,
#' USA, pp. 267–281. Springer Verlag.
#'
#' Arlot, S., Brault, V., Baudry, J., Maugis, C., and Michel, B. (2016).
#' capushe: CAlibrating Penalities Using Slope HEuristics. R package version 1.1.1.
#'
#' Biernacki, C., G. Celeux, and G. Govaert (2000). Assessing a mixture model for
#' clustering with the integrated classification likelihood. \emph{IEEE Transactions
#' on Pattern Analysis and Machine Intelligence} 22.
#'
#' Bozdogan, H. (1994). Mixture-model cluster analysis using model selection criteria
#' and a new informational measure of complexity. In \emph{Proceedings of the First US/Japan
#' Conference on the Frontiers of Statistical Modeling: An Informational Approach:
#' Volume 2 Multivariate Statistical Modeling}, pp. 69–113. Dordrecht: Springer Netherlands.
#'
#' Robinson, M.D., and Oshlack, A. (2010). A scaling normalization method for differential
#' expression analysis of RNA-seq data. \emph{Genome Biology} 11, R25.
#'
#' Schwarz, G. (1978). Estimating the dimension of a model. \emph{The Annals of Statistics}
#' 6.
#'
#' Silva, A. et al. (2019). A multivariate Poisson-log normal mixture model
#' for clustering transcriptome sequencing data. \emph{BMC Bioinformatics} 20.
#' \href{https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2916-0}{Link}
#'
#' @export
#' @importFrom edgeR calcNormFactors
#' @importFrom mclust unmap
#' @importFrom mclust map
#' @importFrom mvtnorm rmvnorm
#' @importFrom utils tail
#'

dataPath <- "/Users/as/Library/CloudStorage/GoogleDrive-anjali@alumni.uoguelph.ca/My Drive/UGuelph/Analysis_AnjaliLaptop/AlgorithmsAndCode/MPLN_Factor_analyzers/For_Anjali_3July2023/Sanjeena_July2023/"
dataset <- as.matrix(read.csv(file = paste0(dataPath, "Yfinal.csv"), header = TRUE, row.names = 1))
dim(dataset) #  712 20502
class(dataset)
typeof(dataset)

subestData <- dataset[1:20, 1:100]
dataset<-subestData

gmin <- 2
gmax <- 3
pmin <- 1
pmax <- 2
modelName <- "CCC"


PMPLNFA <- function(dataset, gmin, gmax, pmin, pmax, modelName) {
}

clustersize <- 2
pSize <-2


PMPLNFAind <- function(dataset, clustersize, pSize){


  # Initialize variables
    lambdanew <- psinew <-list()

    d <- ncol(dataset)
    N <- nrow(dataset)

    mu <- psi<- lambda <- sigmaVar <- isigma <- list()
    m <- S <- P <- Q <- list()

    # Other intermediate terms
    Sk <- array(0, c(d, d, clustersize))
    start <- GX <- dGX <- zS <- list()

    # Normalization factors
    normFactors <- edgeR::calcNormFactors(t(dataset), method = "TMM")
    libMat <- matrix(normFactors, nrow = N, ncol = d, byrow = F)

    kmeansOut <- stats::kmeans(log(dataset + 1), centers = gmax,
                               nstart = 100)$cluster

    # Initialize z
    z <- mclust::unmap(kmeansOut)
    # View(z)
    piG <- colSums(z) / N
    ng <- colSums(z)

    # Initialize parameters
    for (g in 1:clustersize) {

      obs <- which(z[, g] == 1)
        if(length(obs) > 1) {
          mu[[g]] <- colMeans(log(dataset[obs, ] + 1 / d))
          sigmaVar[[g]] <- var(log(dataset[obs, ] + 1 / d))
          # isigma issue
          # tried replacing this with Woodbury and moving it to line
          # 71 but still gives singularity error
          isigma[[g]] <- tryCatch(solve(sigmaVar[[g]]), error = function(err) NA)
          if(all(is.na(isigma[[g]]))) {
            isigma[[g]] <- diag(ncol(dataset[obs, ])) # if error with inverse
          }
        } else if (length(obs) == 1) { # if only one observation in a cluster
          mu[[g]] <- log(dataset[obs, ] + 1 / d)
          sigmaVar[[g]] <- diag(ncol(dataset))
          isigma[[g]] <- tryCatch(solve(sigmaVar[[g]]), error = function(err) NA)
          if(all(is.na(isigma[[g]]))) {
            isigma[[g]] <- diag(ncol(dataset))
          }
        }


      temp <- eigen(sigmaVar[[g]])
      lambda[[g]] <- matrix(NA, ncol = pSize, nrow = d)
      for (q in 1:pSize) {
        lambda[[g]][, q] <- temp$vectors[, q] * sqrt(temp$values[q])
      }
      psi[[g]] <- diag(sigmaVar[[g]] -
                                   lambda[[g]] %*%
                                   t(lambda[[g]])) * diag(d)
    }

    for (g in 1:clustersize) {
      start[[g]] <- log(dataset + 1) ###Starting value for M
      m[[g]] <- log(dataset + 1)
      S[[g]] <- list()
      for (i in 1:N) {
        S[[g]][[i]] <- diag(d) * 0.000000001
      }
    }

    it <- 1
    aloglik <- loglik <- NULL
    checks <- aloglik[c(1, 2, 3)] <- 0
    itMax <- 100

    # Begin clusterig
    while (checks == 0) {


      for (g in 1:clustersize) {
        GX[[g]] <- dGX[[g]] <- zS[[g]] <- list()
        z[is.nan(z)] <- 0
        for (i in 1:N) {
          #print(i)
          dGX[[g]][[i]] <- diag(exp(log(libMat[i, ]) + start[[g]][i, ] +
                                      0.5 * diag(S[[g]][[i]])), d) + isigma[[g]]
          S[[g]][[i]] <- solve(dGX[[g]][[i]])
          zS[[g]][[i]] <- z[i,g] * S[[g]][[i]]
          GX[[g]][[i]] <- dataset[i, ] - exp(start[[g]][i, ] + log(libMat[i, ]) +
                                            0.5 * diag(S[[g]][[i]])) -
                                            (isigma[[g]]) %*% (start[[g]][i, ] -
                                            mu[[g]])
          m[[g]][i, ] <- start[[g]][i, ] + S[[g]][[i]] %*% GX[[g]][[i]]
        }

        start[[g]] <- m[[g]]

        # Updating mu
        mu[[g]] <- colSums(z[, g] * m[[g]]) / sum(z[, g])

        # Updating Sample covariance
        muMatrix <- matrix(rep(mu[[g]], N), nrow = N, byrow = TRUE)
        res <- m[[g]] - muMatrix
        temp <- cov.wt(res, wt = z[,g], center = FALSE, method = "ML")
        Sk[,,g] <- temp$cov
      }


      if (it < 10) repmax <- 10 else repmax <- 1
      for (rep in 1:repmax) {
        lambdaOld <- lambda
        psiOld <-psi
        updates <- modelUpdates(modelName = modelName,
                                zS = zS,
                                ng = ng,
                                z = z,
                                lambda = lambda,
                                isigma = isigma,
                                clustersize = clustersize,
                                pmax.var = pSize,
                                Sk = Sk,
                                psi = psi)


        sigmaVar <- updates$sigmaVar
        psi <- updates$psi
        lambda <- updates$lambda
        par <- updates$par


        for (g in 1:clustersize) {
          ## check diag dimension
          isigma[[g]] <- solve(psi[[g]]) -
            (solve(psi[[g]]) %*% lambda[[g]] %*%
            solve(diag(dim(bigTheta[[g]])[1]) +
            (t(lambda[[g]]) %*% solve(psi[[g]]) %*%
            lambda[[g]])) %*% t(lambda[[g]]) %*%
            solve(psi[[g]]))
        }
      }


      piG <- colSums(z) / N
      ng <- colSums(z)
      # libMat<-matrix(normFactors,ncol=d,nrow=N, byrow=T) ###Matrix containing normalization factor
      ### Some useful functions
      funFive <- function(x, y = isigma[[g]]) {
        sum(diag(x %*% y))
      }

      Ffunction <- matrix(NA, ncol = clustersize, nrow = N)

      for (g in 1:clustersize) {
        two <- rowSums(exp(m[[g]] + log(libMat) + 0.5 *
                             matrix(unlist(lapply(S[[g]], diag)), ncol = d, byrow = TRUE)))
        five <- 0.5 * unlist(lapply(S[[g]], funFive))
        six <- 0.5 * log(unlist(lapply(S[[g]], det)))
        Ffunction[, g] <- piG[g] * exp(rowSums(m[[g]] * dataset) - two -
                                  rowSums(lfactorial(dataset)) +
                                  rowSums(log(libMat) * dataset) -
                                  0.5 * mahalanobis(m[[g]], center = mu[[g]], cov = isigma[[g]], inverted = TRUE) -
                                  five + six - 0.5 * log(det(sigmaVar[[g]])) - d / 2)
      }


      loglik[it] <- sum(log(rowSums(Ffunction)))
      z <- Ffunction / rowSums(Ffunction)
      if (it <= 5) {
        z[z == "NaN"] <- 0
      }
      #print(loglik)
      #plot(loglik,type="l")


      if (it > 5) {
        #Aitkaine's stopping criterion
        if ((loglik[it - 1] - loglik[it - 2]) == 0) {
          checks <- 1
        } else {
          a <- (loglik[it] - loglik[it-1]) / (loglik[it - 1] - loglik[it - 2])
          addTo <- (1 / (1 - a) * (loglik[it] - loglik[it-1]))
          # }
          aloglik[it] <- loglik[it - 1] + addTo
          if (abs(aloglik[it] - loglik[it - 1]) < 0.001) {
            checks <- 1
          } else {
            checks <- checks
          }
        }
      }


      # print(it)
      it <- it + 1
      if (it == itMax) {
        checks <- 1}
    }
    # plot(loglik,type="l")

    # par<-G*(d*pmax-0.5*pmax*(pmax-1))+G*d
    ##par from covariance only has the covariance parameters so now we need to add the parameters for the mean and pi
    BIC <- 2 * loglik[it - 1] - (par + (clustersize - 1) + clustersize * d) * log(N)
    mapz <- matrix(0, ncol = clustersize, nrow = N)
    for (g in 1:clustersize) {
      mapz[which(mclust::map(z) == g), g] <- 1
    }
    forICL <- function(g) {if (sum(mapz[, g]) == 0){0} else
      if (sum(mapz[, g]) > 0){sum(log(z[which(mapz[, g] == 1), g]))}}
    ICL <- BIC - 2 * sum(sapply(1:clustersize, forICL))
    true <- NA

    modelList <- list()
    modelList[[1]]<- piG
    modelList[[2]]<- mu
    modelList[[3]]<- sigmaVar
    modelList[[4]]<- lambda
    modelList[[5]]<- psi
    modelList[[6]]<- z
    modelList[[7]]<- loglik
    modelList[[8]]<- kmeansOut
    modelList[[9]]<- true
    modelList[[10]]<- BIC
    modelList[[11]]<- ICL
    modelList[[12]]<- modelName
    modelList[[13]]<- clustersize
    modelList[[14]]<- pmax
    names(modelList)<-c("piG", "mu", "sigmaVar", "lambda", "psi", "z", "loglik", "kmeans", "true", "BIC", "ICL", "modelName", "G", "p")


    class(modelList) <- "mixMPLNFA"
    ## calc best model
    return(modelList)
  }


modelUpdates <- function(modelName,
                      zS,
                      ng,
                      z,
                      lambda,
                      isigma,
                      clustersize,
                      pmax.var,
                      Sk,
                      psi) {

  betaVar <- bigTheta <- sigmaVar <- list()
  Sk <<- Sk
  zS <<- zS
  ng <<- ng
  z <<- z
  psi <<- psi
  lambda <<- lambda
  if (substr(modelName, 1, 3) == "UUU") {
    for (g in 1:clustersize) {
      betaVar[[g]] <- t(lambda[[g]]) %*% isigma[[g]]
      bigTheta[[g]] <- diag(pmax.var) - betaVar[[g]] %*%
        lambda[[g]] + betaVar[[g]] %*% Sk[, , g] %*%
        t(betaVar[[g]])
    }
    betaVar <<- betaVar
    bigTheta <<- bigTheta
    for (g in 1:clustersize) {
      lambdanew[[g]] <- funLambdag(g)
    }
    lambdanew <<- lambdanew
    for (g in 1:clustersize) {
      psinew[[g]] <- fun_psi_uuu(g)
    }
    for (g in 1:clustersize) {
      sigmaVar[[g]] <- (lambdanew[[g]] %*% t(lambdanew[[g]]) +
                           psinew[[g]])
      lambda[[g]] <- lambdanew[[g]]
      psi[[g]] <- psinew[[g]]
    }
    par <- clustersize * (d * pmax.var - 0.5 * pmax.var * (pmax.var -
                                                   1)) + clustersize * d
  }
  if (substr(modelName, 1, 3) == "UUC") {
    for (g in 1:clustersize) {
      betaVar[[g]] <- t(lambda[[g]]) %*% isigma[[g]]
      bigTheta[[g]] <- diag(pmax.var) - betaVar[[g]] %*%
        lambda[[g]] + betaVar[[g]] %*% Sk[, , g] %*%
        t(betaVar[[g]])
    }
    betaVar <<- betaVar
    bigTheta <<- bigTheta
    for (g in 1:clustersize) {
      lambdanew[[g]] <- funLambdag(g)
    }
    lambdanew <<- lambdanew
    for (g in 1:clustersize) {
      psinew[[g]] <- fun_psi_uuc(g) * diag(d)
    }
    psi <- psinew
    for (g in 1:clustersize) {
      sigmaVar[[g]] <- (lambdanew[[g]] %*% t(lambdanew[[g]]) +
                           psinew[[g]])
      lambda[[g]] <- lambdanew[[g]]
    }
    par <- clustersize * (d * pmax.var - 0.5 * pmax.var * (pmax.var -
                                                   1)) + clustersize
  }
  if (substr(modelName, 1, 3) == "UCU") {
    for (g in 1:clustersize) {
      betaVar[[g]] <- t(lambda[[g]]) %*% isigma[[g]]
      bigTheta[[g]] <- diag(pmax.var) - betaVar[[g]] %*%
        lambda[[g]] + betaVar[[g]] %*% Sk[, , g] %*%
        t(betaVar[[g]])
    }
    betaVar <<- betaVar
    bigTheta <<- bigTheta
    for (g in 1:clustersize) {
      lambdanew[[g]] <- funLambdag(g)
    }
    lambdanew <<- lambdanew
    psinew[[1]] <- rowSums(sapply(1:clustersize, fun_psi_ucu)) * diag(d)
    for (g in 1:clustersize) {
      sigmaVar[[g]] <- (lambdanew[[g]] %*% t(lambdanew[[g]]) +
                           psinew[[1]])
      lambda[[g]] <- lambdanew[[g]]
      psi[[g]] <- psinew[[1]]
    }
    par <- clustersize * (d * pmax.var - 0.5 * pmax.var * (pmax.var -
                                                   1)) + d
  }
  if (substr(modelName, 1, 3) == "UCC") {
    for (g in 1:clustersize) {
      betaVar[[g]] <- t(lambda[[g]]) %*% isigma[[g]]
      bigTheta[[g]] <- diag(pmax.var) - betaVar[[g]] %*%
        lambda[[g]] + betaVar[[g]] %*% Sk[, , g] %*%
        t(betaVar[[g]])
    }
    betaVar <<- betaVar
    bigTheta <<- bigTheta
    for (g in 1:clustersize) {
      lambdanew[[g]] <- funLambdag(g)
    }
    lambdanew <<- lambdanew
    psinew[[1]] <- sum(sapply(1:clustersize, funPsiUCC)) * diag(d)
    for (g in 1:clustersize) {
      sigmaVar[[g]] <- (lambdanew[[g]] %*% t(lambdanew[[g]]) +
                           psinew[[1]])
      lambda[[g]] <- lambdanew[[g]]
      psi[[g]] <- psinew[[1]]
    }
    par <- clustersize * (d * pmax.var - 0.5 * pmax.var * (pmax.var -
                                                   1)) + 1
  }
  if (substr(modelName, 1, 3) == "CUU") {
    for (g in 1:clustersize) {
      betaVar[[g]] <- t(lambda[[1]]) %*% isigma[[g]]
      bigTheta[[g]] <- diag(pmax.var) - betaVar[[g]] %*%
        lambda[[1]] + betaVar[[g]] %*% Sk[, , g] %*%
        t(betaVar[[g]])
    }
    betaVar <<- betaVar
    bigTheta <<- bigTheta
    for_lam_1 <- matrix(rowSums(sapply(1:clustersize, funLambdaCUU1)),
                        d, pmax.var)
    for_lam_2 <- list()
    for (k in 1:d) {
      k <<- k
      if (pmax.var > 1)
        for_lam_2[[k]] <- solve(matrix(rowSums(sapply(1:clustersize,
                                                      funLambdaCUU2)), pmax.var, pmax.var))
      else for_lam_2[[k]] <- solve(matrix(sum(sapply(1:clustersize,
                                                     funLambdaCUU2)), pmax.var, pmax.var))
    }
    lambdanew[[1]] <- matrix(NA, d, pmax.var)
    for (k in 1:d) {
      lambdanew[[1]][k, ] <- for_lam_1[k, ] %*% for_lam_2[[k]]
    }
    lambdanew <<- lambdanew
    for (g in 1:clustersize) {
      psinew[[g]] <- fun_psi_cuu(g) * diag(d)
      sigmaVar[[g]] <- (lambdanew[[1]] %*% t(lambdanew[[1]]) +
                           psinew[[g]])
      lambda[[g]] <- lambdanew[[1]]
    }
    psi <- psinew
    par <- (d * pmax.var - 0.5 * pmax.var * (pmax.var - 1)) +
      clustersize * d
  }
  if (substr(modelName, 1, 3) == "CUC") {
    for (g in 1:clustersize) {
      betaVar[[g]] <- t(lambda[[1]]) %*% isigma[[g]]
      bigTheta[[g]] <- diag(pmax.var) - betaVar[[g]] %*%
        lambda[[1]] + betaVar[[g]] %*% Sk[, , g] %*%
        t(betaVar[[g]])
    }
    betaVar <<- betaVar
    bigTheta <<- bigTheta
    if (pmax.var > 1) {
      lambdanew[[1]] <- matrix(rowSums(sapply(1:clustersize, funLambdaCUC1)),
                               d, pmax.var) %*% solve(matrix(rowSums(sapply(1:clustersize,
                                                                            funLambdaCUC2)), pmax.var, pmax.var))
    }
    else lambdanew[[1]] <- matrix(rowSums(sapply(1:clustersize, funLambdaCUC1)),
                                  d, pmax.var) %*% solve(matrix(sum(sapply(1:clustersize, funLambdaCUC2)),
                                                                pmax.var, pmax.var))
    lambda[[1]] <- lambdanew[[1]]
    lambdanew <<- lambdanew
    for (g in 1:clustersize) {
      psinew[[g]] <- fun_psi_cuc(g) * diag(d)
      sigmaVar[[g]] <- (lambdanew[[1]] %*% t(lambdanew[[1]]) +
                           psinew[[g]])
      lambda[[g]] <- lambdanew[[1]]
      psi[[g]] <- psinew[[g]]
    }
    par <- (d * pmax.var - 0.5 * pmax.var * (pmax.var - 1)) +
      clustersize
  }
  if (substr(modelName, 1, 3) == "CCU") {
    betaVar[[1]] <- t(lambda[[1]]) %*% isigma[[1]]
    Sg_av <- matrix(rowSums(sapply(1:clustersize, funSgav)), d, d)
    bigTheta[[1]] <- diag(pmax.var) - betaVar[[1]] %*% lambda[[1]] +
      betaVar[[1]] %*% Sg_av %*% t(betaVar[[1]])
    bigTheta_av <- bigTheta[[1]]
    betaVar <<- betaVar
    bigTheta_av <<- bigTheta_av
    Sg_av <<- Sg_av
    lambdanew[[1]] <- funLambdaCCUnCCC()
    lambdanew <<- lambdanew
    psinew[[1]] <- fun_psi_ccu() * diag(d)
    bigtheta_old <- bigTheta[[1]]
    for (g in 1:clustersize) {
      sigmaVar[[g]] <- (lambdanew[[1]] %*% t(lambdanew[[1]]) +
                           psinew[[1]])
      lambda[[g]] <- lambdanew[[1]]
      psi[[g]] <- psinew[[1]]
      bigTheta[[g]] <- bigtheta_old
    }
    bigTheta <<- bigTheta
    par <- (d * pmax.var - 0.5 * pmax.var * (pmax.var - 1)) +
      d
  }
  if (substr(modelName, 1, 3) == "CCC") {
    betaVar[[1]] <- t(lambda[[1]]) %*% isigma[[1]]
    Sgav <- matrix(rowSums(sapply(1:clustersize, funSgav)), d, d)
    bigTheta[[1]] <- diag(pmax.var) - betaVar[[1]] %*% lambda[[1]] +
      betaVar[[1]] %*% Sgav %*% t(betaVar[[1]])
    bigThetaav <- bigTheta[[1]]
    betaVar <<- betaVar
    bigThetaav <<- bigThetaav
    Sgav <<- Sgav
    lambdanew[[1]] <- funLambdaCCUnCCC()
    lambdanew <<- lambdanew
    psinew[[1]] <- fun_psi_ccc() * diag(d)
    bigthetaOld <- bigTheta[[1]]
    for (g in 1:clustersize) {
      sigmaVar[[g]] <- (lambdanew[[1]] %*% t(lambdanew[[1]]) +
                           psinew[[1]])
      lambda[[g]] <- lambdanew[[1]]
      psi[[g]] <- psinew[[1]]
      bigTheta[[g]] <- bigthetaOld
    }
    bigTheta <<- bigTheta
    par <- (d * pmax.var - 0.5 * pmax.var * (pmax.var - 1)) +
      1
  }
  return(list(sigmaVar = sigmaVar, psi = psi, lambda = lambda,
              par = par, isigma = isigma))
}



# [END]
