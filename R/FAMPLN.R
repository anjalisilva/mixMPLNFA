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
#'    observations and d is the dimensionality. The program will cluster
#'    n observations into gmin:gmax groups. If rowSums are zero,
#'    these rows should be removed first.
#' @param membership A numeric vector of length nrow(dataset) containing the
#'    cluster membership of each observation. If not available,
#'    leave as "none".
#' @param gmin A positive integer specifying the minimum number of components
#'    to be considered in the clustering run.
#' @param gmax A positive integer, >= gmin, specifying the maximum number of
#'    components to be considered in the clustering run.
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


PMPLNFA <- function(dataset,
                    gmin = 1,
                    gmax = 2,
                    pmin = 1,
                    pmax = 2,
                    modelName = "CCC",
                    trueLabels,
                    normalize = "Yes") {

  # Keeping track of time
  ptm <- proc.time()

  # Performing checks of user input
  if (typeof(dataset) != "double" & typeof(dataset) != "integer") {
    stop("Dataset type needs to be integer.")
  }

  if (any((dataset %% 1 == 0) == FALSE)) {
    stop("Dataset should be a matrix of counts.")
  }

  if (is.matrix(dataset) != TRUE) {
    stop("Dataset needs to be a matrix.")
  }

  if (any(colSums(dataset) <= 0)) {
    stop("Column sums cannot be less than or equal to 0. Double check dataset.")
  }

  dimensionality <- ncol(dataset)
  nObservations <- nrow(dataset)

  if(is.numeric(gmin) != TRUE || is.numeric(gmax) != TRUE) {
    stop("Class of gmin and gmin should be numeric.")
  }

  if(is.numeric(pmin) != TRUE || is.numeric(pmax) != TRUE) {
    stop("Class of pmin and pmin should be numeric.")
  }

  if (gmax < gmin) {
    stop("gmax cannot be less than gmin.")
  }

  if (pmax < pmin) {
    stop("pmax cannot be less than pmin.")
  }

  if (pmax >= dimensionality) {
    stop("pmax has to be less than d (dimensionality)")
  }

  if (gmax > nObservations) {
    stop("gmax cannot be larger than nrow(dataset).")
  }

  if (gmax > nObservations) {
    stop("gmax cannot be larger than nrow(dataset).")
  }

  # To add:
  # Initialization methods
  # Remove rows with only zeros, if present


  if (is.character(normalize) != TRUE) {
    stop("normalize should be a string of class character specifying
      if normalization should be performed.")
  }

  if (normalize != "Yes" && normalize != "No") {
    stop("normalize should be a string indicating Yes or No, specifying
      if normalization should be performed.")
  }

  # Calculating normalization factors
  if(normalize == "Yes") {
    normFactors <- log(as.vector(edgeR::calcNormFactors(as.matrix(dataset),
                                                        method = "TMM")))
    normFactors <- matrix(normFactors, nrow = nObservations,
                          ncol = dimensionality, byrow = T)
  } else if(normalize == "No") {
    normFactors <- matrix(1, nrow = N, ncol = d, byrow = T)
  } else{
    stop("normalize should be 'Yes' or 'No' ")
  }


  clusterResults <- list() # to save cluster output
  for (gmodel in seq_along(1:(gmax - gmin + 1))) {
    for (pmodel in seq_along(1:(pmax - pmin + 1))) {

      if(length(1:(gmax - gmin + 1)) == gmax) {
        clustersize <- gmodel
      } else if(length(1:(gmax - gmin + 1)) < gmax) {
        clustersize <- seq(gmin, gmax, 1)[gmodel]
      }

      # iterating through q
      clusterResults[[gmodel]] < - list()
      if(length(1:(pmax - pmin + 1)) == pmax) {
        pSize <- pmodel
      } else if(length(1:(pmax - pmin + 1)) < pmax) {
        pSize <- seq(pmin, pmax, 1)[pmodel]
      }


      cat("\n Running for g =", clustersize, "and q = ", pSize)
      clusterResults[[gmodel]][[pSize]] <- PMPLNFAind(dataset = dataset,
                                             clustersize = clustersize,
                                             pSize = pSize,
                                             modelName = modelName,
                                             normFactors = normFactors)

    }

  }
    names(clusterResults) <- paste0(rep("G=", length(seq(gmin, gmax, 1))), seq(gmin, gmax, 1))



    BIC <- ICL <- AIC <- AIC3 <- Djump <- DDSE <- nParameters <- logLikelihood <- vector()

    for(g in seq_along(1:(gmax - gmin + 1))) {
      # save the final log-likelihood
      logLikelihood[g] <- unlist(tail(parallelRun[[g]]$allResults$logLikelihood,
                                      n = 1))

      if(length(1:(gmax - gmin + 1)) == gmax) {
        clustersize <- g
      } else if(length(1:(gmax - gmin + 1)) < gmax) {
        clustersize <- seq(gmin, gmax, 1)[g]
      }

    }

    return(invisible(NULL))
  }





PMPLNFAind <- function(dataset,
                       clustersize,
                       pSize,
                       modelName,
                       normFactors) {


  # Initialize variables
    d <- ncol(dataset)
    N <- nrow(dataset)

    mu <- psi<- lambda <- sigmaVar <- isigma <- list()
    m <- S <- P <- Q <- lambdanew <- list()

    # Other intermediate terms
    Sk <- array(0, c(d, d, clustersize))
    start <- GX <- dGX <- zS <- list()

    # Normalization factors
    # normFactors <- edgeR::calcNormFactors(t(dataset), method = "TMM")
    # libMatFull <- libMat <- matrix(normFactors, nrow = N, ncol = d, byrow = F)
    # libMatFull <- libMat <- matrix(1, nrow = N, ncol = d, byrow = T)
    libMatFull <- libMat <- normFactors

    kmeansOut <- stats::kmeans(x = log(dataset + 1),
                               centers = clustersize,
                               nstart = 100)$cluster

    # Initialize z
    z <- mclust::unmap(classification = kmeansOut)
    # View(z)
    piG <- colSums(z) / N
    ng <- colSums(z)

    # Initialize parameters
    for (g in 1:clustersize) {

      obs <- which(z[, g] == 1)
        if(length(obs) > 1) {
          mu[[g]] <- colMeans(log(dataset[obs, ] + 1 / d))
          sigmaVar[[g]] <- var(log(dataset[obs, ] + 1 / d))

          # solve issue anticipation
          # isigma[[g]]<-solve(sigma.var[[g]])
          # SD comment: tried replacing this with Woodbury and moving it to line
          # 71 but still gives singularity error
          isigma[[g]] <- tryCatch(solve(sigmaVar[[g]]), error = function(err) NA)
          if(all(is.na(isigma[[g]]))) {
            isigma[[g]] <- diag(ncol(dataset[obs, ])) # if error with inverse
          }

        } else if (length(obs) == 1) { # if only one observation in a cluster
          mu[[g]] <- log(dataset[obs, ] + 1 / d)
          sigmaVar[[g]] <- diag(ncol(dataset))

          # solve issue anticipation
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
      print(it)

      for (g in 1:clustersize) {
        GX[[g]] <- dGX[[g]] <- zS[[g]] <- list()
        z[is.nan(z)] <- 0
        for (i in 1:N) {
          #print(i)
          dGX[[g]][[i]] <- diag(exp(log(libMat[i, ]) + start[[g]][i, ] +
                                      0.5 * diag(S[[g]][[i]])), d) + isigma[[g]]

          # solve issue anticipation
          # S[[g]][[i]] <- solve(dGX[[g]][[i]])
          # solve if singular
          S[[g]][[i]] <- tryCatch(solve(dGX[[g]][[i]]), error = function(err) NA)
          if(all(is.na(S[[g]][[i]]))) {
            S[[g]][[i]] <- diag(ncol(dGX[[g]][[i]]))
          }

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
        temp <- stats::cov.wt(x = res, wt = z[, g],
                              center = FALSE,
                              method = "ML")
        Sk[,,g] <- temp$cov
      }


      if (it < 10) {
        repmax <- 10
      } else {
        repmax <- 1
      }

      for (rep in 1:repmax) {
        cat("rep = ", rep, "\n")
        lambdaOld <- lambda
        psiOld <- psi
        updates <- modelUpdates(modelName = modelName,
                                zS = zS,
                                ng = ng,
                                z = z,
                                lambda = lambda,
                                isigma = isigma,
                                clustersize = clustersize,
                                pmaxVar = pSize,
                                Sk = Sk,
                                psi = psi,
                                d = d,
                                N = N)


        sigmaVar <- updates$sigmaVar
        psi <- updates$psi
        lambda <- updates$lambda
        par <- updates$par
        bigTheta <- updates$bigTheta
        isigma <- updates$isigma


        for (g in 1:clustersize) {

          # solve issue anticipation
          # isigma[[g]] <- solve(psi[[g]]) -
          #   (solve(psi[[g]]) %*% lambda[[g]] %*%
          #      solve(diag(dim(bigTheta[[g]])[1]) +
          #              (t(lambda[[g]]) %*% solve(psi[[g]]) %*%
          #                 lambda[[g]])) %*% t(lambda[[g]]) %*%
          #      solve(psi[[g]]))

          solvePsig <- tryCatch(solve(psi[[g]]), error = function(err) NA)
          if(all(is.na(solvePsig))) {
            solvePsig <- diag(ncol(psi[[g]]))
          }

          isigma[[g]] <- solvePsig -
            (solvePsig %*% lambda[[g]] %*%
            solve(diag(dim(bigTheta[[g]])[1]) +
            (t(lambda[[g]]) %*% solvePsig %*%
            lambda[[g]])) %*% t(lambda[[g]]) %*%
              solvePsig)
        }
      }


      piG <- colSums(z) / N
      ng <- colSums(z)
      # libMat<-matrix(normFactors,ncol=d,nrow=N, byrow=T) ###Matrix containing normalization factor
      ### Some useful functions
      funFive <- function(x, y, g) {
        ySpecific = y[[g]]
        funFiveReturn <- sum(diag(x %*% ySpecific))
        return(funFiveReturn)
      }

      Ffunction <- matrix(NA, ncol = clustersize, nrow = N)

      for (g in 1:clustersize) {
        two <- rowSums(exp(m[[g]] + log(libMatFull) + 0.5 *
                             matrix(unlist(lapply(S[[g]], diag)), ncol = d, byrow = TRUE)))
        five <- 0.5 * unlist(lapply(S[[g]], funFive, y = isigma, g = g))

        # detecting if six gives -Inf which happens when
        # dimensionality is large, e.g., 500
        six <- 0.5 * log(unlist(lapply(S[[g]], det)))
        if(any(is.infinite(six)) == TRUE) { # detect if -Inf
          six <- 0.5 * log(unlist(lapply(testing, det)) + 1e-10)
        }

        Ffunction[, g] <- piG[g] * exp(rowSums(m[[g]] * dataset) -
                                       two -
                                       rowSums(lfactorial(dataset)) +
                                       rowSums(log(libMatFull) * dataset) -
                                       0.5 * mahalanobis(m[[g]], center = mu[[g]], cov = isigma[[g]], inverted = TRUE) -
                                       five +
                                       six -
                                       0.5 * log(det(sigmaVar[[g]])) -
                                       d / 2)
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


      print(it)
      it <- it + 1
      if (it == itMax) {
        checks <- 1
      }
    }
    # plot(loglik,type="l")

    # par<-G*(d*pmax-0.5*pmax*(pmax-1))+G*d
    ##par from covariance only has the covariance parameters so now we need to add the parameters for the mean and pi
    BIC <- 2 * loglik[it - 1] - (par + (clustersize - 1) + clustersize * d) * log(N)
    mapz <- matrix(0, ncol = clustersize, nrow = N)
    for (g in 1:clustersize) {
      mapz[which(mclust::map(z) == g), g] <- 1
    }
    forICL <- function(g, mapz, z) {
      if (sum(mapz[, g]) == 0) {
        returnValICL <- 0
      } else if (sum(mapz[, g]) > 0) {
        returnValICL <- sum(log(z[which(mapz[, g] == 1), g]))
      }
      return(returnValICL)
    }
    ICL <- BIC - 2 * sum(sapply(1:clustersize, forICL,
                         mapz = mapz, z = z))
    true <- NA

    modelList <- list()
    modelList[[1]] <- piG
    modelList[[2]] <- mu
    modelList[[3]] <- sigmaVar
    modelList[[4]] <- lambda
    modelList[[5]] <- psi
    modelList[[6]] <- z
    modelList[[7]] <- loglik
    modelList[[8]] <- kmeansOut
    modelList[[9]] <- BIC
    modelList[[10]] <- ICL
    modelList[[11]] <- modelName
    modelList[[12]] <- clustersize
    modelList[[13]] <- pSize
    names(modelList) <-c("piG", "mu", "sigmaVar",
                        "lambda", "psi", "z", "loglik",
                        "kmeans", "BIC", "ICL",
                        "modelName", "G", "p")

    class(modelList) <- "mixMPLNFA"
    return(modelList)

  }

# [END]
