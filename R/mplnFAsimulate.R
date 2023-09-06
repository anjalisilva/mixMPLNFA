#' Generating Data Using Mixture of MPLN Factor Analyzers
#'
#' This function simulates data from a mixture of multivariate
#' Poisson-log normal factor analyzers family (PMPLNFA). In the
#' PMPLNFA framework restrictions are introduced to the model
#' parameters with the aim of obtaining parsimonious models, which
#' are sufficiently flexible for clustering purposes. Since the
#' largest contribution of free parameters is through the covariance
#' matrices, it is the focus for introduction of parsimony.
#'
#' @param nObservations A positive integer indicating the number
#'    of observations for the dataset.
#' @param dimensionality A positive integer indicating the dimensionality for the
#'    dataset.
#' @param mixingProportions A numeric vector that length equal to the number of total
#'    components, indicating the proportion of each component. Vector content should
#'    sum to 1.
#' @param mu A matrix of size (dimensionality x number of components), indicating
#'    the mean for each component. See example.
#' @param sigma A matrix of size ((dimensionality * number of components) x
#'    dimensionality), indicating the covariance matrix for each component.
#'    See example.
#' @param produceImage A character string indicating whether or not to
#'    produce an image. Options "Yes" or "No". Image will be produced as
#'    'Pairs plot of log-transformed data.png" in the current working
#'    directory.
#' @param ImageName A character string indicating name for image, if
#'    produceImage is set to "Yes". Default is "TwoComponents".
#'
#' @return Returns an S3 object of class mplnDataGenerator with results.
#' \itemize{
#'   \item dataset - Simulated dataset.
#'   \item trueMembership -A numeric vector indicating the membership of
#'      each observation.
#'   \item probaPost - A matrix indicating the posterior probability that
#'      each observation belong to the component/cluster.
#'   \item truenormfactors - A numeric vector indicating the true
#'      normalization factors used for adjusting the library sizes.
#'   \item observations - Number of observations in the simulated dataset.
#'   \item dimensionality - Dimensionality of the simulated dataset.
#'   \item mixingProportions - A numeric vector indicating the mixing
#'      proportion of each component.
#'   \item mu - True mean used for the simulated dataset.
#'   \item sigma - True covariances used for the simulated dataset.
#'}
#'
#' @examples
#' set.seed(100)
#' pfactors <- 3 # number of true latent factors
#' dimensionality <- 10 # dimensionality of observed data
#' G <- 2 # number of groups/clusters
#' nObservations <- 1000 ### sample size or number of observations
#' totData <- 100 # total number of datasets
#'
#' lambda_2 <- lambda_1 <- matrix(runif(pfactors * dimensionality, -1, 1), nrow = dimensionality)
#' psi_2 <- psi_1 <- diag(dimensionality) * runif(1)
#' sigma_1<-lambda_1%*%t(lambda_1)+psi_1
#' sigma_2<-lambda_2%*%t(lambda_2)+psi_2
#'
#' mu1 <- c(6, 3, 3, 6 ,3)
#' mu2 <- c(5, 3, 5, 3, 5)
#' mixingProportions <- c(0.32, 0.68)
#'
#'
#'
#' par[[1]]<-NULL # for pi
#' mu <- list() # for mu
#' Lambda <- list() # for Lambda
#' Psi <-list() # for Psi
#' par[[5]]<-list() # for Sigma
#' names(par)<-c("pi","mu","Lambda","Psi","Sigma")
#'
#'
#' mu <-list(mu_1,mu_2)
#' Lambda <- list(lambda_1,lambda_2)
#' Psi <- list(psi_1,psi_2)
#' Sigma <- list(sigma_1,sigma_2)
#'
#' @author Anjali Silva, \email{anjali@alumni.uoguelph.ca}
#'
#' @references
#' Aitchison, J. and C. H. Ho (1989). The multivariate Poisson-log normal distribution.
#' \emph{Biometrika} 76.
#'
#' Silva, A. et al. (2019). A multivariate Poisson-log normal mixture model
#' for clustering transcriptome sequencing data. \emph{BMC Bioinformatics} 20.
#' \href{https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2916-0}{Link}
#'
#' @export
#' @import stats
#' @importFrom mvtnorm rmvnorm
#' @importFrom edgeR calcNormFactors
#' @importFrom grDevices dev.off
#' @importFrom grDevices png
#' @importFrom graphics pairs
#'
#'

mplnFADataGenerator <- function(nObservations = 1000,
                                dimensionality = 10,
                                mixingProportions = c(0.32, 0.68),
                                mu = list(c(6, 3, 3, 6 ,3),
                                          c(5, 3, 5, 3, 5)),
                                Lambda = list(matrix(runif(pfactors * dimensionality, -1, 1), nrow = dimensionality),
                                              matrix(runif(pfactors * dimensionality, -1, 1), nrow = dimensionality)),
                                Psi = list(diag(dimensionality) * runif(1),
                                           diag(dimensionality) * runif(1)),
                                trueClusters = 2,
                                pfactors = 3,
                                modelName = "CCC") {
  for (run in 1:totData) {

    if(modelName == "CCC") {
      dat<-list() # Stores information about all run
      Ymatrix <- Xmatrix <- matrix(0, ncol = dimensionality, nrow = nObservations)
      Umatrix <- rmvnorm(nObservations, mean = rep(0, pfactors), sigma = diag(pfactors))
      zmatrix <- t(rmultinom(n = nObservations, size = 1, prob = mixingProportions))

      for (i in 1:nObservations) {
        grp <- which(zmatrix[i,] == 1)
        Xmatrix[i,] <- rmvnorm(1, mean = mu[[grp]] + Lambda[[grp]] %*% Umatrix[i,],
                            sigma = Psi[[grp]])
        for (j in 1:dimensionality) {
          Ymatrix[i, j] <- rpois(1, exp(Xmatrix[i, j]))
        }
      }

    }

    if(modelName == "UCC") {

    }

    if(modelName == "UUU") {

    }

    # sigma_1 <- lambda_1%*%t(lambda_1)+psi_1
    # sigma_2 <- lambda_2%*%t(lambda_2)+psi_2

    dat[[run]] <- list()
    dat[[run]][[1]] <- list(mu,
                            sigma,
                            Lambda,
                            Psi)
    dat[[run]][[2]] <- Ymatrix
    dat[[run]][[3]] <- Xmatrix
    dat[[run]][[4]] <- zmatrix
    dat[[run]][[5]] <- Umatrix

  }
  class(dat) <- "mplnFADataGenerator"
  return(dat)
}
