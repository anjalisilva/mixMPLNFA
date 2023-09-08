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
#' @param numDatasets A positive integer indicating the number
#'    of datasets to be generated. Default value is 10.
#' @param nObservations A positive integer indicating the number
#'    of observations for the dataset. Or the sample size. Default
#'    value is 1000.
#' @param dimensionality A positive integer indicating the
#'    dimensionality for the dataset. Default value is 10.
#' @param dimensionality A positive integer indicating the
#'    dimensionality for the dataset. Default value is 10.
#' @param trueClusters A positive integer indicating the number
#'    of total components or clusters. Default value is 2.
#' @param pfactors A positive integer indicating the number
#'    of total latent factors. Default value is 3.
#' @param modelName A character string indicating the model name to generate
#'     covariance matrix, Sigma. Since the largest contribution of free parameters
#'     is through the covariance matrices, it is the focus for introduction
#'     of parsimony. The constraints can be imposed on Lambda (loading matrix)
#'     and Psi (error variance and isotropic) which are used to generate Sigma.
#'     The 'C' stands for constrained and 'U' stands for unconstrained. The order
#'     goes as loading matrix (Lambda), error variance (Psi) and isotropic (Psi),
#'     respectively. Example, if the loading matrix (Lambda), error variance (Psi)
#'     and isotropic are all constrained, then select 'CCC'. Options are one of
#'     "UUU", "UUC", "UCU", "UCC", "CUU", "CUC", "CCU", "CCC".
#' @param mu A list of length equal to the value provided to 'trueClusters' argument
#'     and each element should have length equal to the value provided to
#'     'dimensionality' argument. See example or default value.
#' @param Lambda A list of length 'trueClusters' with each list element having
#'     a matrix with rows equal to value provided to 'dimensionality' argument
#'     and columns equal to value provided to 'pfactors' argument. The values
#'     should be provided as per the model used in modelName argument. See example
#'     or default values.
#' @param Psi Psi should be a list of length 'trueClusters' with each list
#'     element having a matrix with rows equal to value provided to
#'     'dimensionality' argument and columns equal to value provided to
#'     'dimensionality' argument. The values should be provided as per the
#'     model used in modelName argument. See example or default values.
#'
#' @return Returns an S3 object of class mplnFADataGenerator with results that
#'    equal the length of numDatasets argument. For each dataset generated, the
#'    following values are provided.
#' \itemize{
#'   \item input - User provided parameters and values used to generate the dataset.
#'   \item dataset - Simulated dataset of counts that has dimensions nObservations x
#'      dimensionality.
#'   \item Xmatrix - The matrix of X values that has dimensions nObservations x
#'      dimensionality.
#'   \item zmatrix - The matrix of membership indicator variables indicating which
#'      cluster/component each observation belongs to.
#'   \item trueMembership - A vector of length nObservations giving the true
#'      membership of all observations.
#'   \item Umatrix - The matrix containing latent factors for each observation. This
#'      has dimensions nObservations x pfactors.
#'}
#'
#' @examples
#'
#' # Example 1: Generate 10 datasets from CCC model
#' # Here, Lambda (loading matrix) and Psi (error variance and
#' # isotropic) are all constrained and hence CCC
#'
#' set.seed(100)
#' numDatasets <- 10 # total number of datasets to be generated
#' pfactors <- 3 # number of true latent factors
#' dimensionality <- 10 # dimensionality of observed data
#' trueClusters <- 2 # number of groups/clusters
#' mixingProportions <- c(0.32, 0.68) # mixing proportions for 2 clusters
#' nObservations <- 1000 # sample size or number of observations
#'
#' # set parameter values
#' mu <- list(c(6, 3, 3, 6, 3, 6, 3, 3, 6 ,3),
#'            c(5, 3, 5, 3, 5, 5, 3, 5, 3, 5))
#'
#' Lambda <- list(matrix(runif(pfactors * dimensionality, -1, 1), nrow = dimensionality),
#'                matrix(runif(pfactors * dimensionality, -1, 1), nrow = dimensionality))
#'
#' Psi <- list(diag(dimensionality) * runif(1),
#'             diag(dimensionality) * runif(1))
#'
#' # generate datasets
#' simDataCCC <- mixMPLNFA::mplnFADataGenerator(numDatasets = numDatasets,
#'                                   nObservations = nObservations,
#'                                   dimensionality = dimensionality,
#'                                   mixingProportions = mixingProportions,
#'                                   trueClusters = trueClusters,
#'                                   pfactors = pfactors,
#'                                   modelName = "CCC",
#'                                   mu = mu,
#'                                   Lambda = Lambda,
#'                                   Psi = Psi)
#'
#' # length of results
#' length(simDataCCC)
#'
#' # names of results
#' names(simDataCCC)
#'
#' # access dataset 1
#' simDataCCC$`dataset=1`$dataset
#'
#' # access input used to generate dataset 1
#' simDataCCC$`dataset=1`$input
#'
#' # to generate Sigma values if need for dataset 1
#' Sigma <- list()
#' for (gvalue in 1:simDataCCC$`dataset=1`$input$trueClusters) {
#'  Lambda <- simDataCCC$`dataset=1`$input$Lambda[[gvalue]]
#'  Psi <- simDataCCC$`dataset=1`$input$Psi[[gvalue]]
#'  Sigma[[gvalue]] <- Lambda %*% t(Lambda) + Psi
#' }
#'
#' Sigma[[1]] # Sigma for C1
#' Sigma[[2]] # Sigma for C2
#'
#' # access dataset 2
#' simDataCCC$`dataset=2`$dataset
#'
#'
#' # Example 2
#' # Here, Lambda (loading matrix) is unconstrained and Psi
#' (error variance and isotropic) are all constrained and
#' hence UCC model is used
#'
#' set.seed(100)
#' numDatasets <- 10 # total number of datasets to be generated
#' pfactors <- 2 # number of true latent factors
#' dimensionality <- 8 # dimensionality of observed data
#' trueClusters <- 4 # number of groups/clusters
#' mixingProportions <- c(0.11, 0.43, 0.24, 0.22) # mixing proportions for 4 clusters
#' nObservations <- 1000 # sample size or number of observations
#'
#' # set parameter values
#' mu <- list(c(6, 3, 3, 6, 3, 6, 3, 3),
#'            c(5, 3, 5, 3, 5, 3, 3, 5),
#'            c(4, 2, 6, 4, 2, 6, 4, 4),
#'            c(1, 3, 5, 1, 3, 5, 3, 5))
#'
#' Lambda <- list(matrix(runif(pfactors * dimensionality, -1, 1), nrow = dimensionality),
#'                matrix(runif(pfactors * dimensionality, -1, 1), nrow = dimensionality),
#'                matrix(runif(pfactors * dimensionality, -1, 1), nrow = dimensionality),
#'                matrix(runif(pfactors * dimensionality, -1, 1), nrow = dimensionality))
#'
#' Psi <- list(diag(dimensionality) * runif(1),
#'             diag(dimensionality) * runif(1),
#'             diag(dimensionality) * runif(1),
#'             diag(dimensionality) * runif(1))
#'
#' # generate datasets
#' simDataUCC <- mixMPLNFA::mplnFADataGenerator(numDatasets = numDatasets,
#'                                   nObservations = nObservations,
#'                                   dimensionality = dimensionality,
#'                                   mixingProportions = mixingProportions,
#'                                   trueClusters = trueClusters,
#'                                   pfactors = pfactors,
#'                                   modelName = "UCC",
#'                                   mu = mu,
#'                                   Lambda = Lambda,
#'                                   Psi = Psi)
#'
#' # access dataset 1
#' simDataUCC$`dataset=1`$dataset
#'
#' # access input used to generate dataset 1
#' simDataUCC$`dataset=1`$input
#'
#'
#' # Example 3
#' # Here, Lambda (loading matrix) is unconstrained and Psi
#' (error variance and isotropic) are all unconstrained and
#' hence UUU model is used
#'
#' set.seed(100)
#' numDatasets <- 10 # total number of datasets to be generated
#' pfactors <- 4 # number of true latent factors
#' dimensionality <- 10 # dimensionality of observed data
#' trueClusters <- 3 # number of groups/clusters
#' mixingProportions <- c(0.23, 0.44, 0.33) # mixing proportions for 4 clusters
#' nObservations <- 1000 # sample size or number of observations
#'
#' # set parameter values
#' mu <- list(c(4, 6, 4, 2, 2, 4, 6, 4, 6, 2),
#'            c(5, 5, 3, 3, 7, 5, 3, 3, 7, 7),
#'            c(2, 4, 4, 7, 2, 4, 7, 2, 7, 4))
#'
#' Lambda <- list(matrix(runif(pfactors * dimensionality, -1, 1), nrow = dimensionality),
#'                matrix(runif(pfactors * dimensionality, -1, 1), nrow = dimensionality),
#'                matrix(runif(pfactors * dimensionality, -1, 1), nrow = dimensionality))
#'
#' Psi <- list(diag(dimensionality) * runif(dimensionality),
#'             diag(dimensionality) * runif(dimensionality),
#'             diag(dimensionality) * runif(dimensionality))
#'
#' # generate datasets
#' simDataUUU <- mixMPLNFA::mplnFADataGenerator(numDatasets = numDatasets,
#'                                   nObservations = nObservations,
#'                                   dimensionality = dimensionality,
#'                                   mixingProportions = mixingProportions,
#'                                   trueClusters = trueClusters,
#'                                   pfactors = pfactors,
#'                                   modelName = "UUU",
#'                                   mu = mu,
#'                                   Lambda = Lambda,
#'                                   Psi = Psi)
#'
#' # access dataset 1
#' simDataUUU$`dataset=1`$dataset
#'
#' # access input used to generate dataset 1
#' simDataUUU$`dataset=1`$input
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
#'
mplnFADataGenerator <- function(numDatasets = 10,
                                nObservations = 1000,
                                dimensionality = 10,
                                mixingProportions = c(0.32, 0.68),
                                trueClusters = 2,
                                pfactors = 3,
                                modelName = "CCC",
                                mu = list(c(6, 3, 3, 6, 3, 6, 3, 3, 6 ,3),
                                          c(5, 3, 5, 3, 5, 5, 3, 5, 3, 5)),
                                Lambda = list(matrix(runif(3 * 10, -1, 1), nrow = 10),
                                              matrix(runif(3 * 10, -1, 1), nrow = 10)),
                                Psi = list(diag(10) * runif(1),
                                           diag(10) * runif(1))) {

  # Checking user input
  if(is.numeric(nObservations) != TRUE) {
    stop("nObservations should be of class numeric.")
  }

  if(is.numeric(dimensionality) != TRUE) {
    stop("dimensionality should be of class numeric.")
  }

  if(is.numeric(pfactors) != TRUE) {
    stop("pfactors should be of class numeric.")
  }

  if(dimensionality <= pfactors ) {
    stop("Dimensionality cannot be less than or equal to latent factors.
         Latent factors should be less than dimensionality.")
  }


  if(is.numeric(trueClusters) != TRUE) {
    stop("trueClusters should be of class numeric.")
  }

  if(trueClusters > nObservations) {
    stop("trueClusters argument value cannot be greater than the
         number of observations.")
  }

  if(is.numeric(numDatasets) != TRUE) {
    stop("numDatasets should be of class numeric.")
  }

  if(is.numeric(mixingProportions) != TRUE) {
    stop("mixingProportions should be a vector of class numeric.")
  }

  if (sum(mixingProportions) != 1) {
    stop("mixingProportions should be a vector that sum to 1.")
  }

  if (length(mixingProportions) != trueClusters) {
    stop("mixingProportions should be a vector that has the length trueClusters.")
  }

  if (is.character(modelName) != TRUE) {
    stop("modelName should be a vector of class character.")
  }

  if(length(modelName) != 1L){
    stop("Only one model name can be used for modelName argument.")
  }

  if(! modelName %in% c("UUU", "UUC", "UCU", "UCC", "CUU", "CUC", "CCU", "CCC")) {
    stop("Model name can only be one of UUU, UUC, UCU, UCC, CUU, CUC, CCU, or CCC")
  }

  if (length(mu[[1]]) != dimensionality) {
    stop("mu should be a list of length equal to the value provided to
         'trueClusters' argument and each element should have length equal
         to the value provided to 'dimensionality' argument.")
  }

  if (length(mu) != trueClusters) {
    stop("mu should be a list of length equal to the value provided to
         'trueClusters' argument.")
  }

  if (ncol(Lambda[[1]]) != pfactors) {
    stop("Lambda should be a list of length 'trueClusters' with each
          list element having a matrix with rows equal to value provided to
          'dimensionality' argument and columns equal to value provided to
          'pfactors' argument.")
  }

  if (nrow(Lambda[[1]]) != dimensionality) {
    stop("Lambda should be a list of length 'trueClusters' with each
          list element having a matrix with rows equal to value provided to
          'dimensionality' argument and columns equal to value provided to
          'pfactors' argument.")
    }

  if (length(Lambda) != trueClusters) {
    stop("Lambda should be a list of length 'trueClusters'.")
  }

  if (ncol(Psi[[1]]) != dimensionality) {
    stop("Psi should be a list of length 'trueClusters' with each
          list element having a matrix with rows equal to value provided to
          'dimensionality' argument and columns equal to value provided to
          'dimensionality' argument.")
  }

  if (nrow(Psi[[1]]) != dimensionality) {
    stop("Psi should be a list of length 'trueClusters' with each
          list element having a matrix with rows equal to value provided to
          'dimensionality' argument and columns equal to value provided to
          'dimensionality' argument.")
  }

  if (length(Psi) != trueClusters) {
    stop("Psi should be a list of length 'trueClusters'.")
  }

  # Stores information about all runs
  dat <- list()

  for (run in 1:numDatasets) {
    set.seed(run) # ensure seed alters for each dataset

    Ymatrix <- Xmatrix <- matrix(0, ncol = dimensionality, nrow = nObservations)
    Umatrix <- mvtnorm::rmvnorm(n = nObservations,
                                mean = rep(0, pfactors),
                                sigma = diag(pfactors))
    zmatrix <- t(stats::rmultinom(n = nObservations, size = 1,
                                  prob = mixingProportions))

    for (i in 1:nObservations) {
      grp <- which(zmatrix[i, ] == 1)
      Xmatrix[i,] <- mvtnorm::rmvnorm(n = 1,
                                      mean = mu[[grp]] + Lambda[[grp]] %*% Umatrix[i, ],
                                      sigma = Psi[[grp]])
      for (j in 1:dimensionality) {
        Ymatrix[i, j] <- stats::rpois(n = 1, lambda = exp(Xmatrix[i, j]))
      }
    }


    dat[[run]] <- list()
    dat[[run]][[1]] <- list(mu = mu,
                            Lambda = Lambda,
                            Psi = Psi,
                            trueClusters = trueClusters,
                            pfactors = pfactors,
                            nObservations = nObservations,
                            dimensionality= dimensionality,
                            mixingProportions = mixingProportions)
    dat[[run]][[2]] <- Ymatrix
    dat[[run]][[3]] <- Xmatrix
    dat[[run]][[4]] <- zmatrix
    dat[[run]][[5]] <- mclust::map(zmatrix)
    dat[[run]][[6]] <- Umatrix
    names(dat[[run]]) <- c("input", "dataset", "Xmatrix", "zmatrix", "trueMembership", "Umatrix")
  }
  names(dat) <- paste0(rep("dataset=", length(seq(1, numDatasets, 1))), seq(1, numDatasets, 1))

  class(dat) <- "mplnFADataGenerator"
  return(dat)
}

# [END]
