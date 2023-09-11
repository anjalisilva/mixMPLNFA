#' Clustering via Parsimonious Mixtures of MPLN Factor Analyzers Family
#'
#' Performs simultaneous clustering and factor analysis using
#' parsimonious mixtures of multivariate Poisson-log
#' normal factor analyzers family (PMPLNFA) via
#' variational Gaussian approximations. Model selection can
#' be done using AIC, AIC3, BIC and ICL. In the PMPLNFA framework
#' restrictions are introduced to the model parameters with the
#' aim of obtaining parsimonious models, which are sufficiently
#' flexible for clustering purposes. Since the largest contribution
#' of free parameters is through the covariance matrices, it is
#' the focus for introduction of parsimony here. The constraints can
#' be imposed on Lambda (loading matrix) and Psi (error variance and
#' isotropic) which are used to generate Sigma as per the factor
#' analysis model by Spearman, 1904 and mixture of factor analyzers
#' model  by Ghahramani et al., 1996. This function
#' simultaneously performs factor analysis and cluster analysis, by
#' assuming that the discrete observed data (counts) have been generated
#' by a factor analyzer model with continuous latent variables.
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
#'    leave as "none". Default is "none".
#' @param gmin A positive integer specifying the minimum number of components
#'    to be considered in the clustering run. Default is 1.
#' @param gmax A positive integer, >= gmin, specifying the maximum number of
#'    components to be considered in the clustering run. Default is 3.
#' @param pmin A positive integer specifying the minimum number of latent
#'    factors to be considered in the clustering run. Default is 1.
#' @param pmax A positive integer, >= pmin, specifying the maximum number of
#'    latent factors to be considered in the clustering run. Default is 2.
#' @param initMethod A character vector indicating the initialization method
#'    to be used. Default is "kmeans". Options are currently "kmeans" only.
#' @param modelNames A character string indicating the model names to test for
#'     covariance matrix, Sigma. Since the largest contribution of free parameters
#'     is through the covariance matrices, it is the focus for introduction
#'     of parsimony. The constraints can be imposed on Lambda (loading matrix)
#'     and Psi (error variance and isotropic) which are used to generate Sigma.
#'     The 'C' stands for constrained and 'U' stands for unconstrained. The order
#'     goes as loading matrix (Lambda), error variance (Psi) and isotropic (Psi),
#'     respectively. Example, if the loading matrix (Lambda), error variance (Psi)
#'     and isotropic are all constrained, then select "CCC". Options are "CCC",
#'     "UUU", and "UCC". Default is "CCC".
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
#'   \item pmin - Minimum number of latent factors (p) considered in the clustering
#'      run.
#'   \item pmax - Maximum number of latent factors (p) considered in the clustering
#'      run.
#'   \item allResults - A list with all results. To access results use the format
#'      of 'g, p, model', respectively. See examples.
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
#' # Example 1: Cluster a UCC datset
#' # Here, Lambda (loading matrix) is unconstrained and Psi
#' # (error variance and isotropic) are all constrained and
#' # hence UCC model is used
#'
#' set.seed(100)
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
#' simDataUCC <- mixMPLNFA::mplnFADataGenerator(numDatasets = 1,
#'                                              nObservations = nObservations,
#'                                              dimensionality = dimensionality,
#'                                              mixingProportions = mixingProportions,
#'                                              trueClusters = trueClusters,
#'                                              pfactors = pfactors,
#'                                              modelName = "UCC",
#'                                              mu = mu,
#'                                              Lambda = Lambda,
#'                                              Psi = Psi)
#'
#' # Clustering
#' MPLNFAEx1 <- mixMPLNFA::MPLNFAClust(
#'   dataset = simDataUCC$`dataset=1`$dataset,
#'   membership = simDataUCC$`dataset=1`$trueMembership,
#'   gmin = 3,
#'   gmax = 5,
#'   pmin = 1,
#'   pmax = 3,
#'   modelNames = c("UUU", "UUC", "UCU", "UCC", "CUU", "CUC", "CCU", "CCC"),
#'   normalize = "Yes")
#'
#' # To see BIC results
#' MPLNFAEx1$BICresults
#' MPLNFAEx1$BICresults$BICmodelselected # "G=4,p=2,model=UCC"
#'
#' # Compare with true labels
#' table(MPLNFAEx1$BICresults$BICmodelSelectedLabels,
#'       simDataUCC$`dataset=1`$trueMembership)
#'  #     1   2   3   4
#'  # 1   0   0   1 197
#'  # 2   0   0 222   0
#'  # 3   0 466   0   0
#'  # 4 114   0   0   0
#'
#'
#' # Access all results for g = 4, p = 2, model = "UCC"
#' # UCC is mentioned in fourth place for input string of modelNames argument
#' MPLNFAEx1$allResults[[2]][[2]][[4]]
#'
#'
#' # Example 2
#' # First generate a dataset from CCC model
#' # Here, Lambda (loading matrix) and Psi (error variance and
#' # isotropic) are all constrained and hence CCC
#'
#' set.seed(100)
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
#' # generate  a dataset
#' simDataCCC <- mixMPLNFA::mplnFADataGenerator(numDatasets = 2,
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
#' dim(simDataCCC$`dataset=2`$dataset) # a dataset of size 1000 by 10
#'
#' # Clustering
#' MPLNFAEx2 <- mixMPLNFA::MPLNFAClust(
#'                      dataset = simDataCCC$`dataset=1`$dataset,
#'                      membership = simDataCCC$`dataset=1`$trueMembership,
#'                      gmin = 1,
#'                      gmax = 3,
#'                      pmin = 2,
#'                      pmax = 4,
#'                      modelNames = c("UUU", "UUC", "UCU", "UCC", "CUU", "CUC", "CCU", "CCC"),
#'                      normalize = "No")
#'
#' names(MPLNFAEx2) # see all names of outputs
#'
#' # To see BIC results
#' MPLNFAEx2$BICresults
#' MPLNFAEx2$BICresults$BICmodelselected # "G=2,p=3,model=UUC"
#'
#' # Compare with true labels
#' table(MPLNFAEx2$BICresults$BICmodelSelectedLabels,
#'       simDataCCC$`dataset=1`$trueMembership)
#' #     1   2
#' # 1   0 659
#' # 2 341   0
#'
#' # Access all results for g = 2, p = 2, model = "CCC"
#' # CCC is mentioned in last place (8th) for input string of modelNames argument
#' MPLNFAEx2$allResults[[2]][[1]][[8]]
#'
#'
#' # Example 3
#' # Here, Lambda (loading matrix) is unconstrained and Psi
#' (error variance and isotropic) are all unconstrained and
#' hence UUU model is used
#'
#' set.seed(100)
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
#' simDataUUU <- mixMPLNFA::mplnFADataGenerator(numDatasets = 1,
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
#' # Clustering
#' MPLNFAEx3 <- mixMPLNFA::MPLNFAClust(
#'                      dataset = simDataUUU$`dataset=1`$dataset,
#'                      membership = simDataUUU$`dataset=1`$trueMembership,
#'                      gmin = 2,
#'                      gmax = 4,
#'                      pmin = 3,
#'                      pmax = 5,
#'                      modelNames = c("UUU", "UUC", "UCU", "UCC", "CUU", "CUC", "CCU", "CCC"),
#'                      normalize = "No")
#'
#' names(MPLNFAEx3) # see all names of outputs
#'
#' # To see BIC results
#' MPLNFAEx3$BICresults
#' MPLNFAEx3$BICresults$BICmodelselected # "G=3,p=4,model=UUU"
#'
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
#' Silva, A., Qin, X., S. J. Rothstein, P. D. McNicholas, and S. Subedi (2023).
#' Finite mixtures of matrix variate Poisson-log normal distributions
#' for three-way count data, \emph{Bioinformatics} 39(5).
#' \href{https://doi.org/10.1093/bioinformatics/btad167}{Link}
#'
#' Subedi, S., R.P. Browne (2020). A family of parsimonious mixtures of
#' multivariate Poisson-lognormal distributions for clustering multivariate
#'  count data. \emph{Stat} 9:e310. \href{https://doi.org/10.1002/sta4.310}{Link}
#'
#' @export
#' @importFrom edgeR calcNormFactors
#' @importFrom mclust unmap
#' @importFrom mclust map
#' @importFrom mvtnorm rmvnorm
#' @importFrom utils tail
#'
MPLNFAClust <- function(dataset,
                    membership = "none",
                    gmin = 1,
                    gmax = 3,
                    pmin = 1,
                    pmax = 2,
                    initMethod = "kmeans",
                    modelNames = "CCC",
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

  if(all(membership != "none") && is.numeric(membership) != TRUE) {
    stop("membership should be a numeric vector containing the
      cluster membership. Otherwise, leave as 'none'.")
  }

  if(all(membership != "none") &&
     all((diff(sort(unique(membership))) == 1) != TRUE) ) {
    stop("Cluster memberships in the membership vector
      are missing a cluster, e.g. 1, 3, 4, 5, 6 is missing cluster 2.")
  }

  if(all(membership != "none") && length(membership) != nObservations) {
    stop("membership should be a numeric vector, where length(membership)
      should equal the number of observations. Otherwise, leave as 'none'.")
  }

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

  if(! all(modelNames %in% c("UUU", "UUC", "UCU", "UCC", "CUU", "CUC", "CCU", "CCC"))) {
    stop("Model name can only be of UUU, UUC, UCU, UCC, CUU, CUC, CCU, or CCC")
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
    normFactors <- as.vector(edgeR::calcNormFactors(as.matrix(dataset),
                                                        method = "TMM"))
    normFactors <- matrix(normFactors, nrow = nObservations,
                          ncol = dimensionality, byrow = T)
  } else if(normalize == "No") {
    normFactors <- matrix(1, nrow = nObservations, ncol = dimensionality, byrow = T)
  } else{
    stop("normalize should be 'Yes' or 'No' ")
  }

  # to save cluster output
  clusterResults <- cluslabels <- list()
  BIC <- ICL <- AIC <- AIC3 <-
    nParameters <- logLikelihood <- nameVector <- vector()

  # for saving model selection values
  numberVec <- (1:((gmax - gmin + 1) * (pmax - pmin + 1) * length(modelNames)))
  numberArray <- array(numberVec,
                  dim = c((gmax - gmin + 1), (pmax - pmin + 1), length(modelNames)))
  rownames(numberArray) <- paste0(rep("G=", length(seq(gmin, gmax, 1))), seq(gmin, gmax, 1))
  colnames(numberArray) <- paste0(rep("p=", length(seq(pmin, pmax, 1))), seq(pmin, pmax, 1))

  for (gmodel in 1:(gmax - gmin + 1)) {
    # iterating through p
    clusterResults[[gmodel]] <- vector("list", length = (pmax - pmin + 1))

    for (pmodel in 1:(pmax - pmin + 1)) {

      for (famodel in seq_along(modelNames)) {

        # cat("\n gmodel", gmodel)
        # cat("\n pmodel", pmodel)
        # cat("\n famodel", famodel)

        if(length(1:(gmax - gmin + 1)) == gmax) {
          clustersize <- gmodel
        } else if(length(1:(gmax - gmin + 1)) < gmax) {
          clustersize <- seq(gmin, gmax, 1)[gmodel]
        }


        if(length(1:(pmax - pmin + 1)) == pmax) {
          pSize <- pmodel
        } else if(length(1:(pmax - pmin + 1)) < pmax) {
          pSize <- seq(pmin, pmax, 1)[pmodel]
        }

        # print statement to user
        cat("\n Running for g =", clustersize, "q =",
            pSize, "and model =", modelNames[famodel])

        clusterResults[[gmodel]][[pmodel]][[famodel]] <-
                                               PMPLNFAind(
                                               dataset = dataset,
                                               clustersize = clustersize,
                                               pSize = pSize,
                                               modelName = modelNames[famodel],
                                               normFactors = normFactors)

        inserNum <- numberArray[gmodel, pmodel, famodel]
        BIC[inserNum] <- clusterResults[[gmodel]][[pmodel]][[famodel]]$BIC
        cluslabels[[inserNum]] <- clusterResults[[gmodel]][[pmodel]][[famodel]]$clusterlabels
        ICL[inserNum] <- clusterResults[[gmodel]][[pmodel]][[famodel]]$ICL
        AIC[inserNum] <- clusterResults[[gmodel]][[pmodel]][[famodel]]$AIC
        AIC3[inserNum] <- clusterResults[[gmodel]][[pmodel]][[famodel]]$AIC3
        logLikelihood[inserNum] <- unlist(tail(
          clusterResults[[gmodel]][[pmodel]][[famodel]]$loglik, n = 1))
        nParameters[inserNum] <- clusterResults[[gmodel]][[pmodel]][[famodel]]$kTotal
        nameVector[inserNum] <- paste0("G=", clustersize,",p=", pSize, ",model=", modelNames[famodel])
      }

      clusterResults[[gmodel]][[pmodel]][[famodel]]
    }
    # names(clusterResults[[gmodel]]) <- paste0(rep("p=", length(seq(pmin, pmax, 1))),
    #                                          seq(pmin, pmax, 1))
  }



  # Add g level names to clusterResults
  # names(clusterResults) <- paste0(rep("G=", length(seq(gmin, gmax, 1))),
  #                                 seq(gmin, gmax, 1))

  names(logLikelihood) <- names(nParameters) <- names(BIC) <-
    names(ICL) <- names(AIC) <- names(AIC3) <-
    names(cluslabels) <- nameVector

  # select best model
  BICbestmodel <- names(BIC)[grep(max(BIC, na.rm = TRUE), BIC)]
  BICmodel <- list(allBICvalues = BIC,
                   BICmodelselected = BICbestmodel,
                   BICmodelSelectedLabels = cluslabels[[grep(max(BIC, na.rm = TRUE), BIC)]])

  ICLbestmodel <- names(ICL)[grep(max(ICL, na.rm = TRUE), ICL)]
  ICLmodel <- list(allICLvalues = ICL,
                   ICLmodelselected = ICLbestmodel,
                   ICLmodelSelectedLabels = cluslabels[[grep(max(ICL, na.rm = TRUE), ICL)]])

  AICbestmodel <- names(AIC)[grep(max(AIC, na.rm = TRUE), AIC)]
  AICmodel <- list(allAICvalues = AIC,
                   AICmodelselected = AICbestmodel,
                   AICmodelSelectedLabels = cluslabels[[grep(max(AIC, na.rm = TRUE), AIC)]])

  AIC3bestmodel <- names(AIC3)[grep(max(AIC3, na.rm = TRUE), AIC3)]
  AIC3model <- list(allAIC3values = AIC3,
                   AIC3modelselected = AIC3bestmodel,
                   AIC3modelSelectedLabels = cluslabels[[grep(max(AIC3, na.rm = TRUE), AIC3)]])

  final <- proc.time() - ptm

  RESULTS <- list(dataset = dataset,
                  dimensionality = dimensionality,
                  normalizationFactors = normFactors[1, ],
                  gmin = gmin,
                  gmax = gmax,
                  pmin = pmin,
                  pmax = pmax,
                  modelNames = modelNames,
                  initalizationMethod = "kmeans",
                  allResults = clusterResults,
                  logLikelihood = logLikelihood,
                  numbParameters = nParameters,
                  trueLabels = membership,
                  ICLresults = ICLmodel,
                  BICresults = BICmodel,
                  AICresults = AICmodel,
                  AIC3results = AIC3model,
                  totalTime = final)

    class(RESULTS) <- "PMPLNFA"
    return(RESULTS)
  }




#' @importFrom mclust unmap
#' @importFrom mclust map
PMPLNFAind <- function(dataset,
                       clustersize,
                       pSize,
                       modelName,
                       normFactors) {


    # Initialize variables
    dimensionality <- ncol(dataset)
    nObservations <- nrow(dataset)

    mu <- psi<- lambda <- sigmaVar <- isigma <- list()
    m <- S <- P <- Q <- lambdanew <- list()

    # Other intermediate terms
    Sk <- array(0, c(dimensionality, dimensionality, clustersize))
    start <- GX <- dGX <- zS <- list()

    # Normalization factors
    # normFactors <- edgeR::calcNormFactors(t(dataset), method = "TMM")
    # libMatFull <- libMat <- matrix(normFactors, nrow = nObservations, ncol = dimensionality, byrow = F)
    # libMatFull <- libMat <- matrix(1, nrow = nObservations, ncol = dimensionality, byrow = T)
    libMatFull <- libMat <- normFactors

    kmeansOut <- stats::kmeans(x = log(dataset + 1),
                               centers = clustersize,
                               nstart = 100)$cluster

    # Initialize z
    z <- mclust::unmap(classification = kmeansOut)
    # View(z)
    piG <- colSums(z) / nObservations
    ng <- colSums(z)

    # Initialize parameters
    for (g in 1:clustersize) {

      obs <- which(z[, g] == 1)
        if(length(obs) > 1) {
          mu[[g]] <- colMeans(log(dataset[obs, ] + 1 / dimensionality))
          sigmaVar[[g]] <- stats::var(log(dataset[obs, ] + 1 / dimensionality))

          # solve issue anticipation
          # isigma[[g]]<-solve(sigma.var[[g]])
          # SD comment: tried replacing this with Woodbury and moving it to line
          # 71 but still gives singularity error
          isigma[[g]] <- tryCatch(solve(sigmaVar[[g]]), error = function(err) NA)
          if(all(is.na(isigma[[g]]))) {
            isigma[[g]] <- diag(ncol(dataset[obs, ])) # if error with inverse
          }

        } else if (length(obs) == 1) { # if only one observation in a cluster
          mu[[g]] <- log(dataset[obs, ] + 1 / dimensionality)
          sigmaVar[[g]] <- diag(ncol(dataset))

          # solve issue anticipation
          isigma[[g]] <- tryCatch(solve(sigmaVar[[g]]), error = function(err) NA)
          if(all(is.na(isigma[[g]]))) {
            isigma[[g]] <- diag(ncol(dataset))
          }
        }


      temp <- eigen(sigmaVar[[g]])
      lambda[[g]] <- matrix(NA, ncol = pSize, nrow = dimensionality)
      for (q in 1:pSize) {
        lambda[[g]][, q] <- temp$vectors[, q] * sqrt(temp$values[q])
      }
      psi[[g]] <- diag(sigmaVar[[g]] -
                                   lambda[[g]] %*%
                                   t(lambda[[g]])) * diag(dimensionality)
    }

    for (g in 1:clustersize) {
      start[[g]] <- log(dataset + 1) ###Starting value for M
      m[[g]] <- log(dataset + 1)
      S[[g]] <- list()
      for (i in 1:nObservations) {
        S[[g]][[i]] <- diag(dimensionality) * 0.000000001
      }
    }

    it <- 1
    aloglik <- loglik <- NULL
    checks <- aloglik[c(1, 2, 3)] <- 0
    itMax <- 100

    # Begin clustering
    while (checks == 0) {
      # cat("\n it = ", it)

      for (g in 1:clustersize) {
        GX[[g]] <- dGX[[g]] <- zS[[g]] <- list()
        z[is.nan(z)] <- 0
        for (i in 1:nObservations) {
          #print(i)
          dGX[[g]][[i]] <- diag(exp(log(libMat[i, ]) + start[[g]][i, ] +
                                      0.5 * diag(S[[g]][[i]])), dimensionality) + isigma[[g]]

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
        muMatrix <- matrix(rep(mu[[g]], nObservations), nrow = nObservations, byrow = TRUE)
        res <- m[[g]] - muMatrix
        temp <- stats::cov.wt(x = res,
                              wt = z[, g],
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
        # cat("rep = ", rep, "\n")
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
                                dimensionality = dimensionality,
                                nObservations = nObservations)


        sigmaVar <- updates$sigmaVar
        psi <- updates$psi
        lambda <- updates$lambda
        paras <- updates$paras
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


      piG <- colSums(z) / nObservations
      ng <- colSums(z)
      # libMat<-matrix(normFactors,ncol=dimensionality,nrow=nObservations, byrow=T) ###Matrix containing normalization factor
      ### Some useful functions
      funFive <- function(x, y, g) {
        ySpecific = y[[g]]
        funFiveReturn <- sum(diag(x %*% ySpecific))
        return(funFiveReturn)
      }

      Ffunction <- matrix(NA, ncol = clustersize, nrow = nObservations)

      for (g in 1:clustersize) {
        two <- rowSums(exp(m[[g]] + log(libMatFull) + 0.5 *
                             matrix(unlist(lapply(S[[g]], diag)), ncol = dimensionality, byrow = TRUE)))
        five <- 0.5 * unlist(lapply(S[[g]], funFive, y = isigma, g = g))

        # detecting if six gives -Inf which happens when
        # dimensionality is large, e.g., 500
        six <- 0.5 * log(unlist(lapply(S[[g]], det)))
        if(any(is.infinite(six)) == TRUE) { # detect if -Inf
          six <- 0.5 * log(unlist(lapply(S[[g]], det)) + 1e-10)
        }

        Ffunction[, g] <- piG[g] * exp(rowSums(m[[g]] * dataset) -
                                       two -
                                       rowSums(lfactorial(dataset)) +
                                       rowSums(log(libMatFull) * dataset) -
                                       0.5 * stats::mahalanobis(m[[g]], center = mu[[g]], cov = isigma[[g]], inverted = TRUE) -
                                       five +
                                       six -
                                       0.5 * log(det(sigmaVar[[g]])) -
                                       dimensionality / 2)
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
        checks <- 1
      }
    }
    # plot(loglik,type="l")


    # paras from covariance only has the covariance parameters so now we need
    # to add the parameters for the mean and pi
    kTotal <- paras + (clustersize - 1) + (clustersize * dimensionality)

    # Inf criteria
    AIC <- 2 * loglik[it - 1] - (2 * kTotal)
    AIC3 <- 2 * loglik[it - 1] - (3 * kTotal)
    BIC <- 2 * loglik[it - 1] - (kTotal * log(nObservations))
    mapz <- matrix(0, ncol = clustersize, nrow = nObservations)
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

    modelList <- list(
      piG = piG,
      mu = mu,
      sigmaVar = sigmaVar,
      lambda = lambda,
      psi = psi,
      z = z,
      loglik = loglik,
      kmeansOut = kmeansOut,
      BIC = BIC,
      ICL = ICL,
      AIC = AIC,
      AIC3 = AIC3,
      modelName = modelName,
      clustersize = clustersize,
      pSize = pSize,
      kTotal = kTotal,
      clusterlabels = mclust::map(z))

    class(modelList) <- "mixMPLNFA"
    return(modelList)

  }

# [END]
