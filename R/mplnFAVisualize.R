#' Visualize Clustered Results Via Line Plots
#'
#' A function to visualize clustering results via line plots.
#' Each cluster will have its own plot. Data is log-transformed
#' prior to visualizing. Values for each sample are connected
#' by dashed lines to illustrate the trends (log counts). The
#' yellow line shows the mean value (log counts) for each cluster.
#'
#' @param dataset A dataset of class matrix and type integer such that
#'    rows correspond to observations and columns correspond to variables.
#' @param clusterMembershipVector A numeric vector of length nrow(dataset)
#'    containing the cluster membership of each observation. If not provided,
#'    all observations will be treated as belonging to one cluster. Default is NA.
#' @param LinePlotColours Character string indicating if the line plots
#'    should be multicoloured or monotone, in black. Options are
#'    'multicolour' or 'black'. Default is 'black'.
#' @param printPlot Logical indicating if plot(s) should be saved in local
#'    directory. Default TRUE. Options TRUE or FALSE.
#' @param fileName Unique character string indicating the name for the plot
#'    being generated. Default is Plot_date, where date is obtained from
#'    date().
#' @param format Character string indicating the format of the image to
#'    be produced. Default 'pdf'. Options 'pdf' or 'png'.
#'
#' @return Plotting function provides the possibility for line plots.
#'
#' @examples
#' # Example 1
#' # Here, Lambda (loading matrix) is unconstrained and Psi
#' # (error variance and isotropic) are all unconstrained and
#' # hence UUU model is used
#'
#' set.seed(100)
#' kfactors <- 4 # number of true latent factors
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
#' Lambda <- list(matrix(runif(kfactors * dimensionality, -1, 1), nrow = dimensionality),
#'                matrix(runif(kfactors * dimensionality, -1, 1), nrow = dimensionality),
#'                matrix(runif(kfactors * dimensionality, -1, 1), nrow = dimensionality))
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
#'                                   kfactors = kfactors,
#'                                   modelName = "UUU",
#'                                   mu = mu,
#'                                   Lambda = Lambda,
#'                                   Psi = Psi)
#'
#' # Clustering
#' MPLNFAEx1 <- mixMPLNFA::MPLNFAClust(
#'                      dataset = simDataUUU$`dataset=1`$dataset,
#'                      membership = simDataUUU$`dataset=1`$trueMembership,
#'                      gmin = 2,
#'                      gmax = 4,
#'                      kmin = 3,
#'                      kmax = 5,
#'                      modelNames = c("UUU", "UUC", "UCU", "UCC", "CUU", "CUC", "CCU", "CCC"),
#'                      normalize = "No")
#'
#'  # Visualize data using line plot
#'  # Use navigation buttons to see previous plots
#'  MPLNFABlack <- mixMPLNFA::mplnFAVisLine(dataset = simDataUUU$dataset,
#'                                          clusterMembershipVector =
#'                                          MPLNFAEx1$BICresults$BICmodelSelectedLabels,
#'                                          fileName = 'Example1',
#'                                          printPlot = FALSE)
#'
#'
#'  # Visualize data using line plot with multicolours
#'  # Use navigation buttons to see previous plots
#'  MPLNLineColor <- mixMPLNFA::mplnFAVisLine(dataset = simDataUUU$dataset,
#'                                          clusterMembershipVector =
#'                                          MPLNFAEx1$BICresults$BICmodelSelectedLabels,
#'                                          fileName = 'Example1MultiColor',
#'                                          LinePlotColours = "multicolour",
#'                                          printPlot = FALSE)
#'
#'
#'
#' @author Anjali Silva, \email{anjali@alumni.uoguelph.ca}
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
#' Ghahramani, Z., G. E. Hinton, et al. (1996). The EM algorithm for mixtures of
#' factor analyzers. Technical report, Technical Report CRG-TR-96-1, University
#' of Toronto.
#'
#' Ghahramani, Z. and Beal, M. (1999). Variational inference for bayesian
#' mixtures of factor analysers. \emph{Advances in neural information processing
#' systems} 12.
#'
#' Robinson, M.D., and Oshlack, A. (2010). A scaling normalization method for differential
#' expression analysis of RNA-seq data. \emph{Genome Biology} 11, R25.
#'
#' Schwarz, G. (1978). Estimating the dimension of a model. \emph{The Annals of Statistics}
#' 6.
#'
#' Spearman, C. (1904). The proof and measurement of association between two things.
#' \emph{The American Journal of Psychology}, 15(1).
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
#' count data. \emph{Stat} 9:e310. \href{https://doi.org/10.1002/sta4.310}{Link}
#'
#' @export
#' @import graphics
#' @importFrom grDevices png
#' @importFrom grDevices pdf
#' @importFrom grDevices dev.off
#' @importFrom RColorBrewer brewer.pal.info
#' @importFrom RColorBrewer brewer.pal
mplnFAVisLine <- function(dataset,
                              clusterMembershipVector = NA,
                              fileName = paste0('Plot_',date()),
                              LinePlotColours = "black",
                              printPlot = TRUE,
                              format = 'pdf') {

  # Checking user input
  if (typeof(dataset) != "double" & typeof(dataset) != "integer") {
    stop("\n Dataset type needs to be integer")
  }

  if (is.matrix(dataset) != TRUE) {
    stop("\n Dataset needs to be a matrix")
  }

  if (is.logical(clusterMembershipVector) == TRUE) {
    cat("\n clusterMembershipVector is not provided.")
    clusterMembershipVector <- rep(1, nrow(dataset))

  } else if (is.numeric(clusterMembershipVector) == TRUE) {
    if (nrow(dataset) != length(clusterMembershipVector)) {
      stop("\n length(clusterMembershipVector) should match
          nrow(dataset)")
    }
  }

  # Obtaining path to save images
  pathNow <- getwd()

  # Saving cluster membership for each observation
  DataPlusLabs <- cbind(dataset, clusterMembershipVector)
  ordervector <- anothervector <- list()

  # Divide observations into each cluster based on membership
  for (i in 1:max(clusterMembershipVector)) {
    ordervector[[i]] <- which(DataPlusLabs[,
                                           ncol(dataset) + 1] == i)
    # divide observations as an integer based on cluster membership
    anothervector[[i]] <- rep(i,
                              length(which(DataPlusLabs[,
                                                        ncol(dataset) + 1] == i)))
  }

  # Setting the colours
  if(max(clusterMembershipVector) > 17) {
    qualColPals <- RColorBrewer::brewer.pal.info[brewer.pal.info$category == 'qual', ]
    coloursBarPlot <- unlist(mapply(RColorBrewer::brewer.pal,
                                    qualColPals$maxcolors,
                                    rownames(qualColPals)))
  } else {
    coloursBarPlot <- c('#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6',
                        '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324',
                        '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1',
                        '#000075', '#808080')
  }


  # empty plots
  linePlots <- NULL

  # Line Plots
  if (LinePlotColours == "multicolour") {
    linePlots <- list()
    for(cluster in unique(clusterMembershipVector)) {

      # Save how many observations below to each cluster size,
      # given by 'cluster'
      if (length(which(DataPlusLabs[, ncol(dataset) + 1] == cluster)) == 1) {
        toPlot2 <- as.matrix(DataPlusLabs[which(DataPlusLabs[,
                                                             ncol(dataset) + 1] == cluster), c(1:ncol(dataset))],
                             ncol = ncol(dataset))
        rownames(toPlot2) <- names(which(DataPlusLabs[, ncol(dataset) + 1] == cluster))
      } else if (length(which(DataPlusLabs[, ncol(dataset) + 1] == cluster)) > 1) {
        toPlot2 <- as.matrix(DataPlusLabs[which(DataPlusLabs[,
                                                             ncol(dataset) + 1] == cluster), c(1:ncol(dataset))],
                             ncol = ncol(dataset))
      }

      # Save column mean in last row
      toplot1 <- rbind(log(toPlot2 + 1), colMeans(log(toPlot2 + 1)))
      # If discontinunity is needed between samples (e.g. for 6 samples)
      # toplot1_space=cbind(toplot1[,c(1:3)],rep(NA,nrow(toPlot2)+1),
      # toplot1[,c(4:6)])


      if (printPlot == TRUE) {
        if (format == 'png') {
          grDevices::png(paste0(pathNow, "/LinePlot_Cluster", cluster,
                                "_", fileName, ".png"))
        } else {
          grDevices::pdf(paste0(pathNow, "/LinePlot_Cluster", cluster,
                                "_", fileName, ".pdf"))
        }

        linePlotMultiCol(dataset = dataset,
                         toplot1 = toplot1,
                         toPlot2 = toPlot2,
                         coloursBarPlot = coloursBarPlot,
                         cluster = cluster)
        grDevices::dev.off()
      }

      linePlots[[cluster]] <- linePlotMultiCol(dataset = dataset,
                                               toplot1 = toplot1,
                                               toPlot2 = toPlot2,
                                               coloursBarPlot = coloursBarPlot,
                                               cluster = cluster)
    }
  } else if (LinePlotColours == "black") {
    linePlots <- list()
    for(cluster in unique(clusterMembershipVector)) {

      # Save how many observations below to each cluster size,
      # given by 'cluster'
      if (length(which(DataPlusLabs[, ncol(dataset) + 1] == cluster)) == 1) {
        toPlot2 <- t(as.matrix(DataPlusLabs[which(DataPlusLabs[,
                                                               ncol(dataset) + 1] == cluster), c(1:ncol(dataset))],
                               ncol = ncol(dataset)))
        rownames(toPlot2) <- names(which(DataPlusLabs[, ncol(dataset) + 1] == cluster))
      } else if (length(which(DataPlusLabs[, ncol(dataset) + 1] == cluster)) > 1) {
        toPlot2 <- as.matrix(DataPlusLabs[which(DataPlusLabs[,
                                                             ncol(dataset) + 1] == cluster), c(1:ncol(dataset))],
                             ncol = ncol(dataset))
      }

      # Save column mean in last row
      toplot1 <- rbind(log(toPlot2 + 1), colMeans(log(toPlot2 + 1)))
      # If discontinunity is needed between samples (e.g. for 6 samples)
      # toplot1_space=cbind(toplot1[,c(1:3)],rep(NA,nrow(toPlot2)+1),
      # toplot1[,c(4:6)])

      if (printPlot == TRUE) {
        if (format == 'png') {
          grDevices::png(paste0(pathNow, "/LinePlot_Cluster", cluster,
                                "_", fileName, ".png"))
        } else {
          grDevices::pdf(paste0(pathNow, "/LinePlot_Cluster", cluster,
                                "_", fileName, ".pdf"))
        }
        linePlotMonoCol(dataset = dataset,
                        toplot1 = toplot1,
                        toPlot2 = toPlot2,
                        cluster = cluster)
        grDevices::dev.off()
      }
      linePlots[[cluster]] <- linePlotMonoCol(dataset = dataset,
                                              toplot1 = toplot1,
                                              toPlot2 = toPlot2,
                                              cluster = cluster)
    }
  }
  return(linePlots)
}

# Helper function for line plot
linePlotMultiCol <- function(dataset,
                             toplot1,
                             toPlot2,
                             coloursBarPlot,
                             cluster) {
  linePlotMultiCol <- graphics::matplot(t(toplot1), type = "l", pch = 1,
                                        col = c(rep(coloursBarPlot[cluster], nrow(toPlot2)), 7),
                                        xlab = "Samples", ylab = "Expression (log counts)", cex = 1,
                                        lty = c(rep(2, nrow(toPlot2)), 1),
                                        lwd = c(rep(3, nrow(toPlot2)), 4),
                                        xaxt = "n", xlim = c(1, ncol(toplot1)),
                                        main = paste("Cluster ", cluster))
  linePlotMultiCol <- linePlotMultiCol + axis(1, at = c(1:ncol(dataset)), labels = colnames(dataset))
  return(linePlotMultiCol)
}

# Helper function for line plot
linePlotMonoCol <- function(dataset,
                            toplot1,
                            toPlot2,
                            cluster) {
  linePlotMonoCol <- graphics::matplot(t(toplot1), type = "l", pch = 1,
                                       col = c(rep(1, nrow(toPlot2)), 7),
                                       xlab = "Samples", ylab = "Expression (log counts)", cex = 1,
                                       lty = c(rep(2, nrow(toPlot2)), 1),
                                       lwd = c(rep(3, nrow(toPlot2)), 4),
                                       xaxt = "n", xlim = c(1, ncol(toplot1)),
                                       main = paste("Cluster ", cluster))
  linePlotMonoCol <- linePlotMonoCol + axis(1, at = c(1:ncol(dataset)), labels = colnames(dataset))
  return(linePlotMonoCol)
}




# [END]
