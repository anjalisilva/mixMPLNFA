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
#' trueMu1 <- c(6.5, 6, 6, 6, 6, 6, 6, 6, 6)
#' trueMu2 <- c(2, 2.5, 2, 2, 2, 2, 2, 2, 2)
#' trueMu3 <- c(1, 1, 1, 1, 1, 1, 1, 1, 1)
#'
#' trueSigma1 <- diag(length(trueMu1)) * 2
#' trueSigma2 <- diag(length(trueMu1)) * 0.05
#' trueSigma3 <- diag(length(trueMu1)) * 0.01
#'
#' # Generating simulated data with 3 clusters
#' sampleData <- MPLNClust::mplnDataGenerator(nObservations = 2000,
#'                      dimensionality = length(trueMu1),
#'                      mixingProportions = c(0.2, 0.3, 0.5),
#'                      mu = rbind(trueMu1, trueMu2, trueMu3),
#'                      sigma = rbind(trueSigma1, trueSigma2, trueSigma3),
#'                      produceImage = "No")
#'
#'  dim(sampleData$dataset) # a dataset of size 2000 by 9
#'
#' # Clustering
#' MPLNFAResults <- mixMPLNFA::MPLNFAClust(
#'                      dataset = sampleData$dataset,
#'                      membership = sampleData$trueMembership,
#'                      gmin = 1,
#'                      gmax = 4,
#'                      pmin = 1,
#'                      pmax = 1,
#'                      modelNames = c("CCU", "UUU"),
#'                      normalize = "Yes")
#'
#'  # Visualize data using line plot
#'  # Use navigation buttons to see previous plots
#'  MPLNFABlack <- mixMPLNFA::mplnVisualizeLine(dataset = sampleData$dataset,
#'                                          clusterMembershipVector =
#'                                          MPLNFAResults$BICresults$BICmodelSelectedLabels,
#'                                          fileName = 'Example1',
#'                                          printPlot = FALSE)
#'
#'
#'  # Visualize data using line plot with multicolours
#'  # Use navigation buttons to see previous plots
#'  MPLNLineColor <- mixMPLNFA::mplnVisualizeLine(dataset = sampleData$dataset,
#'                                          clusterMembershipVector =
#'                                          MPLNFAResults$BICresults$BICmodelSelectedLabels,
#'                                          fileName = 'Example1MultiColor',
#'                                          LinePlotColours = "multicolour",
#'                                          printPlot = FALSE)
#'
#'  # Example 2
#' trueMu1 <- c(6.5, 6, 6, 6, 6, 6, 6, 6, 6)
#' trueMu2 <- c(1, 2.5, 2, 2, 2, 2, 1, 1, 1)
#'
#' trueSigma1 <- diag(length(trueMu1)) * 2.5
#' trueSigma2 <- diag(length(trueMu1)) * 0.5
#'
#' # Generating simulated data with 2 clusters
#' sampleData2 <- MPLNClust::mplnDataGenerator(nObservations = 2000,
#'                      dimensionality = length(trueMu1),
#'                      mixingProportions = c(0.6, 0.4),
#'                      mu = rbind(trueMu1, trueMu2),
#'                      sigma = rbind(trueSigma1, trueSigma2),
#'                      produceImage = "No")
#'
#'  dim(sampleData2$dataset) # a dataset of size 2000 by 9
#'
#' # Clustering
#' MPLNFAResults2 <- mixMPLNFA::MPLNFAClust(
#'                      dataset = sampleData2$dataset,
#'                      membership = sampleData2$trueMembership,
#'                      gmin = 1,
#'                      gmax = 3,
#'                      pmin = 1,
#'                      pmax = 2,
#'                      modelNames = c("CCU"),
#'                      normalize = "Yes")
#'
#'  # Visualize data using line plot with multicolours
#'  # Use navigation buttons to see previous plots
#'  MPLNLineColor <- mixMPLNFA::mplnVisualizeLine(dataset = sampleData2$dataset,
#'                                          clusterMembershipVector =
#'                                          MPLNFAResults2$BICresults$BICmodelSelectedLabels,
#'                                          fileName = 'TwoClusterModel',
#'                                          LinePlotColours = "multicolour",
#'                                          printPlot = FALSE)
#'
#'
#' @author Anjali Silva, \email{anjali@alumni.uoguelph.ca}
#'
#' @export
#' @import graphics
#' @importFrom grDevices png
#' @importFrom grDevices pdf
#' @importFrom grDevices dev.off
#' @importFrom RColorBrewer brewer.pal.info
#' @importFrom RColorBrewer brewer.pal
mplnVisualizeLine <- function(dataset,
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