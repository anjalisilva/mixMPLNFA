# "Checking Variational Gaussian Approximations Approach"
# library(mixMPLNFA)

test_that("Checking clustering results", {

  trueMu1 <- c(6.5, 6, 6, 6, 6, 6)
  trueMu2 <- c(2, 2.5, 2, 2, 2, 2)

  trueSigma1 <- diag(6) * 2
  trueSigma2 <- diag(6)

  # Generating simulated data
  sampleData <- MPLNClust::mplnDataGenerator(nObservations = 100,
                                             dimensionality = 6,
                                             mixingProportions = c(0.79, 0.21),
                                             mu = rbind(trueMu1, trueMu2),
                                             sigma = rbind(trueSigma1, trueSigma2),
                                             produceImage = "No",
                                             ImageName = "TwoComponents")


  # Clustering simulated count data
  clusResultsUUU <- mixMPLNFA::MPLNFAClust(dataset = sampleData$dataset,
                                  gmin = 1,
                                  gmax = 3,
                                  pmin = 1,
                                  pmax = 2,
                                  modelName = "UUU",
                                  normalize = "Yes")

  expect_type(clusResultsUUU, "list")
  expect_length(clusResultsUUU, 18)
  expect_s3_class(clusResultsUUU, "PMPLNFA")
  expect_identical(clusResultsUUU$gmin, 1)
  expect_identical(clusResultsUUU$gmax, 3)
  expect_identical(clusResultsUUU$pmin, 1)
  expect_identical(clusResultsUUU$pmax, 2)
  expect_named(clusResultsUUU, c("dataset", "dimensionality", "normalizationFactors",
                                 "gmin",  "gmax", "pmin",
                                 "pmax", "modelNames", "initalizationMethod", "allResults",
                                 "logLikelihood", "numbParameters", "trueLabels",
                                 "ICLresults", "BICresults", "AICresults",
                                 "AIC3results", "totalTime"))
  expect_output(str(clusResultsUUU), "List of 17")
  expect_identical(clusResultsUUU$BICresults$BICmodelselected, "G=2,p=1,model=UUU")
})

# "Checking for invalid user input"
test_that("Data clustering error upon invalid user input", {

  # dataset provided as character
  testthat::expect_error(mixMPLNFA::MPLNFAClust(dataset = "sampleData$dataset",
                                            gmin = 1,
                                            gmax = 3,
                                            pmin = 1,
                                            pmax = 2,
                                            modelName = "UUU",
                                            normalize = "Yes"))

  # dataset provided as logical
  testthat::expect_error(mixMPLNFA::MPLNFAClust(dataset = FALSE,
                                            gmin = 1,
                                            gmax = 3,
                                            pmin = 1,
                                            pmax = 2,
                                            modelName = "UUU",
                                            normalize = "Yes"))

  # Incorrect size for p as has to be p < d
  testthat::expect_error(mixMPLNFA::MPLNFAClust(dataset = sampleData$dataset,
                                            gmin = 1,
                                            gmax = 3,
                                            pmin = 1,
                                            pmax = ncol(sampleData$dataset) + 1,
                                            modelName = "UUU",
                                            normalize = "Yes"))

  # Incorrect g as gmax cannot be larger than d
  testthat::expect_error(mixMPLNFA::MPLNFAClust(dataset = sampleData$dataset,
                                            gmin = 1,
                                            gmax = nrow(sampleData$dataset) + 1,
                                            pmin = 1,
                                            pmax = 2,
                                            modelName = "UUU",
                                            normalize = "Yes"))


  # Incorrect input type for gmin
  testthat::expect_error(mixMPLNFA::MPLNFAClust(dataset = sampleData$dataset,
                                            gmin = "1",
                                            gmax = 3,
                                            pmin = 1,
                                            pmax = 2,
                                            modelName = "UUU",
                                            normalize = "Yes"))


  # Incorrect input for gmin and gmax
  testthat::expect_error(mixMPLNFA::MPLNFAClust(dataset = sampleData$dataset,
                                            gmin = 5,
                                            gmax = 2,
                                            pmin = 1,
                                            pmax = 2,
                                            modelName = "UUU",
                                            normalize = "Yes"))

  # Incorrect input for pmin and pmax
  testthat::expect_error(mixMPLNFA::MPLNFAClust(dataset = sampleData$dataset,
                                            gmin = 1,
                                            gmax = 2,
                                            pmin = 3,
                                            pmax = 1,
                                            modelName = "UUU",
                                            normalize = "Yes"))

  # Incorrect input type for normalize
  testthat::expect_error(mixMPLNFA::MPLNFAClust(dataset = sampleData$dataset,
                                            gmin = 1,
                                            gmax = 2,
                                            pmin = 3,
                                            pmax = 1,
                                            modelName = "UUU",
                                            normalize = NA))


  # Incorrect input type for modelName
  testthat::expect_error(mixMPLNFA::MPLNFAClust(dataset = sampleData$dataset,
                                            gmin = 1,
                                            gmax = 2,
                                            pmin = 3,
                                            pmax = 1,
                                            modelName = "UUU",
                                            normalize = NA))

  # Incorrect input type for membership
  testthat::expect_error(mixMPLNFA::MPLNFAClust(dataset = sampleData$dataset,
                                            membership = "",
                                            gmin = 1,
                                            gmax = 2,
                                            pmin = 1,
                                            pmax = 2,
                                            modelName = "UUU",
                                            normalize = NA))
})
# [END]
