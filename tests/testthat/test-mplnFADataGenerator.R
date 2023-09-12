# "Checking Data Generation Via mplnFADataGenerator"
# library(mixMPLNFA)

test_that("Checking data generation results", {

  set.seed(100)
  numDatasets <- 10 # total number of datasets to be generated
  kfactors <- 4 # number of true latent factors
  dimensionality <- 10 # dimensionality of observed data
  trueClusters <- 3 # number of groups/clusters
  mixingProportions <- c(0.23, 0.44, 0.33) # mixing proportions for 4 clusters
  nObservations <- 1000 # sample size or number of observations

  # set parameter values
  mu <- list(c(4, 6, 4, 2, 2, 4, 6, 4, 6, 2),
             c(5, 5, 3, 3, 7, 5, 3, 3, 7, 7),
             c(2, 4, 4, 7, 2, 4, 7, 2, 7, 4))

  Lambda <- list(matrix(runif(kfactors * dimensionality, -1, 1), nrow = dimensionality),
                 matrix(runif(kfactors * dimensionality, -1, 1), nrow = dimensionality),
                 matrix(runif(kfactors * dimensionality, -1, 1), nrow = dimensionality))

  Psi <- list(diag(dimensionality) * runif(dimensionality),
              diag(dimensionality) * runif(dimensionality),
              diag(dimensionality) * runif(dimensionality))

  # generate datasets
  simDataUUU <- mixMPLNFA::mplnFADataGenerator(numDatasets = numDatasets,
                                               nObservations = nObservations,
                                               dimensionality = dimensionality,
                                               mixingProportions = mixingProportions,
                                               trueClusters = trueClusters,
                                               kfactors = kfactors,
                                               modelName = "UUU",
                                               mu = mu,
                                               Lambda = Lambda,
                                               Psi = Psi)

  expect_type(simDataUUU, "list")
  expect_length(simDataUUU, 10)
  expect_s3_class(simDataUUU, "mplnFADataGenerator")
  expect_identical(simDataUUU$`dataset=1`$input$trueClusters, 3)
  expect_identical(simDataUUU$`dataset=1`$input$kfactors, 4)
  expect_identical(simDataUUU$`dataset=1`$input$nObservations, 1000)
  expect_identical(simDataUUU$`dataset=1`$input$dimensionality, 10)
  expect_named(simDataUUU$`dataset=1`$input, c("mu", "Lambda", "Psi", "trueClusters",
                                               "kfactors", "nObservations", "dimensionality",
                                               "mixingProportions"))
  expect_output(str(simDataUUU), "List of 10")
})

# "Checking for invalid user input"
test_that("Data generation error upon invalid user input", {

  set.seed(100)
  numDatasets <- 10 # total number of datasets to be generated
  kfactors <- 4 # number of true latent factors
  dimensionality <- 10 # dimensionality of observed data
  trueClusters <- 3 # number of groups/clusters
  mixingProportions <- c(0.23, 0.44, 0.33) # mixing proportions for 4 clusters
  nObservations <- 1000 # sample size or number of observations

  # set parameter values
  mu <- list(c(4, 6, 4, 2, 2, 4, 6, 4, 6, 2),
             c(5, 5, 3, 3, 7, 5, 3, 3, 7, 7),
             c(2, 4, 4, 7, 2, 4, 7, 2, 7, 4))

  Lambda <- list(matrix(runif(kfactors * dimensionality, -1, 1), nrow = dimensionality),
                 matrix(runif(kfactors * dimensionality, -1, 1), nrow = dimensionality),
                 matrix(runif(kfactors * dimensionality, -1, 1), nrow = dimensionality))

  Psi <- list(diag(dimensionality) * runif(dimensionality),
              diag(dimensionality) * runif(dimensionality),
              diag(dimensionality) * runif(dimensionality))


  # number of dataset provided as character
  testthat::expect_error(mixMPLNFA::mplnFADataGenerator(numDatasets = "numDatasets",
                                                        nObservations = nObservations,
                                                        dimensionality = dimensionality,
                                                        mixingProportions = mixingProportions,
                                                        trueClusters = trueClusters,
                                                        kfactors = kfactors,
                                                        modelName = "UUU",
                                                        mu = mu,
                                                        Lambda = Lambda,
                                                        Psi = Psi))

  # nObservations provided as logical
  testthat::expect_error(mixMPLNFA::mplnFADataGenerator(numDatasets = numDatasets,
                                                        nObservations = FALSE,
                                                        dimensionality = dimensionality,
                                                        mixingProportions = mixingProportions,
                                                        trueClusters = trueClusters,
                                                        kfactors = kfactors,
                                                        modelName = "UUU",
                                                        mu = mu,
                                                        Lambda = Lambda,
                                                        Psi = Psi))

  # Incorrect size for p as has to be p < d
  testthat::expect_error(mixMPLNFA::mplnFADataGenerator(numDatasets = numDatasets,
                                                        nObservations = nObservations,
                                                        dimensionality = 10,
                                                        mixingProportions = mixingProportions,
                                                        trueClusters = trueClusters,
                                                        kfactors = 11,
                                                        modelName = "UUU",
                                                        mu = mu,
                                                        Lambda = Lambda,
                                                        Psi = Psi))

  # Incorrect trueClusters as trueClusters cannot be larger than nObservations
  testthat::expect_error(mixMPLNFA::mplnFADataGenerator(numDatasets = numDatasets,
                                                        nObservations = 100,
                                                        dimensionality = dimensionality,
                                                        mixingProportions = mixingProportions,
                                                        trueClusters = 200,
                                                        kfactors = kfactors,
                                                        modelName = "UUU",
                                                        mu = mu,
                                                        Lambda = Lambda,
                                                        Psi = Psi))


  # Incorrect modelName
  testthat::expect_error(mixMPLNFA::mplnFADataGenerator(numDatasets = numDatasets,
                                                        nObservations = nObservations,
                                                        dimensionality = dimensionality,
                                                        mixingProportions = mixingProportions,
                                                        trueClusters = trueClusters,
                                                        kfactors = kfactors,
                                                        modelName = "AAA",
                                                        mu = mu,
                                                        Lambda = Lambda,
                                                        Psi = Psi))


})
# [END]
