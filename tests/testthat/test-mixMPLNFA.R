context("Checking Variational Gaussian Approximations Approach")
library(mixMPLNFA)

test_that("Checking clustering results", {

  testData1 <- dat[[2]][[2]]

  # Clustering simulated matrix variate count data
  # "UUU", "UUC", "UCU", "UCC", "CUU", "CUC", "CCU", and "CCC".
  clusResultsUUU <- mixMPLNFA::MPLNFA(dataset = testData1,
                                  gmin = 1,
                                  gmax = 3,
                                  pmin = 1,
                                  pmax = 2,
                                  modelName = c("UUU"),
                                  normalize = "Yes")

  expect_type(clusResultsUUU, "list")
  expect_length(clusResultsUUU, 17)
  expect_s3_class(clusResultsUUU, "mvplnVGA")
  expect_identical(clusResultsUUU$nUnits, 1000L)
  expect_identical(clusResultsUUU$nVariables, 3L)
  expect_identical(clusResultsUUU$nOccassions, 2L)
  expect_named(clusResultsUUU, c("dataset", "nUnits",
                                  "nVariables", "nOccassions",
                                  "normFactors", "gmin", "gmax",
                                  "initalizationMethod", "allResults",
                                  "loglikelihood", "nParameters",
                                  "trueLabels", "ICLAll",
                                  "BICAll", "AICAll",
                                  "AIC3All", "totalTime"))
  expect_output(str(clusResultsUUU), "List of 17")
  expect_vector(clusResultsUUU$trueLabels, ptype = double(), size = 1000)
  expect_identical(clusResultsUUU$BICAll$BICmodelselected, 1)
})

context("Checking for invalid user input")
test_that("Data clustering error upon invalid user input", {

  # dataset provided as character
  testthat::expect_error(mixMPLNFA::MPLNFA(dataset = "dataset",
                                            gmin = 1,
                                            gmax = 3,
                                            pmin = 1,
                                            pmax = 2,
                                            modelName = c("UUU"),
                                            normalize = "Yes"))

  # dataset provided as logical
  testthat::expect_error(mixMPLNFA::MPLNFA(dataset = FALSE,
                                            gmin = 1,
                                            gmax = 3,
                                            pmin = 1,
                                            pmax = 2,
                                            modelName = c("UUU"),
                                            normalize = "Yes"))

  # Incorrect size for p as has to be p < d
  testthat::expect_error(mixMPLNFA::MPLNFA(dataset = testData1,
                                            gmin = 1,
                                            gmax = 3,
                                            pmin = 1,
                                            pmax = ncol(testData1) + 1,
                                            modelName = c("UUU"),
                                            normalize = "Yes"))

  # Incorrect g as gmax cannot be larger than d
  testthat::expect_error(mixMPLNFA::MPLNFA(dataset = testData1,
                                            gmin = 1,
                                            gmax = ncol(testData1) + 1,
                                            pmin = 1,
                                            pmax = 2,
                                            modelName = c("UUU"),
                                            normalize = "Yes"))


  # Incorrect input type for gmin
  testthat::expect_error(mixMPLNFA::MPLNFA(dataset = testData1,
                                            gmin = "1",
                                            gmax = 3,
                                            pmin = 1,
                                            pmax = 2,
                                            modelName = c("UUU"),
                                            normalize = "Yes"))


  # Incorrect input for gmin and gmax
  testthat::expect_error(mixMPLNFA::MPLNFA(dataset = testData1,
                                            gmin = 5,
                                            gmax = 2,
                                            pmin = 1,
                                            pmax = 2,
                                            modelName = c("UUU"),
                                            normalize = "Yes"))

  # Incorrect input for pmin and pmax
  testthat::expect_error(mixMPLNFA::MPLNFA(dataset = testData1,
                                            gmin = 1,
                                            gmax = 2,
                                            pmin = 3,
                                            pmax = 1,
                                            modelName = c("UUU"),
                                            normalize = "Yes"))

  # Incorrect input type for normalize
  testthat::expect_error(mixMPLNFA::MPLNFA(dataset = testData1,
                                            gmin = 1,
                                            gmax = 2,
                                            pmin = 3,
                                            pmax = 1,
                                            modelName = c("UUU"),
                                            normalize = NA))


  # Incorrect input type for modelName
  testthat::expect_error(mixMPLNFA::MPLNFA(dataset = testData1,
                                            gmin = 1,
                                            gmax = 2,
                                            pmin = 3,
                                            pmax = 1,
                                            modelName = c("UUU"),
                                            normalize = NA))

  # Incorrect input type for membership
  testthat::expect_error(mixMPLNFA::MPLNFA(dataset = testData1,
                                            membership = "",
                                            gmin = 1,
                                            gmax = 2,
                                            pmin = 1,
                                            pmax = 2,
                                            modelName = c("UUU"),
                                            normalize = NA))
})
# [END]
