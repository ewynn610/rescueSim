blankParams <- RescueParams()
fullParams <- RescueParams(
  logLibFacVar = 1, logLibMean = 1, logLibSD = 1,
  exprsMean = c(1, 1, 1), dispersion = c(1, 1, 1),
  sampleFacVarMean = 1, sampleFacVarSD = 1,
  subjectFacVarMean = 1, subjectFacVarSD = 1,
  nSubjsPerGroup = 1, twoGroupDesign = F,
  nTimepoints = 1, maxCellsPerSamp = 1, minCellsPerSamp = 1,
  propDE = 0, deLogFC = 0
)

test_that("getRescueParam working", {
  expect_equal(
    lapply(slotNames(fullParams), getRescueParam, paramObj = fullParams),
    list(1, FALSE, 1, 1, 1, 1, 1, 1, c(1, 1, 1), c(1, 1, 1), 1, 1, 1, 1, 0, 0)
  )
})

test_that("updateRescueParams working", {
  paramLS <- list(
    logLibFacVar = 1, logLibMean = 1, logLibSD = 1,
    exprsMean = c(1, 1, 1), dispersion = c(1, 1, 1),
    sampleFacVarMean = 1, sampleFacVarSD = 1,
    subjectFacVarMean = 1, subjectFacVarSD = 1,
    nSubjsPerGroup = 1, twoGroupDesign = F,
    nTimepoints = 1, maxCellsPerSamp = 1, minCellsPerSamp = 1,
    propDE = 0, deLogFC = 0
  )
  expect_equal(updateRescueParams(blankParams, paramLS), fullParams)
})

test_that(".checkSingleValueParams working", {
  singleValueSlots <- c(
    "nTimepoints", "nSubjsPerGroup",
    "logLibFacVar", "logLibMean", "logLibSD",
    "sampleFacVarMean",
    "sampleFacVarSD",
    "subjectFacVarMean",
    "subjectFacVarSD",
    "propDE",
    "twoGroupDesign"
  )
  singleValList <- c(rep(list(c(1, 1)), 10), rep(list(c(F, F)), 1))
  names(singleValList) <- singleValueSlots


  lapply(singleValueSlots, function(x) {
    expect_error(
      updateRescueParams(fullParams, singleValList),
      paste("Parameter", x, "should contain a single value")
    )
  })
})

test_that(".checkExprsDispLength working", {
  expect_error(
    updateRescueParams(
      fullParams,
      list(
        exprsMean = c(1, 1),
        dispersion = c(1, 1, 1, 1)
      )
    ),
    "Parameters exprsMean and dispersion must be vectors of equal length"
  )
})


test_that(".checkCellsPerSampLength working", {
  wrongParams <- list(
    list(
      nSubjsPerGroup = 2,
      twoGroupDesign = T,
      nTimepoints = c(2),
      maxCellsPerSamp = c(1, 1, 1),
      minCellsPerSamp = c(1, 1, 1)
    )
  )

  rightParams <- list(
    list(
      nSubjsPerGroup = 2,
      twoGroupDesign = T,
      nTimepoints = c(2),
      maxCellsPerSamp = rep(1, 8),
      minCellsPerSamp = rep(1, 8)
    )
  )
  lapply(wrongParams, function(ls) {
    lapply(c("maxCellsPerSamp", "minCellsPerSamp"), function(x) {
      expect_error(
        updateRescueParams(fullParams, ls),
        paste0(
          "Parameter ", x,
          " should contain a single value (",
          x, " for all samples), a vector with length equal to the number of conditions (group/timepoint combinations),or a vector with length equal to the number of total samples."
        ),
        fixed = T
      )
    })
  })

  lapply(rightParams, function(ls) {
    lapply(c("maxCellsPerSamp", "minCellsPerSamp"), function(x) {
      expect_no_error(
        updateRescueParams(fullParams, ls)
      )
    })
  })
})
