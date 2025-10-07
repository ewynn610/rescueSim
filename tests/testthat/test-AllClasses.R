blankParams <- RescueSimParams()
fullParams <- RescueSimParams(
  logLibFacVar = 1, logLibMean = 1, logLibSD = 1,
  exprsMean = c(1, 1, 1), dispersion = c(1, 1, 1),
  sampleFacVarMean = 1, sampleFacVarSD = 1,
  subjectFacVarMean = 1, subjectFacVarSD = 1,
  nSubjsPerGroup = 1, twoGroupDesign = F,
  nTimepoints = 1, maxCellsPerSamp = 1, minCellsPerSamp = 1,
  propDE = 0, deLog2FC = 0
)

test_that("getRescueSimParam working", {
  expect_equal(
    lapply(slotNames(fullParams), getRescueSimParam, paramObj = fullParams),
    list(1, FALSE, 1, 1, 1, 1, 1, 1, numeric(0), c(1, 1, 1), c(1, 1, 1), 1, 1,numeric(0), 1, 1,numeric(0), 0, 0)
  )
})

test_that("updateRescueSimParams working", {
  paramLS <- list(
    logLibFacVar = 1, logLibMean = 1, logLibSD = 1,
    exprsMean = c(1, 1, 1), dispersion = c(1, 1, 1),
    sampleFacVarMean = 1, sampleFacVarSD = 1,
    subjectFacVarMean = 1, subjectFacVarSD = 1,
    nSubjsPerGroup = 1, twoGroupDesign = F,
    nTimepoints = 1, maxCellsPerSamp = 1, minCellsPerSamp = 1,
    propDE = 0, deLog2FC = 0
  )
  expect_equal(updateRescueSimParams(blankParams, paramLS), fullParams)
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
      updateRescueSimParams(fullParams, singleValList),
      paste("Parameter", x, "should contain a single value")
    )
  })
})

test_that(".checkGeneParamLengths working", {
  expect_error(
    updateRescueSimParams(
      fullParams,
      list(
        exprsMean = c(1, 1),
        dispersion = c(1, 1, 1, 1)
      )
    ),
    "Parameters exprsMean and dispersion must be vectors of equal length"
  )
    expect_error(
        updateRescueSimParams(
            fullParams,
            list(
                twoGroupDesign=F,
                nTimepoints=3,
                exprsMean = c(1, 1),
                dispersion = c(1,1),
                deLog2FC=list(time1 = c(1,1,1), time2=c(1,1,1))
            )
        ),
        "The lengths of the following Log2FC vectors do not match exprsMean: time1, time2"
    )
})

test_that("deLog2FC: numeric values pass validation", {
    expect_silent(updateRescueSimParams(fullParams, list(deLog2FC = 1)))
    expect_silent(updateRescueSimParams(fullParams, list(deLog2FC = rnorm(3))))
})

test_that("deLog2FC: list with missing nTimepoints triggers error", {
    broken <- fullParams
    methods::slot(broken, "nTimepoints") <- numeric(0)
    expect_error(updateRescueSimParams(broken, list(deLog2FC = list(time1 = rnorm(3)))),
                 regexp = "must be specified")
})

test_that("deLog2FC: one timepoint and one group triggers error", {
    expect_error(updateRescueSimParams(fullParams, list(
        nTimepoints = 1,
        twoGroupDesign = FALSE,
        deLog2FC = list(time1 = rnorm(3))
    )), regexp = "DE cannot be defined")
})

test_that("deLog2FC: one timepoint, two groups - group1 valid, group0 invalid", {
    expect_silent(updateRescueSimParams(fullParams, list(
        nTimepoints = 1,
        twoGroupDesign = TRUE,
        deLog2FC = list(group1 = rnorm(3))
    )))
    expect_error(updateRescueSimParams(fullParams, list(
        nTimepoints = 1,
        twoGroupDesign = TRUE,
        deLog2FC = list(group0 = rnorm(3))
    )), regexp = "group0.*reference level")
})

test_that("deLog2FC: multi-timepoint, one group - time0 invalid", {
    expect_error(updateRescueSimParams(fullParams, list(
        nTimepoints = 2,
        twoGroupDesign = FALSE,
        deLog2FC = list(time0 = rnorm(3))
    )), regexp = "time0.*reference level")
})

test_that("deLog2FC: multi-timepoint, two groups - valid names only", {
    good <- list(
        time1_group0 = rnorm(3),
        time1_group1 = rnorm(3),
        time2_group0 = rnorm(3),
        time2_group1 = rnorm(3)
    )
    expect_silent(updateRescueSimParams(fullParams, list(
        nTimepoints = 3,
        twoGroupDesign = TRUE,
        deLog2FC = good
    )))

    bad <- good
    bad$time0_group0 <- rnorm(3)
    expect_error(updateRescueSimParams(fullParams, list(
        nTimepoints = 3,
        twoGroupDesign = TRUE,
        deLog2FC = bad
    )), regexp = "time0_group0.*reference level")
})

test_that("deLog2FC: unnamed list triggers error", {
    unnamed <- list(rnorm(3), rnorm(3))
    expect_error(updateRescueSimParams(fullParams, list(
        nTimepoints = 2,
        twoGroupDesign = FALSE,
        deLog2FC = unnamed
    )), regexp = "must be named")
})

test_that("deLog2FC: unexpected and missing names are caught", {
    expect_error(updateRescueSimParams(fullParams, list(
        nTimepoints = 2,
        twoGroupDesign = FALSE,
        deLog2FC = list(foo = rnorm(3))
    )), regexp = "Unexpected names")

    expect_error(updateRescueSimParams(fullParams, list(
        nTimepoints = 3,
        twoGroupDesign = FALSE,
        deLog2FC = list(time1 = rnorm(3))  # missing time2
    )), regexp = "Missing expected names")
})

test_that("deLog2FC: non-numeric elements trigger error", {
    expect_error(updateRescueSimParams(fullParams, list(
        nTimepoints = 2,
        twoGroupDesign = FALSE,
        deLog2FC = list(time1 = "not numeric")
    )), regexp = "not numeric vectors")
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
        updateRescueSimParams(fullParams, ls),
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
        updateRescueSimParams(fullParams, ls)
      )
    })
  })
})
