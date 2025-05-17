#' The RescueParams Class
#'
#' Class for holding parameters used to simulate data using \code{RESCUE}
#'
#' @slot nTimepoints Number of timepoints (i.e. samples) per subject.
#' Holds a single numeric value >0 representing the number of timepoints for
#' all subjects.
#' @slot twoGroupDesign Logical value indicating whether to simulate two groups
#' (ex. treatment and control group).
#' @slot nSubjsPerGroup Number of subjects per group (if using two group design)
#'  or number of total subjects (if only single group). Holds a single
#'  numeric value >0 indicating the number of subjects per group.
#' @slot maxCellsPerSamp Maximum parameter used when drawing number of cells per
#'  sample from a discrete uniform distribution. Holds a single numeric value
#'   >0 indicating the maximum number of cells per sample, or a vector with
#'   length equal to the number of conditions (group/timepoint combinations)
#'   with each value representing the maximum for a single condition, or a
#'   vector with length equal to the number of total samples with each value
#'   representing the maximum for a single sample.
#' @slot minCellsPerSamp Minimum parameter used when drawing number of cells per
#'  sample from a discrete uniform distribution. Holds a single numeric value >0
#'  indicating the minimum number of cells per sample, or a vector with length
#'  equal to the number of conditions (group/timepoint combinations) with each
#'  value representing the minimum for a single condition, or a vector with
#'  length equal to the number of total samples with each value representing
#'  the minimum for a single sample.
#' @slot logLibMean Mean library size (log scale) parameter. Used to draw
#'  library size from a log-normal distribution. Holds a single value >0
#'  indicating the mean library size on a log scale.
#' @slot logLibSD Library size standard deviation (log scale) parameter.
#' Used to draw library size from a log-normal distribution. Holds a
#' single value >=0 indicating the standard deviation of library size on
#' a log scale.
#' @slot logLibFacVar Variance used for drawing sample-level multiplicative
#' factors which give different library size distributions for each sample.
#' Larger values result in larger variation in library size distributions by
#' sample. Holds a single value >=0 with 0 indicating no difference in the
#' library size distributions by sample.
#' @slot exprsMean Gene-specific mean expression value representing the average
#' expression of each gene in the dataset. Holds a vector of numeric values >=0
#' with length equal to desired number of genes in the simulated data where each
#' value indicates the average expression for a single gene.
#' @slot dispersion Gene-specific dispersion value representing the variation in
#' expression for each gene in the dataset. Holds a vector of numeric values >0
#' with length equal to desired number of genes in the simulated data where each
#' value indicates the dispersion for a single gene.
#' @slot sampleFacVarMean Mean used for drawing variance (log-scale) of
#' sample-level multiplicative factors. Larger values result in
#' more between-sample variation.
#' @slot sampleFacVarSD Standard deviation used for drawing variance (log-scale)
#' of sample-level multiplicative factors. Larger values result in more
#' variation in amount of between-sample variation across different genes.
#' Must be a value >=0
#' @slot subjectFacVarMean Mean used for drawing variance (log-scale) of
#' subject-level multiplicative factors. Larger values result in larger
#' between-subject variation.
#' @slot subjectFacVarSD Standard deviation used for drawing variance
#' (log-scale) of subject-level multiplicative factors. Larger values result in
#' more variation in amount of between-subject variation across different genes.
#' Must be a value >=0.
#' @slot propDE Proportion of genes differentially expressed between
#' timepoints/groups. Must be a numeric value between 0 and 1.
#' @slot deLogFC Fold change values used for differentially expressed genes.
#' A single numeric value >=0 indicating the positive log-fold change for
#' differentially expressed genes (log2 fold change values will be simulated
#' evenly to be positive or negative) or a vector of values (positive of
#' negative) to draw log2 fold change values from.
#'
#' @author Elizabeth Wynn
#'
#'
#' @export


setClass("RescueParams",
  slots = c(
    nTimepoints = "numeric",
    twoGroupDesign = "logical",
    nSubjsPerGroup = "numeric",
    maxCellsPerSamp = "numeric",
    minCellsPerSamp = "numeric",
    logLibMean = "numeric",
    logLibSD = "numeric",
    logLibFacVar = "numeric",
    exprsMean = "numeric",
    dispersion = "numeric",
    sampleFacVarMean = "numeric",
    sampleFacVarSD = "numeric",
    subjectFacVarMean = "numeric",
    subjectFacVarSD = "numeric",
    propDE = "numeric",
    deLogFC = "numeric"
  )
)


#' Generate RescueParams object
#'
#' @return An object of class \code{RescueParams}
#'
#' @param ... Any parameter name followed by initial valus
#'
#' @author Elizabeth Wynn
#'
#' @examples
#' ## Create params object
#' myParams <- RescueParams()
#'
#' ## Create parameter object with nTimepoints and nSubjsPerGroup pre-set
#' myParams <- RescueParams(nTimepoints = 2, nSubjsPerGroup = 5)
#'
#' @export

RescueParams <- function(...) {
  obj <- methods::new("RescueParams", ...)
  methods::validObject(obj)
  obj
}

#' Update RescueParams object
#'
#' Manually update parameters in RescueParams object
#'
#' @param paramObj Object of class \code{\link{RescueParams-class}}
#' @param paramValues List of parameter values with list names are the parameter
#'  names
#'
#' @return Object of class \code{\link{RescueParams-class}} with specified parameters updated
#'
#' @examples
#' ## Create Parameter object
#' myParams <- RescueParams()
#'
#' ## Update nTimepoints and nSubjsPerGroup
#' myParams <- updateRescueParams(
#'   paramObj = myParams,
#'   paramValues = list(
#'     nTimepoints = 3,
#'     nSubjsPerGroup = 10
#'   )
#' )
#' ## Check Values
#' getRescueParam(myParams, "nTimepoints")
#' getRescueParam(myParams, "nSubjsPerGroup")
#'
#' @export

updateRescueParams <- function(paramObj, paramValues) {
  for (name in names(paramValues)) {
    methods::slot(paramObj, name) <- paramValues[[name]]
  }
  methods::validObject(paramObj)
  paramObj
}

#' Extract Parameters
#'
#' Extract individual parameters from RescueParams object
#'
#' @param paramObj Object of class \code{\link{RescueParams-class}}
#' @param paramName String containing parameter you would like to extract
#'
#' @return Specified parameter value
#'
#' @examples
#' ## Create Parameter object
#' myParams <- RescueParams(nTimepoints = 2)
#'
#' ## Check Values
#' getRescueParam(myParams, "nTimepoints")
#'
#' @export

getRescueParam <- function(paramObj, paramName) {
  methods::slot(paramObj, paramName)
}


## Helper Functions

## Check param validity
rescueParamsValidity <- function(object) {
  ## Make sure values have the right lengths
  param_lengths <- vapply(
    methods::slotNames(object), function(slot) {
      length(slot(object, slot))
    },
    numeric(1)
  )

  ## Single Value Parameters
  singleValueError <- .checkSingleValueParams(param_lengths)

  ## ExprsMean and dispersion
  exprsDispError <- .checkExprsDispLength(param_lengths)

  ## maxCellsPerSamp, minCellsPerSamp
  cellsPerSampError <- .checkCellsPerSampLength(
    object, param_lengths
  )

  ## Make sure specified parameters are not negative
  paramValuesError <- .checkParamValues(object, param_lengths)

  error_string <- c(
    singleValueError, exprsDispError, cellsPerSampError
  )


  names(error_string) <- NULL
  if (length(error_string) == 0) {
    return(TRUE)
  } else {
    return(error_string)
  }
}

setValidity("RescueParams", rescueParamsValidity)




.checkSingleValueParams <- function(param_lengths) {
  singleValueSlots <- c(
    "nTimepoints", "nSubjsPerGroup",
    "logLibFacVar", "logLibMean", "logLibSD",
    "sampleFacVarMean",
    "sampleFacVarSD",
    "subjectFacVarMean",
    "subjectFacVarSD",
    "twoGroupDesign",
    "propDE"
  )
  singleValueIndicator <- param_lengths[singleValueSlots] > 1
  error_string <- vector()
  if (any(singleValueIndicator)) {
    singleValueError <- vapply(
      names(singleValueIndicator[singleValueIndicator]),
      function(x) {
        paste(
          "Parameter", x,
          "should contain a single value"
        )
      },
      character(1)
    )
  } else {
    singleValueError <- NULL
  }
  return(singleValueError)
}

.checkExprsDispLength <- function(param_lengths) {
  if (param_lengths["exprsMean"] != param_lengths["dispersion"] &
    (param_lengths["exprsMean"] != 0 & param_lengths["dispersion"] != 0)) {
    exprsDispError <-
      "Parameters exprsMean and dispersion must be vectors of equal length"
  } else {
    exprsDispError <- NULL
  }
  return(exprsDispError)
}


.checkCellsPerSampLength <- function(object, param_lengths) {
  cellsPerSampIndicator <- param_lengths[c(
    "maxCellsPerSamp",
    "minCellsPerSamp"
  )] > 0
  cellsPerSampErr <- NULL
  if (any(cellsPerSampIndicator) &
    param_lengths["twoGroupDesign"] == 1) {
    nsubjs <- ifelse(getRescueParam(object, "twoGroupDesign"),
        2 * getRescueParam(object, "nSubjsPerGroup"),
        getRescueParam(object, "nSubjsPerGroup")
      )

    nsamps <- nsubjs * getRescueParam(object, "nTimepoints")

    nconditions <- ifelse(getRescueParam(object, "twoGroupDesign"),
                          2 * getRescueParam(object, "nTimepoint"),
                          getRescueParam(object, "nTimepoints")
    )

    cellsPerSampIndicator <-
      param_lengths[c("maxCellsPerSamp", "minCellsPerSamp")] %in%
      c(1, 2, nsamps)
    names(cellsPerSampIndicator) <- c("maxCellsPerSamp", "minCellsPerSamp")
    if (any(!cellsPerSampIndicator)) {
      cellsPerSampErr <-
        vapply(
          names(cellsPerSampIndicator[!cellsPerSampIndicator]),
          function(x) {
            paste0(
              "Parameter ", x,
              " should contain a single value (",
              x, " for all samples), a vector with length equal to the number of conditions (group/timepoint combinations),or a vector with length equal to the number of total samples."
            )
          },
          character(1)
        )
    }
  }
  return(cellsPerSampErr)
}


.checkParamValues <- function(object, param_lengths) {
  ## lt = parameters that can't be less than 0
  ## lteq = parameters that can't be less than or equal to 0
  paramsGT0 <- c(
    logLibFacVar = "lt", logLibMean = "lteq",
    logLibSD = "lt", exprsMean = "lteq",
    dispersion = "lteq",
    sampleFacVarSD = "lt",
    subjectFacVarSD = "lt",
    nSubjsPerGroup = "lteq",
    nTimepoints = "lteq",
    maxCellsPerSamp = "lteq",
    minCellsPerSamp = "lteq",
    propDE = "lt",
    deLogFC = "lt"
  )

  LT0Indicator <- vapply(names(paramsGT0), function(x) {
    if (param_lengths[x] != 0) {
      ifelse(paramsGT0[x] == "lt", any(methods::slot(object, x) < 0),
        any(methods::slot(object, x) < 0)
      )
    } else {
      FALSE
    }
  }, logical(1))

  if (any(LT0Indicator)) {
    LT0Error <- vapply(
      names(LT0Indicator[LT0Indicator]), function(x) {
        lteq <- ifelse(paramsGT0[x] == "lt", "less than", "less than or equal to")
        paste("Parameter", x, "should not contain values", lteq, "0")
      },
      character(1)
    )
  } else {
    LT0Error <- NULL
  }
  return(LT0Error)
}

.checkExprsDispLength <- function(param_lengths) {
  if (param_lengths["exprsMean"] != param_lengths["dispersion"] &
    (param_lengths["exprsMean"] != 0 | param_lengths["dispersion"] != 0)) {
    exprsDispError <-
      "Parameters exprsMean and dispersion must be vectors of equal length"
  } else {
    exprsDispError <- NULL
  }
  return(exprsDispError)
}
