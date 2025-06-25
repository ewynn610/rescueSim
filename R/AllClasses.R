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
#' @slot customLibSizes An optional numeric vector of user-provided library sizes. When
#' provided, library sizes for simulation will be sampled (with replacement)
#' from this vector instead of the default log-normal distribution.
#' Log-transformed values are adjusted by a sample specific factor so that the
#' average log library size equals the overall average multiplied by the
#' sample-specific factor.
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
#'
#' Specifies the log2 fold changes for differentially expressed (DE) genes.
#' All values are interpreted as relative to a common baseline condition:
#' \strong{time0 and/or group0}.
#'
#' Acceptable formats include:
#' \itemize{
#'   \item A single numeric value \code{>= 0}, indicating the absolute log2 fold
#'    change between baseline and the final condition. Log2 fold changes will be
#'     randomly assigned as positive or negative across DE genes.
#'
#'   \item A numeric vector of possible log2 fold change values (positive and/or
#'    negative) from which values will be randomly drawn for DE genes.
#'
#'   \item A named list of numeric vectors specifying gene-specific log2 fold
#'   changes for each experimental condition, relative to the baseline
#'   (\code{"time0"}, \code{"group0"}, or \code{"time0_group0"}
#'   depending on design). Each vector must have length equal to the number of
#'   genes. Valid names for list elements depend on the experimental design:
#'   \itemize{
#'     \item For single-group, multi-timepoint designs: \code{"time1"},
#'     \code{"time2"}, etc.
#'     \item For two-group, single-timepoint designs: \code{"group1"}
#'     \item For multi-timepoint, two-group designs: \code{"time1_group0"},
#'     \code{"time1_group1"}, etc.
#'   }
#'   The reference condition — \code{"time0"}, \code{"group0"}, or
#'   \code{"time0_group0"} — must \strong{not} be included in the list, as it is
#'    implicitly treated as having log2FC = 0.
#' }
#'
#' When a single value or vector is provided (i.e., not a list), log2 fold
#' changes are applied in a structured way:
#' \itemize{
#'   \item For two-group, two-timepoint designs: group 0 shows no change over
#'   time, while group 1 exhibits a linear change over time from 0 to the
#'   specified log2FC.
#'   \item For designs with more than two timepoints: a linear trajectory is
#'   simulated, and the log2FC value represents the total change from time 0 to
#'   the final timepoint.
#' }
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
             customLibSizes = "numeric",
             exprsMean = "numeric",
             dispersion = "numeric",
             sampleFacVarMean = "numeric",
             sampleFacVarSD = "numeric",
             subjectFacVarMean = "numeric",
             subjectFacVarSD = "numeric",
             propDE = "numeric",
             deLogFC = "ANY"
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

RescueParams <- function(deLogFC = 0, propDE = 0, ...) {
    obj <- methods::new("RescueParams",deLogFC = deLogFC, propDE = propDE,...)
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
    geneParamLengthError <- .checkGeneParamLengths(param_lengths, methods::slot(object, "deLogFC"))

    deLogFCError<-.checkDeLogFCFormat(methods::slot(object, "deLogFC"),
                                      methods::slot(object, "nTimepoints"),
                                      methods::slot(object, "twoGroupDesign"))

    ## maxCellsPerSamp, minCellsPerSamp
    cellsPerSampError <- .checkCellsPerSampLength(
        object, param_lengths
    )

    ## Make sure specified parameters are not negative
    paramValuesError <- .checkParamValues(object, param_lengths)

    propDEError<-.checkPropDEUpperBound(methods::slot(object, "propDE"))

    error_string <- c(
        singleValueError, deLogFCError,geneParamLengthError, cellsPerSampError,
        propDEError
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
                              2 * getRescueParam(object, "nTimepoints"),
                              getRescueParam(object, "nTimepoints")
        )

        cellsPerSampIndicator <-
            param_lengths[c("maxCellsPerSamp", "minCellsPerSamp")] %in%
            c(1, nconditions, nsamps)
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
        logLibSD = "lt", customLibSizes = "lt",
        exprsMean = "lteq",
        dispersion = "lteq",
        sampleFacVarSD = "lt",
        subjectFacVarSD = "lt",
        nSubjsPerGroup = "lteq",
        nTimepoints = "lteq",
        maxCellsPerSamp = "lteq",
        minCellsPerSamp = "lteq",
        propDE = "lt"
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

.checkDeLogFCFormat <- function(deLogFC, nTimepoints, twoGroupDesign) {

    # CASE 1: single value or numeric vector
    if (is.numeric(deLogFC)) {
        return(NULL)  # valid
    }

    # CASE 2: list of named vectors
    # CASE 3: list of named numeric vectors
    if (is.list(deLogFC)) {
        if (length(nTimepoints) == 0 || length(twoGroupDesign) == 0) {
            return("If deLogFC is a list, both nTimepoints and twoGroupDesign must be specified.")
        }

        if (nTimepoints == 1 && !twoGroupDesign) {
            return("With only one timepoint and one group variable, DE cannot be defined.")
        }

        # Define expected and forbidden names
        if (nTimepoints == 1 && twoGroupDesign) {
            expected_names <- "group1"
            forbidden_name <- "group0"
        } else if (!twoGroupDesign) {
            expected_names <- paste0("time", 1:(nTimepoints - 1))
            forbidden_name <- "time0"
        } else {
            expected_names <- unlist(lapply(1:(nTimepoints - 1), function(t) {
                paste0("time", t, "_group", 0:1)
            }))
            forbidden_name <- "time0_group0"
        }

        provided_names <- names(deLogFC)

        if (is.null(provided_names)) {
            return("If deLogFC is a list, all elements must be named.")
        }

        if (forbidden_name %in% provided_names) {
            return(paste0("'", forbidden_name, "' should not be included in deLogFC; it is the reference level."))
        }

        # Check for unexpected or missing names
        unexpected <- setdiff(provided_names, expected_names)
        missing <- setdiff(expected_names, provided_names)

        if (length(unexpected) > 0) {
            return(paste("Unexpected names in deLogFC list:", paste(unexpected, collapse = ", ")))
        }

        if (length(missing) > 0) {
            return(paste("Missing expected names in deLogFC list:", paste(missing, collapse = ", ")))
        }

        # Check all elements are numeric vectors
        bad_type <- names(deLogFC)[!vapply(deLogFC, is.numeric, logical(1))]
        if (length(bad_type) > 0) {
            return(paste("The following elements of deLogFC are not numeric vectors:",
                         paste(bad_type, collapse = ", ")))
        }

        return(NULL)
    }

    return("deLogFC must be NULL, a single numeric value, a numeric vector, or a named list of numeric vectors.")
}

.checkPropDEUpperBound <- function(propDE) {
    if(length(propDE)==1){
    if (propDE > 1 ) {
        return("propDE must be less than or equal to 1.")
    }
    }
    return(NULL)
}

.checkGeneParamLengths <- function(param_lengths, deLogFC) {
    errors <- character()

    # Check exprsMean and dispersion lengths
    if (param_lengths["exprsMean"] != param_lengths["dispersion"] &&
        (param_lengths["exprsMean"] != 0 || param_lengths["dispersion"] != 0)) {
        errors <- c(errors, "Parameters exprsMean and dispersion must be vectors of equal length")
    }

    # Check deLogFC list element lengths (if it's a list)
    if (is.list(deLogFC)) {
        bad_lengths <- vapply(deLogFC, function(x) length(x) != param_lengths["exprsMean"], logical(1))
        if (any(bad_lengths)) {
            bad_keys <- names(deLogFC)[bad_lengths]
            msg <- paste0("The lengths of the following logFC vectors do not match exprsMean: ",
                          paste(bad_keys, collapse = ", "))
            errors <- c(errors, msg)
        }
    }

    if (length(errors) == 0) return(NULL)
    return(paste(errors, collapse = "\n"))
}
