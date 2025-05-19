#' Estimate RESCUE Parameters
#'
#' Estimate parameters to be used to simulate repeated measures scRNA-Seq data
#' @param sce \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#' containing empirical data with a counts matrix in the \code{counts} slot and
#' cell level data (ex. sample/subject/timepoint labels) in the \code{colData}
#' slot.
#' @param paramObj \code{\link{RescueParams-class}} object with empty slots
#' for parameters which need to be estimated and parameters. If \code{NULL},
#' all parameter values will be estimated.
#' @param sampleVariable String denoting name of sample identifier variable in
#' the \code{colData} of \code{SingleCellExperiment} object. If \code{NULL},
#' data is assumed to contain only one sample parameters
#' \code{logLibFacVar}, \code{sampleFacVarMean}, \code{sampleFacVarSD},
#' \code{subjectFacVarMean}, and \code{subjectFacVarSD} parameters
#' cannot be estimated.
#' @param subjectVariable String denoting name of subject identifier variable in
#'  the \code{colData} of \code{SingleCellExperiment} object. If \code{NULL},
#'  parameters \code{nSubjsPerGroup}, \code{subjectFacVarMean}, and
#'  \code{subjectFacVarSD} parameters cannot be estimated.
#' @param timepointVariable String denoting name of timepoint identifier
#' variable in the \code{colData} of \code{SingleCellExperiment} object.
#' If \code{NULL}, the number of timepoints for the simulated data is set to 2.
#' @param groupVariable String denoting name of group identifier variable in
#' the \code{colData} of \code{SingleCellExperiment} object. If \code{NULL}, a
#' single group design is assumed.
#' @param nonDEGs Vector containing names or row indices of genes to be used to
#' estimate sample and subject factor parameters. Using genes that are
#' differentially expressed across timepoints may lead to inaccurate estimates.
#' @param cellParamsByCondition Logical value indicating whether
#' \code{maxCellsPerSamp} and \code{minCellsPerSamp} should be estimated for
#' each condition (group/timepoint) or overall.
#'
#' @details
#' All parameters in \code{\link{RescueParams-class}} object can be
#' estimated/extracted from empirical data except \code{propDE} and
#' \code{deLogFC} which are both set to 0 (no differential expression) if not
#' set manually. If a paramObj is provided, only parameters with empty slots in
#' the object will be estimated.
#'
#' @author Elizabeth Wynn
#'
#' @examples
#' ## Load the data
#' data("RecAM_sce")
#'
#' ## Subset data to 40 genes
#' RecAM_sce_small <- RecAM_sce[1:40, ]
#'
#' ## Estimate all parameters from the data
#' myParams_estimated <- estRescueParams(
#'   sce = RecAM_sce_small, paramObj = NULL,
#'   sampleVariable = "sampleID",
#'   subjectVariable = "subjectID",
#'   timepointVariable = "time",
#'   groupVariable = NULL, nonDEGs = NULL,
#'   cellParamsByCondition = FALSE
#' )
#'
#' ## Example using Single timepoint data
#'
#' ## Subset the data to a single timepoint
#' RecAM_sce_single_time <- RecAM_sce_small[, RecAM_sce_small$time == "time0"]
#' ## Create a params object with necessary parameters filled
#' myParams <- RescueParams(subjectFacVarMean = -4, subjectFacVarSD = 1, nSubjsPerGroup = 5)
#' ## Estimate remaining parameters from the data
#' myParams_estimated <- estRescueParams(
#'   sce = RecAM_sce_single_time,
#'   paramObj = myParams,
#'   sampleVariable = "sampleID"
#' )
#'
#' @export

estRescueParams <- function(sce, paramObj = NULL,
                            sampleVariable = NULL,
                            subjectVariable = NULL,
                            groupVariable = NULL,
                            timepointVariable = NULL,
                            nonDEGs = NULL,
                            cellParamsByCondition = FALSE) {
  ## Check sce is good
  checkmate::assertClass(sce, "SingleCellExperiment")

  ## Check params is good
    if(!is.null(paramObj)){
        checkmate::assertClass(paramObj, "RescueParams")
    }

  ## Check sample var and subject and group variable var are in sce
  colDat_var <- c(
    "sampleVariable" = sampleVariable, "subjectVariable" = subjectVariable,
    "groupVariable" = groupVariable, "timepointVariable" = timepointVariable
  )
  test <- lapply(seq_len(length(colDat_var)), function(idx) {
    varName <- names(colDat_var)[idx]
    varVal <- colDat_var[idx]
    if (!is.character(varVal)) {
      stop(paste(varName, "must be NULL or a character string"))
    } else {
      if (!(varVal %in% names(colData(sce)))) {
        stop(paste(varName, "must be the name of a variable stored in colData of the sce object"))
      } else {
        return(invisible(TRUE))
      }
    }
  })



  ## Check nonde genes are in sce
  .checkNonDEGs(nonDEGs, sce)

  ## Check cellParamsByCondition
  checkmate::assertLogical(cellParamsByCondition)

  ## Make params object if not provided
  if (is.null(paramObj)) paramObj <- RescueParams()

  ## Find which slots need to be estimated
  est_indicator <- vapply(
    methods::slotNames(paramObj), function(slot) {
      length(getRescueParam(paramObj, slot)) == 0
    },
    logical(1)
  )

  ## Go through estimation process if any need to be estimated
  if (any(est_indicator)) {
    if (length(getRescueParam(paramObj, "propDE")) == 0) {
      paramObj <- updateRescueParams(paramObj, list(propDE = 0))
      est_indicator["propDE"]=F
    }
    if (length(getRescueParam(paramObj, "deLogFC")) == 0) {
      paramObj <- updateRescueParams(paramObj, list(deLogFC = 0))
      est_indicator["deLogFC"]=F

    }




    ## Check if data needs to be normalized
    norm_data_slots <- c(
      "dispersion", "exprsMean", "sampleFacVarMean",
      "sampleFacVarSD", "subjectFacVarMean", "subjectFacVarSD"
    )
    norm_data_indicator <- any(est_indicator[norm_data_slots])

    ## Check if batch parameters need to be estimated
    batch_var_slots <- c("sampleFacVarSD", "subjectFacVarMean", "subjectFacVarSD")
    batch_var_indicator <- any(est_indicator[batch_var_slots])

    ## Check if library mean/SD needs to be estimated
    logLibMeanSD_slots <- c("logLibMean", "logLibSD")
    logLibMeanSD_indicator <- any(est_indicator[logLibMeanSD_slots])

    ## Check if minCellsPerSamp/maxCellsPerSamp needs to be estimated
    cellsPerSamp_slots <- c("minCellsPerSamp", "maxCellsPerSamp")
    cellsPerSamp_indicator <- any(est_indicator[cellsPerSamp_slots])

    ## Normalize data if needed
    if (norm_data_indicator) {
      sce <- scuttle::computePooledFactors(sce)
      SingleCellExperiment::normcounts(sce) <-
        scuttle::normalizeCounts(sce, log = FALSE)
    }

    ## Estimate needed parameters
    if (logLibMeanSD_indicator) {
      logLibMeanSD <- .estLibMeanSD(sce, sampleVariable)
    } else {
      logLibMeanSD <- NULL
    }

    params_est_first <- c(
      "exprsMean", "dispersion",
      "twoGroupDesign", "nTimepoints"
    )
    est_indicator_fil <- est_indicator[params_est_first]

    paramObj <- .getParams(param_names=names(est_indicator_fil[est_indicator_fil]),
                                    paramObj=paramObj, sce=sce, sampleVariable=sampleVariable,
                                    subjectVariable=subjectVariable, groupVariable=groupVariable,
                                    logLibMeanSD=logLibMeanSD, timepointVariable=timepointVariable)

    twoGroupDesign <- getRescueParam(paramObj, "twoGroupDesign")
    nTimepoints <- getRescueParam(paramObj, "nTimepoints")

    if (batch_var_indicator) {
      if (!est_indicator[["dispersion"]]) {
        phi <- .estPhi(sce, sampleVariable)
      } else {
        phi <- getRescueParam(paramObj, "dispersion")
      }
      if (!est_indicator[["exprsMean"]]) {
        means <- .estMeanExpr(sce, sampleVariable)
      } else {
        means <- getRescueParam(paramObj, "exprsMean")
      }

      batchVarParams <- .estBatchVarParams(
        sce, sampleVariable,
        subjectVariable, nonDEGs,
        phi, means
      )
    }

    est_indicator[params_est_first] <- FALSE

    if (cellsPerSamp_indicator) {
      cellsPerSamp <- .estNCellParams(
        sce, sampleVariable,
        cellParamsByCondition, groupVariable,
        timepointVariable, nTimepoints,
        twoGroupDesign
      )
    } else {
      cellsPerSamp <- NULL
    }



    paramObj <- .getParams(
      names(est_indicator[est_indicator]),
      paramObj, sce, sampleVariable, subjectVariable,
      groupVariable, logLibMeanSD, batchVarParams,
      twoGroupDesign, cellsPerSamp,
      timepointVariable
    )


    return(paramObj)
  } else {
    warning("All parameters in paramObj already filled. No parameters estimated.")
    return(paramObj)
  }
}

## Helper Functions

.checkNonDEGs <- function(nonDEGs, sce) {
  checkmate::assert(
    checkmate::checkNull(nonDEGs),
    checkmate::checkCharacter(nonDEGs),
    checkmate::checkNumeric(nonDEGs)
  )
  if (is.numeric(nonDEGs)) {
    if (any(nonDEGs) <= 0) {
      stop("Index values for nonDEGs must be positive")
    }
    if (max(nonDEGs) > nrow(sce)) {
      stop("Index value(s) provided in nonDEGs greater than number of genes in sce")
    }
  } else if (is.character(nonDEGs)) {
    if (any(!(nonDEGs %in% rownames(sce)))) {
      stop("Gene names provided in nonDEGs not present in sce")
    }
  }
}

.getParams <- function(param_names, paramObj, sce = NULL, sampleVariable = NULL,
                       subjectVariable = NULL,
                       groupVariable = NULL, logLibMeanSD = NULL,
                       batchVarParams = NULL, twoGroupDesign = NULL,
                       cellsPerSamp = NULL, timepointVariable = NULL) {
  param_list <- lapply(param_names, function(type) {
    switch(type,
      logLibFacVar = .estLibFactorVar(sce, sampleVariable),
      logLibMean = logLibMeanSD[["mu"]],
      logLibSD = logLibMeanSD[["sd"]],
      exprsMean = .estMeanExpr(sce, sampleVariable),
      dispersion <- .estPhi(sce, sampleVariable),
      sampleFacVarMean = batchVarParams[["samp_mean"]],
      sampleFacVarSD = batchVarParams[["samp_sd"]],
      subjectFacVarMean = batchVarParams[["subj_mean"]],
      subjectFacVarSD = batchVarParams[["subj_sd"]],
      nSubjsPerGroup = .estnSubjsPerGroup(
        sce, subjectVariable,
        groupVariable,
        twoGroupDesign
      ),
      twoGroupDesign = .estTwoGroupDesignParam(sce, groupVariable),
      nTimepoints = .estNTimepoints(sce, timepointVariable),
      maxCellsPerSamp = cellsPerSamp[["maxCells"]],
      minCellsPerSamp = cellsPerSamp[["minCells"]]
    )
  })
  names(param_list) <- param_names
  updateRescueParams(paramObj, param_list)
}


## Function to estimate library factor variances
.estLibFactorVar <- function(sce, sampleVariable) {
  ## make dataframe of library sizes and sample the cell came from
  lib_sizes <- data.frame(
    log_lib_size = log(Matrix::colSums(SingleCellExperiment::counts(sce))),
    sample_id = sce[[sampleVariable]]
  )

  ## Find multiplicative factors showing how average sample lib. size differs

  lib_means <- stats::aggregate(lib_sizes$log_lib_size,
    list(sample_id = lib_sizes$sample_id),
    FUN = mean
  )
  overall_mean <- mean(lib_means$x)
  lib_size_facs <- lib_means$x / overall_mean

  ## Take variance of library size factors
  stats::var(lib_size_facs)
}

## Get mean/sd log library size
.estLibMeanSD <- function(sce, sampleVariable) {
  ## Log lib size by sample
  lib_sizes <- data.frame(
    log_lib_size = log(Matrix::colSums(SingleCellExperiment::counts(sce)))
  )
  if (!is.null(sampleVariable)) {
    lib_sizes$sample_id <- sce[[sampleVariable]]
    ## Get mean and sd log lib size for each sample
    lib_means <- stats::aggregate(lib_sizes$log_lib_size,
      list(sample_id = lib_sizes$sample_id),
      FUN = mean
    )
    lib_sds <- stats::aggregate(lib_sizes$log_lib_size,
      list(sample_id = lib_sizes$sample_id),
      FUN = stats::sd
    )

    ## Take average mean and sd across samples
    lib_mu <- mean(lib_means$x)
    lib_sd <- mean(lib_sds$x)
  } else {
    ## Take average mean and sd across all cells
    lib_mu <- mean(lib_sizes$log_lib_size)
    lib_sd <- sd(lib_sizes$log_lib_size)
  }


  return(c(mu = lib_mu, sd = lib_sd))
}

## Estimate mean (normalized) expression
.estMeanExpr <- function(sce, sampleVariable) {
  if (!is.null(sampleVariable)) {
    samp_means=suppressMessages(dplyr::bind_cols(lapply(unique(sce[[sampleVariable]]), function(x){
      rowMeans(SingleCellExperiment::normcounts(sce)[,sce[[sampleVariable]]==x])
    })))

    gene_means=rowMeans(samp_means)
    names(gene_means)=names(sce)
  } else {
    gene_means <- Matrix::rowMeans(SingleCellExperiment::normcounts(sce))
  }
  return(gene_means)
}


## Estimate dispersion parameter for each gene
.estPhi <- function(sce, sampleVariable) {
  if (!is.null(sampleVariable)) {
    samp_tab <- table(sce[[sampleVariable]])
    dispersionEstSample <- names(which.max(samp_tab))
    sce <- sce[, sce[[sampleVariable]] == dispersionEstSample]
  }

  design_mat <- stats::model.matrix(~1, data = data.frame(colnames(sce)))
  dge <- scran::convertTo(sce, "edgeR")
  dge <- edgeR::estimateDisp(dge, design_mat,
    prior.df = 0, min.row.sum = 0,
    trend.method = "none"
  )
  phi <- dge$tagwise.dispersion
  names(phi) <- rownames(sce)
  phi
}


## Estimate mean/variance for batch effect parameters
.estBatchVarParams <- function(sce, sampleVariable, subjectVariable, nonDEGs, phi,
                               means, perc_0_keep = .6) {
  if (!is.null(nonDEGs)) {
    sce <- sce[nonDEGs]
    phi <- phi[nonDEGs]
    means <- means[nonDEGs]
  }
  perc0 <- Matrix::rowSums(SingleCellExperiment::counts(sce) == 0) /
    ncol(SingleCellExperiment::counts(sce))
  idx_keep <- perc0 < perc_0_keep
  sce <- sce[idx_keep, ]
  phi <- phi[idx_keep]
  means <- means[idx_keep]
  batch_vars <- .calcBatchVars(sce, sampleVariable, subjectVariable)
  if (is.null(subjectVariable)) batch_vars <- data.frame(samp_var = batch_vars)
  mean_cells <- ifelse(is.null(subjectVariable), mean(table(sce[[sampleVariable]])),
    mean(table(sce[[subjectVariable]]))
  )
  null_vars <- (phi * means^2 + means) / (mean_cells * means^2)
  types <- if (is.null(subjectVariable)) c("samp") else c("samp", "subj")
  batchParams <- lapply(types, function(type) {
    name <- paste0(type, "_var")
    df_fil <- data.frame(
      tot_var = batch_vars[, name], means = means, phi = phi,
      null_vars = null_vars
    )
    df_fil$corrected <- df_fil$tot_var - df_fil$null_vars
    n <- ifelse(type == "samp", length(unique(sce[[sampleVariable]])),
      length(unique(sce[[subjectVariable]]))
    )
    var_tot <- stats::var(df_fil$tot_var)
    var_null <- stats::var(df_fil$null_vars)
    mean_tot <- mean(df_fil$tot_var^2)
    mean_tot_2 <- mean(df_fil$tot_var^2)
    alpha <- var_tot - 2 * mean_tot_2 / (n - 1)
    alpha <- (var_tot - 2 * mean_tot^2 / (n - 1)) / ((n - 1 + 2) / (n - 1))
    my_mean <- mean(df_fil$corrected)
    my_var <- alpha - var_null
    my_lvar <- log(my_var / (my_mean^2) + 1)
    my_lmean <- log(my_mean) - my_lvar / 2
    my_lsd <- sqrt(my_lvar)
    list(my_lmean = my_lmean, my_lsd = my_lsd)
  })
  names(batchParams) <- types
  if (is.null(subjectVariable)) {
    batch_vars <- list(samp = batchParams$samp$batch_vars)
    ret <- list(
      samp_mean = batchParams$samp$my_lmean, samp_sd = batchParams$samp$my_lsd
    )
  } else {
    batch_vars <- list(samp = batchParams$samp$batch_vars, subj = batchParams$subj$batch_vars)
    ret <- list(
      samp_mean = batchParams$samp$my_lmean, samp_sd = batchParams$samp$my_lsd,
      subj_mean = batchParams$subj$my_lmean, subj_sd = batchParams$subj$my_lsd
    )
  }
}

.calcBatchVars <- function(sce, sampleVariable, subjectVariable) {
  if (!is.null(subjectVariable)) {
    samp_means <- stats::aggregate(
      as.matrix(Matrix::t(SingleCellExperiment::normcounts(sce))),
      list(
        samples = sce[[sampleVariable]],
        subjects = sce[[subjectVariable]]
      ),
      mean
    )
    subj_means <- stats::aggregate(
      samp_means[, -seq_len(2)],
      list(subjects = samp_means$subjects), mean
    )
    overall_mean <- colMeans(samp_means[, -seq_len(2)])

    subj_batch_eff <- t(subj_means[, -1]) / overall_mean

    rownames(subj_means) <- subj_means$subjects
    subj_means2 <- subj_means[paste0(samp_means$subjects), ]
    samp_batch_eff <- as.matrix(t(samp_means[, -seq_len(2)] / subj_means2[, -1]))
    samp_batch_eff[is.nan(samp_batch_eff)] <- 1

    batch_vars <- data.frame(
      samp_var = apply(samp_batch_eff, 1, var),
      subj_var = apply(subj_batch_eff, 1, var)
    )
  } else {
    samp_means <- stats::aggregate(
      as.matrix(Matrix::t(SingleCellExperiment::normcounts(sce))),
      list(samples = sce[[sampleVariable]]), mean
    )
    overall_mean <- colMeans(samp_means[, -1])


    samp_batch_eff <- t(samp_means[, -1]) / overall_mean
    batch_vars <- data.frame(samp_var = apply(samp_batch_eff, 1, var))
  }
}

.estTwoGroupDesignParam <- function(sce, groupVariable) {
  if (is.null(groupVariable)) {
    return(FALSE)
  } else if (length(unique(sce[[groupVariable]])) == 1) {
    return(FALSE)
  } else if (length(unique(sce[[groupVariable]])) == 2) {
    return(TRUE)
  } else {
    warning("More than two groups present in user provided data.
                Simulating more than two groups is not currently supported.
                Setting simulation parameter to simulate two groups.")
    return(TRUE)
  }
}

.estNTimepoints <- function(sce, timepointVariable) {
  if (is.null(timepointVariable)) {
    return(2)
  } else {
    return(length(unique(sce[[timepointVariable]])))
  }
}

.estnSubjsPerGroup <- function(sce, subjectVariable, groupVariable, twoGroupDesign) {
  nGroups <- ifelse(twoGroupDesign, 2, 1)
  if (is.null(groupVariable)) {
    nSubjsPerGroup <- length(table(sce[[subjectVariable]]))
    return(nSubjsPerGroup)
  }
  subj_group <- table(sce[[subjectVariable]], sce[[groupVariable]])
  nSubjsPerGroup <- apply(subj_group, 2, function(x) sum(x != 0))
  nGroupsPerSubj <- apply(subj_group, 1, function(x) sum(x != 0))
  if (any(nGroupsPerSubj) > 1) stop("More than one group assigned to
                                       subject(s) in user supplied data")
  if (length(unique(nSubjsPerGroup)) == 1) {
    return(nSubjsPerGroup[[1]])
  } else if (length(nSubjsPerGroup) != nGroups) {
    stop("Different number of subjects in each group in user supplied data
        and number of groups in data does not match twoGroupDesign parameter")
  } else {
    names(nSubjsPerGroup) <- paste0("group", seq_len(nGroups))
    return(nSubjsPerGroup)
  }
}

.estNCellParams <- function(sce, sampleVariable, cellParamsByCondition,
                            groupVariable, timepointVariable,
                            nTimepoints, twoGroupDesign) {

  if (!cellParamsByCondition) {
    if (is.null(sampleVariable)) {
      minCells <- maxCells <- ncol(sce)
    } else {
      num_cell_tab <- table(sce[[sampleVariable]])
      minCells <- round(quantile(num_cell_tab, .1))
      maxCells <- round(quantile(num_cell_tab, .9))
    }
  } else {
    nGroups <- ifelse(twoGroupDesign, 2, 1)
    if (is.null(groupVariable) & is.null(timepointVariable)) {
      stop("To estimate number of cell parameters by condition,
                  groupVariable and/or timepointVariable must be provided")
    } else if (is.null(groupVariable)) {
      groupVariable <- sce[[timepointVariable]]
      num_cell_tab <- table(sce[[sampleVariable]], groupVariable)
      minCells <- apply(num_cell_tab, 2, function(x) round(quantile(x[x != 0], .1)))
      maxCells <- apply(num_cell_tab, 2, function(x) round(quantile(x[x != 0], .9)))
      if (twoGroupDesign) {
        minCells <- rep(minCells, 2)
        maxCells <- rep(maxCells, 2)
      }
    } else if (is.null(timepointVariable)) {
      groupVariable <- sce[[groupVariable]]
      num_cell_tab <- table(sce[[sampleVariable]], groupVariable)
      minCells <- apply(num_cell_tab, 2, function(x) round(quantile(x[x != 0], .1)))
      maxCells <- apply(num_cell_tab, 2, function(x) round(quantile(x[x != 0], .9)))
      if (nTimepoints > 1) {
        minCells <- rep(minCells, each = nTimepoints)
        maxCells <- rep(maxCells, each = nTimepoints)
      }
    } else {
      groupVariable <- paste(sce[[groupVariable]],
        sce[[timepointVariable]],
        sep = "_"
      )
      #groupVariable <- sce[[groupVariable]]
      num_cell_tab <- table(sce[[sampleVariable]], groupVariable)
      minCells <- apply(num_cell_tab, 2, function(x) round(quantile(x[x != 0], .1)))
      maxCells <- apply(num_cell_tab, 2, function(x) round(quantile(x[x != 0], .9)))
    }
    group_df <- expand.grid(
      group = paste0(
        "group",
        seq_len(nGroups)
      ),
      timepoint = paste0("timepoint", seq_len(nTimepoints))
    )
    names(minCells) <- names(maxCells) <- paste(group_df$group,
      group_df$timepoint,
      sep = "_"
    )
  }
  list(minCells = minCells, maxCells = maxCells)
}
