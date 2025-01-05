#' Simulate Repeated Measures scRNA-seq Data
#'
#' Simulate Repeated Measures scRNA-seq Data using the RESCUE method
#'
#' @param paramObj \code{\link{RescueParams-class}} object with no empty slots.
#'
#' @return \code{\link[SingleCellExperiment]{SingleCellExperiment}} object with
#' the data in the following slots
#'
#' \describe{
#' \item{\code{counts}}{Matrix of raw counts with genes represented by rows and
#' cells represented by columns.}
#'     \item{\code{colData}}{
#'         \describe{
#'             \item{sampleID}{Sample identifier}
#'             \item{subjectID}{Subject identifier}
#'             \item{time}{Timepoint identifier}
#'             \item{group}{Group identifier}
#'         }
#'     }
#'     \item{\code{rowData}}{
#'         \describe{
#'             \item{deLogFC}{Log2 fold change between first and last timepoints
#'             or groups, or, if multi-timepoint and multi-group data, Log2 fold
#'              change between the last timepoint and group vs. other
#'              timepoint/group combinations}
#'         }
#'     }
#'
#' }
#'
#' @examples
#' # Read in parameter object
#' data("RecAM_params")
#'
#' # Simulate data
#' simDat=simRescueData(RecAM_params)
#'
#' # Examine
#' library(SingleCellExperiment)
#' counts(simDat)[1:5, 1:5]
#'
#' @export

simRescueData <- function(paramObj) {
  ## Check that it is a RescueParams class
  checkmate::assertClass(paramObj, "RescueParams")

  ## Make sure no parameter is equal to zero
  check_empty <- vapply(
    methods::slotNames(paramObj), function(slot) {
      length(getRescueParam(paramObj, slot)) == 0
    },
    logical(1)
  )
  if (any(check_empty)) {
    stop("One or more parameter slots are empty. All parameter slots must be filled")
  }

  ## Generate basic meta data
  colDat <- .getBasicMeta(paramObj)

  ## Get Library Sizes
  libSizeFacs <- .drawLibFacs(
    getRescueParam(paramObj, "logLibFacVar"),
    colDat$sampleID
  )
  libSizes <- .drawLibSizes(
    libSizeFacs,
    getRescueParam(paramObj, "logLibMean"),
    getRescueParam(paramObj, "logLibSD")
  )

  ## Simulate batch factors
  exprsMean <- getRescueParam(paramObj, "exprsMean")
  nGenes <- length(exprsMean)
  batchVars <- .drawBatchVars(
    nGenes,
    getRescueParam(paramObj, "sampleFacVarMean"),
    getRescueParam(paramObj, "sampleFacVarSD"),
    getRescueParam(paramObj, "subjectFacVarMean"),
    getRescueParam(paramObj, "subjectFacVarSD")
  )

  dir_batch <- .drawBatchFacs(
    batchVars, colDat$sampleID, colDat$subjectID,
    getRescueParam(paramObj, "nSubjsPerGroup")
  )



  de <- .setDE(
    getRescueParam(paramObj, "deLogFC"),
    getRescueParam(paramObj, "propDE"),
    nGenes,
    colDat
  )

  ## Simulate True Expression values
  trueExprs <- .simTrueExprs(getRescueParam(paramObj, "dispersion"),
    exprsMean = exprsMean,
    de = de$de, a = dir_batch$a, b = dir_batch$b
  )

  counts <- .simCounts(trueExprs, libSizes)


  ## Create object
  colDat$time <- paste0("time", colDat$time)
  colDat$group <- paste0("group", colDat$group)
  rownames(colDat) <- paste0("cell", seq_len(nrow(colDat)))

  rowDat <- data.frame(deLogFC = de$deFac, row.names = rownames(counts))


  rownames(colDat) <- colnames(counts)

  sce <- SingleCellExperiment::SingleCellExperiment(list(counts = counts),
    colData = colDat,
    rowData = rowDat
  )
  ## Order by subjectID and sampleID
  sce=sce[,order(colDat[,"subjectID"],colDat[,"sampleID"] )]

  sce

}




.getBasicMeta <- function(paramObj) {
  ## Extract paramObj
  nSubjsPerGroup <- getRescueParam(paramObj, "nSubjsPerGroup")
  nTimepoints <- getRescueParam(paramObj, "nTimepoints")
  nGroups <- ifelse(getRescueParam(paramObj, "twoGroupDesign"), 2, 1)
  maxCellsPerSamp <- getRescueParam(paramObj, "maxCellsPerSamp")
  minCellsPerSamp <- getRescueParam(paramObj, "minCellsPerSamp")


  ## Get nsamps
  nsamps <- nSubjsPerGroup * nTimepoints * nGroups

  ## Generate sample ID
  sampMeta <- data.frame(sampleID = paste0("sample", seq_len(nsamps)))


  ## Generate subject ID
  sampMeta$subjectID <- rep(paste0("subject", seq_len(nSubjsPerGroup * nGroups)), nTimepoints)

  ## Generat time
  sampMeta$time <- rep(seq.int(0, nTimepoints - 1), each = nSubjsPerGroup * nGroups)

  ## Generate Group
  sampMeta$group <- rep(rep(seq.int(0, nGroups - 1), each = nSubjsPerGroup), nTimepoints)


  conditions <- paste0("group", sampMeta$group, "_time", sampMeta$time)

  if (length(minCellsPerSamp) == 1) {
    sampMeta$minCellsPerSamp <- minCellsPerSamp
    sampMeta$maxCellsPerSamp <- maxCellsPerSamp
  } else if (length(minCellsPerSamp) == length(unique(conditions))) {
    names(minCellsPerSamp) <- unique(conditions)
    names(maxCellsPerSamp) <- unique(conditions)
    sampMeta$minCellsPerSamp <- minCellsPerSamp[conditions]
    sampMeta$maxCellsPerSamp <- maxCellsPerSamp[conditions]
  } else if (length(minCellsPerSamp) == nrow(sampMeta)) {
    names(minCellsPerSamp) <- unique(sampMeta$sampleID)
    names(maxCellsPerSamp) <- unique(sampMeta$sampleID)
    sampMeta$minCellsPerSamp <- minCellsPerSamp[sampMeta$sampleID]
    sampMeta$maxCellsPerSamp <- maxCellsPerSamp[sampMeta$sampleID]
  }
  nCells <- apply(sampMeta, 1, function(x) {
    minCells <- x["minCellsPerSamp"]
    maxCells <- x["maxCellsPerSamp"]
    if (minCells == maxCells) {
      return(minCells)
    } else {
      return(sample(seq.int(minCells, maxCells), 1))
    }
  })

  idx <- unlist(lapply(1:nrow(sampMeta), function(x) {
    rep(x, nCells[x])
  }))

  sampMeta[idx, c("sampleID", "subjectID", "time", "group")]
}


.drawLibFacs <- function(logLibFacVar, sampleID) {
  ## Find alpha
  n <- length(unique(sampleID))

  if (logLibFacVar == 0) {
    ## If variance is 0, all batch factors are 1
    batch_facs <- rep(1, n)
  } else {
    alpha <- .calcAlpha(logLibFacVar, n)

    ## For every sampleID, draw library size factor
    batch_facs <- gtools::rdirichlet(1, rep(alpha, n)) * n
  }

  names(batch_facs) <- unique(sampleID)
  fac_vec <- batch_facs[sampleID]
}

.drawLibSizes <- function(lib_facs, mu, sd) {
  lib_sizes <- stats::rlnorm(length(lib_facs), mu * lib_facs, sd)
}


.drawBatchVars <- function(n_genes,
                           sampleFacVarMean, sampleFacVarSD,
                           subjectFacVarMean, subjectFacVarSD) {
  samp_var <- stats::rlnorm(n_genes, sampleFacVarMean, sampleFacVarSD)
  subj_var <- stats::rlnorm(n_genes, subjectFacVarMean, subjectFacVarSD)
  ret <- data.frame(samp_var = samp_var, subj_var = subj_var)
}

############ Check about having uneqaul samps. per subj
.drawBatchFacs <- function(batchVars, sampleID, subjectID,
                           nSubjsPerGroup) {
  samp_subj_df <- data.frame(
    sampleID = unique(sampleID),
    subjectID = subjectID[!duplicated(sampleID)]
  )
  nTimepoints <- mean(table(samp_subj_df$subjectID))

  n_subjs <- length(unique(subjectID))

  alpha_subj <- vapply(batchVars$subj_var, .calcAlpha,
    n = n_subjs,
    FUN.VALUE = numeric(1)
  )
  alpha_subj <- ifelse(alpha_subj <= 0, min(alpha_subj[alpha_subj > 0]), alpha_subj)
  alpha_samp <- vapply(batchVars$samp_var, .calcAlpha,
    n = nTimepoints,
    FUN.VALUE = numeric(1)
  )
  alpha_samp <- ifelse(alpha_samp <= 0, min(alpha_samp[alpha_samp > 0]), alpha_samp)
  a <- vapply(alpha_subj, function(alpha) {
    x <- gtools::rdirichlet(1, rep(alpha, n_subjs)) * n_subjs
    if (any(is.na(x))) {
      x <- rep(0, n_subjs)
      x[sample(seq_len(n_subjs), 1)] <- 1
      x <- x * n_subjs
    }
    return(x)
  }, FUN.VALUE = numeric(n_subjs))

  b <- matrix(nrow = length(alpha_samp), ncol = nrow(samp_subj_df))
  colnames(b) <- samp_subj_df$sampleID
  for (i in seq_len(length(unique(subjectID)))) {
    subject <- unique(subjectID)[i]
    subj_sampleID <- samp_subj_df$sampleID[samp_subj_df$subjectID == subject]
    nSampsPerSubj <- length(subj_sampleID)
    b_temp <- vapply(alpha_samp, function(alpha) {
      x <- gtools::rdirichlet(1, rep(alpha, nSampsPerSubj)) * nSampsPerSubj
      if (any(is.na(x))) {
        x <- rep(0, nSubjsPerGroup)
        x[sample(seq_len(nSubjsPerGroup), 1)] <- 1
        x <- x * nSubjsPerGroup
      }
      return(x)
    }, FUN.VALUE = numeric(nSampsPerSubj))
    b[, subj_sampleID] <- t(b_temp)
  }



  a <- t(a)
  colnames(a) <- unique(subjectID)
  a <- a[, subjectID]

  b <- b[, sampleID]

  return(list(a = a, b = b))
}

## Calculate alpha for Dirichlet distribution
.calcAlpha <- function(v, n) {
  (n - 1 - v) / (n * v)
}

.setDE <- function(deLogFC, propDE, nGenes, colDat) {
  time <- colDat$time
  group <- colDat$group
  nGroups <- length(unique(group))
  nTimepoints <- length(unique(time))

  if (nGroups == 1) group <- group + 1
  if (nTimepoints == 1) {
    time <- time + 1
    ## Just changing this so that don't get 0 in denominator below
    nTimepoints <- 2
  }

  de_genes <- sample(c(TRUE, FALSE), nGenes,
    replace = TRUE,
    prob = c(propDE, 1 - propDE)
  )


  if (length(deLogFC) == 1) {
    deLogFC <- c(-deLogFC, deLogFC)
  }

  de_facs <- rep(0, nGenes)
  de_facs[de_genes] <- sample(deLogFC, sum(de_genes), replace = TRUE)


  de <- matrix((time / (nTimepoints - 1)) * group,
    nrow = nGenes,
    ncol = length(time), byrow = TRUE
  )

  de <- de * de_facs
  de <- 2^de
  return(list(de = de, deFacs = de_facs))
}

## Changed to round(ncells) so it will work for really big numbers
.simTrueExprs=function (dispersion, exprsMean, de, a, b)
{
  n_cells <- ncol(de)
  x_new <- a * b * de * exprsMean
  ngenes <- nrow(x_new)
  ncells <- ncol(x_new)
  mat <- matrix(stats::rgamma(ngenes * round(ncells), shape = 1/dispersion,
                              scale = x_new * dispersion), nrow = ngenes, ncol = ncells)
  rownames(mat) <- names(exprsMean)
  idx_invar <- which(dispersion == 0)
  for (idx in idx_invar) {
    mat[idx, ] <- exprsMean[idx]
  }
  mat
}

.simCounts <- function(trueExprs, libSizes) {
  ncells <- ncol(trueExprs)
  ngenes <- nrow(trueExprs)

  ## Adjusting for library size
  trueExprsCorrected <- t(libSizes * t(trueExprs / colSums(trueExprs)))

  ## For each cell, draw count from Poisson with mean trueExprs
  counts <- matrix(stats::rpois(round(ncells) * ngenes, lambda = trueExprsCorrected),
    nrow = ngenes, ncol = ncells
  )
  if (is.null(rownames(trueExprsCorrected))) {
    row.names(counts) <- paste0("gene_", seq_len(nrow(counts)))
  } else {
    row.names(counts) <- row.names(trueExprsCorrected)
  }
  colnames(counts) <- paste0("cell_", seq_len(ncol(counts)))
  counts
}
