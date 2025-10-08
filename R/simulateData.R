#' Simulate Repeated Measures scRNA-seq Data
#'
#' Simulate Repeated Measures scRNA-seq Data using the RESCUE method
#'
#' @param paramObj \code{\link{RescueSimParams-class}} object with no empty slots.
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
#' \item{\code{rowData}}{
#'   \describe{
#'     \item{deLog2FC}{A \code{DataFrame} containing log2 fold change information for each gene.
#'     Each column corresponds to a non-reference experimental condition
#'     (e.g., \code{"time1"}, \code{"group1"}, \code{"time1_group1"}, etc.), and values represent
#'     gene-level log2 fold changes relative to the baseline condition
#'     (\code{"time0"}, \code{"group0"}, or \code{"time0_group0"} depending on the design).
#'     The reference condition itself is not included as a column.}
#'   }
#' }
#'
#' }
#'
#' @examples
#' # Read in data
#'  data("RecAM_sce")
#'
#'  # Calculate sim parameters for first 50 genes
#' RecAM_sce <- RecAM_sce[1:50,]
#' RecAM_params<-estRescueSimParams(RecAM_sce, sampleVariable = "sampleID",
#' subjectVariable = "subjectID", timepointVariable = "time")
#'
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
    ## Check that it is a RescueSimParams class
    checkmate::assertClass(paramObj, "RescueSimParams")

    ## Make sure no parameter (besides customLibSizes) is equal to zero
    check_empty <- vapply(
        methods::slotNames(paramObj)[methods::slotNames(paramObj)!="customLibSizes"], function(slot) {
            length(getRescueSimParam(paramObj, slot)) == 0
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
        getRescueSimParam(paramObj, "logLibFacVar"),
        colDat$sampleID
    )
    libSizes <- .drawLibSizes(
        libSizeFacs,
        getRescueSimParam(paramObj, "logLibMean"),
        getRescueSimParam(paramObj, "logLibSD"),
        getRescueSimParam(paramObj, "customLibSizes")
    )

    ## Simulate batch factors
    exprsMean <- getRescueSimParam(paramObj, "exprsMean")
    nGenes <- length(exprsMean)
    batchVars <- .drawBatchVars(
        nGenes,
        getRescueSimParam(paramObj, "sampleFacVarMean"),
        getRescueSimParam(paramObj, "sampleFacVarSD"),
        getRescueSimParam(paramObj, "subjectFacVarMean"),
        getRescueSimParam(paramObj, "subjectFacVarSD"),
        nTimepoints=getRescueSimParam(paramObj, "nTimepoints")
    )

    dir_batch <- .drawBatchFacs(
        batchVars, colDat$sampleID, colDat$subjectID,
        getRescueSimParam(paramObj, "nSubjsPerGroup")
    )



    de <- .setDE(
        deLog2FC =getRescueSimParam(paramObj, "deLog2FC"),
        propDE = getRescueSimParam(paramObj, "propDE"),
        nGenes = nGenes,
        colDat = colDat
    )

    ## Simulate True Expression values
    trueExprs <- .simTrueExprs(getRescueSimParam(paramObj, "dispersion"),
                               exprsMean = exprsMean,
                               de = de$de, a = dir_batch$a, b = dir_batch$b
    )

    counts <- .simCounts(trueExprs, libSizes)


    ## Create object
    colDat$time <- paste0("time", colDat$time)
    colDat$group <- paste0("group", colDat$group)
    rownames(colDat) <- paste0("cell", seq_len(nrow(colDat)))

    if(!is.null(de$deFacs)){
        rowDat <- data.frame(de$deFacs, row.names = rownames(counts))
        colnames(rowDat) <- paste0("deLog2FC.", names(de$deFacs))
    }else rowDat=NULL


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
    nSubjsPerGroup <- getRescueSimParam(paramObj, "nSubjsPerGroup")
    nTimepoints <- getRescueSimParam(paramObj, "nTimepoints")
    nGroups <- ifelse(getRescueSimParam(paramObj, "twoGroupDesign"), 2, 1)
    maxCellsPerSamp <- getRescueSimParam(paramObj, "maxCellsPerSamp")
    minCellsPerSamp <- getRescueSimParam(paramObj, "minCellsPerSamp")


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
    } else if (length(minCellsPerSamp) == length(unique(conditions))) {
        names(minCellsPerSamp) <- unique(conditions)
        sampMeta$minCellsPerSamp <- minCellsPerSamp[conditions]
    } else if (length(minCellsPerSamp) == nrow(sampMeta)) {
        names(minCellsPerSamp) <- unique(sampMeta$sampleID)
        sampMeta$minCellsPerSamp <- minCellsPerSamp[sampMeta$sampleID]
    }

    if (length(maxCellsPerSamp) == 1) {
        sampMeta$maxCellsPerSamp <- maxCellsPerSamp
    } else if (length(maxCellsPerSamp) == length(unique(conditions))) {
        names(maxCellsPerSamp) <- unique(conditions)
        sampMeta$maxCellsPerSamp <- maxCellsPerSamp[conditions]
    } else if (length(maxCellsPerSamp) == nrow(sampMeta)) {
        names(maxCellsPerSamp) <- unique(sampMeta$sampleID)
        sampMeta$maxCellsPerSamp <- maxCellsPerSamp[sampMeta$sampleID]
    }

    bad_idx <- which(sampMeta$minCellsPerSamp > sampMeta$maxCellsPerSamp)
    if (length(bad_idx) > 0) {
        warning(sprintf(
            "minCellsPerSamp was greater than maxCellsPerSamp for %d sample(s): %s. Values were switched.",
            length(bad_idx),
            paste(sampMeta$sampleID[bad_idx], collapse = ", ")
        ))
    }

    nCells <- apply(sampMeta, 1, function(x) {
        minCells <- as.numeric(x["minCellsPerSamp"])
        maxCells <- as.numeric(x["maxCellsPerSamp"])
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

.drawLibSizes <- function(lib_facs, mu, sd, customLibSizes) {
    if(length(customLibSizes)>0){
        customLogLibSizes<-log(customLibSizes)
        mu=mean(customLogLibSizes)
        customMean=mean(customLogLibSizes)
        rawLogLibSizes<-sample(customLogLibSizes, length(lib_facs), replace = T)

        ## Adjusting so that there is a location shift for sample-specific mu
        ## sample specific mu will equal mu*lib_fac
        ## lib_fac is sample specific factor
        adjLogLibSizes <- rawLogLibSizes + mu * (lib_facs - 1)
        lib_sizes=as.numeric(exp(adjLogLibSizes))
    }else{
        lib_sizes <- stats::rlnorm(length(lib_facs), mu * lib_facs, sd)
    }
    return(lib_sizes)
}


.drawBatchVars <- function(n_genes,
                           sampleFacVarMean, sampleFacVarSD,
                           subjectFacVarMean, subjectFacVarSD, nTimepoints) {
    if(nTimepoints==1){
        samp_var=NA
    }else{
        samp_var <- stats::rlnorm(n_genes, sampleFacVarMean, sampleFacVarSD)
    }
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

    a <- vapply(alpha_subj, function(alpha) {
        x <- gtools::rdirichlet(1, rep(alpha, n_subjs)) * n_subjs
        if (any(is.na(x))) {
            x <- rep(0, n_subjs)
            x[sample(seq_len(n_subjs), 1)] <- 1
            x <- x * n_subjs
        }
        return(x)
    }, FUN.VALUE = numeric(n_subjs))

    ## If only one timepoint, then all values of b should be 1
    if(nTimepoints==1){
        b <- matrix(1,nrow = nrow(batchVars), ncol = nrow(samp_subj_df))
        colnames(b) <- samp_subj_df$sampleID
    }else{
        alpha_samp <- vapply(batchVars$samp_var, .calcAlpha,
                             n = nTimepoints,
                             FUN.VALUE = numeric(1)
        )
        alpha_samp <- ifelse(alpha_samp <= 0, min(alpha_samp[alpha_samp > 0]), alpha_samp)

        b <- matrix(nrow = length(alpha_samp), ncol = nrow(samp_subj_df))
        colnames(b) <- samp_subj_df$sampleID
        for (i in seq_len(length(unique(subjectID)))) {
            subject <- unique(subjectID)[i]
            subj_sampleID <- samp_subj_df$sampleID[samp_subj_df$subjectID == subject]
            nSampsPerSubj <- length(subj_sampleID)
            b_temp <- vapply(alpha_samp, function(alpha) {
                x <- gtools::rdirichlet(1, rep(alpha, nSampsPerSubj)) * nSampsPerSubj
                if (any(is.na(x))) {
                    x <- rep(0, nSampsPerSubj)
                    x[sample(seq_len(nSampsPerSubj), 1)] <- 1
                    x <- x * nSampsPerSubj
                }
                return(x)
            }, FUN.VALUE = numeric(nSampsPerSubj))
            b[, subj_sampleID] <- t(b_temp)
        }
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

.setDE <- function(deLog2FC, propDE, nGenes, colDat) {
    time <- colDat$time
    group <- colDat$group
    nGroups <- length(unique(group))
    nTimepoints <- length(unique(time))
    nCells <- nrow(colDat)

    # CASE 1: Numeric input (scalar or vector)
    if (is.numeric(deLog2FC)) {

        # Handle case where there's no DE to model
        if (nTimepoints == 1 && nGroups == 1) {
            warning("No differential expression modeled: only one timepoint and one group provided.")
            de <- matrix(1, nrow = nGenes, ncol = nCells)
            return(list(de = de, deFacs = NULL))
        }

        if (nGroups == 1) group <- group + 1
        if (nTimepoints == 1) time <- time + 1

        # For design factor scaling
        denom <- max(nTimepoints - 1, 1)
        design_factor <- (time / denom) * group

        # Sample DE genes
        de_genes <- sample(c(TRUE, FALSE), nGenes,
                           replace = TRUE,
                           prob = c(propDE, 1 - propDE))

        # Symmetric Log2FCs if scalar
        if (length(deLog2FC) == 1) {
            deLog2FC <- c(-deLog2FC, deLog2FC)
        }

        # Assign Log2FCs to DE genes
        de_facs <- rep(0, nGenes)
        de_facs[de_genes] <- sample(deLog2FC, sum(de_genes), replace = TRUE)

        de_log <- matrix(design_factor, nrow = nGenes, ncol = nCells, byrow = TRUE)
        de_log <- de_log * de_facs

        # Generate condition names
        colNames <- character(nCells)
        if (nTimepoints == 1 && nGroups == 2) {
            colNames <- paste0("group", group)
            condition_names <- "group1"
        } else if (nGroups == 1 && nTimepoints > 1) {
            colNames <- paste0("time", time)
            condition_names <- paste0("time", 1:(nTimepoints - 1))
        } else {
            colNames <- paste0("time", time, "_group", group)
            condition_names <- unlist(lapply(0:(nTimepoints - 1), function(t) {
                paste0("time", t, "_group", 0:1)
            }))
            ## Remove reference
            condition_names<-condition_names[condition_names!="time0_group0"]
        }

        # Build deLog2FC list
        deLog2FC_list <- lapply(condition_names, function(cname) {
            idx <- which(colNames == cname)
            if (length(idx) == 0) {
                rep(0, nGenes)
            } else {
                rowMeans(de_log[, idx, drop = FALSE])
            }
        })
        names(deLog2FC_list) <- condition_names

        de <- 2^de_log
        return(list(de = de, deFacs = deLog2FC_list))
    }

    # CASE 2: List of vectors with gene-specific Log2FCs
    if (is.list(deLog2FC)) {
        colNames <- character(nCells)

        if (nTimepoints == 1 && nGroups == 2) {
            colNames <- paste0("group", group)
        } else if (nGroups == 1 && nTimepoints > 1) {
            colNames <- paste0("time", time)
        } else if (nGroups == 2 && nTimepoints > 1) {
            colNames <- paste0("time", time, "_group", group)
        }

        # Build DE matrix from log2FCs
        de_log <- matrix(0, nrow = nGenes, ncol = nCells)
        for (cond_name in names(deLog2FC)) {
            de_log[, colNames == cond_name] <- deLog2FC[[cond_name]]
        }

        de <- 2^de_log
        return(list(de = de, deFacs = deLog2FC))
    }

    stop("Invalid format for deLog2FC.")
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
