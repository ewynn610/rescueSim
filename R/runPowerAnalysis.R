#' Run Power Simulations for rescueSim Data
#'
#' Simulate scRNA-seq data under different experimental settings and compute
#' statistical power to detect differential expression (DE) using a user-specified DE function.
#'
#' @param baseParams A \code{RescueSimParams} object containing
#' baseline simulation parameters.
#' @param scenarios A data.frame specifying simulation settings to vary.
#' Column names must match slots in the \code{RescueSimParams} object.
#' Each row defines one scenario.
#' @param deFunction A function that takes a \code{SingleCellExperiment} object
#' simulated using \code{\link{simRescueData}} and returns a data.frame with
#' columns \code{gene} and \code{padj}.
#' @param nSim Integer specifying the number of simulations to run per scenario.
#'  Default is 1.
#' @param padjThresh Significance threshold to call DE (e.g., FDR < 0.05).
#' Default is 0.05.
#' @param returnFDR Logical. Whether to calculate and return the false discovery rate (FDR) along with power.
#' Default is TRUE.
#' @param conditions A character vector of length 2 specifying the conditions to
#'  compare (e.g., \code{c("time1", "time3_group1")}). If \code{NULL} (default),
#'  the comparison is inferred based on available conditions in \code{rowData},
#'  comparing the baseline condition (e.g., \code{"time0"}, \code{"group0"}, or
#'  \code{"time0_group0"}) to the final condition (e.g., last timepoint, or
#'  group 1 at the final timepoint).
#' @param saveSimPath Optional path to directory where simulated SCE objects
#' should be saved as .rds files. If NULL (default), simulated data is not saved.
#' @param saveDePath Optional path to directory where DE result data.frames
#' should be saved as .rds files. If NULL (default), results are not saved.
#' @param verbose Logical. Whether to print progress messages. Default is TRUE.
#' @param ... Additional arguments passed to \code{deFunction()}.
#'
#' @details
#' For each scenario, this function simulates \code{nSim} datasets based on the,
#' supplied settings, applies the user-supplied DE method, and
#' calculates power as the proportion of true DE genes correctly identified.
#' DE status is determined by assessing equality in the  simulated log2 fold-change values
#' between conditions. False discovery rate (FDR) is optionally
#' calculated as the proportion of genes called DE that are not truly DE.
#'
#' @return A data.frame with one row per simulation replicate, including:
#' \itemize{
#'   \item Power: proportion of truly DE genes detected
#'   \item FDR: false discovery rate (if \code{returnFDR = TRUE})
#'   \item Scenario settings and simulation number
#'   \item The reference and comparison condition used
#' }
#'
#' Optionally, simulated datasets and DE results can be saved to disk.
#'
#' @seealso \code{\link{isDEbetweenConditions}}, \code{\link{simRescueData}},
#' \code{\link{updateRescueSimParams}}, \code{\link{RescueSimParams}}
#'
#' @export

runRescueSimPower <- function(baseParams, scenarios,
                              deFunction,
                              nSim = 1,
                              padjThresh = 0.05,
                              returnFDR = FALSE,
                              conditions=NULL,
                              saveSimPath = NULL,
                              saveDePath = NULL,
                              verbose = TRUE,
                              ...) {
    ## Check if conditions are formated correctly
    if (!is.null(conditions)) {
        if (!is.character(conditions) || length(conditions) != 2) {
            stop("`conditions` must be a character vector of length 2.")
        }
    }

    ## Check if scenario names are correct:
    valid_params <- methods::slotNames(baseParams)

    invalid <- setdiff(colnames(scenarios), valid_params)
    if (length(invalid) > 0) {
        stop(sprintf(
            "The following scenario parameters are not valid RescueSimParams fields: %s",
            paste(invalid, collapse = ", ")
        ))
    }


    results <- list()
    save_sce <- !is.null(saveSimPath)
    save_de <- !is.null(saveDePath)

    if (save_sce && !dir.exists(saveSimPath)) dir.create(saveSimPath, recursive = TRUE)
    if (save_de && !dir.exists(saveDePath)) dir.create(saveDePath, recursive = TRUE)

    for (i in seq_len(nrow(scenarios))) {
        scenario <- scenarios[i, ]
        if (verbose) message("Running scenario ", i, " of ", nrow(scenarios))

        paramList <- lapply(scenario, function(x) if (is.list(x)) x[[1]] else x)

        updatedParams <- tryCatch({
            updateRescueSimParams(baseParams, paramList)
        }, error = function(e) {
            warning(sprintf("Scenario %d failed: %s", i, e$message))
            return(NULL)
        })

        if (is.null(updatedParams)) next

        for (sim in seq_len(nSim)) {
            if (verbose) message("  Sim ", sim)

            sce <- simRescueData(updatedParams)

            rd_cols <- colnames(SummarizedExperiment::rowData(sce))
            expected_conds <- gsub("^deLog2FC\\.","" ,rd_cols)

            if(!is.null(conditions)){

                missing_cols <- setdiff(conditions, expected_conds)


                if (length(missing_cols) > 0) {
                    stop(
                        sprintf(
                            "The specified condition(s) %s not found in rowData columns. Expected conditions: %s",
                            paste(conditions, collapse = ", "),
                            paste(expected_conds, collapse = ", ")
                        )
                    )
                }
                cond1=conditions[1]
                cond2=conditions[2]
            }else{
                if (all(grepl("^time\\d+$", expected_conds))) {
                    cond1 <- "time0"
                    cond2 <- paste0("time", max(as.integer(gsub("time", "", expected_conds))))
                } else if (all(grepl("^group\\d+$", expected_conds))) {
                    cond1 <- "group0"
                    cond2 <- "group1"
                } else if (all(grepl("^time\\d+_group\\d+$", expected_conds))) {
                    cond1 <- "time0_group0"
                    split_cond <- do.call(rbind, strcapture("time(\\d+)_group(\\d+)",
                                                            expected_conds,
                                                            list(time = integer(), group = integer())))
                    max_time <- max(split_cond$time)
                    cond2 <- paste0("time", max_time, "_group1")
                } else {
                    stop("Could not determine comparison conditions from rowData.")
                }
                myConditions=c(cond1, cond2)

            }



            ## Run DE
            de_res <- tryCatch({
                deFunction(sce, ...)
            }, error = function(e) {
                warning(sprintf("DE failed for scenario %d sim %d: %s", i, sim, e$message))
                return(NULL)
            })

            if (is.null(de_res)) next

            if (!all(c("gene", "padj") %in% colnames(de_res))) {
                stop("deFunction must return a data.frame with columns: gene, padj")
            }

            if (anyNA(de_res$padj)) {
                warning("Some padj values are NA and will be removed before calculating power.")
                de_res <- de_res[!is.na(de_res$padj), , drop = FALSE]
            }

            ## Get power etc.
            is_de <- isDEbetweenConditions(sce, myConditions[1], myConditions[2])


            matched <- intersect(de_res$gene, names(is_de))
            called_de <- de_res$padj[match(matched, de_res$gene)] < padjThresh
            truth <- is_de[matched]

            power <- mean(truth[called_de])

            fdr <- if (returnFDR) {
                n_fp <- sum(!truth & called_de)
                n_called <- sum(called_de)
                if (n_called == 0) 0 else n_fp / n_called
            }

            # Save .rds files if desired
            sce_file <- NULL
            de_file <- NULL

            if (save_sce) {
                sce_file <- file.path(saveSimPath, sprintf("sce_scenario%d_sim%d.rds", i, sim))
                saveRDS(sce, sce_file)
            }

            if (save_de) {
                de_file <- file.path(saveDePath, sprintf("de_scenario%d_sim%d.rds", i, sim))
                saveRDS(de_res, de_file)
            }

            result_row <- cbind(
                as.data.frame(scenario),
                sim = sim,
                condition1 = cond1,
                condition2 = cond2,
                power = power)

            if(returnFDR){
                result_row$fdr<-fdr
            }

            results[[length(results) + 1]] <- result_row


            rm(sce, de_res)
            gc()
        }
    }

    out <- do.call(rbind, results)
    rownames(out) <- NULL
    return(out)
}

#' Identify Genes Simulated to Be Differentially Expressed Between Conditions
#'
#' Returns a logical vector indicating whether each gene was simulated to be
#' differentially expressed between two conditions in a
#' \code{SingleCellExperiment} generated by the \pkg{rescueSim} package.
#'
#' @details
#' All log2 fold changes are defined relative to a baseline or reference
#' condition—typically \code{"time0"}, \code{"group0"}, or
#' \code{"time0_group0"}—and stored in the \code{rowData} with column names of
#' the form \code{"deLog2FC.<condition>"}. Reference conditions are not included
#' in the \code{rowData}.
#'
#' If one of the specified conditions is a reference, the function checks whether
#' the corresponding log2 fold change column is non-zero. If both are
#' non-reference, it checks whether the two log2 fold change columns differ.
#'
#' This function does not perform a statistical test; it simply reports whether
#' each gene was \emph{simulated} to be differentially expressed between the two
#' conditions based on the stored simulation parameters.
#'
#' @param sce A \code{SingleCellExperiment} object simulated using
#' \code{rescueSim::simRescueData()}, with log2 fold change annotations stored in
#' the \code{rowData}.
#' @param condition1 A character string specifying the first condition. This can
#' be a baseline reference (e.g. \code{"time0"}, \code{"group0"}, or
#' \code{"time0_group0"}) or a condition present in \code{rowData} (without the
#' \code{"deLog2FC."} prefix).
#' @param condition2 A character string specifying the second condition,
#' interpreted the same way as \code{condition1}.
#'
#' @return A named logical vector indicating whether each gene was simulated to
#' be differentially expressed between \code{condition1} and \code{condition2}.
#' Names correspond to gene names.
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
#' simDat=simRescueDataa(RecAM_params)
#'
#' is_de <- isDEbetweenConditions(simDat, "time1", "time0")
#'
#' @export



isDEbetweenConditions <- function(sce, condition1, condition2) {
    rd <- SummarizedExperiment::rowData(sce)
    rd_names <- colnames(rd)

    # Helper to check if a condition is the reference
    is_reference <- function(cond) {
        cond %in% c("time0", "group0", "time0_group0")
    }

    # Extract actual column names (they include a prefix like 'deLog2FC.')
    full_cond1 <- paste0("deLog2FC.", condition1)
    full_cond2 <- paste0("deLog2FC.", condition2)

    cond1_in_rd <- full_cond1 %in% rd_names
    cond2_in_rd <- full_cond2 %in% rd_names

    if (!cond1_in_rd && !is_reference(condition1)) {
        stop(paste0("Condition '", condition1, "' not found in rowData and is not a recognized reference condition."))
    }

    if (!cond2_in_rd && !is_reference(condition2)) {
        stop(paste0("Condition '", condition2, "' not found in rowData and is not a recognized reference condition."))
    }

    # Determine DE status
    is_de <- NULL
    if (!cond1_in_rd && cond2_in_rd) {
        is_de <- rd[[full_cond2]] != 0
    } else if (cond1_in_rd && !cond2_in_rd) {
        is_de <- rd[[full_cond1]] != 0
    } else if (cond1_in_rd && cond2_in_rd) {
        is_de <- rd[[full_cond1]] != rd[[full_cond2]]
    } else {
        stop("Both conditions appear to be reference groups or invalid. Cannot determine DE status.")
    }

    # Add gene names
    names(is_de) <- rownames(sce)
    return(is_de)
}

