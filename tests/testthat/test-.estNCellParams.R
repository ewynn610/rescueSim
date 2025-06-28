library(SingleCellExperiment)

# Helper to create dummy SCE
create_dummy_sce <- function(samples, groups = NULL, timepoints = NULL) {
    n_cells <- length(samples)
    sce <- SingleCellExperiment(assays = list(counts = matrix(1, nrow = 10, ncol = n_cells)))
    colData(sce)$sampleID <- samples
    if (!is.null(groups)) colData(sce)$group <- groups
    if (!is.null(timepoints)) colData(sce)$time <- timepoints
    sce
}

test_that(".estNCellParams works when no cellParamsByCondition", {
    sce <- create_dummy_sce(samples = rep("A", 10))
    out <- rescueSim:::.estNCellParams(sce, sampleVariable = "sampleID",
                           cellParamsByCondition = FALSE,
                           groupVariable = NULL, timepointVariable = NULL,
                           nTimepoints = 1, twoGroupDesign = FALSE)
    expect_named(out, c("minCells", "maxCells"))
    expect_true(is.numeric(out$minCells))
    expect_true(is.numeric(out$maxCells))
})

test_that(".estNCellParams errors on nTimepoints mismatch", {
    sce <- create_dummy_sce(
        samples = rep(c("A", "B"), each = 5),
        timepoints = c(rep("t0", 5), rep("t1", 5))
    )
    expect_error(
        rescueSim:::.estNCellParams(sce, sampleVariable = "sampleID",
                        cellParamsByCondition = TRUE,
                        groupVariable = NULL, timepointVariable = "time",
                        nTimepoints = 3, twoGroupDesign = FALSE),
        "Number of unique timepoints in data \\(2\\) does not match nTimepoints \\(3\\)"
    )
})

test_that(".estNCellParams errors on twoGroupDesign mismatch (too many groups)", {
    sce <- create_dummy_sce(
        samples = rep(c("A", "B", "C"), each = 3),
        groups = rep(c("g0", "g1", "g2"), each = 3)
    )
    expect_error(
        rescueSim:::.estNCellParams(sce, sampleVariable = "sampleID",
                        cellParamsByCondition = TRUE,
                        groupVariable = "group", timepointVariable = NULL,
                        nTimepoints = 1, twoGroupDesign = TRUE),
        "twoGroupDesign = TRUE but data has 3 unique groups"
    )
})

test_that(".estNCellParams errors on twoGroupDesign mismatch (too few groups)", {
    sce <- create_dummy_sce(
        samples = rep(c("A", "B"), each = 3),
        groups = rep("g0", 6)
    )
    expect_error(
        rescueSim:::.estNCellParams(sce, sampleVariable = "sampleID",
                        cellParamsByCondition = TRUE,
                        groupVariable = "group", timepointVariable = NULL,
                        nTimepoints = 1, twoGroupDesign = TRUE),
        "twoGroupDesign = TRUE but data has 1 unique groups"
    )
})

test_that(".estNCellParams produces names for min/maxCells", {
    sce <- create_dummy_sce(
        samples = rep(c("A", "B", "C"), each = 3),
        groups = c(rep(c("g0"), 9)),
        timepoints = rep(c("t0", "t1", "t2"), each = 3)
    )
    out <- rescueSim:::.estNCellParams(sce, sampleVariable = "sampleID",
                           cellParamsByCondition = TRUE,
                           groupVariable = "group", timepointVariable = "time",
                           nTimepoints = 3, twoGroupDesign = F)
    expect_named(out, c("minCells", "maxCells"))
    expect_true(all(grepl("time", names(out$minCells))))
})
