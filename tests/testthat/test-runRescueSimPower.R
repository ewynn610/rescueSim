## Load data
data("RecAM_sce")
RecAM_sce <- RecAM_sce[1:50,]
baseParams<-estRescueSimParams(RecAM_sce, sampleVariable = "sampleID",
                              subjectVariable = "subjectID",
                              timepointVariable = "time")

test_that("runRescueSimPower returns valid output structure", {
    # Example minimal scenario
    scenarios <- data.frame(
        minCellsPerSamp = 100,
        maxCellsPerSamp=100,
        nSubjsPerGroup = 3
    )

    # Example trivial DE function
    fake_de_fun <- function(sce, ...) {
        data.frame(
            gene = rownames(sce),
            padj = runif(nrow(sce))
        )
    }

    # Run
    out <- runRescueSimPower(
        baseParams = baseParams,
        scenarios = scenarios,
        deFunction = fake_de_fun,
        nSim = 1,
        verbose = FALSE
    )

    expect_s3_class(out, "data.frame")
    expect_true(all(c("minCellsPerSamp", "maxCellsPerSamp", "nSubjsPerGroup",
                      "sim", "condition1", "condition2", "power") %in%
                        colnames(out)))
})

test_that("runRescueSimPower errors if invalid scenario parameter", {
    scenarios <- data.frame(
        badParam = 1
    )

    expect_error(
        runRescueSimPower(baseParams, scenarios, function(sce) data.frame(gene=character(), padj=numeric())),
        "not valid RescueSimParams fields"
    )
})

test_that("runRescueSimPower errors if conditions badly specified", {
    scenarios <- data.frame(minCellsPerSamp = 100)

    expect_error(
        runRescueSimPower(baseParams, scenarios, function(sce)
            data.frame(gene=character(), padj=numeric()),
                          conditions = "only_one_condition"),
        "`conditions` must be a character vector of length 2"
    )
})

test_that("runRescueSimPower handles NA padj values and issues warning", {
    scenarios <- data.frame(minCellsPerSamp = 100)

    de_fun_na <- function(sce, ...) {
        data.frame(
            gene = rownames(sce),
            padj = c(NA, runif(nrow(sce) - 1))
        )
    }

    expect_warning(
        out <- runRescueSimPower(baseParams, scenarios, de_fun_na,
                                 nSim = 1, verbose = FALSE),
        "Some padj values are NA"
    )

    expect_s3_class(out, "data.frame")
})

test_that("runRescueSimPower optionally computes FDR", {
    scenarios <- data.frame(minCellsPerSamp = 100)

    fake_de_fun <- function(sce, ...) {
        data.frame(
            gene = rownames(sce),
            padj = runif(nrow(sce))
        )
    }

    out <- runRescueSimPower(
        baseParams, scenarios, fake_de_fun, nSim = 1,
        returnFDR = TRUE, verbose = FALSE
    )

    expect_true("fdr" %in% colnames(out))
})

test_that("runRescueSimPower saves files correctly", {
    # Create temp directories
    tmp_sce_dir <- tempfile("sce_dir_")
    tmp_de_dir <- tempfile("de_dir_")
    dir.create(tmp_sce_dir)
    dir.create(tmp_de_dir)

    # Minimal scenario
    scenarios <- data.frame(minCellsPerSamp = 10, maxCellsPerSamp = 10, nSubjsPerGroup = 2)

    # Minimal deFunction that returns required columns
    de_fun <- function(sce) {
        data.frame(
            gene = rownames(sce),
            padj = runif(nrow(sce))
        )
    }

    # Run the function
    out <- runRescueSimPower(
        baseParams = baseParams,
        scenarios = scenarios,
        deFunction = de_fun,
        nSim = 1,
        saveSimPath = tmp_sce_dir,
        saveDePath = tmp_de_dir,
        verbose = FALSE
    )

    # Check files exist
    sce_files <- list.files(tmp_sce_dir, pattern = "\\.rds$", full.names = TRUE)
    de_files <- list.files(tmp_de_dir, pattern = "\\.rds$", full.names = TRUE)

    expect_true(length(sce_files) == 1)
    expect_true(length(de_files) == 1)

    expect_true(file.exists(sce_files[1]))
    expect_true(file.exists(de_files[1]))

    # Optional cleanup
    unlink(tmp_sce_dir, recursive = TRUE)
    unlink(tmp_de_dir, recursive = TRUE)
})
