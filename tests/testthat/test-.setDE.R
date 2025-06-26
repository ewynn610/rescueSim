test_that(".setDE: scalar deLogFC, 1 group, 2 timepoints", {
    colDat <- data.frame(group = rep(0, 4), time = rep(c(0, 1), each = 2))
    result <- rescueSim:::.setDE(deLogFC = 1, propDE = 1, nGenes = 3, colDat = colDat)

    expect_equal(dim(result$de), c(3, 4))
    expect_true(all(result$de > 0))
    expect_length(result$deFacs, 1)
    expect_length(result$deFacs[[1]], 3)
    expect_true(all(result$deFacs[[1]] %in% c(1, -1)))  # since only +/-1 possible
})

test_that(".setDE: vector deLogFC, propDE < 1", {
    colDat <- data.frame(group = rep(0, 4), time = rep(c(0, 1), each = 2))
    set.seed(24)
    result <- rescueSim:::.setDE(deLogFC = c(-1, 1), propDE = 0.5, nGenes = 4, colDat = colDat)

    expect_equal(dim(result$de), c(4, 4))
    expect_true(any(result$deFacs[[1]] == 0))
    expect_true(all(result$deFacs[[1]] %in%c(0,1,-1)))
})

test_that(".setDE: list deLogFC, one timepoint, two groups", {
    colDat <- data.frame(group = rep(0:1, each = 2), time = rep(0, 4))
    custom_deLogFC <- list(group1 = c(1, 0, 0))

    result <- rescueSim:::.setDE(deLogFC = custom_deLogFC, propDE = numeric(0), nGenes = 3, colDat = colDat)

    expect_equal(dim(result$de), c(3, 4))
    expect_equal(result$de[1, colDat$group == 1], rep(2, 2))
    expect_equal(result$de[2:3, ], matrix(1, nrow=2, ncol=4))
})

test_that(".setDE: list deLogFC, two timepoints, one group", {
    colDat <- data.frame(group = rep(0, 4), time = rep(0:1, each = 2))
    geneNames <- c("geneA", "geneB", "geneC")
    custom_deLogFC <- list(time1 = setNames(c(0, 1, 0), geneNames))

    result <-rescueSim:::.setDE(deLogFC = custom_deLogFC, propDE = numeric(0), nGenes = 3, colDat = colDat)

    expect_equal(dim(result$de), c(3, 4))
    expect_equal(result$de[2, colDat$time == 1], rep(2, 2))
    expect_equal(result$de[c(1,3), ], matrix(1, nrow=2, ncol=4))
})

test_that(".setDE: list deLogFC, two timepoints and two groups", {
    colDat <- expand.grid(time = 0:1, group = 0:1)
    geneNames <- c("geneA", "geneB", "geneC")
    custom_deLogFC <- list(
        time1_group0 = setNames(c(1, 0, 0), geneNames),
        time1_group1 = setNames(c(0, 1, 0), geneNames),
        time0_group1 = setNames(c(0, 0, 1), geneNames)
    )

    result <- rescueSim:::.setDE(deLogFC = custom_deLogFC, propDE = numeric(0), nGenes = 3, colDat = colDat)
    colNames <- paste0("time", colDat$time, "_group", colDat$group)

    expect_equal(result$de[1, colNames == "time1_group0"], rep(2, sum(colNames == "time1_group0")))
    expect_equal(result$de[2, colNames == "time1_group1"], rep(2, sum(colNames == "time1_group1")))
    expect_equal(result$de[3, colNames == "time0_group1"], rep(2, sum(colNames == "time0_group1")))
})

test_that(".setDE: returns all 1s if propDE = 0", {
    colDat <- data.frame(group = rep(0, 4), time = rep(c(0, 1), each = 2))
    result <-rescueSim:::.setDE(deLogFC = 1, propDE = 0, nGenes = 3, colDat = colDat)

    expect_true(all(result$de == 1))
    expect_true(all(result$deFacs[[1]] == 0))
})

test_that(".setDE: list deLogFC, 3 timepoints, 1 group", {
    # Simulate colDat for 3 timepoints, 1 group
    colDat <- data.frame(
        group = rep(0, 6),
        time = rep(0:2, each = 2)
    )

    geneNames <- c("geneA", "geneB", "geneC")

    # Provide logFC for time1 and time2 vs. time0 (baseline)
    custom_deLogFC <- list(
        time1 = setNames(c(1, 0, 0), geneNames),  # geneA upregulated at time1
        time2 = setNames(c(0, 1, 0), geneNames)   # geneB upregulated at time2
    )

    result <- rescueSim:::.setDE(
        deLogFC = custom_deLogFC,
        propDE = numeric(0),
        nGenes = 3,
        colDat = colDat
    )

    # Expect correct shape
    expect_equal(dim(result$de), c(3, 6))

    # Column labels expected inside function: time0, time1, time2
    expect_equal(result$de[1, colDat$time == 1], rep(2, 2))  # geneA at time1
    expect_equal(result$de[2, colDat$time == 2], rep(2, 2))  # geneB at time2
    expect_equal(result$de[3, ], rep(1, 6))                  # geneC unchanged
})

test_that(".setDE: numeric deLogFC gives linear change over 3 timepoints, 1 group", {
    set.seed(12)
    colDat <- data.frame(
        group = rep(0, 6),
        time = rep(0:2, each = 2)
    )

    result <- rescueSim:::.setDE(deLogFC = 1, propDE = 1, nGenes = 1, colDat = colDat)

    # Reverse-engineer the design factor:
    # For group=0, function internally adds 1, so group=1
    # time / (nTimepoints - 1) * group = time / 2
    expected_logFC <- c(0, 0.5, 1)  # for time0, time1, time2
    expected_FC <- 2^expected_logFC

    actual <- result$de[1, ]

    expect_equal(actual[colDat$time == 0], rep(expected_FC[1], 2))
    expect_equal(actual[colDat$time == 1], rep(expected_FC[2], 2))
    expect_equal(actual[colDat$time == 2], rep(expected_FC[3], 2))
})

test_that(".setDE: numeric deLogFC gives linear interaction over 3 timepoints, 2 groups", {
    set.seed(24)
    colDat <- expand.grid(
        time = 0:2,
        group = 0:1
    )
    colDat <- colDat[rep(1:nrow(colDat), each = 2), ]  # 2 replicates per condition

    result <- rescueSim:::.setDE(deLogFC = 1, propDE = 1, nGenes = 1, colDat = colDat)

    # Design factor is (time / 2) * group
    expected_logFC <- outer(0:2, 0:1, function(t, g) (t / 2) * g)*-1
    expected_FC <- 2^expected_logFC

    for (t in 0:2) {
        for (g in 0:1) {
            idx <- which(colDat$time == t & colDat$group == g)
            expected <- expected_FC[t + 1, g + 1]
            expect_equal(result$de[1, idx], rep(expected, length(idx)))
        }
    }
})

