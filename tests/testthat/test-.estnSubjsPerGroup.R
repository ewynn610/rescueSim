test_that(".estnSubjsPerGroup works without groupVariable (single group)", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        colData = data.frame(
            subjectID = rep(paste0("S", 1:3), each = 10)
        )
    )
    nSubjs <- rescueSim:::.estnSubjsPerGroup(sce, subjectVariable = "subjectID", groupVariable = NULL, twoGroupDesign = TRUE)
    expect_equal(nSubjs, 3)
})

test_that(".estnSubjsPerGroup works with two groups and balanced subjects", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        colData = data.frame(
            subjectID = rep(paste0("S", 1:4), each = 5),
            group = rep(rep(c("A", "B"), each = 10), 1)
        )
    )
    nSubjs <- rescueSim:::.estnSubjsPerGroup(sce, subjectVariable = "subjectID", groupVariable = "group", twoGroupDesign = TRUE)
    expect_equal(nSubjs, 2)
})

test_that(".estnSubjsPerGroup returns named vector for unbalanced subjects", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        colData = data.frame(
            subjectID = c(rep("S1", 5), rep("S2", 5), rep("S3", 5), rep("S4", 5), rep("S5", 5)),
            group = c(rep("A", 10), rep("B", 15))
        )
    )
    nSubjs <- rescueSim:::.estnSubjsPerGroup(sce, subjectVariable = "subjectID", groupVariable = "group", twoGroupDesign = TRUE)
    expect_named(nSubjs)
    expect_equal(length(nSubjs), 2)
    expect_equal(sum(nSubjs), 5) # 2 in A, 3 in B
})

test_that(".estnSubjsPerGroup throws error if subject assigned to multiple groups", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        colData = data.frame(
            subjectID = rep(c("S1", "S2"), each = 5),
            group = c(rep("A", 5), rep("B", 3), rep("A", 2))
        )
    )
    expect_error(
        rescueSim:::.estnSubjsPerGroup(sce, subjectVariable = "subjectID", groupVariable = "group", twoGroupDesign = TRUE),
        "Each subject must belong to only one group"
    )
})

test_that(".estnSubjsPerGroup throws error if group count mismatch twoGroupDesign", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        colData = data.frame(
            subjectID = c(rep("s1", 3), rep("s2", 5), rep("s3", 5), rep("s4", 5)),
            group = c(rep("A", 3), rep("B", 5), rep("C", 10))
        )
    )
    expect_error(
        rescueSim:::.estnSubjsPerGroup(sce, subjectVariable = "subjectID", groupVariable = "group", twoGroupDesign = TRUE),
        "Number of groups in data does not match twoGroupDesign setting"
    )
})
