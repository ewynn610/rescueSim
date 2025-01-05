#' Example Repeated Measures scRNA-Seq Data
#'
#' Sample \code{\link[SingleCellExperiment]{SingleCellExperiment}} object with
#' gene expression data for recruited airspace macrophage cells from
#' bronchoalveolar lavage samples collected from healthy adults.
#' The data has been subset to contain gene expression for 1940 genes from
#' 976 cells. Data was collected from five subjects at two timepoints per
#' subject. Genes included in the dataset were assessed to be invariant
#' across timepoints.
#'
#'
#' .
#'
#' @format \code{\link[SingleCellExperiment]{SingleCellExperiment}} object with
#' the following data in object slots:
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
#' }
#' @usage data(RecAM_sce)
"RecAM_sce"
