###############################################################################

#' Setup
#'
#' @description Setup a SingleCellExperiment by adding necessary metadata.
#'
#' @param sce A \code{SingleCellExperiment}.
#' @param time_col The name of the column containing pseudotimes.
#' @param features Genes to use. If \code{NULL} (the default), all genes in the
#' selected assay are used.
#' @param assay The assay to use. Default is \code{"counts"}.
#' @param cores The number of cores to use for parallel computations. Default
#' is \code{1L}.
#'
#' @returns The same \code{SingleCellExperiment} as was passed as input, but
#' modified to include a named metadata entry \code{dyn_corr}. This entry is a
#' list recording \code{time_col}, \code{features}, \code{assay}, and
#' \code{cores}.
#'
#' @export
setup <- function(sce,
                  time_col,
                  features = NULL,
                  assay = "counts",
                  cores = 1L) {
    assert_that(
        methods::is(sce, "SingleCellExperiment"),
        time_col %in% colnames(SummarizedExperiment::colData(sce)),
        is.null(features) || is.character(features),
        is.character(assay) && assay %in% SummarizedExperiment::assayNames(sce),
        is.numeric(cores) && cores >= 1L
    )

    if (is.null(features)) {
        features <- rownames(SummarizedExperiment::assay(sce, assay))
    }

    metadata(sce)$dyn_corr <- list(
        time_col = time_col,
        features = features,
        assay = assay,
        cores = as.integer(cores)
    )
    return(sce)
}

###############################################################################
