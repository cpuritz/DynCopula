###############################################################################

#' Setup
#'
#' @description Setup
#'
#' @param sce A \code{SingleCellExperiment}.
#' @param tcol The name of the column containing pseudotimes.
#' @param features Genes to use. If \code{NULL} (the default), all genes in the
#' selected assay are used.
#' @param assay The assay to use. Default is \code{"counts"}.
#' @param cores The number of cores to use for parallel computations. Default
#' is \code{1L}.
#'
#' @returns The same \code{SingleCellExperiment} as was passed as input, but
#' modified to include a named metadata entry \code{dyn_corr_info}. This entry
#' is a list recording \code{tcol}, \code{features}, \code{assay}, and
#' \code{cores}.
#'
#' @export
setup <- function(sce,
                  tcol,
                  features = NULL,
                  assay = "counts",
                  cores = 1L) {
    assertthat::assert_that(
        methods::is(sce, "SingleCellExperiment"),
        tcol %in% colnames(SummarizedExperiment::colData(sce)),
        is.null(features) || is.character(features),
        is.character(assay) && assay %in% SummarizedExperiment::assayNames(sce),
        is.numeric(cores) && cores >= 1L
    )

    if (is.null(features)) {
        features <- rownames(SummarizedExperiment::assay(sce, assay))
    }

    S4Vectors::metadata(sce)$dyn_corr_info <- list(
        tcol = tcol,
        features = features,
        assay = assay,
        cores = as.integer(cores)
    )
    return(sce)
}

###############################################################################
