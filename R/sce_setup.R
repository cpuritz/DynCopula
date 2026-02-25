###############################################################################

#' Setup a SingleCellExperiment for analysis
#'
#' @description Setup a \code{SingleCellExperiment} for analysis by adding
#' necessary metadata.
#'
#' @param sce A \code{SingleCellExperiment}.
#' @param assay Which assay to use. Either \code{"counts"} or \code{"logcounts"}.
#' @param time_col The name of the \code{colData} column containing pseudotimes.
#' @param features Which genes to model. Default is all genes.
#' @param cores The number of cores to use for parallel computations. Default
#' is \code{1L}.
#'
#' @returns The same \code{SingleCellExperiment} as was passed as input, but
#' modified to include a named metadata entry \code{dyn_corr} that is used to
#' store all information related to this analysis.
#'
#' @export
setup <- function(sce,
                  assay = c("counts", "logcounts"),
                  time_col,
                  features = rownames(sce),
                  cores = 1L) {
    assert_that(
        methods::is(sce, "SingleCellExperiment"),
        time_col %in% colnames(SummarizedExperiment::colData(sce)),
        is.null(features) || is.character(features),
        is.numeric(cores) && cores >= 1L,
        all(features %in% rownames(sce))
    )

    assay <- match.arg(assay)
    assert_that(assay %in% SummarizedExperiment::assayNames(sce))

    metadata(sce)$dyn_corr <- list(
        time_col = time_col,
        features = features,
        assay = assay,
        cores = as.integer(cores)
    )
    return(sce)
}

###############################################################################
