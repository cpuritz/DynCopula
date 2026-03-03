###############################################################################

#' Fit a dynamic Gaussian copula to scRNA-seq data
#'
#' @description Fit a dynamic Gaussian copula to scRNA-seq data using local
#' likelihood.
#'
#' @param sce A \code{SingleCellExperiment}.
#' @param lambda Vector of smoothing parameters to test. Default is
#' \code{10^(seq(-5, 5, length.out = 7))}.
#' @param K Dimension of the spline basis matrix. Default is \code{30}.
#' @param nfold Number of folds for cross-validation. Default is \code{5}.
#' @param control A \code{list} of control parameters for optimization.
#'
#' @details Cross-validation is used to select the value of \code{lambda}.
#' Optimization is performed using L-BFGS. The \code{control} argument is a list
#' that supplies control parameters for optimization. The following parameters
#' can be supplied:
#' \itemize{
#'   \item \code{max_itr} Maximum number of iterations. Default is \code{100}.
#'   \item \code{history_size} History size. Default is \code{30}.
#'   \item \code{tolerance_grad} Termination tolerance for gradient. Default is
#'   \code{1e-7}.
#'   \item \code{tolerance_change} Termination tolerance for log-likelihood.
#'   Default is \code{1e-9}.
#'   \item \code{precision}: Floating point precision for calculations. Either
#'   \code{"float32"} or \code{"float64"}. Default is \code{"float64"}.
#'   \item \code{boundary}: Boundary extension for spline knot boundaries.
#'   Default is \code{0.05}.
#' }
#' Any control parameters not specified are replaced by their default values.
#'
#' @returns The same \code{SingleCellExperiment} as was passed as input, but
#' with the metadata entry \code{dyn_corr} updated to include the following
#' elements:
#' \itemize{
#'   \item \code{rho}: Matrix of estimated correlation coefficients.
#'   \item \code{eta}: Matrix of estimated calibrations coefficients.
#'   \item \code{cv}: Cross-validation results.
#' }
#'
#' @export
fit_dyn_corr <- function(sce,
                         lambda = 10^(seq(-5, 5, length.out = 7)),
                         K = 30,
                         nfold = 5,
                         control = list()) {
    assert_that(
        methods::is(sce, "SingleCellExperiment"),
        "dyn_corr" %in% names(metadata(sce))
    )

    dyn_corr <- metadata(sce)$dyn_corr
    assert_that("FX" %in% names(dyn_corr))
    assay <- dyn_corr$assay
    pseudotimes <- sce[[dyn_corr$time_col]]

    if (assay == "counts") {
        # Construct jittered pseudo-observations
        FX <- dyn_corr$FXm + (dyn_corr$FX - dyn_corr$FXm) * dyn_corr$V
    } else {
        # Already jittered for logcounts
        FX <- dyn_corr$FX
    }

    # Estimate copula parameters
    design <- data.frame(
        "time" = pseudotimes
    )
    res <- fit_dyn_gc(
        FX = FX,
        design = design,
        lambda = lambda,
        K = K,
        nfold = nfold,
        cores = dyn_corr$cores,
        control = control
    )

    # Convert numeric labels to gene names
    gene_names <- sapply(colnames(res$rho), function(x) {
        x2 <- unlist(strsplit(x, split = "rho"))[2]
        ix <- as.numeric(unlist(strsplit(x2, split = '_')))
        return(paste(dyn_corr$features[ix], collapse = '_'))
    })
    colnames(res$rho) <- gene_names

    # Save results in metadata
    metadata(sce)$dyn_corr$eta <- res$eta
    metadata(sce)$dyn_corr$rho <- res$rho
    metadata(sce)$dyn_corr$cv <- res$cv

    return(sce)
}

###############################################################################
