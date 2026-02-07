###############################################################################

#' Fit a dynamic Gaussian copula to scRNA-seq data
#'
#' @description Fit a dynamic Gaussian copula to scRNA-seq data using local
#' likelihood.
#'
#' @param sce A \code{SingleCellExperiment}.
#' @param lambda Vector of smoothing parameters. Default is
#' \code{10^(seq(-5, 5, length.out = 7))}.
#' @param df Vector of degrees of freedom. Default is \code{c(10, 50, 100)}.
#' @param nfold Number of folds for cross-validation. Default is \code{5}.
#' @param control A \code{list} of control parameters for optimization.
#'
#' @details This function can only be run after
#' \link[DynCopula]{generate_metacells} has been run.
#'
#' Cross-validation (CV) is used to select the values of \code{lambda} and
#' \code{df}.
#'
#' Optimization is performed using L-BFGS. The \code{control} argument
#' is a list that supplies control parameters for optimization. The following
#' parameters can be supplied:
#' \itemize{
#'   \item \code{max_outer} Maximum number of outer iterations. Default is
#'   \code{1}.
#'   \item \code{max_itr} Maximum number of inner iterations. Default is
#'   \code{100}.
#'   \item \code{history_size} History size. Default is \code{30}.
#'   \item \code{tolerance_grad} Termination tolerance for gradient. Default is
#'   \code{1e-7}.
#'   \item \code{tolerance_change} Termination tolerance for log-likelihood.
#'   Default is \code{1e-9}.
#' }
#'
#' @returns The same \code{SingleCellExperiment} as was passed as input, but
#' with the metadata entry \code{dyn_corr} updated to include the following
#' elements:
#' \itemize{
#'   \item \code{rho}: Matrix of estimated correlation coefficients.
#'   \item \code{eta}: Matrix of estimated calibrations coefficients.
#'   \item \code{lambda}: The smoothing parameter selected via CV.
#'   \item \code{df}: The degrees of freedom selected CV.
#' }
#'
#' @export
fit_dyn_corr <- function(sce,
                         lambda = 10^(seq(-5, 5, length.out = 7)),
                         df = c(10, 50, 100),
                         nfold = 5,
                         control = list()) {
    assert_that(
        methods::is(sce, "SingleCellExperiment"),
        "dyn_corr" %in% names(metadata(sce)),
        "metacell_sce" %in% names(metadata(sce)$dyn_corr)
    )

    sce_mc <- metadata(sce)$dyn_corr$metacell_sce
    dyn_corr <- metadata(sce_mc)$dyn_corr
    assert_that("margins" %in% names(dyn_corr))

    pseudotimes <- sce_mc[[dyn_corr$time_col]]

    # Construct jittered pseudo-observations
    FX <- dyn_corr$FXm + (dyn_corr$FX - dyn_corr$FXm) * dyn_corr$V

    # Estimate copula parameters
    res <- fit_spline_gaussian(
        FX = FX,
        x = pseudotimes,
        lambda = lambda,
        df = df,
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
    metadata(sce)$dyn_corr$lambda <- res$lambda
    metadata(sce)$dyn_corr$df <- res$df

    return(sce)
}

###############################################################################
