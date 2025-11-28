###############################################################################

#' Fit a dynamic Gaussian copula to scRNA-seq data
#'
#' @description Fit a dynamic Gaussian copula to scRNA-seq data using local
#' likelihood.
#'
#' @param sce A \code{SingleCellExperiment}.
#' @param bandwidth Kernel bandwidth. Must be between \code{0} and \code{1}. If
#' a vector of bandwidths is passed, then cross validation is used to choose
#' the optimal bandwidth.
#' @param control A \code{list} of control parameters for optimization.
#' @param ncv If multiple bandwidths are passed, leave-one-out cross validation
#' is used to choose the optimal bandwidth. This argument specifies the number
#' of pseudotime values to use for cross validation. If \code{NULL} (the
#' default), all pseudotime values are used. If only a single bandwidth is
#' passed, this parameter has no effect.
#'
#' @details This function can only be run after
#' \link[DynCopula]{generate_metacells} has been run.
#'
#' Optimization is performed using L-BFGS. The \code{control} argument
#' is a list that supplies control parameters for optimization. The following
#' parameters can be supplied:
#' \itemize{
#'   \item \code{max_epoch} Maximum number of epochs. Default is \code{1}.
#'   \item \code{max_itr} Maximum number of internal iterations. Default isjls
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
#'   \item \code{rho}: Matrix of estimated copula parameters.
#'   \item \code{bandwidth}: The kernel bandwidth used to estimate copula
#'   parameters.
#' }
#'
#' @export
fit_dyn_corr <- function(sce,
                         bandwidth,
                         control = list(),
                         ncv = NULL) {
    assert_that(
        methods::is(sce, "SingleCellExperiment"),
        "dyn_corr" %in% names(metadata(sce)),
        "metacell_sce" %in% names(metadata(sce)$dyn_corr),
        is.null(ncv) || is.numeric(ncv)
    )

    sce_mc <- metadata(sce)$dyn_corr$metacell_sce
    dyn_corr <- metadata(sce_mc)$dyn_corr
    assert_that("margins" %in% names(dyn_corr))

    X <- SummarizedExperiment::assay(sce_mc, dyn_corr$assay)
    X <- Matrix::t(X[dyn_corr$features, ])
    pseudotimes <- sce_mc[[dyn_corr$time_col]]
    t0 <- sort(pseudotimes)

    FX <- dyn_corr$FX
    FXm <- dyn_corr$FXm
    V <- dyn_corr$V
    NX <- stats::qnorm(FXm + (FX - FXm) * V)

    if (length(bandwidth) == 1) {
        message("Estimating correlation coefficients")
        res <- fit_dynamic_gaussian(
            NX = NX,
            x = pseudotimes,
            x0 = t0,
            h = bandwidth,
            control = control,
            cores = dyn_corr$cores
        )
    } else {
        if (is.null(ncv)) {
            ncv <- length(pseudotimes)
        } else {
            ncv <- as.integer(ncv)
            assertthat::assert_that(ncv > 0 && ncv <= length(pseudotimes))
        }

        message("Performing cross validation to select bandwidth")
        cv <- bandwidth_select_cv(
            NX = NX,
            x = pseudotimes,
            bandwidths = bandwidth,
            xind = ncv,
            control = control,
            cores = dyn_corr$cores
        )
        h_opt <- cv$bandwidth[which.max(cv$ll)]

        message("Estimating correlation coefficients")
        res <- fit_dynamic_gaussian(
            NX = NX,
            x = pseudotimes,
            x0 = t0,
            h = h_opt,
            control = control,
            cores = dyn_corr$cores
        )
    }

    # Convert numeric labels to gene names
    gene_names <- sapply(colnames(res$rho), function(x) {
        x2 <- unlist(strsplit(x, split = "rho"))[2]
        ix <- as.numeric(unlist(strsplit(x2, split = '_')))
        return(paste(dyn_corr$features[ix], collapse = '_'))
    })
    colnames(res$rho) <- gene_names

    # Save results in metadata
    metadata(sce)$dyn_corr$rho <- res$rho
    metadata(sce)$dyn_corr$bandwidth <- res$h

    return(sce)
}

###############################################################################
