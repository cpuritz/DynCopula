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
#' The argument \code{ncv} specifies the number of pseudotime values to use for
#' leave-one-out cross validation (LOOCV). If \code{ncv} is less than the
#' number of metacells (each of which is assigned a unique pseudotime value),
#' then only a subset of the pseudotime values are used. This reduces run time
#' but is only an estimate of full LOOCV.
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
        is.numeric(bandwidth) && bandwidth > 0 && bandwidth < 1,
        is.list(control),
        is.null(ncv) || (is.numeric(ncv) && ncv > 1)
    )

    sce_mc <- metadata(sce)$dyn_corr$metacell_sce
    dyn_corr <- metadata(sce_mc)$dyn_corr
    assert_that("margins" %in% names(dyn_corr))

    X <- SummarizedExperiment::assay(sce_mc, dyn_corr$assay)
    X <- Matrix::t(X[dyn_corr$features, ])
    pseudotimes <- sce_mc[[dyn_corr$time_col]]
    t0 <- sort(pseudotimes)

    # Construct jittered pseudo-observations
    FX <- dyn_corr$FXm + (dyn_corr$FX - dyn_corr$FXm) * dyn_corr$V

    if (length(bandwidth) == 1) {
        message("Estimating correlation coefficients")
        res <- fit_dynamic_gaussian(
            FX = FX,
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
            if (ncv > length(pseudotimes)) {
                ncv <- length(pseudotimes)
            }
        }

        message("Performing cross validation to select bandwidth")
        cv <- bandwidth_select(
            FX = FX,
            x = pseudotimes,
            bandwidths = bandwidth,
            xind = ncv,
            control = control,
            cores = dyn_corr$cores
        )
        h_opt <- cv$bandwidth[which.max(cv$cv)]

        message("Estimating correlation coefficients using optimal bandwidth")
        res <- fit_dynamic_gaussian(
            FX = FX,
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
