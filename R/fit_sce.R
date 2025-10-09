###############################################################################

#' Time-varying gene correlations
#'
#' @description Estimate a time-varying gene-gene correlation matrix.
#'
#' @param sce A \code{SingleCellExperiment}.
#' @param t0 A vector of pseudotimes to estimate copula parameters at.
#' @param h Kernel bandwidth. Must satisfy \code{0 < h < 1}.
#' @param control A \code{list} of control parameters for optimization.
#'
#' @details Optimization is performed using L-BFGS. The \code{control} argument
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
#' modified to include a named metadata entry \code{dyn_corr}. This entry is a
#' list with the following components:
#' \itemize{
#'   \item \code{rho}: Matrix of estimated pairwise correlation coefficients.
#'   \item \code{t0}: The input argument \code{t0}.
#'   \item \code{h}: The input argument \code{h}.
#' }
#'
#' @export
fit_dyn_corr <- function(sce,
                         t0,
                         h,
                         control = list()) {
    assert_that(
        methods::is(sce, "SingleCellExperiment"),
        "dyn_corr_info" %in% names(metadata(sce)),
        "margins" %in% names(metadata(sce)),
        is.numeric(t0)
    )

    if ("metacell_sce" %in% names(metadata(sce))) {
        sce_comp <- metadata(sce)$metacell_sce
    } else {
        sce_comp <- sce
    }

    info <- metadata(sce_comp)$dyn_corr_info
    X <- SummarizedExperiment::assay(sce_comp, info$assay)
    X <- Matrix::t(X[info$features, ])
    pseudotimes <- SummarizedExperiment::colData(sce_comp)[[info$time_col]]
    pobs <- .pseudo_obs(sce_comp)

    message("Estimating correlation coefficients")
    res <- fit_dynamic_gaussian(
        FX = pobs$FX,
        FXm = pobs$FXm,
        x = pseudotimes,
        x0 = t0,
        h = h,
        control = control,
        cores = info$cores
    )

    # Map numeric labels to gene names
    colnames(res$rho) <- sapply(colnames(res$rho), function(x) {
        x2 <- unlist(strsplit(x, split = "rho"))[2]
        ix <- as.numeric(unlist(strsplit(x2, split = '_')))
        return(paste(info$features[ix], collapse = '_'))
    })

    res <- res[c("rho", "x0", "h")]
    names(res)[names(res) == "x0"] <- "t0"
    metadata(sce)$dyn_corr <- res
    return(sce)
}

###############################################################################

#' Time-varying gene correlations
#'
#' @description Estimate a time-varying gene-gene correlation matrix.
#'
#' @param sce A \code{SingleCellExperiment}.
#' @param bandwidths A vector of kernel bandwidths.
#' @param control A \code{list} of control parameters for optimization.
#' @param return_all Whether all models should be returned, or just the best.
#' Default is \code{FALSE}.
#'
#' @details The optimal kernel bandwidth is selected via AIC. This requires
#' parameter estimation at every pseudotime value for every bandwidth and thus
#' may take a while to run.
#'
#' See \link[DynCopula]{fit_dyn_corr} for optimization details.
#'
#' @returns The same \code{SingleCellExperiment} as was passed as input, but
#' modified to include a named metadata entry \code{dyn_corr}. This entry is a
#' list with the following components:
#' \itemize{
#'   \item \code{rho}: Matrix of estimated pairwise correlation coefficients.
#'   If \code{return_all = TRUE}, this will be a list of matrices, one for
#'   each bandwidth.
#'   \item \code{t0}: The time points coefficients were estimated at.
#'   \item \code{aic}: Vector of AIC values.
#'   \item \code{h}: Optimal kernel bandwidth.
#'   \item \code{bandwidths}: The input argument \code{bandwidths}.
#' }
#'
#' @export
fit_dyn_corr_sel <- function(sce,
                             bandwidths,
                             control = list(),
                             return_all = FALSE) {
    assert_that(
        methods::is(sce, "SingleCellExperiment"),
        "dyn_corr_info" %in% names(metadata(sce)),
        "margins" %in% names(metadata(sce)),
        is.numeric(bandwidths),
        is.logical(return_all)
    )

    if ("metacell_sce" %in% names(metadata(sce))) {
        sce_comp <- metadata(sce)$metacell_sce
    } else {
        sce_comp <- sce
    }

    info <- metadata(sce_comp)$dyn_corr_info
    X <- SummarizedExperiment::assay(sce_comp, info$assay)
    X <- Matrix::t(X[info$features, ])
    pseudotimes <- SummarizedExperiment::colData(sce_comp)[[info$time_col]]
    pobs <- .pseudo_obs(sce_comp)

    res <- bandwidth_select(
        FX = pobs$FX,
        FXm = pobs$FXm,
        x = pseudotimes,
        bandwidths = bandwidths,
        control = control,
        cores = info$cores,
        return_all = return_all
    )

    # Convert column names of coefficient matrices to gene names
    cnames <- ifelse(return_all, colnames(res$rho[[1]]), colnames(res$rho))
    cnames <- sapply(cnames, function(x) {
        x2 <- unlist(strsplit(x, split = "rho"))[2]
        ix <- as.numeric(unlist(strsplit(x2, split = '_')))
        return(paste(info$features[ix], collapse = '_'))
    })
    if (return_all) {
        for (i in seq_along(res$rho)) {
            colnames(res$rho[[i]]) <- cnames
        }
    } else {
        colnames(res$rho) <- cnames
    }

    res <- res[c("rho", "x0", "aic", "h", "bandwidths")]
    names(res)[names(res) == "x0"] <- "t0"
    metadata(sce)$dyn_corr <- res
    return(sce)
}

###############################################################################

#' Get results
#'
#' @description Get results from running either \link[DynCopula]{fit_dyn_corr}
#' or \link[DynCopula]{fit_dyn_corr_sel}.
#'
#' @param sce A \code{SingleCellExperiment}.
#'
#' @returns A list.
#'
#' @export
get_results <- function(sce) {
    return(metadata(sce)$dyn_corr)
}

###############################################################################
