###############################################################################

#' Time-varying gene correlations
#'
#' @description Estimate a time-varying gene-gene correlation matrix.
#'
#' @param sce A \code{SingleCellExperiment}.
#' @param t0 A vector of pseudotimes to estimate copula parameters at. If
#' \code{NULL} (the default), parameters are estimated at all pseudotimes.
#' @param h Kernel bandwidth. Must satisfy \code{0 < h < 1}.
#' @param control A \code{list} of control parameters for optimization.
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
#'   \item \code{rho}: Matrix of estimated pairwise correlation coefficients.
#'   \item \code{eta}: Matrix of estimated coefficients in the unconstrained
#'   space.
#'   \item \code{t0}: The input argument \code{t0}.
#'   \item \code{h}: The input argument \code{h}.
#' }
#'
#' @export
fit_dyn_corr <- function(sce,
                         t0 = NULL,
                         h,
                         control = list()) {
    assert_that(
        methods::is(sce, "SingleCellExperiment"),
        "dyn_corr" %in% names(metadata(sce)),
        "metacell_sce" %in% names(metadata(sce)$dyn_corr),
        is.null(t0) || is.numeric(t0)
    )

    sce_mc <- metadata(sce)$dyn_corr$metacell_sce
    dyn_corr <- metadata(sce_mc)$dyn_corr
    assert_that("margins" %in% names(dyn_corr))

    X <- SummarizedExperiment::assay(sce_mc, dyn_corr$assay)
    X <- Matrix::t(X[dyn_corr$features, ])
    pseudotimes <- sce_mc[[dyn_corr$time_col]]

    if (is.null(t0)) {
        t0 <- sort(unique(pseudotimes))
    }

    FX <- dyn_corr$FX
    FXm <- dyn_corr$FXm
    V <- dyn_corr$V
    NX <- stats::qnorm(FXm + (FX - FXm) * V)

    message("Estimating correlation coefficients")
    res <- fit_dynamic_gaussian(
        NX = NX,
        x = pseudotimes,
        x0 = t0,
        h = h,
        control = control,
        cores = dyn_corr$cores
    )

    # Interpolate unconstrained coefficients to original pseudotime values
    eta <- res$eta
    t_old <- res$x0
    t_new <- sce[[dyn_corr$time_col]]
    nc <- dim(eta)[2]
    eta_int <- matrix(nrow = length(t_new), ncol = nc, dimnames = dimnames(eta))
    for (j in nc) {
        eta_int[, j] <- stats::approx(
            x = t_old,
            y = eta[, j],
            xout = t_new,
            rule = 2
        )$y
    }

    # Interpolated correlation coefficients
    rho_int <- t(apply(eta_int, 1, function(v) {
        copula::P2p(vec2cor(v))
    }))
    if (nc == 1) {
        rho_int <- t(rho_int)
    }

    # Map numeric labels to gene names
    gene_names <- sapply(colnames(res$rho), function(x) {
        x2 <- unlist(strsplit(x, split = "rho"))[2]
        ix <- as.numeric(unlist(strsplit(x2, split = '_')))
        return(paste(dyn_corr$features[ix], collapse = '_'))
    })
    colnames(res$eta) <- colnames(res$rho) <- gene_names
    colnames(eta_int) <- colnames(rho_int) <- gene_names

    metadata(sce)$dyn_corr$x <- res$x
    metadata(sce)$dyn_corr$x0 <- res$x0
    metadata(sce)$dyn_corr$h <- res$h
    metadata(sce)$dyn_corr$eta <- res$eta
    metadata(sce)$dyn_corr$rho <- res$rho
    metadata(sce)$dyn_corr$eta_int <- eta_int
    metadata(sce)$dyn_corr$rho_int <- rho_int

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
#'
#' @details This function can only be run after
#' \link[DynCopula]{generate_metacells} has been run.
#'
#' The optimal kernel bandwidth is selected via AIC. This requires
#' parameter estimation at every pseudotime value for every bandwidth and thus
#' may take a while to run.
#'
#' See \link[DynCopula]{fit_dyn_corr} for optimization details.
#'
#' @returns The same \code{SingleCellExperiment} as was passed as input, but
#' with the metadata entry \code{dyn_corr} updated to include the following
#' elements:
#' \itemize{
#'   \item \code{rho}: Matrix of estimated pairwise correlation coefficients.
#'   \code{eta}: Matrix of estimated coefficients in the unconstrained space.
#'   \item \code{t0}: The time points coefficients were estimated at.
#'   \item \code{aic}: Vector of AIC values.
#'   \item \code{h}: Optimal kernel bandwidth.
#'   \item \code{bandwidths}: The input argument \code{bandwidths}.
#' }
#'
#' @export
fit_dyn_corr_sel <- function(sce,
                             bandwidths,
                             control = list()) {
    assert_that(
        methods::is(sce, "SingleCellExperiment"),
        "dyn_corr" %in% names(metadata(sce)),
        "metacell_sce" %in% names(metadata(sce)$dyn_corr),
        is.numeric(bandwidths)
    )

    sce_comp <- metadata(sce)$dyn_corr$metacell_sce
    dyn_corr <- metadata(sce_comp)$dyn_corr
    X <- SummarizedExperiment::assay(sce_comp, dyn_corr$assay)
    X <- Matrix::t(X[dyn_corr$features, ])
    pseudotimes <- sce_comp[[dyn_corr$time_col]]

    FX <- dyn_corr$FX
    FXm <- dyn_corr$FXm
    V <- dyn_corr$V
    NX <- stats::qnorm(FXm + (FX - FXm) * V)

    res <- bandwidth_select(
        NX = NX,
        x = pseudotimes,
        bandwidths = bandwidths,
        control = control,
        cores = dyn_corr$cores,
        return_all = FALSE
    )

    # Interpolate unconstrained coefficients to original pseudotime values
    eta <- res$eta
    t_old <- res$x0
    t_new <- sce[[dyn_corr$time_col]]
    nc <- dim(eta)[2]
    eta_int <- matrix(nrow = length(t_new), ncol = nc, dimnames = dimnames(eta))
    for (j in nc) {
        eta_int[, j] <- stats::approx(
            x = t_old,
            y = eta[, j],
            xout = t_new,
            rule = 2
        )$y
    }

    # Interpolated correlation coefficients
    rho_int <- t(apply(eta_int, 1, function(v) {
        copula::P2p(vec2cor(v))
    }))
    if (nc == 1) {
        rho_int <- t(rho_int)
    }

    # Map numeric labels to gene names
    gene_names <- sapply(colnames(res$rho), function(x) {
        x2 <- unlist(strsplit(x, split = "rho"))[2]
        ix <- as.numeric(unlist(strsplit(x2, split = '_')))
        return(paste(dyn_corr$features[ix], collapse = '_'))
    })
    colnames(res$eta) <- colnames(res$rho) <- gene_names
    colnames(eta_int) <- colnames(rho_int) <- gene_names

    metadata(sce)$dyn_corr$x <- res$x
    metadata(sce)$dyn_corr$x0 <- res$x0
    metadata(sce)$dyn_corr$h <- res$h
    metadata(sce)$dyn_corr$eta <- res$eta
    metadata(sce)$dyn_corr$rho <- res$rho
    metadata(sce)$dyn_corr$eta_int <- eta_int
    metadata(sce)$dyn_corr$rho_int <- rho_int
    metadata(sce)$dyn_corr$aic <- res$aic
    metadata(sce)$dyn_corr$bandwidths <- res$bandwidths

    return(sce)
}

###############################################################################
