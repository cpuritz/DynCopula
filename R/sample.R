###############################################################################

#' Sample cells
#'
#' @description Sample cells from the fitted joint distribution given new
#' covariate values.
#'
#' @param sce A \code{SingleCellExperiment}.
#' @param design A \code{data.frame} specifying the new design matrix.
#'
#' @returns A \code{SingleCellExperiment} with simulated counts.
#'
#' @export
sample_cells <- function(sce,
                         design) {
    assert_that(
        methods::is(sce, "SingleCellExperiment"),
        "copula_fit" %in% names(metadata(sce)),
        is.data.frame(design)
    )

    if (is.null(metadata(sce)$copula_fit$margins)) {
        stop("Cells can't be sampled since the margins were not saved. ",
             "Rerun 'fit_margins' with 'save = TRUE'.")
    }

    copula_fit <- metadata(sce)$copula_fit
    object <- copula_fit$fit

    # Predict correlation matrices
    pred <- stats::predict(
        object = copula_fit$fit,
        design = design,
        type = "response_matrix"
    )

    # Sample from copula
    d <- object$dim
    U <- do.call(rbind, lapply(pred, function(R) {
        cop <- copula::normalCopula(
            param = copula::P2p(R),
            dim = d,
            dispstr = "un"
        )
        return(copula::rCopula(1L, cop))
    }))

    # Extract counts matrix
    features <- copula_fit$features
    counts <- SummarizedExperiment::assay(sce, "counts")
    counts <- Matrix::t(counts[features, ])

    # Convert margins
    design_old <- as.data.frame(SummarizedExperiment::colData(sce))
    counts_sim <- lapply(seq_along(features), function(i) {
        # Marginal model
        mfun <- copula_fit$margins[[i]]

        # Original model data
        ddata <- cbind(design_old, data.frame(countsforgene = counts[, i]))

        # Get model parameters at new covariates
        par_pred <- gamlss::predictAll(
            object = mfun,
            newdata = design,
            type = "response",
            data = ddata
        )

        # Load the correct quantile function
        qfun <- getExportedValue("gamlss.dist", paste0("q", mfun$family[1]))

        # Transform margins
        V <- do.call(qfun, c(list(p = U[, i]), par_pred))
        return(V)
    })

    # Construct new SingleCellExperiment
    counts_sim <- t(do.call(cbind, counts_sim))
    assay <- list(counts = methods::as(counts_sim, "dgCMatrix"))
    sce_sim <- SingleCellExperiment::SingleCellExperiment(assay)
    rownames(sce_sim) <- features
    colnames(sce_sim) <- paste0("sim", seq_len(dim(design)[1]))
    SummarizedExperiment::colData(sce_sim) <- methods::as(design, "DataFrame")

    return(sce_sim)
}

###############################################################################
