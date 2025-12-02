###############################################################################

#' Sample cells
#'
#' @description Sample cells at specified pseudotimes.
#'
#' @param sce A \code{SingleCellExperiment}.
#' @param times A vector of pseudotimes to sample at.
#'
#' @returns A \code{SingleCellExperiment}.
#'
#' @export
sample_cells <- function(sce,
                         times) {
    assert_that(
        methods::is(sce, "SingleCellExperiment"),
        "dyn_corr" %in% names(metadata(sce)),
        all(c("margins", "rho") %in% names(metadata(sce)$dyn_corr)),
        is.numeric(times)
    )

    dyn_corr <- metadata(sce)$dyn_corr
    ngenes <- length(dyn_corr$features)
    counts <- SummarizedExperiment::assay(sce, dyn_corr$assay)
    counts <- Matrix::t(counts[dyn_corr$features, ])
    pseudotimes <- SummarizedExperiment::colData(sce)[[dyn_corr$time_col]]

    # Out of the time points at which correlation coefficients were estimated,
    # use the closest one to the specified time
    t0 <- dyn_corr$x0
    min_ix <- sapply(times, function(x) { which.min(abs(x - t0)) })

    message("Sampling copula")
    rho_split <- apply(
        X = dyn_corr$rho,
        MARGIN = 1,
        FUN = c,
        simplify = FALSE
    )
    rho_slices <- rho_split[min_ix]
    progressr::with_progress({
        pbar <- progressr::progressor(along = seq_len(ngenes))
        U <- lapply(
            X = rho_slices,
            FUN = function(r) {
                copula <- copula::normalCopula(
                    param = unlist(r),
                    dim = ngenes,
                    dispstr = "un"
                )
                V <- copula::rCopula(1L, copula)
                pbar()
                return(V)
            }
        )
    })
    U <- do.call(rbind, U)

    # Convert to specified margins
    message("Converting margins")
    new_times <- stats::setNames(data.frame(times), dyn_corr$time_col)
    progressr::with_progress({
        pbar <- progressr::progressor(along = seq_len(ngenes))
        counts_sim <- lapply(seq_len(ngenes), function(i) {
            mdat <- data.frame(counts[, i], pseudotimes)
            names(mdat) <- c("x", dyn_corr$time_col)
            mfun <- dyn_corr$margins[[i]]
            par_pred <- gamlss::predictAll(
                object = mfun,
                newdata = new_times,
                type = "response",
                data = mdat
            )
            # Load the correct quantile function
            qfun <- getExportedValue("gamlss.dist", paste0("q", mfun$family[1]))
            V <- do.call(qfun, c(list(p = U[, i]), par_pred))
            pbar()
            return(V)
        })
    })
    counts_sim <- Matrix::t(do.call(cbind, counts_sim))
    counts_sim <- methods::as(counts_sim, "dgCMatrix")

    # Create a SingleCellExperiment with the simulated counts
    sce_sim <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = counts_sim)
    )
    colnames(sce_sim) <- paste0("sim", seq_along(times))
    rownames(sce_sim) <- dyn_corr$features
    sce_sim$pseudotime <- times
    return(sce_sim)
}

###############################################################################
