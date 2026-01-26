###############################################################################

#' Sample cells
#'
#' @description Sample cells at specified pseudotimes.
#'
#' @param sce A \code{SingleCellExperiment}.
#' @param new_times A vector of pseudotimes to sample at.
#'
#' @returns A \code{SingleCellExperiment}.
#'
#' @export
sample_cells <- function(sce, new_times) {
    assert_that(
        methods::is(sce, "SingleCellExperiment"),
        "dyn_corr" %in% names(metadata(sce)),
        all(c("margins", "rho") %in% names(metadata(sce)$dyn_corr)),
        is.numeric(new_times)
    )

    dyn_corr <- metadata(sce)$dyn_corr
    genes <- dyn_corr$features
    time_col <- dyn_corr$time_col

    counts <- SummarizedExperiment::assay(sce, dyn_corr$assay)
    counts <- Matrix::t(counts[genes, ])

    # Out of the time points at which correlation coefficients were estimated,
    # use the closest one to the specified time.
    metacell_pseudotimes <- metadata(sce)$dyn_corr$metacell_sce[[time_col]]
    min_ix <- sapply(new_times, function(x) {
        which.min(abs(x - metacell_pseudotimes))
    })

    message("Sampling copula")
    # Each step is very quick, but we generally expect new_times to be large. So
    # a progress bar is still needed, but we'll only update the progress bar
    # every 100 steps to avoid it flashing due to rapid updates.
    nsteps <- ceiling(length(min_ix) / 100)
    progressr::with_progress({
        pbar <- progressr::progressor(steps = nsteps)
        U <- lapply(seq_along(min_ix), function(i) {
            r <- dyn_corr$rho[min_ix[i], , drop = TRUE]
            copula <- copula::normalCopula(
                param = r,
                dim = length(genes),
                dispstr = "un"
            )
            V <- copula::rCopula(1, copula)
            if (i %% 100 == 0) {
                pbar()
            }
            return(V)
        })
    })
    U <- do.call(rbind, U)

    message("Converting margins")
    new_times_df <- stats::setNames(data.frame(new_times), dyn_corr$time_col)
    pseudotimes <- sce[[time_col]]
    progressr::with_progress({
        pbar <- progressr::progressor(steps = length(genes))
        counts_sim <- lapply(seq_along(genes), function(i) {
            mdat <- data.frame(counts[, i], pseudotimes)
            names(mdat) <- c("x", dyn_corr$time_col)
            mfun <- dyn_corr$margins[[i]]

            # Get model parameters at new pseudotimes
            par_pred <- gamlss::predictAll(
                object = mfun,
                newdata = new_times_df,
                type = "response",
                data = mdat
            )

            # Load the correct quantile function
            qfun <- getExportedValue("gamlss.dist", paste0("q", mfun$family[1]))
            # Transform margins
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
    colnames(sce_sim) <- paste0("sim", seq_along(new_times))
    rownames(sce_sim) <- genes
    sce_sim$pseudotime <- new_times
    return(sce_sim)
}

###############################################################################
