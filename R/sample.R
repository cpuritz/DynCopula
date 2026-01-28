###############################################################################

#' Sample cells along trajectory
#'
#' @description Sample cells along the trajectory at specified pseudotimes.
#'
#' @param sce A \code{SingleCellExperiment}.
#' @param new_times A vector of pseudotimes to sample at.
#' @param interpolation Interpolation method, either linear or spline.
#'
#' @returns A \code{SingleCellExperiment}.
#'
#' @details To sample from the copula at the new times, calibration coefficients
#' are interpolated at new times using either linear or spline
#' interpolation. The copula is then sampled from at each new time, and the
#' margins are converted to those fit to the original counts.
#'
#' @export
sample_cells <- function(sce,
                         new_times,
                         interpolation = c("linear", "spline")) {
    assert_that(
        methods::is(sce, "SingleCellExperiment"),
        "dyn_corr" %in% names(metadata(sce)),
        all(c("margins", "rho") %in% names(metadata(sce)$dyn_corr)),
        is.numeric(new_times)
    )
    interpolation <- match.arg(interpolation)

    dyn_corr <- metadata(sce)$dyn_corr
    genes <- dyn_corr$features
    time_col <- dyn_corr$time_col

    counts <- SummarizedExperiment::assay(sce, dyn_corr$assay)
    counts <- Matrix::t(counts[genes, ])

    # Interpolate calibration coefficients at new times
    metacell_pseudotimes <- dyn_corr$metacell_sce[[time_col]]
    if (interpolation == "spline") {
        interp_fun <- function(y) {
            spline_fit <- splines::interpSpline(metacell_pseudotimes, y)
            stats::predict(spline_fit, new_times)$y
        }
    } else {
        interp_fun <- function(y) {
            stats::approx(
                x = metacell_pseudotimes,
                y = y,
                xout = new_times,
                method = "linear"
            )$y
        }
    }
    eta_int <- apply(dyn_corr$eta, 2, interp_fun)

    message("Sampling copula")
    # Each step is very quick, but we generally expect new_times to be large. So
    # a progress bar is still needed, but we'll only update the progress bar
    # every 100 steps to avoid it flashing due to rapid updates.
    dstep <- 100
    nsteps <- ceiling(length(metacell_pseudotimes) / dstep)
    progressr::with_progress({
        pbar <- progressr::progressor(steps = nsteps)
        U <- lapply(seq_along(metacell_pseudotimes), function(i) {
            # Convert interpolated calibration coefficients to correlation
            # coefficients
            R <- vec2cor(eta_int[i, ])
            copula <- copula::normalCopula(
                param = copula::P2p(R),
                dim = length(genes),
                dispstr = "un"
            )
            V <- copula::rCopula(1, copula)
            if (i %% dstep == 0) {
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

    # Create a SingleCellExperiment with the simulated counts
    counts_sim <- Matrix::t(do.call(cbind, counts_sim))
    counts_sim <- methods::as(counts_sim, "dgCMatrix")
    sce_sim <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = counts_sim)
    )
    colnames(sce_sim) <- paste0("sim", seq_along(new_times))
    rownames(sce_sim) <- genes
    sce_sim$pseudotime <- new_times

    return(sce_sim)
}

###############################################################################
