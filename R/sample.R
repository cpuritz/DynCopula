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
#' are interpolated at new times using either linear or spline interpolation.
#'
#' @export
sample_cells <- function(sce,
                         new_times,
                         interpolation = c("linear", "spline")) {
    assert_that(
        methods::is(sce, "SingleCellExperiment"),
        "copula_fit" %in% names(metadata(sce)),
        all(c("margins", "rho") %in% names(metadata(sce)$copula_fit)),
        is.numeric(new_times)
    )
    interpolation <- match.arg(interpolation)

    if (!metadata(sce)$copula_fit$full_margins) {
        stop("Cells can't be sampled since the margins were not saved. ",
             "Rerun 'fit_margins' with 'save = TRUE'.")
    }

    copula_fit <- metadata(sce)$copula_fit
    genes <- copula_fit$features
    time_col <- copula_fit$time_col

    if (copula_fit$assay == "logcounts") {
        stop("NOT IMPLEMENTED YET")
    }

    counts <- SummarizedExperiment::assay(sce, copula_fit$assay)
    counts <- Matrix::t(counts[genes, ])

    # Interpolate calibration coefficients at new times
    pseudotimes <- sce[[time_col]]
    if (interpolation == "spline") {
        interp_fun <- function(y) {
            spline_fit <- splines::interpSpline(pseudotimes, y)
            stats::predict(spline_fit, new_times)$y
        }
    } else {
        interp_fun <- function(y) {
            stats::approx(
                x = pseudotimes,
                y = y,
                xout = new_times,
                method = "linear"
            )$y
        }
    }
    eta_int <- apply(copula_fit$eta, 2, interp_fun)

    message("Sampling copula")
    # Each step is very quick, but we generally expect new_times to be large.
    # So a progress bar is still needed, but we'll only update the progress
    # bar every 100 steps to avoid it flashing due to rapid updates.
    dstep <- 100
    nsteps <- ceiling(length(new_times) / dstep)
    progressr::with_progress({
        pbar <- progressr::progressor(steps = nsteps)
        U <- lapply(seq_along(new_times), function(i) {
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
    new_times_df <- stats::setNames(data.frame(new_times), time_col)
    progressr::with_progress({
        pbar <- progressr::progressor(steps = length(genes))
        counts_sim <- lapply(seq_along(genes), function(i) {
            mdat <- data.frame(counts[, i], pseudotimes)
            names(mdat) <- c("x", time_col)
            mfun <- copula_fit$margins[[i]]

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
