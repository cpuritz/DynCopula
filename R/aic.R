###############################################################################

#' Model AIC
#'
#' @description Compute the model AIC. Requires that coefficients were computed
#' at all time points.
#'
#' @param res A \code{list} or a \code{SingleCellExperiment}.
#' @param cores Number of cores to use. Default is \code{1}. Only used when
#' \code{res} is a \code{list}.
#' @param ... Other arguments.
#'
#' @details If \code{res} is a \code{list}, it must be the output of
#' \link[DynCopula]{fit_dynamic_gaussian}. In this case, correlation
#' coefficients must have been estimated at all time points (i.e.
#' \link[DynCopula]{fit_dynamic_gaussian} must have been run with
#' \code{x0 = x}). If \code{res} is a \code{SingleCellExperiment}, is must be
#' the output of either \link[DynCopula]{fit_dyn_corr} or
#' \link[DynCopula]{fit_dyn_corr_sel}.
#'
#' @returns If \code{res} is a \code{list}, returns the model AIC. If \code{res}
#' is a \code{SingleCellExperiment}, returns the input
#' \code{SingleCellExperiment} but with the metadata entry \code{dyn_corr}
#' updated to include a new entry \code{aic} storing the model AIC.
#'
#' @export
model_aic <- function(res, ...) {
    UseMethod("model_aic")
}

###############################################################################

#' @describeIn model_aic Method for \code{list}.
#' @export
model_aic.list <- function(res, cores = 1L, ...) {
    assert_that(
        all(c("x", "x0", "h", "eta", "NX") %in% names(res)),
        all(res$x == res$x0),
        cores >= 1L
    )

    cl <- parallel::makeCluster(as.integer(cores))
    future::plan(future::cluster, workers = cl)
    on.exit({ future::plan(future::sequential); parallel::stopCluster(cl) },
            add = TRUE)

    progressr::with_progress({
        pbar <- progressr::progressor(along = res$x)
        aic <- future.apply::future_lapply(
            X = seq_along(res$x),
            FUN = function(i) {
                y <- py_load("dynamic_gaussian")$model_aic(
                    eta_i = res$eta[i, ],
                    x = res$x,
                    NX = res$NX,
                    h = res$h,
                    i = i - 1
                )
                pbar()
                return(y)
            },
            future.seed = TRUE,
            future.globals = c("res", "pbar")
        )
    })
    return(sum(unlist(aic)))
}

###############################################################################

#' @describeIn model_aic Method for \code{SingleCellExperiment}.
#' @export
model_aic.SingleCellExperiment <- function(res, ...) {
    assert_that(
        "dyn_corr" %in% names(metadata(res)),
        "margins" %in% names(metadata(res)$dyn_corr)
    )

    message("Computing model AIC")
    dyn_corr <- metadata(res)$dyn_corr
    x <- res[[dyn_corr$time_col]]
    cores <- dyn_corr$cores

    FX <- dyn_corr$FX
    FXm <- dyn_corr$FXm
    V <- dyn_corr$V
    NX <- stats::qnorm(FXm + (FX - FXm) * V)

    dat <- list(
        x = x,
        x0 = x,
        h = dyn_corr$h,
        eta = dyn_corr$eta_int,
        NX = NX
    )
    metadata(res)$dyn_corr$aic <- model_aic(dat, cores)
    return(res)
}

###############################################################################
