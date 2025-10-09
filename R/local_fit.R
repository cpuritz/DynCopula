###############################################################################

#' Fit a dynamic Gaussian copula model
#'
#' @description Fit a time-varying Gaussian copula to a time series.
#'
#' @param FX Matrix of pseudo-observations at time points.
#' @param FXm Matrix of left limits of pseudo-observations at time points. If
#' \code{NULL}, data is assumed to be continuous. Otherwise, data is assumed
#' to be count-valued.
#' @param x Vector of time points corresponding to \code{FX} and \code{FXm}.
#' Must be sorted and have no duplicates.
#' @param x0 Time points to estimate copula parameters at.
#' @param h Kernel bandwidth. Must satisfy \code{0 < h < 1}.
#' @param control A \code{list} of control parameters for optimization.
#' @param cores Number of cores to use for parallel optimization. Parallelized
#' over \code{x0}. Default is \code{1}.
#'
#' @details Optimization is performed using L-BFGS. The \code{control} argument
#' is a list that supplies control parameters for optimization. The following
#' parameters can be supplied:
#' \itemize{
#'   \item \code{max_epoch} Maximum number of epochs. Default is \code{1}.
#'   \item \code{max_itr} Maximum number of internal iterations. Default is
#'   \code{100}.
#'   \item \code{history_size} History size. Default is \code{30}.
#'   \item \code{tolerance_grad} Termination tolerance for gradient. Default is
#'   \code{1e-7}.
#'   \item \code{tolerance_change} Termination tolerance for log-likelihood.
#'   Default is \code{1e-9}.
#' }
#'
#' @return A list with the following components:
#' \itemize{
#'   \item \code{x}: The input argument \code{x}.
#'   \item \code{x0}: The input argument \code{x0}.
#'   \item \code{h}: The input argument \code{h}.
#'   \item \code{NX}: The normal-transformed pseudo-observations.
#'   \item \code{eta}: Matrix of estimated coefficients in the unconstrained
#'   space.
#'   \item \code{rho}: Matrix of estimated pairwise correlation coefficients.
#' }
#'
#' @export
fit_dynamic_gaussian <- function(FX,
                                 FXm = NULL,
                                 x,
                                 x0,
                                 h,
                                 control = list(),
                                 cores = 1L) {
    # Basic checks
    assert_that(
        is.vector(x, mode = "numeric"),
        is.numeric(FX) && is.matrix(FX),
        is.null(FXm) || (is.numeric(FXm) && is.matrix(FXm)),
        dim(FX)[1] == length(x),
        is.null(FXm) || all(dim(FX) == dim(FXm)),
        is.numeric(h) && h > 0 && h < 1,
        is.numeric(cores) && cores >= 1,
        is.list(control),
        !anyDuplicated(x) && !is.unsorted(x),
        is.numeric(x0)
    )
    cores <- as.integer(cores)

    # Default control parameters
    defaults <- list(
        max_epoch = 1L,
        max_itr = 100L,
        history_size = 30L,
        tolerance_grad = 1e-7,
        tolerance_change = 1e-9
    )
    control <- utils::modifyList(defaults, control)
    assert_that(all(names(control) %in% names(defaults)))
    control$max_itr <- as.integer(control$max_itr)
    control$patience <- as.integer(control$patience)

    # Verify control parameters
    assert_that(
        all(sapply(control, is.numeric)),
        control$max_epoch >= 1L,
        control$max_itr >= 1L,
        control$history_size >= 1L,
        control$tolerance_grad > 0,
        control$tolerance_change > 0
    )

    # If data is count-valued, jitter to produce continuous data
    if (!is.null(FXm)) {
        FX <- .jitter(FX, FXm)
    }

    # Scale times to [0, 1]
    min_x <- x[1]
    dx <- x[length(x)] - min_x
    x <- (x - min_x) / dx
    x0 <- (x0 - min_x) / dx

    # Convert to standard normal margins
    NX <- stats::qnorm(FX)

    # Set up futures plan
    cl <- parallel::makeCluster(cores)
    future::plan(future::cluster, workers = cl)
    on.exit({ future::plan(future::sequential); parallel::stopCluster(cl) },
            add = TRUE)

    # Parallelized with progress bar
    progressr::with_progress({
        pbar <- progressr::progressor(along = x0)
        eta_est <- future.apply::future_lapply(
            X = x0,
            FUN = function(t0) {
                y <- py_load("dynamic_gaussian")$fit_gaussian(
                    par0 = .init_par(t0, x, NX, h),
                    x = x,
                    NX = NX,
                    h = h,
                    control = control,
                    x0 = t0
                )
                pbar()
                return(y)
            },
            future.seed = TRUE,
            future.globals = c("x", "NX", "h", "control", "pbar")
        )
    })

    # Estimated eta matrix
    Hhat <- do.call(rbind, eta_est)
    # Estimated correlation matrix
    Rhat <- t(apply(Hhat, 1, function(v) {
        copula::P2p(vec2cor(v))
    }))

    # Needed to ensure consistent shape of Rhat across all dimensions
    d <- dim(FX)[2]
    if (d == 2) {
        Rhat <- t(Rhat)
    }

    # Add numbered eta/rho labels
    ix_lab <- apply(utils::combn(seq_len(d), 2), 2, function(x) {
        paste(x, collapse = '_')
    })
    colnames(Hhat) <- paste0("eta", ix_lab)
    colnames(Rhat) <- paste0("rho", ix_lab)

    # Rescale times back to original scale
    x <- x * dx + min_x
    x0 <- x0 * dx + min_x

    return(list(
        NX = NX,
        x = x,
        x0 = x0,
        h = h,
        eta = data.frame(Hhat),
        rho = data.frame(Rhat)
    ))
}

###############################################################################

#' Bandwidth selection for a dynamic Gaussian copula model
#'
#' @description Select the optimal bandwidth for fitting a time-varying Gaussian
#' copula.
#'
#' @param FX Matrix of pseudo-observations at time points.
#' @param FXm Matrix of left limits of pseudo-observations at time points. If
#' \code{NULL}, data is assumed to be continuous. Otherwise, data is assumed
#' to be count-valued.
#' @param x Vector of time points corresponding to \code{FX} and \code{FXm}.
#' Must be sorted and have no duplicates.
#' @param bandwidths Vector of kernel bandwidths to test.
#' @param control A \code{list} of control parameters for optimization.
#' @param cores Number of cores to use for parallel optimization. Default is
#' \code{1}.
#' @param return_all Whether all models should be returned, or just the best.
#' Default is \code{FALSE}.
#'
#' @details See \link[DynCopula]{fit_dynamic_gaussian} for details.

#' @return A list with the following components:
#' \itemize{
#'   \item \code{x}: The input argument \code{x}.
#'   \item \code{x0}: The time points at which coefficients were estimated.
#'   \item \code{h}: The optimal bandwidth selected.
#'   \item \code{eta}: Matrix of coefficients in the unconstrained space
#'   estimated using optimal bandwidth. If \code{return_all = TRUE}, this will
#'   be a list of matrices, one for each bandwidth.
#'   \item \code{rho}: Matrix of pairwise correlation coefficients estimated
#'   using optimal bandwidth. If \code{return_all = TRUE}, this will
#'   be a list of matrices, one for each bandwidth.
#'   \item \code{bandwidths}: The input argument \code{bandwidths}.
#'   \item \code{aic}: Vector of AIC values at each bandwidth.
#' }
#'
#' @export
bandwidth_select <- function(FX,
                             FXm = NULL,
                             x,
                             bandwidths,
                             control = list(),
                             cores = 1L,
                             return_all = FALSE) {
    assert_that(
        is.vector(x, mode = "numeric"),
        is.numeric(FX) && is.matrix(FX),
        is.null(FXm) || (is.numeric(FXm) && is.matrix(FXm)),
        dim(FX)[1] == length(x),
        is.null(FXm) || all(dim(FX) == dim(FXm)),
        is.numeric(bandwidths),
        !anyDuplicated(x),
        is.logical(return_all)
    )

    # If data, is count-valued, jitter to produce continuous data
    if (!is.null(FXm)) {
        FX <- .jitter(FX, FXm)
    }

    if (return_all) {
        res_list <- list()
        for (i in seq_along(bandwidths)) {
            message("Testing bandwidth = ", bandwidths[i])
            res <- fit_dynamic_gaussian(
                FX = FX,
                x = x,
                x0 = x,
                h = bandwidths[i],
                control = control,
                cores = cores
            )
            aic[i] <- model_aic(res, cores)
            res_list[[i]] <- res
        }
        aic <- sapply(res_list, '[[', "aic")
        return(list(
            x = x,
            x0 = x,
            h = bandwidths[which.min(aic)],
            eta = lapply(res_list, '[[', "eta"),
            rho = lapply(res_list, '[[', "rho"),
            bandwidths = bandwidths,
            aic = aic
        ))
    } else {
        # Select model with lowest AIC
        aic <- numeric(length = length(bandwidths))
        best_aic <- Inf
        best_res <- NULL
        for (i in seq_along(bandwidths)) {
            message("Testing bandwidth = ", bandwidths[i])
            res <- fit_dynamic_gaussian(
                FX = FX,
                x = x,
                x0 = x,
                h = bandwidths[i],
                control = control,
                cores = cores
            )
            aic[i] <- model_aic(res, cores)

            if (aic[i] <= best_aic) {
                best_aic <- aic[i]
                best_res <- res
            }
        }
        best_res$aic <- aic
        best_res$bandwidths <- bandwidths
        return(best_res)
    }
}

###############################################################################

#' Model AIC
#'
#' @description Compute the model AIC. Requires that coefficients were computed
#' at all time points.
#'
#' @param res Output of \link[DynCopula]{fit_dynamic_gaussian}.
#' @param cores Number of cores to use. Default is \code{1}.
#'
#' @details Correlation coefficients must have been estimated at all time
#' points. That is, \link[DynCopula]{fit_dynamic_gaussian} must have been run
#' with \code{x0 = x}.
#'
#' @returns The AIC of the model.
#'
#' @export
model_aic <- function(res,
                      cores = 1L) {
    assert_that(
        is.list(res) && all(c("x", "x0", "h", "eta", "NX") %in% names(res)),
        all(res$x == res$x0),
        cores >= 1L
    )
    cores <- as.integer(cores)

    eta <- as.matrix(res$eta)
    h <- res$h
    x <- res$x
    NX <- res$NX

    cl <- parallel::makeCluster(cores)
    future::plan(future::cluster, workers = cl)
    on.exit({ future::plan(future::sequential); parallel::stopCluster(cl) },
            add = TRUE)

    progressr::with_progress({
        pbar <- progressr::progressor(along = x)
        aic <- future.apply::future_lapply(
            X = seq_along(x),
            FUN = function(i) {
                y <- py_load("dynamic_gaussian")$model_aic(
                    eta_i = eta[i, ],
                    x = x,
                    NX = NX,
                    h = h,
                    i = i - 1
                )
                pbar()
                return(y)
            },
            future.seed = TRUE,
            future.globals = c("eta", "x", "NX", "h", "pbar")
        )
    })
    return(sum(unlist(aic)))
}

###############################################################################
