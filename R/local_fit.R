###############################################################################

#' Fit a dynamic Gaussian copula model
#'
#' @description Fit a time-varying Gaussian copula to a time series.
#'
#' @param NX Matrix of normal-transformed pseudo-observations at time points.
#' @param x Vector of time points corresponding to \code{NX}.
#' Must be sorted and have no duplicates.
#' @param x0 Time points to estimate copula parameters at.
#' @param h Kernel bandwidth. Must satisfy \code{0 < h < 1}.
#' @param degree Degree of local polynomial approximation. Default is \code{0}.
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
#'   \item \code{NX}: The input argument \code{NX}.
#'   \item \code{eta}: Matrix of estimated coefficients in the unconstrained
#'   space.
#'   \item \code{rho}: Matrix of estimated pairwise correlation coefficients.
#' }
#'
#' @export
fit_dynamic_gaussian <- function(NX,
                                 x,
                                 x0,
                                 h,
                                 degree = 0L,
                                 control = list(),
                                 cores = 1L) {
    # Basic checks
    assert_that(
        is.vector(x, mode = "numeric"),
        is.numeric(NX) && is.matrix(NX),
        dim(NX)[1] == length(x),
        is.numeric(h) && h > 0 && h < 1,
        is.numeric(cores) && cores >= 1,
        is.list(control),
        !anyDuplicated(x) && !is.unsorted(x),
        is.numeric(x0),
        is.numeric(degree) && degree >= 0
    )
    cores <- as.integer(cores)
    degree <- as.integer(degree)

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

    # Verify control parameters
    assert_that(
        all(sapply(control, is.numeric)),
        control$max_epoch >= 1L,
        control$max_itr >= 1L,
        control$history_size >= 1L,
        control$tolerance_grad > 0,
        control$tolerance_change > 0
    )

    # Scale times to [0, 1]
    min_x <- x[1]
    dx <- x[length(x)] - min_x
    x <- (x - min_x) / dx
    x0 <- (x0 - min_x) / dx

    if (cores > 1L) {
        # Set up futures plan
        cl <- parallel::makeCluster(cores)
        future::plan(future::cluster, workers = cl)
        on.exit({ future::plan(future::sequential); parallel::stopCluster(cl) },
                add = TRUE)

        # This environment stores the fit_gaussian function once loaded from the
        # Python module to avoid having to keep reloading it. The module is
        # loaded one on each worker.
        .fit_env <- new.env(parent = emptyenv())
        fit_gaussian <- function(par0, x, NX, h, control, x0) {
            # Only load the module once per worker
            if (!exists("fit_fun", envir = .fit_env, inherits = FALSE)) {
                .fit_env$fit_fun <- reticulate::import_from_path(
                    module = "dynamic_gaussian",
                    path = system.file("python", package = "DynCopula"),
                    delay_load = FALSE
                )$fit_gaussian
            }
            .fit_env$fit_fun(
                par0 = par0,
                x = x,
                NX = NX,
                h = h,
                control = control,
                x0 = x0,
                degree = degree
            )
        }

        # Parallelized with progress bar
        progressr::with_progress({
            pbar <- progressr::progressor(along = x0)
            eta_est <- future.apply::future_lapply(
                X = x0,
                FUN = function(t0) {
                    y <- fit_gaussian(
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
                future.globals = TRUE
            )
        })
    } else {
        # Load Python module
        fit_fun <- reticulate::import_from_path(
            module = "dynamic_gaussian",
            path = system.file("python", package = "DynCopula"),
            delay_load = FALSE
        )$fit_gaussian

        progressr::with_progress({
            pbar <- progressr::progressor(along = x0)
            eta_est <- lapply(
                X = x0,
                FUN = function(t0) {
                    y <- fit_fun(
                        par0 = .init_par(t0, x, NX, h),
                        x = x,
                        NX = NX,
                        h = h,
                        control = control,
                        x0 = t0,
                        degree = degree
                    )
                    pbar()
                    return(y)
                }
            )
        })
    }

    # Estimated eta matrix
    Hhat <- do.call(rbind, eta_est)
    # Estimated correlation matrix
    Rhat <- t(apply(Hhat, 1, function(v) {
        copula::P2p(vec2cor(v))
    }))

    # Needed to ensure consistent shape of Rhat across all dimensions
    d <- dim(NX)[2]
    if (d == 2) {
        Rhat <- t(Rhat)
    }

    # Add numbered eta/rho labels
    colnames(Hhat) <- paste0("eta", seq_len(choose(d, 2)))
    ix_lab <- apply(utils::combn(seq_len(d), 2), 2, function(x) {
        paste(x, collapse = '_')
    })
    colnames(Rhat) <- paste0("rho", ix_lab)

    # Rescale times back to original scale
    x <- x * dx + min_x
    x0 <- x0 * dx + min_x

    return(list(
        x = x,
        x0 = x0,
        h = h,
        NX = NX,
        eta = Hhat,
        rho = Rhat
    ))
}

###############################################################################

#' Bandwidth selection for a dynamic Gaussian copula model
#'
#' @description Select the optimal kernel bandwidth using leave-one-out cross
#' validation (LOOCV).
#'
#' @param NX Matrix of normal-transformed pseudo-observations at time points.
#' @param x Vector of time points corresponding to \code{NX}.
#' Must be sorted and have no duplicates.
#' @param bandwidths Vector of kernel bandwidths to test.
#' @param xind Number of points for LOOCV.
#' @param degree Degree of local polynomial approximation. Default is \code{0}.
#' @param control A \code{list} of control parameters for optimization.
#' @param cores Number of cores to use. Parallelized over \code{bandwidths}.
#' Default is \code{1}.
#'
#' @details See \link[DynCopula]{fit_dynamic_gaussian} for details on parameter
#' estimation.

#' @return A \code{data.frame} specifying the LOOCV log-likelihood at each
#' bandwidth value.
#'
#' @export
bandwidth_select_cv <- function(NX,
                                x,
                                bandwidths,
                                xind,
                                degree = 0L,
                                control = list(),
                                cores = 1L) {
    assert_that(
        is.vector(x, mode = "numeric"),
        is.numeric(NX) && is.matrix(NX),
        dim(NX)[1] == length(x),
        is.numeric(bandwidths) && all(bandwidths > 0) && all(bandwidths < 1),
        !anyDuplicated(x),
        is.numeric(xind) && xind > 1,
        is.numeric(degree) && degree >= 0
    )

    # Use xind equally spaced covariate values
    xind <- as.integer(xind)
    xind <- unique(floor(seq(1, length(x), length.out = xind)))

    # Compute log-likelihood of eta at ith observation under a Gaussian copula
    # model
    loglik <- function(NX_i, eta) {
        copula_log_dens <- mvtnorm::dmvnorm(
            x = NX_i,
            sigma = vec2cor(eta),
            log = TRUE,
            checkSymmetry = FALSE
        )
        margin_log_dens <- sum(stats::dnorm(NX_i, log = TRUE))
        return(copula_log_dens - margin_log_dens)
    }

    ll_all <- numeric(length = length(bandwidths))
    if (cores > 1L) {
        # Set up futures plan
        cl <- parallel::makeCluster(cores)
        future::plan(future::cluster, workers = cl)
        on.exit({ future::plan(future::sequential); parallel::stopCluster(cl) },
                add = TRUE)

        for (i in seq_along(bandwidths)) {
            h <- bandwidths[i]
            message("Testing bandwidth = ", h)

            # Parallelized with progress bar
            progressr::with_progress({
                pbar <- progressr::progressor(along = xind)
                ll_h <- future.apply::future_lapply(
                    X = xind,
                    FUN = function(ix) {
                        # Estimate copula parameters when leaving out
                        # observation at ix
                        eta <- fit_dynamic_gaussian(
                            NX = NX[-ix, ],
                            x = x[-ix],
                            x0 = x[ix],
                            h = h,
                            degree = degree,
                            control = control
                        )$eta
                        # Log-likelihood of eta at observation ix
                        ll <- loglik(NX[ix, ], as.vector(eta))
                        pbar()
                        return(ll)
                    },
                    future.seed = TRUE,
                    future.globals = c("x", "NX", "h", "control")
                )
            })
            ll_all[i] <- sum(unlist(ll_h))
        }
    } else {
        for (i in seq_along(bandwidths)) {
            h <- bandwidths[i]
            message("Testing bandwidth = ", h)

            progressr::with_progress({
                pbar <- progressr::progressor(along = xind)
                ll_h <- lapply(
                    X = xind,
                    FUN = function(ix) {
                        # Estimate copula parameters when leaving out
                        # observation at ix
                        eta <- fit_dynamic_gaussian(
                            NX = NX[-ix, ],
                            x = x[-ix],
                            x0 = x[ix],
                            h = h,
                            degree = degree,
                            control = control
                        )$eta
                        # Log-likelihood of eta at observation ix
                        ll <- loglik(NX[ix, ], as.vector(eta))
                        pbar()
                        return(ll)
                    }
                )
            })
            ll_all[i] <- sum(unlist(ll_h))
        }
    }
    return(data.frame(bandwidth = bandwidths, ll = ll_all))
}

###############################################################################
