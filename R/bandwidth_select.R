###############################################################################

#' Bandwidth selection for a dynamic Gaussian copula model
#'
#' @description Select the optimal kernel bandwidth for a dynamic Gaussian
#' copula model using leave-one-out cross validation (LOOCV).
#'
#' @param FX Matrix of pseudo-observations at covariate values.
#' @param x Vector of covariate values corresponding to \code{FX}. Must be
#' sorted and have no duplicates.
#' @param bandwidths Vector of kernel bandwidths to test.
#' @param xind Number of covariate values to use for LOOCV. Default is
#' \code{length(x)}.
#' @param degree Degree of local polynomial approximation. Default is \code{0}.
#' @param control A \code{list} of control parameters for optimization.
#' @param cores Number of cores to use. Parallelized over \code{bandwidths}.
#' Default is \code{1}.
#'
#' @details See \link[DynCopula]{fit_dynamic_gaussian} for details on parameter
#' estimation.
#'
#' The argument \code{xind} specifies the number of covariate values to use for
#' LOOCV. Full LOOCV corresponds to \code{xind = length(x)}. If
#' \code{xind < length(x)}, then only a subset of the covariates are used to
#' reduce the run time. The indices of the chosen covariates are equally spaced
#' along \code{seq_along(x)}. This is only an estimate to full LOOCV and may be
#' quite inaccurate for \code{xind << length(x)}.
#'
#' The optimal bandwidth is the one that maximizes the cross-validated
#' likelihood criterion.
#'
#' @return A \code{data.frame} specifying the cross-validated likelihood
#' criterion at each bandwidth value.
#'
#' @export
bandwidth_select <- function(FX,
                             x,
                             bandwidths,
                             xind = length(x),
                             degree = 0L,
                             control = list(),
                             cores = 1L) {
    assert_that(
        is.vector(x, mode = "numeric"),
        is.matrix(FX),
        is.numeric(FX),
        dim(FX)[1] == length(x),
        is.numeric(bandwidths) && all(bandwidths > 0) && all(bandwidths < 1),
        !anyDuplicated(x),
        is.numeric(xind) && xind > 1,
        is.numeric(degree) && degree >= 0
    )

    # Use xind equally spaced covariate values
    xind <- as.integer(xind)
    if (xind == length(x)) {
        xind <- seq_along(x)
    } else {
        xind <- unique(floor(seq(1, length(x), length.out = xind)))
    }

    # Log-likelihood of eta given X under a Gaussian copula model
    loglik <- function(X, eta) {
        NX <- stats::qnorm(X)
        copula_ll <- mvtnorm::dmvnorm(
            x = NX,
            sigma = vec2cor(eta),
            log = TRUE,
            checkSymmetry = FALSE
        )
        margin_ll <- sum(stats::dnorm(NX, log = TRUE))
        return(copula_ll - margin_ll)
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
                        # Estimate calibration coefficients when leaving out
                        # observation at ix
                        eta <- fit_dynamic_gaussian(
                            FX = FX[-ix, ],
                            x = x[-ix],
                            x0 = x[ix],
                            h = h,
                            degree = degree,
                            control = control
                        )$eta
                        # Log-likelihood of estimated calibration coefficients
                        # at observation ix
                        ll <- loglik(FX[ix, ], as.vector(eta))
                        pbar()
                        return(ll)
                    },
                    future.seed = TRUE,
                    future.globals = TRUE
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
                        # Estimate calibration coefficients when leaving out
                        # observation at ix
                        eta <- fit_dynamic_gaussian(
                            FX = FX[-ix, ],
                            x = x[-ix],
                            x0 = x[ix],
                            h = h,
                            degree = degree,
                            control = control
                        )$eta
                        # Log-likelihood of estimated calibration coefficients
                        # at observation ix
                        ll <- loglik(FX[ix, ], as.vector(eta))
                        pbar()
                        return(ll)
                    }
                )
            })
            ll_all[i] <- sum(unlist(ll_h))
        }
    }
    return(data.frame(bandwidth = bandwidths, cv = ll_all))
}

###############################################################################
