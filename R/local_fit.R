###############################################################################

#' Local Fit Gaussian
#'
#' @description Fit a dynamic Gaussian copula to time series data.
#'
#' @param FX Matrix of pseudo-observations F_i(X_i) at time points.
#' @param FXm Matrix of left limits of distribution functions F_i(X_i - 1) at
#' time points. If \code{NULL}, it is assumed that all margins are continuous.
#' @param x Vector of time points corresponding to \code{FX} and \code{FXm}.
#' @param x0 Time points to estimate copula parameters at.
#' @param band Kernel bandwidth.
#' @param scale Scale factor used for unconstrained parametrization of
#' correlation matrices.
#' @param R0 Fraction of neighbors to use to estimate the initial correlation
#' matrix. Default is \code{0.05}. Must be between \code{0} (exclusive) and
#' \code{1} (inclusive).
#' @param optMethod Optimization method. Either \code{"SGD"} or \code{"L-BFGS"}.
#' \code{"SGD"} implements stochastic gradient descent with gradients computed
#' using automatic differentiation. \code{"L-BFGS"} implements the quasi-Newton
#' limited memory BFGS with gradients estimated numerically.
#' @param control A \code{list} of control parameters for optimization.
#' @param cores Number of cores to use for parallel optimization. Parallelized
#' over \code{x0}. Default is \code{1}.
#' @param cores2 If using \code{"L-BFGS"}, each optimization step can be
#' parallelized, in addition to parallelization over \code{x0}. The total
#' number of cores used is \code{cores * cores2}. Default is \code{1}.
#'
#' @details The \code{control} argument is a list that supplies control
#' parameters for optimization. For SGD, the following parameters can be
#' supplied:
#' \itemize{
#'   \item \code{maxit} Maximum number of iterations. Default is \code{100}.
#'   \item \code{lr} Learning rate for RMSProp. Default is \code{1e-2}.
#'   \item \code{reltol} Relative convergence tolerance.
#'   \item \code{patience} Optimization stops if the relative log-likelihood
#'   has not decreased by a factor of \code{reltol} within the last
#'   \code{patience} iterations. Default is \code{5}.
#'   \item \code{weight_decay} Weight decay for RMSProp. Default is \code{0}.
#' }
#' For L-BFGS, the available control parameters and default values are the same
#' as those of the \code{\link[stats]{optim}} function from the \strong{stats}
#' package.
#'
#' @return \itemize{
#'   \item \code{x}: The input argument \code{x}.
#'   \item \code{x0}: The input argument \code{x0}.
#'   \item \code{eta}: Matrix of estimated coefficients in the unconstrained
#'   space.
#'   \item \code{rho}: Matrix of estimated pairwise correlation coefficients.
#'   \item \code{scale}: The input argument \code{scale}.
#'   \item \code{convergence}: Convergence codes for each coefficient.
#'   \code{0} indicates successful completion. \code{1} indicates that iteration
#'   limit had been reached.
#'   \item \code{loss}: Loss history for each coefficient. Only available if
#'   \code{optMethod} is \code{"SGD"}.
#' }
#'
#' @export
local_fit_gaussian <- function(FX,
                               FXm = NULL,
                               x,
                               x0,
                               band,
                               scale,
                               R0 = 0.05,
                               optMethod = c("L-BFGS", "SGD"),
                               control = list(),
                               cores = 1L,
                               cores2 = 1L) {
    optMethod <- match.arg(optMethod)

    # Basic checks
    assertthat::assert_that(
        dim(FX)[1] == length(x),
        is.null(FXm) || all(dim(FX) == dim(FXm)),
        min(x0) >= min(x),
        max(x0) <= max(x),
        is.numeric(band) && band > 0,
        is.numeric(scale) && scale > 0,
        is.numeric(R0) && R0 > 0 && R0 <= 1,
        is.character(optMethod),
        is.numeric(cores) && cores >= 1,
        is.numeric(cores2) && cores2 >= 1
    )

    cores <- as.integer(cores)
    cores2 <- as.integer(cores2)
    if (cores2 > 1L && optMethod == "SGD") {
        message("'cores2' is not used for SGD optimization.")
    }

    # Control parameters
    assertthat::assert_that(all(sapply(control, is.numeric)))
    if (!"maxit" %in% names(control)) {
        control <- c(list(maxit = 100L), control)
    } else {
        control$maxit <- as.integer(control$maxit)
        assertthat::assert_that(control$maxit >= 1)
    }
    if (optMethod == "SGD") {
        if (!"lr" %in% names(control)) {
            control <- c(list(lr = 1e-2), control)
        }
        if (!"reltol" %in% names(control)) {
            control <- c(list(reltol = 1e-5), control)
        }
        if (!"patience" %in% names(control)) {
            control <- c(list(patience = 5L), control)
        } else {
            control$patience <- as.integer(control$patience)
        }
        if (!"weight_decay" %in% names(control)) {
            control <- c(list(weight_decay = 0), control)
        }
        assertthat::assert_that(
            control$lr > 0,
            control$reltol > 0,
            control$patience >= 1L,
            control$weight_decay >= 0
        )
    }

    # Sort covariate values and scale to [0, 1]
    ord <- order(x)
    x <- x[ord]
    FX <- FX[ord, ]
    if (!is.null(FXm)) {
        FXm <- FXm[ord, ]
    }
    min_x <- x[1]
    max_x <- x[length(x)]
    x <- (x - min_x) / (max_x - min_x)
    x0 <- (x0 - min_x) / (max_x - min_x)

    args <- as.list(environment())
    if (is.null(FXm)) {
        fun <- .local_fit_cts_gaussian
    } else {
        fun <- .local_fit_discrete_gaussian
    }
    res <- do.call(fun, args[names(formals(fun))])
    eta_vals <- res[["eta_vals"]]

    # Estimated eta matrix
    Hhat <- do.call(rbind, eta_vals)
    # Estimated correlation matrix
    Rhat <- t(apply(Hhat, 1, function(v) {
        copula::P2p(vec2cor(v, scale = scale))
    }))

    ix_lab <- apply(utils::combn(seq(dim(FX)[2]), 2), 2, function(x) {
        paste(x, collapse = '')
    })
    colnames(Hhat) <- paste0("eta", ix_lab)
    colnames(Rhat) <- paste0("rho", ix_lab)

    res <- c(res[setdiff(names(res), "eta_vals")],
             list(x = x * (max_x - min_x) + min_x,
                  x0 = x0 * (max_x - min_x) + min_x,
                  eta = data.frame(Hhat),
                  rho = data.frame(Rhat),
                  scale = scale))
    return(res)
}

###############################################################################

#' Discrete Local Fit Gaussian
#'
#' @inheritParams local_fit_gaussian
#'
#' @return Estimated coefficients and information about optimization
.local_fit_discrete_gaussian <- function(FX,
                                         FXm,
                                         x,
                                         x0,
                                         band,
                                         scale,
                                         R0,
                                         optMethod,
                                         control,
                                         cores,
                                         cores2) {
    # Set up futures plan
    future::plan("multisession", workers = cores)

    # Convert to standard normal margins
    NX <- stats::qnorm(FX)
    NXm <- stats::qnorm(FXm)

    lbfgs_optim <- function(x0i) {
        # Use points nearby to estimate initial correlation matrix
        dx <- x - x0i
        max_thr <- sort(abs(dx))[min(length(dx), dim(NX)[2])]
        thr <- max(stats::quantile(abs(dx), R0), max_thr)
        NX_loc <- NX[which(abs(dx) <= thr), ]

        eta0 <- cor2vec(stats::cor(NX_loc, method = "pearson"), scale)
        par0 <- c(eta0, rep(0, length(eta0)))

        # Negative local log likelihood function
        loglik <- function(par) {
            par <- as.numeric(par)

            # Epanechnikov kernel weights
            wgt <- 3/(4 * band) * pmax(1 - (dx / band)^2, 0)

            # Only compute densities when weight is nonzero
            pos_ix <- which(wgt > 0)
            P <- sapply(pos_ix, function(j) {
                # Reconstruct correlation matrix
                npar <- as.integer(length(par) / 2)
                eta0 <- par[1:npar]
                eta1 <- par[(npar + 1):(2 * npar)]
                R <- vec2cor(eta0 + eta1 * dx[j], scale = scale)

                if (min(eigen(R)$values) <= 0) {
                    # Matrix is theoretically PD but may not be
                    # numerically PD
                    return(0)
                } else {
                    # Numerical evaluation of Gaussian CDF
                    return(tryCatch(
                        TruncatedNormal::mvNcdf(
                            l = NXm[j, ],
                            u = NX[j, ],
                            Sig = R,
                            n = 1e3
                        )$prob,
                        error = function(e) {
                            # Error often thrown for poorly conditioned
                            # matrices, return probability of 0
                            0
                        }
                    ))
                }
            })
            if (any(P <= 0)) {
                # Probabilities are theoretically nonnegative but may
                # be numerically negative. Return large number instead
                # of infinity.
                return(1e20)
            } else {
                return(-sum(wgt[pos_ix] * log(P)))
            }
        }

        # Create inner cluster
        cl <- parallel::makeCluster(cores2)
        on.exit(parallel::stopCluster(cl))
        parallel::clusterExport(
            cl = cl,
            varlist = c("vec2cor", "scale", "band", "NX", "NXm")
        )

        # Parallelized L-BFGS optimization
        res <- optimParallel::optimParallel(
            par = par0,
            fn = loglik,
            parallel = list(cl = cl),
            control = control
        )
        pbar()
        return(list(
            par = res$par[1:length(eta0)],
            convergence = res$convergence
        ))
    }

    sgd_safe <- function(x0i) {
        # Use points nearby to estimate initial correlation matrix
        dx <- x - x0i
        max_thr <- sort(abs(dx))[min(length(dx), dim(NX)[2])]
        thr <- max(stats::quantile(abs(dx), R0), max_thr)
        NX_loc <- NX[which(abs(dx) <= thr), ]

        eta0 <- cor2vec(stats::cor(NX_loc, method = "pearson"), scale)
        par0 <- c(eta0, rep(0, length(eta0)))

        fit_discrete <- py_load()$fit_discrete
        res <- fit_discrete(par0, dx, NX, NXm, scale, band, control)
        loss <- res$hist

        if (res$convergence == 2) {
            message("Failed prematurely")
            par1 <- res$par[1:length(eta0)]
            n_itr <- max(control$maxit - length(res$hist), 1L)
            message("Running ", n_itr, " additional iterations")

            loglik_lbfgs <- function(par) {
                par <- as.numeric(par)

                # Epanechnikov kernel weights
                wgt <- 3/(4 * band) * pmax(1 - (dx / band)^2, 0)

                # Only compute densities when weight is nonzero
                pos_ix <- which(wgt > 0)
                P <- sapply(pos_ix, function(j) {
                    # Reconstruct correlation matrix
                    npar <- as.integer(length(par) / 2)
                    eta0 <- par[1:npar]
                    eta1 <- par[(npar + 1):(2 * npar)]
                    R <- vec2cor(eta0 + eta1 * dx[j], scale = scale)

                    if (min(eigen(R)$values) <= 0) {
                        # Matrix is theoretically PD but may not be
                        # numerically PD
                        return(0)
                    } else {
                        # Numerical evaluation of Gaussian CDF
                        return(tryCatch(
                            TruncatedNormal::mvNcdf(
                                l = NXm[j, ],
                                u = NX[j, ],
                                Sig = R,
                                n = 1e3
                            )$prob,
                            error = function(e) {
                                # Error often thrown for poorly conditioned
                                # matrices, return probability of 0
                                0
                            }
                        ))
                    }
                })
                if (any(P <= 0)) {
                    # Probabilities are theoretically nonnegative but may
                    # be numerically negative. Return large number instead
                    # of infinity.
                    return(1e20)
                } else {
                    return(-sum(wgt[pos_ix] * log(P)))
                }
            }
            res <- stats::optim(
                par = par1,
                fn = loglik_lbfgs,
                method = "L-BFGS-B",
                control = list(maxit = n_itr)
            )
        }
        pbar()

        return(list(
            par = res$par[1:length(eta0)],
            convergence = res$convergence,
            loss = loss
        ))
    }

    progressr::with_progress({
        pbar <- progressr::progressor(along = x0)
        if (optMethod == "L-BFGS") {
            opt_res <- future.apply::future_lapply(
                X = x0,
                FUN = lbfgs_optim,
                future.seed = TRUE,
                future.packages = c("parallel", "optimParallel")
            )
            info <- NULL
        } else if (optMethod == "SGD") {
            opt_res <- future.apply::future_lapply(
                X = x0,
                FUN = sgd_safe,
                future.seed = TRUE
            )
            info <- list(loss = lapply(opt_res, '[[', "loss"))
        }
        return(c(list(eta_vals = lapply(opt_res, '[[', "par"),
                      convergence = sapply(opt_res, '[[', "convergence")),
                 info))
    })
}

###############################################################################

#' Continuous Local Fit Gaussian
#'
#' @inheritParams local_fit_gaussian
#'
#' @return Estimated coefficients and information about optimization
.local_fit_cts_gaussian <- function(FX,
                                    x,
                                    x0,
                                    band,
                                    scale,
                                    R0,
                                    optMethod,
                                    control,
                                    cores,
                                    cores2)  {
    # Set up futures plan
    future::plan("multisession", workers = cores)

    # Convert to standard normal margins
    NX <- stats::qnorm(FX)

    lbfgs_optim <- function(x0i) {
        # Use points nearby to estimate initial correlation matrix
        dx <- x - x0i
        max_thr <- sort(abs(dx))[min(length(dx), dim(NX)[2])]
        thr <- max(stats::quantile(abs(dx), R0), max_thr)
        NX_loc <- NX[which(abs(dx) <= thr), ]

        eta0 <- cor2vec(stats::cor(NX_loc, method = "pearson"), scale)
        par0 <- c(eta0, rep(0, length(eta0)))

        # Negative local log likelihood function
        loglik <- function(par) {
            par <- as.numeric(par)

            # Epanechnikov kernel weights
            wgt <- 3/(4 * band) * pmax(1 - (dx / band)^2, 0)

            # Only compute densities when weight is nonzero
            pos_ix <- which(wgt > 0)
            P <- sapply(pos_ix, function(j) {
                # Reconstruct correlation matrix
                npar <- as.integer(length(par) / 2)
                eta0 <- par[1:npar]
                eta1 <- par[(npar + 1):(2 * npar)]
                R <- vec2cor(eta0 + eta1 * dx[j], scale = scale)

                if (min(eigen(R)$values) <= 0) {
                    # Matrix is theoretically PD but may not be
                    # numerically PD
                    return(0)
                } else {
                    return(mvtnorm::dmvnorm(x = NX[j, ], sigma = R))
                }
            })

            if (any(P <= 0)) {
                # Probabilities are theoretically nonnegative but may
                # be numerically negative. Return large number instead
                # of infinity.
                return(1e20)
            } else {
                return(-sum(wgt[pos_ix] * log(P)))
            }
        }

        # Create inner cluster
        cl <- parallel::makeCluster(cores2)
        on.exit(parallel::stopCluster(cl))
        parallel::clusterExport(
            cl = cl,
            varlist = c("vec2cor", "scale", "band", "NX")
        )

        # Parallelized L-BFGS optimization
        opt <- optimParallel::optimParallel(
            par = par0,
            fn = loglik,
            parallel = list(cl = cl),
            control = control
        )
        pbar()
        return(list(par = opt$par[1:length(eta0)],
                    convergence = opt$convergence))
    }

    sgd_optim <- function(x0i) {
        # Use points nearby to estimate initial correlation matrix
        dx <- x - x0i
        max_thr <- sort(abs(dx))[min(length(dx), dim(NX)[2])]
        thr <- max(stats::quantile(abs(dx), R0), max_thr)
        NX_loc <- NX[which(abs(dx) <= thr), ]

        eta0 <- cor2vec(stats::cor(NX_loc, method = "pearson"), scale)
        par0 <- c(eta0, rep(0, length(eta0)))

        fit_continuous <- py_load()$fit_continuous
        res <- fit_continuous(par0, dx, NX, scale, band, control)
        pbar()

        return(list(par = res$par[1:length(eta0)],
                    loss = res$hist,
                    convergence = res$convergence))
    }

    progressr::with_progress({
        pbar <- progressr::progressor(along = x0)
        if (optMethod == "L-BFGS") {
            res <- future.apply::future_lapply(
                X = x0,
                FUN = lbfgs_optim,
                future.seed = TRUE,
                future.packages = c("parallel", "optimParallel")
            )
            info <- NULL
        } else if (optMethod == "SGD") {
            res <- future.apply::future_lapply(
                X = x0,
                FUN = sgd_optim,
                future.seed = TRUE
            )
            info <- list(loss = lapply(res, '[[', "loss"))
        }
        return(c(list(eta_vals = lapply(res, '[[', "par"),
                      convergence = sapply(res, '[[', "convergence")),
                 info))
    })
}

###############################################################################

#' Log-likelihood
#'
#' Compute log-likelihood
#'
#' @param FX Matrix of pseudo-observations F_i(X_i).
#' @param FXm Matrix of left limits of distribution functions F_i(X_i - 1).
#' If \code{NULL}, it is assumed that all margins are continuous.
#' @param R Correlation matrix.
#'
#' @return Log-likelihood
#'
#' @export
loglik <- function(FX, FXm = NULL, R)  {
    NX <- stats::qnorm(FX)
    if (is.null(FXm)) {
        ll <- apply(NX, 1, function(x) { mvtnorm::dmvnorm(x = x, sigma = R) })
    } else {
        NXm <- stats::qnorm(FXm)
        ll <- sapply(seq(dim(FX)[1]), function(i) {
            TruncatedNormal::mvNcdf(
                l = NXm[i, ],
                u = NX[i, ],
                Sig = R,
                n = 1e3
            )$prob
        })
    }
    return(sum(ll))
}

###############################################################################

#' Estimate correlation coefficients
#'
#' Estimate correlation coefficients using smoothing splines fit in the
#' unconstrained space.
#'
#' @param res Output of \code{local_fit_gaussian}.
#' @param x Times at which to estimate coefficients.
#' @param df Number of degrees of freedom for spline fitting. If \code{NULL}
#' (the default), leave-one-out cross validation is used instead.
#'
#' @return Predicted correlation coefficients.
#'
#' @export
predict_rho <- function(res,
                        x = res$x,
                        df = NULL) {
    # Predict eta values using smooth splines
    eta_pred <- apply(res$eta, 2, function(y) {
        if (is.null(df)) {
            fit <- stats::smooth.spline(res$x0, y, cv = TRUE)
        } else {
            fit <- stats::smooth.spline(res$x0, y, df = df)
        }
        return(stats::predict(fit, x)$y)
    })

    # Convert to correlation matrices
    v2c <- function(eta) { vec2cor(eta, scale = res$scale) }
    R_pred <- apply(eta_pred, 1, v2c, simplify = FALSE)

    # Extract time series for each coefficient
    ix <- which(lower.tri(R_pred[[1]]), arr.ind = TRUE)
    rho <- apply(ix, 1, function(v) {
        sapply(R_pred, function(R) { R[v[1], v[2]] })
    })
    return(rho)
}

###############################################################################
